import random
from math import sqrt, log
from osgeo import gdal
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import collections
import geopandas as gpd
from scipy.sparse import csr_matrix
import pickle
import os
import glob


def compute_M(data):
    cols = np.arange(data.size)
    return csr_matrix((cols, (data.ravel(), cols)),
                      shape=(data.max() + 1, data.size))


def get_indices_sparse(data):
    M = compute_M(data)
    return [np.unravel_index(row.data, data.shape) for row in M]


def write_array_disk_universal(array, reference_raster, outPath, scaler=1, out_file_type='.tif', edit_meta=False,
                               noDataValue=-999, dtype=gdal.GDT_Int16, adapt_pixel_size=False, adapted_pixel_size=None):
    """
    outpath should be without the file ending
    """
    print('writing raster: ', outPath)
    nd_mask = np.where(array == noDataValue)
    array = array * scaler
    array[nd_mask] = noDataValue
    bands = array.shape[0]
    mowSumRst = gdal.Open(reference_raster)
    drvMemR = gdal.GetDriverByName('MEM')

    RasterXSize = array.shape[2]
    RasterYSize = array.shape[1]
    mowDoy_ds = drvMemR.Create('', RasterXSize, RasterYSize, bands, dtype)
    if adapt_pixel_size:
        reference_gt = mowSumRst.GetGeoTransform()
        new_gt = (
            reference_gt[0], float(adapted_pixel_size * reference_gt[1]), reference_gt[2], reference_gt[3],
            reference_gt[4],
            float(adapted_pixel_size * reference_gt[5]))
        mowDoy_ds.SetGeoTransform(new_gt)
    else:
        mowDoy_ds.SetGeoTransform(mowSumRst.GetGeoTransform())
    mowDoy_ds.SetProjection(mowSumRst.GetProjection())

    for b in range(bands):
        mowDoy_ds.GetRasterBand(b + 1).WriteArray(array[b, :, :])
        mowDoy_ds.GetRasterBand(b + 1).SetNoDataValue(noDataValue)
    if out_file_type == '.tif':
        driver = gdal.GetDriverByName("GTiff")
    elif out_file_type == '.asc':
        driver = gdal.GetDriverByName("AAIGrid")
    copy_ds = driver.CreateCopy(outPath + out_file_type, mowDoy_ds, 0, ['COMPRESS=LZW'])
    copy_ds = None
    mowSumRst = None


def get_entropy(map_2d, agg_len, return_img=False, return_count=False, return_agr_area=False):

    side_length_y = map_2d.shape[0]
    side_length_x = map_2d.shape[1]
    side_length_block = agg_len
    total_n_pixel_block = agg_len ** 2
    entropy_list = []
    n_unique_list = []
    area_list = []
    for y in range(side_length_block, side_length_y + 1, side_length_block):
        for x in range(side_length_block, side_length_x + 1, side_length_block):
            block = map_2d[(y - side_length_block):y, (x - side_length_block):x]

            unique, counts = np.unique(block, return_counts=True)
            n_pixel_agr = np.where(block > 0, True, False).sum()
            area_ha = (n_pixel_agr * 100) * 0.0001

            area_list.append(area_ha)
            if unique[0] == 0:
                unique = unique[1:]
                counts = counts[1:]

            n_unique = len(unique)
            n_unique_list.append(n_unique)

            prob = counts / counts.sum()
            entropy = -np.sum(prob * np.log(prob))

            if area_ha == 0:
                entropy_list.append(-999)
                continue
            #counts = counts / total_n_pixel_block
            #entropy = -sum(np.array(counts) * np.log(counts))
            entropy_list.append(entropy)
            # -sum( p*log(p) )
    if return_agr_area:
        img_2d = np.reshape(area_list, (int(side_length_y / agg_len), int(side_length_x / agg_len)))
        return area_list, img_2d
    if return_count:
        img_2d = np.reshape(n_unique_list, (int(side_length_y / agg_len), int(side_length_x / agg_len)))
        return n_unique_list, img_2d

    if return_img:
        img_2d = np.reshape(entropy_list, (int(side_length_y / agg_len), int(side_length_x / agg_len)))
        return sum(entropy_list), img_2d
    else:
        return entropy_list


def get_shares_farm(farm_id_arr, crop_id_arr):
    # iterate through farms and crop types; get the share of each crop type in each farm
    unique_crops = np.unique(crop_id_arr)
    shares_list = []
    crop_list = []
    farm_list = []

    sparse_idx = get_indices_sparse(farm_id_arr.astype(int))
    for i, farm_id in enumerate(sparse_idx):
        if farm_id == -999:
            continue
        farm_mask = farm_id
        # farm_mask = np.where(farm_id_arr == farm_id, True, False)
        farm_crop_id_arr = crop_id_arr[farm_mask]
        collect = collections.Counter(farm_crop_id_arr)

        for crop in unique_crops:
            if crop == -999:
                continue
            # Only true if we used 10m Pixels
            crop_share_farm_m2 = collect[crop] * 100
            shares_list.append(crop_share_farm_m2)
            crop_list.append(crop)
            farm_list.append(i)
    return farm_list, crop_list, shares_list


def get_historic_croptypes_artificialdata(unique_field_ids):
    # get the crop types that have not been cultivated on a field in the past
    historic_croptypes = []
    for id in unique_field_ids:
        # create a random number of crops between 1 and 10 that cannot be cultivated
        historic_croptypes.append(np.random.randint(1, 12, size=np.random.randint(1, 10, size=1)[0]))
    historic_croptypes_dict = dict(zip(unique_field_ids, historic_croptypes))
    return historic_croptypes_dict


def get_historic_croptypes(field_id_array, historic_croptypes_array, unique_croptypes):
    # get the crop types that have not been cultivated on a field in the past
    print('creating historic crop type dictionary...')
    taboo_crops_list = []
    unique_id_list = []
    n_bands = historic_croptypes_array.shape[0]
    historic_croptypes_array = np.reshape(historic_croptypes_array, newshape=(historic_croptypes_array.shape[0], historic_croptypes_array.shape[1]*historic_croptypes_array.shape[2]))

    # Aggregate crop types to coarser classes
    #historic_croptypes_array[np.where((historic_croptypes_array==9) | (historic_croptypes_array==10) | (historic_croptypes_array==7) )] = 2
    #historic_croptypes_array[np.where((historic_croptypes_array==14))] = 12

    field_id_array_rav = field_id_array.ravel()

    sparse_idx = get_indices_sparse(field_id_array_rav)

    for i, id in enumerate(sparse_idx):
        if id[0].size == 0:
            continue
        # create a mask for the id
        if id == -999:
            unique_id_list.append(0)
        else:
            unique_id_list.append(float(i))
        mask = id
        # iterate through all bands to get only the majority class per year
        historic_cultivations = []

        for year in range(n_bands):
            historic, counts = np.unique(historic_croptypes_array[year, mask], return_counts=True)
            # In case of a very small field it is possible that the 0 class has the highest count; In this case the second
            # highest value is chosen
            if historic[np.argmax(counts)] == 0:
                historic_cultivations.append(historic[np.argsort(counts)[-2]])
            else:
                historic_cultivations.append(historic[np.argmax(counts)])

        taboo_crops = []

        for x in unique_croptypes:
            if x not in historic_cultivations and x != 0:
                taboo_crops.append(float(x))
        # unique croptypes contains 0 taboo crops cannot contain 0 -> len - 1
        if len(taboo_crops) >= len(unique_croptypes)-1:
            print(i, id, historic_cultivations)
            print('Error; setting no crops as taboo', taboo_crops, unique_croptypes, i)
            #taboo_crops = taboo_crops[:-1]
            taboo_crops = []
        taboo_crops_list.append(taboo_crops)

    historic_croptypes_dict = dict(zip(unique_id_list, taboo_crops_list))
    return historic_croptypes_dict
