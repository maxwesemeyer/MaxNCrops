import random
from math import sqrt, log, ceil
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
import os.path
import rasterio
from rasterio import features
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from itertools import chain
from config import *
import itertools


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
    if verbatim:
        print('writing raster: ', outPath)
    nd_mask = np.where(array == noDataValue)
    array = array * scaler
    array[nd_mask] = noDataValue
    bands = array.shape[0]
    outraster = gdal.Open(reference_raster)
    drvMemR = gdal.GetDriverByName('MEM')

    RasterXSize = array.shape[2]
    RasterYSize = array.shape[1]
    outds = drvMemR.Create('', RasterXSize, RasterYSize, bands, dtype)
    if adapt_pixel_size:
        reference_gt = outraster.GetGeoTransform()
        new_gt = (
            reference_gt[0], float(adapted_pixel_size * reference_gt[1]), reference_gt[2], reference_gt[3],
            reference_gt[4],
            float(adapted_pixel_size * reference_gt[5]))
        outds.SetGeoTransform(new_gt)
    else:
        outds.SetGeoTransform(outraster.GetGeoTransform())
    outds.SetProjection(outraster.GetProjection())

    for b in range(bands):
        outds.GetRasterBand(b + 1).WriteArray(array[b, :, :])
        outds.GetRasterBand(b + 1).SetNoDataValue(noDataValue)
    if out_file_type == '.tif':
        driver = gdal.GetDriverByName("GTiff")
    elif out_file_type == '.asc':
        driver = gdal.GetDriverByName("AAIGrid")
    copy_ds = driver.CreateCopy(outPath + out_file_type, outds, 0, ['COMPRESS=LZW'])
    copy_ds = None
    outraster = None


def get_entropy(map_2d, agg_len, return_type=None):
    # return_type = area, count, Shannon diversity
    side_length_y = map_2d.shape[0]
    side_length_x = map_2d.shape[1]
    side_length_block = agg_len
    #total_n_pixel_block = agg_len ** 2
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
                entropy_list.append(nd_value)
                continue
            #counts = counts / total_n_pixel_block
            #entropy = -sum(np.array(counts) * np.log(counts))
            entropy_list.append(entropy)
            # -sum( p*log(p) )
    if return_type == 'area':
        img_2d = np.reshape(area_list, (int(side_length_y / agg_len), int(side_length_x / agg_len)))
        return area_list, img_2d
    if return_type == 'count':
        img_2d = np.reshape(n_unique_list, (int(side_length_y / agg_len), int(side_length_x / agg_len)))
        return n_unique_list, img_2d
    if return_type == 'Shannon diversity':
        img_2d = np.reshape(entropy_list, (int(side_length_y / agg_len), int(side_length_x / agg_len)))
        return sum(entropy_list), img_2d
    if return_type == '':
        # Returns shannon diversity as well but not as 2d np array
        return entropy_list


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
    historic_crops_list = []
    n_bands = historic_croptypes_array.shape[0]
    historic_croptypes_array = np.reshape(historic_croptypes_array, newshape=(historic_croptypes_array.shape[0], historic_croptypes_array.shape[1]*historic_croptypes_array.shape[2]))

    # Aggregate crop types to coarser classes
    # TODO hard coded crop classes; only needed when working with ID_KTYP_2
    historic_croptypes_array[np.where((historic_croptypes_array==9) | (historic_croptypes_array==10) | (historic_croptypes_array==7))] = 2
    historic_croptypes_array[np.where((historic_croptypes_array==14))] = 12
    #historic_croptypes_array[np.where((historic_croptypes_array==7))] = 6
    historic_croptypes_array[np.where((historic_croptypes_array==255))] = 0
    historic_croptypes_array[np.where((historic_croptypes_array==99))] = 0

    field_id_array_rav = field_id_array.ravel()

    sparse_idx = get_indices_sparse(field_id_array_rav)

    for i, id in enumerate(sparse_idx):
        if id[0].size == 0:
            continue
        # create a mask for the id
        if id == nd_value:
            unique_id_list.append(0)
        else:
            unique_id_list.append(float(i))
        mask = id
        # iterate through all bands to get only the majority class per year
        historic_cultivations = []

        for year in range(n_bands):
            historic, counts = np.unique(historic_croptypes_array[year, mask], return_counts=True)
            try:
                if i == 0:
                    historic_cultivations.append(historic[np.argsort(counts)[0]])
                else:
                    # TODO changed the function a bit
                    # In case of a very small field it is possible that the 0 class has the highest count; In this case the second
                    # highest value is chosen

                    #elif historic[np.argmax(counts)] == 0:
                    #    historic_cultivations.append(historic[np.argsort(counts)[-2]])
                    #else:
                    historic_cultivations.append(historic[np.argmax(counts)])
            except:
                historic_cultivations.append(historic[np.argmax(counts)])

        taboo_crops = []

        for x in unique_croptypes:
            if x not in historic_cultivations and x != 0:
                taboo_crops.append(float(x))
        # unique croptypes contains 0 taboo crops cannot contain 0 -> len - 1
        if len(taboo_crops) >= len(unique_croptypes)-1:
            #print(i, id, historic_cultivations)
            if i > 0:
                # i == 0 refers to the area that surrounds all fields; So in this case the warning can be ignored
                print('no data available; setting no crops as taboo for field', i, taboo_crops, unique_croptypes)
            taboo_crops = []
        taboo_crops_list.append(taboo_crops)
        historic_crops_list.append(historic_cultivations)
    taboo_croptypes_dict = dict(zip(unique_id_list, taboo_crops_list))
    historic_croptypes_dict = dict(zip(unique_id_list, historic_crops_list))
    return taboo_croptypes_dict, historic_croptypes_dict


def longest_sequence(binary_array):
    max_sequence = 0
    current_sequence = 0
    start_index = 0
    max_start_index = 0

    for i, value in enumerate(binary_array):
        if value == 1:
            current_sequence += 1
            if current_sequence == 1:
                start_index = i  # Update start index when a new sequence begins
            if current_sequence > max_sequence:
                max_sequence = current_sequence
                max_start_index = start_index
        else:
            current_sequence = 0

    return max_sequence, max_start_index


def shortest_sequence(binary_array):
    # This function measures the shortest sequence of 0s between 1s
    binary_array_rev = ~binary_array
    # in case there is only one 1, the length of the bin array is returned
    if sum(binary_array) == 1:
        return len(binary_array), 0

    if sum(binary_array_rev) == 0:
        return 0, 0
    if sum(binary_array_rev) == len(binary_array_rev):
        return sum(binary_array_rev), 0
    # this happens if the crop was cultivated in two consecutive years;
    if longest_sequence(~binary_array_rev)[0] > 1:
        return 0, longest_sequence(~binary_array_rev)[1]

    max_sequence = 0
    current_sequence = 0
    start_index = 0

    sequence_lengths = []
    start_indices = []
    for i, value in enumerate(binary_array_rev):
        # we are only interested in the gaps between 1s; this way we avoid measuring the gap starting at index 0
        # e.g. : [0, 1, 0, 0, 1]
        if value == 1 and sum(binary_array[:i]) >= 1:
            current_sequence += 1
            if current_sequence == 1:
                start_index = i  # Update start index when a new sequence begins
            max_sequence = current_sequence

        else:
            if current_sequence >= 1:
                sequence_lengths.append(max_sequence)
                start_indices.append(start_index)
            current_sequence = 0
    # this happens if the sequence of 1s goes until the last item
    if not sequence_lengths:
        return current_sequence, start_index
    argmin = np.argmin(sequence_lengths)
    return sequence_lengths[argmin], start_indices[argmin]
