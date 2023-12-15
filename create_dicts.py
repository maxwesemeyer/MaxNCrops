import matplotlib.pyplot as plt
import numpy as np
import rasterio.mask

from __functions import *


def prepare_data(temp_path, out_path, agg_length=100, crop_type_column=None, farm_id_column=None):

    shp_p = './' + temp_path + '/' + 'iacs.shp'
    crop_arr = gdal.Open('./' + temp_path + '/' + 'IDKTYP.tif').ReadAsArray().astype(int)
    ds_fid = gdal.Open('./' + temp_path + '/' + 'Field_ID.tif')
    gt_fid = ds_fid.GetGeoTransform()
    field_id_arr = ds_fid.ReadAsArray().astype(int)

    farmid_arr = gdal.Open('./' + temp_path + '/' + 'Farm_ID.tif').ReadAsArray().astype(int)

    sparse_idx = get_indices_sparse(field_id_arr.astype(int))

    unique_crops = np.unique(crop_arr).astype(int).tolist()
    if verbatim:
        print('number of unique crop types:', len(unique_crops), 'crops: ', unique_crops, 'including no data value')
    unique_farms = np.unique(farmid_arr)
    if verbatim:
        print('number of farms', len(unique_farms))

    iacs_gp = gpd.read_file(shp_p)
    iacs_gp['area_m2'] = iacs_gp.area
    unique_field_ids = np.unique(iacs_gp['field_id'].tolist()).astype(float).tolist()
    unique_field_ids.insert(0, 0)
    # extract starting values from initial situation
    start_vals = []
    for crop in unique_crops:
        if crop == 0:
            continue
        start_vals_crop = []
        for i, feat in iacs_gp.iterrows():
            if feat[crop_type_column] == crop:
                start_vals_crop.append(1)
            else:
                start_vals_crop.append(0)
        start_vals.append(start_vals_crop)

    field_id_arr_1d = field_id_arr.flatten()
    field_id_arr_1d[np.where(field_id_arr_1d == nd_value)] = 0

    #########################################################################
    if not os.path.isfile('./' + temp_path + '/' + 'shares_iacs.csv'):
        dissolved_shares_farm_crop = iacs_gp.copy().dissolve(by=[farm_id_column, crop_type_column])
        dissolved_shares_farm_crop['area_m2'] = dissolved_shares_farm_crop.area
        dissolved_shares_farm_crop.to_csv('./' + temp_path + '/' + 'shares_iacs.csv')
    #########################################################################
    # we get the crop types for each field that have not been on that field in the past
    # this will be used to makes constraints in the optimization later
    if not os.path.isfile('./' + temp_path + '/' + 'historic_croptypes_dict.pkl'):
        if verbatim:
            print('Using as historic croptypes: ', glob.glob('./input/*.tif')[0])
        ds_hist = gdal.Open(glob.glob('./input/*.tif')[0])
        gt_hist = ds_hist.GetGeoTransform()
        if verbatim:
            print('checking geotransforms of historic croptypes and rasterized data')
            print(gt_hist, gt_fid, gt_hist[0] == gt_fid[0], gt_hist[3] == gt_fid[3])
        if not gt_hist[0] == gt_fid[0] or gt_hist[3] == gt_fid[3]:
            #######################################################################################
            # Path to the reference raster for the extent
            reference_raster_path = './' + temp_path + '/' + 'reference_raster.tif'
            # Open the reference raster to get the extent
            with rasterio.open(reference_raster_path) as ref_raster:
                bounds = ref_raster.bounds

            # Open the raster to be clipped
            with rasterio.open(glob.glob('./input/*.tif')[0]) as src:
                from rasterio.windows import from_bounds
                # Clip the raster to the extent of the reference raster
                window = from_bounds(bounds.left, bounds.bottom, bounds.right, bounds.top, src.transform)

                # Read the data using the window
                out_image = src.read(window=window, boundless=True)
                # Update metadata with new dimensions and transform
                out_meta = src.meta.copy()
                out_meta.update({
                    'driver': 'GTiff',
                    'height': out_image.shape[1],
                    'width': out_image.shape[2],
                    'transform': src.window_transform(window)
                })

                # Save the clipped raster
                output_raster_path = './' + temp_path + '/' + 'hist_croptypes_adapted.tif'
                with rasterio.open(output_raster_path, 'w', **out_meta) as dest:
                    dest.write(out_image)
            #######################################################################################
            ds_hist = gdal.Open('./' + temp_path + '/' + 'hist_croptypes_adapted.tif')
            hist_croptypes = ds_hist.ReadAsArray()
        else:

            hist_croptypes = ds_hist.ReadAsArray()

            if hist_croptypes.shape[1] - field_id_arr.shape[0] > 0:
                # slice array to match the shape of the other rasters
                diff_abs = abs(field_id_arr.shape[0] - hist_croptypes.shape[1])
                new_length = hist_croptypes.shape[1] - diff_abs
                hist_croptypes = hist_croptypes[:, :new_length, :]
            elif hist_croptypes.shape[1] - field_id_arr.shape[0] < 0:
                # add emtpy cols and rows to hist croptypes
                add = field_id_arr.shape[0] - hist_croptypes.shape[1]
                hist_croptypes = np.pad(hist_croptypes, ((0, 0), (0, add), (0, 0)), 'constant', constant_values=(0))
            if hist_croptypes.shape[2] - field_id_arr.shape[1] > 0:
                # slice array to match the shape of the other rasters
                diff_abs = abs(field_id_arr.shape[1] - hist_croptypes.shape[2])
                new_length = hist_croptypes.shape[2] - diff_abs
                hist_croptypes = hist_croptypes[:, :, :new_length]
            elif hist_croptypes.shape[2] - field_id_arr.shape[1] < 0:
                # add emtpy cols and rows to hist croptypes
                add = field_id_arr.shape[1] - hist_croptypes.shape[2]
                hist_croptypes = np.pad(hist_croptypes, ((0, 0), (0, 0), (0, add)), 'constant', constant_values=(0))
            write_array_disk_universal(hist_croptypes, './' + temp_path + '/' + 'reference_raster.tif',
                                       outPath='./' + temp_path + '/' + 'hist_croptypes_adapted',
                                       dtype=gdal.GDT_Int16, scaler=1, adapt_pixel_size=False)
        taboo_croptypes_dict, historic_croptypes_dict = get_historic_croptypes(field_id_array=field_id_arr, historic_croptypes_array=hist_croptypes, unique_croptypes=unique_crops)
        with open('./' + temp_path + '/' + 'historic_croptypes_dict.pkl', 'wb') as f:
            pickle.dump(historic_croptypes_dict, f)
        with open('./' + temp_path + '/' + 'taboo_croptypes_dict.pkl', 'wb') as f:
            pickle.dump(taboo_croptypes_dict, f)

    with open('./' + temp_path + '/' + 'historic_croptypes_dict.pkl', 'rb') as f:
        historic_croptypes_dict = pickle.load(f)
    # append the historic crop types to the shapefile to calculate the crop acreage per farm later
    hist_crops_list = []
    for fid in unique_field_ids:
        hist_crops_list.append(historic_croptypes_dict[fid])
    hist_df = pd.DataFrame(np.stack(hist_crops_list, axis=0))

    hist_df.columns = ['crp_yr_' + str(year) for year in range(len(hist_df.columns))]
    hist_df['field_id'] = unique_field_ids
    iacs_gp = pd.merge(iacs_gp, hist_df, on="field_id")

    n_years = len(historic_croptypes_dict[1])

    # Calculate the crop acreage for each farm and crop for each year
    # Here we have to make the following assumptions:
    #   a) each field belongs to the same farm throughout the entire time period
    #   b) the field geometries do not change in the entire period
    if not os.path.isfile('./' + temp_path + '/' + 'shares_iacs_seq.csv'):
        for year in range(n_years):
            if year == 0:
                dissolved_shares_farm_crop = iacs_gp.copy().dissolve(by=[farm_id_column, 'crp_yr_' + str(year)])
                dissolved_shares_farm_crop['area_m2'] = dissolved_shares_farm_crop.area
                dissolved_shares_farm_crop['year'] = year
                dissolved_shares_farm_crop = dissolved_shares_farm_crop.reset_index()
                dissolved_shares_farm_crop = dissolved_shares_farm_crop[[farm_id_column, 'crp_yr_' + str(year), 'area_m2', 'year']]
                dissolved_shares_farm_crop.columns = [farm_id_column, crop_type_column, 'area_m2', 'year']
            else:
                dissolved_shares_farm_crop_iter = iacs_gp.copy().dissolve(by=[farm_id_column, 'crp_yr_' + str(year)])
                dissolved_shares_farm_crop_iter['area_m2'] = dissolved_shares_farm_crop_iter.area
                dissolved_shares_farm_crop_iter['year'] = year
                dissolved_shares_farm_crop_iter = dissolved_shares_farm_crop_iter.reset_index()
                dissolved_shares_farm_crop_iter = dissolved_shares_farm_crop_iter[
                    [farm_id_column, 'crp_yr_' + str(year), 'area_m2', 'year']]
                dissolved_shares_farm_crop_iter.columns = [farm_id_column, crop_type_column, 'area_m2', 'year']
                dissolved_shares_farm_crop = pd.concat([dissolved_shares_farm_crop, dissolved_shares_farm_crop_iter], ignore_index=True)

        dissolved_shares_farm_crop.to_csv('./' + temp_path + '/' + 'shares_iacs_seq.csv')
    ####################################################################################################################
    # creates a dictionary which states which indices of a one-dimensional array belong to a certain block (landscape)
    count_pixel_per_block = agg_length ** 2
    num_block_col = int(field_id_arr.shape[0] / agg_length)
    num_block_row = int(field_id_arr.shape[1] / agg_length)
    num_blocks_total = int(num_block_row * num_block_col)

    if not os.path.isfile('./' + temp_path + '/' + 'spatial_aggregation_dict.pkl'):
        if verbatim:
            print(num_blocks_total)
        agg_dicts = []
        block_counter = 0
        for block in range(num_blocks_total):

            for i in range(0, agg_length * num_block_row, agg_length):
                if block_counter == num_blocks_total:
                    break
                indices_list = []
                for a in range(agg_length):

                    if block_counter >= num_block_row:
                        start = i + (a * agg_length) * num_block_row + ((count_pixel_per_block * num_block_row) * block)
                        end = (i + agg_length) + (a * agg_length) * num_block_row + (
                                (count_pixel_per_block * num_block_row) * block)
                    elif block_counter < num_block_row:
                        # for the first blocks
                        start = i + (a * agg_length) * num_block_row
                        end = (i + agg_length) + (a * agg_length) * num_block_row

                    indices = list(range(start, end))
                    indices_list.append(indices)

                agg_dicts.append([item for sublist in indices_list for item in sublist])
                block_counter += 1
        del block_counter
        # 10x10 blocks
        spatial_aggregration_dictionary = dict(zip(list(range(num_blocks_total)), agg_dicts))

        with open('./' + temp_path + '/' + 'spatial_aggregation_dict.pkl', 'wb') as f:
            pickle.dump(spatial_aggregration_dictionary, f)
        del indices_list, agg_dicts

    with open('./' + temp_path + '/' + 'spatial_aggregation_dict.pkl', 'rb') as f:
        spatial_aggregration_dictionary = pickle.load(f)

    ####################################################################################################################
    # create farm_field_dict
    # the shape is [farm1 [field1, field2...], farm2 [field1, field2...] ... ]
    #
    if not os.path.isfile('./' + temp_path + '/' + 'farm_field_dict.pkl'):
        out_shares_farms = []
        for ct, farm in enumerate(unique_farms):
            if farm == 0:
                out_shares_farms.append(0)
                continue
            if verbatim:
                print(ct, 'of', len(unique_farms), farm)
            out_shares = []
            iacs_gp_subset = iacs_gp[iacs_gp[farm_id_column] == farm]
            field_ids_iter = iacs_gp_subset['field_id'].tolist()
            for field_id in unique_field_ids:
                if field_id in field_ids_iter:
                    subset_field = iacs_gp_subset[iacs_gp_subset['field_id'] == field_id]
                    out_shares.append(subset_field['area_m2'].values[0])
                else:
                    out_shares.append(0)
            if sum(out_shares) == 0:
                print('error farm ', farm, out_shares_farms)
            out_shares_farms.append(csr_matrix(out_shares))
        farm_field_dict = dict(zip(unique_farms, out_shares_farms))
        with open('./' + temp_path + '/' + 'farm_field_dict.pkl', 'wb') as f:
            pickle.dump(farm_field_dict, f)

        # the shape is [farm1 [field1, field2...], farm2 [field1, field2...] ... ]
        #

    ####################################################################################################################
    # create block shares dict
    # block id: [share of fieldid1 in this block, share of fieldid2 in this block,
    # share of fieldid3 in this block, (=number of decision units)...]
    if not os.path.isfile('./' + temp_path + '/' + 'block_dict.pkl'):
        out_shares = []
        for block_i in range(num_blocks_total):
            if verbatim:
                print(block_i, 'of', num_blocks_total)
            # 10x10 m Pixels -> have to multiply by 10*10 to convert to square meters
            ids = field_id_arr_1d[spatial_aggregration_dictionary[block_i]]
            collect = collections.Counter(ids)
            block_values_list = []
            for id in unique_field_ids:

                if id in collect.keys():
                    block_values_list.append(collect[id])

                else:
                    block_values_list.append(0)
            out_shares.append(csr_matrix(block_values_list))
        block_dict = dict(zip(list(range(num_blocks_total)), out_shares))
        del spatial_aggregration_dictionary, out_shares

        with open('./' + temp_path + '/' + 'block_dict.pkl', 'wb') as f:
            pickle.dump(block_dict, f)

    return field_id_arr, farmid_arr, sparse_idx, unique_crops, unique_field_ids, iacs_gp, unique_farms, num_blocks_total, \
           start_vals