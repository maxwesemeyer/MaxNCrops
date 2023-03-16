from __functions import *


def prepare_data(agg_length=100, crop_type_column=None, farm_id_column=None):

    shp_p = './temp/iacs.shp'
    crop_arr = gdal.Open('./temp/IDKTYP.tif').ReadAsArray().astype(int)
    field_id_arr = gdal.Open('./temp/Field_ID.tif').ReadAsArray().astype(int)

    farmid_arr = gdal.Open('./temp/Farm_ID.tif').ReadAsArray().astype(int)

    sparse_idx = get_indices_sparse(field_id_arr.astype(int))

    unique_crops = np.unique(crop_arr).astype(int).tolist()

    print('number of unique crop types:', len(unique_crops), 'crops: ', unique_crops, 'including no data value')
    unique_farms = np.unique(farmid_arr)
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
    field_id_arr_1d[np.where(field_id_arr_1d == -999)] = 0

    size_of_test_data = field_id_arr_1d.shape[0]

    len_raster = sqrt(size_of_test_data)
    num_block_y = int(len_raster / agg_length)
    num_blocks = int(num_block_y ** 2)

    if not os.path.isfile('./temp/shares_iacs.csv'):
        dissolved_shares_farm_crop = iacs_gp.copy().dissolve(by=[farm_id_column, crop_type_column])
        dissolved_shares_farm_crop['area_m2'] = dissolved_shares_farm_crop.area
        dissolved_shares_farm_crop.to_csv('./temp/shares_iacs.csv')

    ##########################################################
    # get the initial landscape shannon diversity, number of unique crops per landscape and the proportion of
    # agricultural area per landscape and write to disk
    if os.path.isfile('./output/initial_count.tif') and os.path.isfile('./output/initial_entropy.tif') and \
            os.path.isfile('./output/agricultural_area_ha.tif'):
        print('skipping entropy calc')

    else:
        total_entropy, img_ha = get_entropy(crop_arr, agg_len=100, return_agr_area=True)
        total_entropy, img_ct = get_entropy(crop_arr, agg_len=100, return_count=True)
        total_entropy, img_div = get_entropy(crop_arr, agg_len=100, return_img=True)

        write_array_disk_universal(np.expand_dims(img_div, axis=0), './temp/reference_raster.tif',
                                   outPath='./output/initial_entropy',
                                   dtype=gdal.GDT_Int32, noDataValue=0, scaler=1000, adapt_pixel_size=True,
                                   adapted_pixel_size=100)

        write_array_disk_universal(np.expand_dims(img_ct, axis=0), './temp/reference_raster.tif',
                                   outPath='./output/initial_count',
                                   dtype=gdal.GDT_Int32, noDataValue=0, scaler=1, adapt_pixel_size=True,
                                   adapted_pixel_size=100)

        write_array_disk_universal(np.expand_dims(img_ha, axis=0), './temp/reference_raster.tif',
                                   outPath='./output/agricultural_area_ha',
                                   dtype=gdal.GDT_Int32, noDataValue=0, scaler=1, adapt_pixel_size=True,
                                   adapted_pixel_size=100)

    #########################################################################
    # we get the crop types for each field that have not been on that field in the past
    # this will be used to makes constraints in the optimization later
    if not os.path.isfile('./temp/historic_croptypes_dict.pkl'):
        print('Using as historic croptypes: ', glob.glob('./input/*.tif')[0])
        hist_croptypes = gdal.Open(glob.glob('./input/*.tif')[0]).ReadAsArray()
        if hist_croptypes.shape[1] - field_id_arr.shape[0] > 0:
            print('not implemented; Historic croptypes extent does not match')
        elif hist_croptypes.shape[1] - field_id_arr.shape[0] < 0:
            # add emtpy cols and rows to hist croptypes
            add = field_id_arr.shape[0] - hist_croptypes.shape[1]
            hist_croptypes = np.pad(hist_croptypes, ((0, 0), (0, add), (0, 0)), 'constant', constant_values=(0))
        if hist_croptypes.shape[2] - field_id_arr.shape[1] > 0:
            print('not implemented; Historic croptypes extent does not match')
        elif hist_croptypes.shape[2] - field_id_arr.shape[1] < 0:
            add = field_id_arr.shape[1] - hist_croptypes.shape[2]
            hist_croptypes = np.pad(hist_croptypes, ((0, 0), (0, 0), (0, add)), 'constant', constant_values=(0))

        historic_croptypes_dict = get_historic_croptypes(field_id_array=field_id_arr, historic_croptypes_array=hist_croptypes, unique_croptypes=unique_crops)
        with open('./temp/historic_croptypes_dict.pkl', 'wb') as f:
            pickle.dump(historic_croptypes_dict, f)
    ####################################################################################################################
    # creates a dictionary which states which indices of a one-dimensional array belong to a certain block
    if not os.path.isfile('./temp/spatial_aggregation_dict.pkl'):
        print(num_blocks)
        count_pixel_per_block = agg_length ** 2
        agg_dicts = []
        block_counter = 0
        for row in range(num_block_y):
            print(row, num_block_y)
            for i in range(0, agg_length * num_block_y, agg_length):
                if block_counter == num_blocks:
                    break
                indices_list = []
                for a in range(agg_length):

                    if block_counter >= num_block_y:
                        start = i + (a * agg_length) * num_block_y + ((count_pixel_per_block * num_block_y) * row)
                        end = (i + agg_length) + (a * agg_length) * num_block_y + (
                                    (count_pixel_per_block * num_block_y) * row)
                    elif block_counter < num_block_y:
                        # for the first blocks
                        start = i + (a * agg_length) * num_block_y
                        end = (i + agg_length) + (a * agg_length) * num_block_y

                    indices = list(range(start, end))
                    indices_list.append(indices)

                agg_dicts.append([item for sublist in indices_list for item in sublist])
                block_counter += 1
        del block_counter
        # 10x10 blocks
        spatial_aggregration_dictionary = dict(zip(list(range(num_blocks)), agg_dicts))
        with open('./temp/spatial_aggregation_dict.pkl', 'wb') as f:
            pickle.dump(spatial_aggregration_dictionary, f)
        del indices_list, agg_dicts

    with open('./temp/spatial_aggregation_dict.pkl', 'rb') as f:
        spatial_aggregration_dictionary = pickle.load(f)

    ####################################################################################################################
    # create farm_field_dict
    # the shape is [farm1 [field1, field2...], farm2 [field1, field2...] ... ]
    #
    if not os.path.isfile('./temp/farm_field_dict.pkl'):
        out_shares_farms = []
        for ct, farm in enumerate(unique_farms):
            if farm == 0:
                out_shares_farms.append(0)
                continue
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
        with open('./temp/farm_field_dict.pkl', 'wb') as f:
            pickle.dump(farm_field_dict, f)

        with open('./temp/farm_field_dict.pkl', 'rb') as f:
            farm_field_dict = pickle.load(f)
        # the shape is [farm1 [field1, field2...], farm2 [field1, field2...] ... ]
        #

    ####################################################################################################################
    # create block shares dict
    # block id: [share of fieldid1 in this block, share of fieldid2 in this block,
    # share of fieldid3 in this block, (=number of decision units)...]
    if not os.path.isfile('./temp/block_dict.pkl'):

        out_shares = []
        for block_i in range(num_blocks):
            print(block_i, 'of', num_blocks)
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
        block_dict = dict(zip(list(range(num_blocks)), out_shares))
        del spatial_aggregration_dictionary, out_shares

        with open('./temp/block_dict.pkl', 'wb') as f:
            pickle.dump(block_dict, f)

        with open('./temp/block_dict.pkl', 'rb') as f:
            block_dict = pickle.load(f)
        print(len(block_dict))

    return crop_arr, field_id_arr, farmid_arr, sparse_idx, \
           unique_crops, unique_field_ids, iacs_gp, unique_farms, field_id_arr_1d, size_of_test_data, len_raster, \
           num_blocks, num_block_y, start_vals