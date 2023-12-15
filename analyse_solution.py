from __functions import *
import matplotlib.pyplot as plt


def analyse_solution_seq(temp_path, out_path, selected_farm_ids, landscape_size=agg_length):
    rst = rasterio.open('' + temp_path + '/' + 'reference_raster.tif')
    opt = gdal.Open('./' + out_path + '/' + 'maxent_croptypes_' + str(tolerance) + '.tif').ReadAsArray()
    n_years = opt.shape[0]

    ####################################################################################################################
    a, agr_area = get_entropy(opt[0, :, :], landscape_size, return_type='area')
    write_array_disk_universal(np.expand_dims(agr_area, axis=0), './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + temp_path + '/' + str(landscape_size) + 'reference_landscape',
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=landscape_size)

    iacs_orig = gpd.read_file('' + out_path + '/' + 'iacs_opt.shp')
    selected_rows = iacs_orig[iacs_orig[farm_id_column].isin(selected_farm_ids)].copy().dissolve()
    with rasterio.open('./' + temp_path + '/' + str(landscape_size) + 'reference_landscape.tif') as src:
        # Create a mask for the geometries within the bounds of the raster
        mask = geometry_mask(selected_rows.geometry, out_shape=src.shape, transform=src.transform, invert=False,
                             all_touched=True)

    mask_landscapes = np.tile(mask, reps=(n_years, 1, 1))
    #####################################################
    meta = rst.meta.copy()
    meta.update(compress='lzw')
    meta.update(count=n_years)
    iacs = gpd.read_file('./' + out_path + '/iacs_opt.shp')
    with rasterio.open('./' + out_path + '/init_crop_allocation.tif', 'w+', **meta) as out:
        for year in range(n_years):
            out.nodata = 0
            out_arr = out.read(year+1)
            # this is where we create a generator of geom, value pairs to use in rasterizing
            shapes = ((geom, value) for geom, value in zip(iacs.geometry, iacs['crp_yr_' + str(year)]))
            burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
            out.write_band(year+1, burned)
    init = gdal.Open('./' + out_path + '/' + 'init_crop_allocation.tif').ReadAsArray()
    init[np.where((init == 255) | (init == 99), True, False)] = 0
    mask = np.where(opt > 0, True, False)
    #init[~mask] = 0

    print('UNIQUE INIT: ', np.unique(init), 'UNIQUE OPT: ', np.unique(opt))
    #print('check the number of pixels: ', np.where(init > 0, True, False).sum(), np.where(opt > 0, True, False).sum())
    #####################################################
    img_entr_init_list = []
    img_entr_opt_list = []

    img_ct_init_list = []
    img_ct_opt_list = []
    agr_area_list = []

    for year in range(n_years):

        a, img_init_ct = get_entropy(init[year, :, :], landscape_size, return_type='count')
        a, img_opt_ct = get_entropy(opt[year, :, :], landscape_size, return_type='count')

        img_ct_init_list.append(img_init_ct)
        img_ct_opt_list.append(img_opt_ct)

        a, img_entr_init = get_entropy(init[year, :, :], landscape_size, return_type='Shannon diversity')
        a, img_entr_opt = get_entropy(opt[year, :, :], landscape_size, return_type='Shannon diversity')

        img_entr_init_list.append(img_entr_init)
        img_entr_opt_list.append(img_entr_opt)

        a, agr_area = get_entropy(opt[year, :, :], landscape_size, return_type='area')
        agr_area_list.append(agr_area)

    img_init = np.stack(img_ct_init_list, axis=0)
    img_init = img_init.astype(float)
    img_init[mask_landscapes] = np.nan

    img_opt = np.stack(img_ct_opt_list, axis=0)
    img_opt = img_opt.astype(float)
    img_opt[mask_landscapes] = np.nan

    img_entr_init = np.stack(img_entr_init_list, axis=0)
    img_entr_init = img_entr_init.astype(float)
    img_entr_init[mask_landscapes] = np.nan

    img_entr_opt = np.stack(img_entr_opt_list, axis=0)
    img_entr_opt = img_entr_opt.astype(float)
    img_entr_opt[mask_landscapes] = np.nan
    agr_area_repeated = np.stack(agr_area_list, axis=0)

    write_array_disk_universal(img_init, './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + out_path + '/' + str(landscape_size) + 'inital_ct',
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=landscape_size)

    write_array_disk_universal(img_opt,
                               './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + out_path + '/' + str(landscape_size) + 'opt_ct_' + str(tolerance),
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=landscape_size)

    write_array_disk_universal(img_entr_init,
                               './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + out_path + '/' + str(landscape_size) + 'initial_entropy',
                               dtype=gdal.GDT_Float32, noDataValue=nd_value, scaler=1, adapt_pixel_size=True,
                               adapted_pixel_size=landscape_size)

    write_array_disk_universal(img_entr_opt,
                                   './' + temp_path + '/' + 'reference_raster.tif',
                                   outPath='./' + out_path + '/' + str(landscape_size) + 'opt_entropy_' + str(tolerance),
                                   dtype=gdal.GDT_Float32, noDataValue=nd_value, scaler=1, adapt_pixel_size=True,
                                   adapted_pixel_size=landscape_size)
    pd.DataFrame(
        {'entropy_init': img_entr_init.ravel(), 'entropy_opt': img_entr_opt.ravel(), 'initial_ct': img_init.ravel(),
         'opt_ct': img_opt.ravel(), 'agr_area': agr_area_repeated.ravel()}).to_csv(
        './' + out_path + '/' + 'entropy_ct_rav' + str(tolerance) + '.csv')

    img_nan_init = img_entr_init.astype(float)
    img_nan_init[np.where(img_nan_init == nd_value, True, False)] = np.nan

    img_nan_opt = img_entr_opt.astype(float)
    img_nan_opt[np.where(img_nan_opt == nd_value, True, False)] = np.nan

    print('intial entropy avg: ', np.nanmean(img_nan_init.flatten()),
          'optimized avg: ', np.nanmean(img_nan_opt.flatten()),
          'difference abs: ', np.nanmean(img_nan_opt.flatten()) - np.nanmean(img_nan_init.flatten()),
          'Diff in percent:',
          ((np.nanmean(img_nan_opt.flatten()) - np.nanmean(img_nan_init.flatten())) / np.nanmean(
              img_nan_init.flatten())) * 100)


def analyse_solution(temp_path, out_path, selected_farm_ids):
    init = gdal.Open('./' + temp_path + '/' + 'IDKTYP.tif').ReadAsArray()
    opt = gdal.Open('./' + out_path + '/' + 'opt_crop_allocation_' + str(tolerance) + '.tif').ReadAsArray()
    ####################################################################################################################
    a, agr_area = get_entropy(opt, agg_length, return_type='area')
    write_array_disk_universal(np.expand_dims(agr_area, axis=0), './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + temp_path + '/' + 'reference_landscape',
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=agg_length)

    iacs_orig = gpd.read_file('' + out_path + '/' + 'iacs_opt.shp')
    selected_rows = iacs_orig[iacs_orig[farm_id_column].isin(selected_farm_ids)].copy().dissolve()
    with rasterio.open('./' + temp_path + '/' + 'reference_landscape.tif') as src:
        # Create a mask for the geometries within the bounds of the raster
        mask = geometry_mask(selected_rows.geometry, out_shape=src.shape, transform=src.transform, invert=False, all_touched=True)

    ####################################################################################################################
    a, img_init_ct = get_entropy(init, agg_length, return_type='count')
    img_init_ct = img_init_ct.astype(float)
    img_init_ct[mask] = np.nan
    a, img_opt_ct = get_entropy(opt, agg_length, return_type='count')
    img_opt_ct = img_opt_ct.astype(float)
    img_opt_ct[mask] = np.nan
    a, img_entr_init = get_entropy(init, agg_length, return_type='Shannon diversity')
    img_entr_init = img_entr_init.astype(float)
    img_entr_init[mask] = np.nan
    a, img_entr_opt = get_entropy(opt, agg_length, return_type='Shannon diversity')
    img_entr_opt = img_entr_opt.astype(float)
    img_entr_opt[mask] = np.nan
    a, agr_area = get_entropy(opt, agg_length, return_type='area')

    write_array_disk_universal(np.expand_dims(img_init_ct, axis=0), './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + out_path + '/' + 'inital_ct',
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=agg_length)

    write_array_disk_universal(np.expand_dims(img_opt_ct, axis=0),
                               './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + out_path + '/' + 'opt_ct_' + str(tolerance),
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=agg_length)

    write_array_disk_universal(np.expand_dims(img_entr_init, axis=0),
                               './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + out_path + '/' + 'initial_ShanDiv',
                               dtype=gdal.GDT_Float32, noDataValue=nd_value, scaler=1, adapt_pixel_size=True,
                               adapted_pixel_size=agg_length)

    write_array_disk_universal(np.expand_dims(img_entr_opt, axis=0),
                               './' + temp_path + '/' + 'reference_raster.tif',
                               outPath='./' + out_path + '/' + 'opt_ShanDiv_' + str(tolerance),
                               dtype=gdal.GDT_Float32, noDataValue=nd_value, scaler=1, adapt_pixel_size=True,
                               adapted_pixel_size=agg_length)
    pd.DataFrame(
        {'entropy_init': img_entr_init.ravel(), 'entropy_opt': img_entr_opt.ravel(), 'initial_ct': img_init_ct.ravel(),
         'opt_ct': img_opt_ct.ravel(), 'agr_area': agr_area.ravel()}).to_csv('./' + out_path + '/' + 'entropy_ct_rav' + str(tolerance) + '.csv')

    ##########################################################
    # get shares of inital crop shares and save them to a csv file
    img_nan_init = img_entr_init.astype(float)
    img_nan_init[np.where(img_nan_init == nd_value, True, False)] = np.nan

    img_nan_opt = img_entr_opt.astype(float)
    img_nan_opt[np.where(img_nan_opt == nd_value, True, False)] = np.nan

    print('intial Shannon diversity average: ', np.nanmean(img_nan_init.flatten()),
          'optimized Shannon diversity average: ', np.nanmean(img_nan_opt.flatten()),
          'difference absolute: ', np.nanmean(img_nan_opt.flatten())-np.nanmean(img_nan_init.flatten()),
          'difference in percent:',
          ((np.nanmean(img_nan_opt.flatten())-np.nanmean(img_nan_init.flatten()))/np.nanmean(img_nan_init.flatten()))*100)


def get_change_map(temp_path, out_path):
    iacs_gp = gpd.read_file('./' + out_path + '/' + 'iacs_opt.shp')
    print(iacs_gp.columns)
    iacs_gp['crp_chngd'] = iacs_gp['OPT_KTYP'] != iacs_gp[crop_type_column]
    number_of_fields = len(iacs_gp.index)
    number_of_changes = iacs_gp['crp_chngd'].sum()
    print('number of changes in crop rotation', iacs_gp['crp_chngd'].sum(), 'total fields:', number_of_fields)
    print((number_of_changes/number_of_fields)*100, '% of fields changes')
    iacs_gp.to_file('./' + out_path + '/' + 'iacs_opt.shp')


def get_change_map_seq(n_years, temp_path, out_path):
    for year in range(n_years):
        iacs_gp = gpd.read_file('./' + out_path + '/' + 'iacs_opt.shp')
        iacs_gp['crp_chngd_' + str(year)] = iacs_gp['OPT_KTYP_' + str(year)] != iacs_gp['crp_yr_' + str(year)]
        number_of_fields = len(iacs_gp.index)
        number_of_changes = iacs_gp['crp_chngd_' + str(year)].sum()
        print('number of changes in crop rotation', iacs_gp['crp_chngd_' + str(year)].sum(), 'total fields:', number_of_fields)
        print((number_of_changes/number_of_fields)*100, '% of fields changes')
        iacs_gp.to_file('./' + out_path + '/' + 'iacs_opt.shp')


def get_shares_seq(iacs_gp, n_years, temp_path, out_path):
    # this function calculates the area for each croptype for the entire study area and checks if the area of the
    # optimized allocation is equal to the initial allocation with a tolerance value

    for year in range(n_years):
        diss_init = iacs_gp.copy().dissolve(by=['crp_yr_' + str(year)], as_index=False)
        diss_init['area_init'] = diss_init.area * 0.0001

        diss_opt = iacs_gp.copy().dissolve(by=['OPT_KTYP_' + str(year)], as_index=False)
        diss_opt['area_opt'] = diss_opt.area * 0.0001

        merged = pd.merge(diss_init, diss_opt, left_on=[farm_id_column, 'crp_yr_' + str(year)], right_on=[farm_id_column, 'OPT_KTYP_' + str(year)])

        # diss_init = diss_init.drop('geometry')
        merged = merged[['crp_yr_' + str(year) + '_x', 'area_init', 'area_opt']]
        merged.to_csv('./' + out_path + '/shares_bb_iacs' + str(year) + '.csv')

        diss_init = iacs_gp.copy().dissolve(by=[farm_id_column, 'crp_yr_' + str(year)], as_index=False).copy()
        diss_init['area_init'] = diss_init.area * 0.0001
        diss_opt = iacs_gp.copy().dissolve(by=[farm_id_column, 'OPT_KTYP_' + str(year)], as_index=False).copy()
        diss_opt['area_opt'] = diss_opt.area * 0.0001
        merged = pd.merge(diss_init, diss_opt, left_on=[farm_id_column, 'crp_yr_' + str(year)], right_on=[farm_id_column, 'OPT_KTYP_' + str(year)])

        if diversity_type == 'attainable':
            # Check if farm crop acreage constraint was violated; There should be zero violations;
            error = 0
            for i, row in merged.iterrows():
                # print(row['area_opt'], row['area_init'])
                if (row['area_opt'] > row['area_init'] + row['area_init'] * (tolerance/100)) or (
                        row['area_opt'] < row['area_init'] - row['area_init'] * (tolerance/100)):
                    error += 1
            print('errors: ', error)


def get_shares(iacs_gp, temp_path, out_path):
    diss_init = iacs_gp.dissolve(by=[crop_type_column], as_index=False)
    diss_init['area_init'] = diss_init.area * 0.0001

    diss_opt = iacs_gp.dissolve(by=['OPT_KTYP'], as_index=False)
    diss_opt['area_opt'] = diss_opt.area * 0.0001

    merged = pd.merge(diss_init, diss_opt, left_on=[crop_type_column],
                      right_on=['OPT_KTYP'])

    merged = merged[[crop_type_column + '_x', 'area_init', 'area_opt']]
    merged.to_csv('./' + out_path + '/shares_bb_iacs.csv')

    diss_init = iacs_gp.dissolve(by=[farm_id_column, crop_type_column], as_index=False).copy()
    diss_init['area_init'] = diss_init.area * 0.0001

    diss_opt = iacs_gp.dissolve(by=[farm_id_column, 'OPT_KTYP'], as_index=False).copy()
    diss_opt['area_opt'] = diss_opt.area * 0.0001
    merged = pd.merge(diss_init, diss_opt, left_on=[farm_id_column, crop_type_column], right_on=[farm_id_column, 'OPT_KTYP'])

    if diversity_type == 'attainable':
        # Check if farm crop acreage constraint was violated; There should be zero violations;
        error = 0
        for i, row in merged.iterrows():

            # print(row['area_opt'], row['area_init'])
            if (row['area_opt'] > row['area_init'] + row['area_init'] * (tolerance / 100)) or (
                    row['area_opt'] < row['area_init'] - row['area_init'] * (tolerance / 100)):
                error += 1

        print('errors: ', error)

