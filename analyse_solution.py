from __functions import *
import matplotlib.pyplot as plt

def analyse_solution_seq(tolerance=1):
    init = gdal.Open('./temp/IDKTYP.tif').ReadAsArray()
    # TODO init should be the historic crop types; Should be written to file earlier
    init = gdal.Open('./temp/hist_croptypes_adapted.tif').ReadAsArray()
    init[np.where((init == 255) | (init == 99), True, False)] = 0
    opt = gdal.Open('./output/maxent_croptypes_' + str(tolerance) + '.tif').ReadAsArray()
    mask = np.where(opt > 0, True, False)
    init[~mask] = 0

    print('UNIQUE INIT: ', np.unique(init), 'UNIQUE OPT: ', np.unique(opt))
    farm_id = gdal.Open('./temp/Farm_ID.tif').ReadAsArray()
    nd_value = -999
    n_years = opt.shape[0]
    img_entr_init_list = []
    img_entr_opt_list = []

    img_ct_init_list = []
    img_ct_opt_list = []
    agr_area_list = []
    print(opt.shape)
    for year in range(n_years):

        a, img_init = get_entropy(init[year, :, :], 100, return_count=True)
        a, img_opt = get_entropy(opt[year, :, :], 100, return_count=True)
        # TODO calculate init for each year; set 99 and 255 to the same class (0?) to use similar classes to opt
        a, img_entr_init = get_entropy(init[year, :, :], 100, return_img=True)
        a, img_entr_opt = get_entropy(opt[year, :, :], 100, return_img=True)

        img_ct_init_list.append(img_init)
        img_ct_opt_list.append(img_opt)
        img_entr_init_list.append(img_entr_init)
        img_entr_opt_list.append(img_entr_opt)

        a, agr_area = get_entropy(opt[year, :], 100, return_agr_area=True)
        agr_area_list.append(agr_area)

    img_init = np.stack(img_ct_init_list, axis=0)
    img_opt = np.stack(img_ct_opt_list, axis=0)
    img_entr_init = np.stack(img_entr_init_list, axis=0)
    img_entr_opt = np.stack(img_entr_opt_list, axis=0)
    agr_area_repeated = np.stack(agr_area_list, axis=0)

    write_array_disk_universal(img_init, './temp/reference_raster.tif',
                               outPath='./output/inital_ct',
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=100)

    write_array_disk_universal(img_opt,
                               './temp/reference_raster.tif',
                               outPath='./output/opt_ct_' + str(tolerance),
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=100)

    write_array_disk_universal(img_entr_init,
                               './temp/reference_raster.tif',
                               outPath='./output/initial_entropy',
                               dtype=gdal.GDT_Float32, noDataValue=nd_value, scaler=1, adapt_pixel_size=True,
                               adapted_pixel_size=100)

    write_array_disk_universal(img_entr_opt,
                                   './temp/reference_raster.tif',
                                   outPath='./output/opt_entropy_' + str(tolerance),
                                   dtype=gdal.GDT_Float32, noDataValue=nd_value, scaler=1, adapt_pixel_size=True,
                                   adapted_pixel_size=100)

    pd.DataFrame(
        {'entropy_init': img_entr_init.ravel(), 'entropy_opt': img_entr_opt.ravel(), 'initial_ct': img_init.ravel(),
         'opt_ct': img_opt.ravel(), 'agr_area': agr_area_repeated.ravel()}).to_csv(
        './output/entropy_ct_rav' + str(tolerance) + '.csv')

    ##########################################################
    # get shares of inital crop shares and save them to a csv file
    """
    for year in range(n_years):
        shares = get_shares_farm(farm_id.astype(int), init.astype(int)[year, :, :])
        out_dict = {'farm_id': shares[0], 'ID_KTYP': shares[1], 'area_m2': shares[2]}
        pd.DataFrame(out_dict).to_csv('./output/shares_init.csv')
    
        shares = get_shares_farm(farm_id.astype(int).flatten(), opt.astype(int)[year, :, :].flatten())
        out_dict = {'farm_id': shares[0], 'ID_KTYP': shares[1], 'area_m2': shares[2]}
        pd.DataFrame(out_dict).to_csv('./output/shares_opt_' + str(tolerance) + '.csv')
    """
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


def analyse_solution(tolerance=1):
    init = gdal.Open('./temp/IDKTYP.tif').ReadAsArray()
    opt = gdal.Open('./output/maxent_croptypes_' + str(tolerance) + '.tif').ReadAsArray()
    farm_id = gdal.Open('./temp/Farm_ID.tif').ReadAsArray()
    nd_value = -999
    a, img_init = get_entropy(init, 100, return_count=True)
    a, img_opt = get_entropy(opt, 100, return_count=True)

    a, img_entr_init = get_entropy(init, 100, return_img=True)
    a, img_entr_opt = get_entropy(opt, 100, return_img=True)

    a, agr_area = get_entropy(opt, 100, return_agr_area=True)

    write_array_disk_universal(np.expand_dims(img_init, axis=0), './temp/reference_raster.tif',
                               outPath='./output/inital_ct',
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=100)

    write_array_disk_universal(np.expand_dims(img_opt, axis=0),
                               './temp/reference_raster.tif',
                               outPath='./output/opt_ct_' + str(tolerance),
                               dtype=gdal.GDT_Int32, noDataValue=0, scaler=100, adapt_pixel_size=True,
                               adapted_pixel_size=100)

    write_array_disk_universal(np.expand_dims(img_entr_init, axis=0),
                               './temp/reference_raster.tif',
                               outPath='./output/initial_entropy',
                               dtype=gdal.GDT_Float32, noDataValue=nd_value, scaler=1, adapt_pixel_size=True,
                               adapted_pixel_size=100)

    write_array_disk_universal(np.expand_dims(img_entr_opt, axis=0),
                               './temp/reference_raster.tif',
                               outPath='./output/opt_entropy_' + str(tolerance),
                               dtype=gdal.GDT_Float32, noDataValue=nd_value, scaler=1, adapt_pixel_size=True,
                               adapted_pixel_size=100)
    pd.DataFrame(
        {'entropy_init': img_entr_init.ravel(), 'entropy_opt': img_entr_opt.ravel(), 'initial_ct': img_init.ravel(),
         'opt_ct': img_opt.ravel(), 'agr_area': agr_area.ravel()}).to_csv('./output/entropy_ct_rav' + str(tolerance) + '.csv')

    ##########################################################
    # get shares of inital crop shares and save them to a csv file
    shares = get_shares_farm(farm_id.astype(int), init.astype(int))
    out_dict = {'farm_id': shares[0], 'ID_KTYP': shares[1], 'area_m2': shares[2]}
    pd.DataFrame(out_dict).to_csv('./output/shares_init.csv')

    shares = get_shares_farm(farm_id.astype(int).flatten(), opt.astype(int).flatten())
    out_dict = {'farm_id': shares[0], 'ID_KTYP': shares[1], 'area_m2': shares[2]}
    pd.DataFrame(out_dict).to_csv('./output/shares_opt_' + str(tolerance) + '.csv')

    img_nan_init = img_entr_init.astype(float)
    img_nan_init[np.where(img_nan_init == nd_value, True, False)] = np.nan

    img_nan_opt = img_entr_opt.astype(float)
    img_nan_opt[np.where(img_nan_opt == nd_value, True, False)] = np.nan

    print('intial entropy avg: ', np.nanmean(img_nan_init.flatten()),
          'optimized avg: ', np.nanmean(img_nan_opt.flatten()),
          'difference abs: ', np.nanmean(img_nan_opt.flatten())-np.nanmean(img_nan_init.flatten()),
          'Diff in percent:',
          ((np.nanmean(img_nan_opt.flatten())-np.nanmean(img_nan_init.flatten()))/np.nanmean(img_nan_init.flatten()))*100)
