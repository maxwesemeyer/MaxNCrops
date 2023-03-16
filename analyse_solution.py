from __functions import *


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
