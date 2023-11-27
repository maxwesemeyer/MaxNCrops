from __functions import *


def create_reference_raster(gpd_frame, reference_gt):
    # creates a reference raster in 10m resolution that covers the study area and will be used for rasterization
    resolution = 10
    bounds = gpd_frame.total_bounds
    crs = gpd_frame.crs
    print(bounds)
    print(bounds[2] - reference_gt[0], bounds[2] , reference_gt[0])
    print(reference_gt[3] - bounds[1], reference_gt[3] , bounds[1])
    # this way the extent becomes divisible by 100
    range_y = ceil((bounds[2] - bounds[0]) / 1000) * 100
    range_x = ceil((bounds[3] - bounds[1]) / 1000) * 100
    arr = np.random.randint(5, size=(range_x, range_y)).astype('int16')

    transform = from_origin(bounds[0], bounds[3], resolution, resolution)
    new_dataset = rasterio.open('' + temp_path + '/' + 'reference_raster.tif', 'w', driver='GTiff',
                                height=arr.shape[0], width=arr.shape[1],
                                count=1,
                                dtype=rasterio.int16,
                                crs=crs,
                                transform=transform)
    new_dataset.write(arr, 1)
    new_dataset.close()


def rasterize_input_shp(crop_type_column=None, farm_id_column=None):
    #path = os.path.dirname(__file__)
    shp_p = glob.glob('./input/*.shp')[0]
    gt_ref = gdal.Open(glob.glob('./input/*.tif')[0]).GetGeoTransform()
    iacs = gpd.read_file(shp_p)
    iacs['field_id'] = range(1, len(iacs.index) + 1)
    iacs.to_file('./' + temp_path + '/' + 'iacs.shp')

    create_reference_raster(iacs, gt_ref)

    rst = rasterio.open('' + temp_path + '/' + 'reference_raster.tif')
    meta = rst.meta.copy()
    meta.update(compress='lzw')
    with rasterio.open('./' + temp_path + '/' + 'IDKTYP.tif', 'w+', **meta) as out:
        out.nodata = 0
        out_arr = out.read(1)
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom, value) for geom, value in zip(iacs.geometry, iacs[crop_type_column]))
        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)

    with rasterio.open('./' + out_path + '/init_crop_allocation.tif', 'w+', **meta) as out:
        out.nodata = 0
        out_arr = out.read(1)
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom, value) for geom, value in zip(iacs.geometry, iacs[crop_type_column]))
        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)
    with rasterio.open('./' + temp_path + '/' + 'Field_ID.tif', 'w+', **meta) as out:
        out.nodata = 0
        out_arr = out.read(1)
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom, value) for geom, value in zip(iacs.geometry, iacs.field_id))
        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)
    with rasterio.open('./' + temp_path + '/' + 'Farm_ID.tif', 'w+', **meta) as out:
        out.nodata = 0
        out_arr = out.read(1)
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom, value) for geom, value in zip(iacs.geometry, iacs[farm_id_column]))
        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)
