# MaxNCrops
MaxNCrops is a spatially explicit multiscale optimization model to determine attainable and potential crop diversity on the landscape scale (1x1km grid) using data from the Integrated Administration and Control System (IACS). It is based on linear programming (Gurobi) and models the relationships between all farms, their fields and the landscapes they are located in to find crop allocations with maximal diversity on the landscape scale. 
Farm level constraints ensure that the initial crop composition per farm is maintained. 
Field boundaries do not change and it consideres only crops that have been cultivated on a field in the previous years to ensure that plots are feasible for a given crop.

You can download the Project and easily run it with our example data. 
You can also run it with your own IACS data. You only need IACS inlcuding a farm ID and a crop type information per agricultural field as well as a raster with crop sequences. 
## Concept 
Decision are made at the field scale, respecting farm level constraints while maximizing crop diversity on the landscape-scale. 
The figure shows how the field, landscape and farm scale are interrelated; 

![modellingappraoch_cropped](https://github.com/maxwesemeyer/MaxNCrops/assets/49986729/5bd4ff1e-87c0-4892-a9c9-23f9e55ab35a)

The spatial crop allocation changes: 

![github_figure](https://github.com/maxwesemeyer/MaxNCrops/assets/49986729/4a6149ae-0a85-4046-aac7-3aef49e41788)

left: the observed crop allocation in IACS; right: the optimized crop allocation with maximum landscape scale-crop diversity



## Use the requirements.txt to install the necessary packages like this: 

Build a new environment from spec list file

```
conda create --name myenv -f requirements.txt
```

Recreate from environment.yml file

```
conda env create -f environment.yml
```
## Run optimization
In order to run the optimization just run the main_optimize.py file. Adjustable parameters can be found in the config file. These are: 

+ ***agg_length***, controls the size of the landscape pixels; *agg_length* of 100 and 10m pixels equals a km² 
+ ***tolerance***, controls the crop composition tolerance per farm and crop type
+ ***crop_type_column***, in case you use your own data state the name of the crop type column here
+ ***farm_id_column***, in case you use your own data state the name of the farm identifier column here
+ ***diversity_type***, set to attainable to enforce the crop composition constraints or to potential to not enforce 
+ ***verbatim***, set this to True to print information about the data preprocessing;
+ ***nd_value***, the no data value for the project can be set here
+
## Output
A variety of output is generated and saved in the output folder. 
The most important are: 

+ ***init_crop_allocation.tif***, a raster with the initial crop allocation
+ ***opt_crop_allocation.tif***, a raster with the optimized crop allocation
+ ***initial_ShanDiv.tif***, a raster with the initial landscape-scale Shannon crop diversity 
+ ***opt_ShanDiv_.tif***, a raster with the optimized landscape-scale Shannon crop diversity 
+ ***iacs_opt.shp***, a shapefile with a cloumn OPT_KTYP that contains the optimized crop type allocation
