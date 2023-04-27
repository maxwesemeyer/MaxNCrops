# MaxNCrops
MaxNCrops is a spatially explicit multiscale optimization model to determine attainable and potential crop diversity on the landscape scale (1x1km grid) using data from the Integrated Administration and Control System (IACS). It is based on linear programming (Gurobi) and models the relationships between all farms, their fields and the landscapes they are located in to find crop allocations with maximal diversity on the landscape scale. 
Farm level constraints ensure that the initial crop composition per farm is maintained. 
Field boundaries do not change and it consideres only crops that have been cultivated on a field in the previous years to ensure that plots are feasible for a given crop.

You can download the Project and easily run it with our example data. 
You can also run it with your own IACS data. You only need IACS inlcuding a farm ID and a crop type information per agricultural field as well as a raster with crop sequences. 

## You can use the requirements.txt to install the necessary packages like this: 

Build a new environment from spec list file

```
conda create --name myenv -f requirements.txt
```

Recreate from environment.yml file

```
conda env create -f environment.yml
```
