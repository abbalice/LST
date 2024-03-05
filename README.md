# Laplacian Smoothing Tool (LST) - 2D
Repository for the paper: "Modeling tsunami initial conditions due to rapid co-seismic sea floor displacement: efficient numerical integration and a tool to build unit source databases"


The tool consists of three steps, detailed in the Supplementary Materials.

A working directory must be created, where the following folders and files are stored:

1. *EVENT_LOC* folder, containing the input data for the algorithm (coseismic deformation and sea-depth), organized in separated sub-folders. Specifically:
   (a) The folder *bathymetry* contains an *".xyz"* file for the sea-depth values in the region;
   (b) The folder *coseismic_deformation* contains sub-folders, each named *COMPONENT*, according to the specific component for the seafloor deformation. Within each of them, a single *".xyz"* file is stored, which is the result of the Okada (1985) algorithm.
2. The script *Library_2d.py*, where all the functions needed are written   
3. The script *step1.py*, that will create a subfolder *outputs_step1* in *EVENT_LOC*
4. The script *step2.py*, that will create the local database "database2_<EVENT_LOC>" in the working directory
5. The script *step3.py*, that will create a subfolder "initial_conditions" in *EVENT_LOC* where both a netCDF file *Xi0_<EVENT_LOC>_<COMPONENT>.nc* containing the final filtered initial condition according to a certain coseismic deformation and a *.png* plot showing the comparison between the unfiltered and the filtered initial condition will be stored.
6. The scripts *launch_mercalli.sh* and *test_launch.sh*. The first will launch automatically the latter on the nodes of the cluster Mercalli @INGV.
   
    
To run the first step, the command required is the following:

> python3 step1.py --event_loc EVENT_LOC

The second step, which is the construction of the local database, can be accomplished by these command lines:

> chmod +x launch_mercalli.sh
>  
> nohup ./launch_mercalli.sh &

If the tool is tested on another EVENT_LOC, line number # in *launch_mercalli.sh* must be changed accordingly.

N.B. To improve the efficiency of the database construction, it could be possible to comment lines number # and # in *launch_mercalli.sh* and to change manually the *iint* variable at line # and # in *test_launch.sh*, such that smaller chunks are run. For the application of the paper, we run chunks of 20 elements each. This will be fixed soon. 

To run the third step, the command is the following:

> nohup python3 step3.py --workdir $(pwd) --event_loc EVENT_LOC --component COMPONENT &

For the application shown in the manuscript, the *EVENT_LOC* folder is named *Kurils* and can be found in this repository.

# Laplacian Smoothing Tool (LST) - 1D

The equivalent 1D case is also provided. 
A working directory must be created, where the following folders and files are stored:

1. *EVENT_LOC* folder, containing the input data for the algorithm (coseismic deformation and sea-depth), organized in separated sub-folders. Specifically:
   - The folder *bathymetry* contains an *.xyz* file for the sea-depth values in the region;
   * The folder *coseismic_deformation* contains sub-folders, each named *COMPONENT*, according to the specific component for the seafloor deformation. Within each of them, a single *.xyz* file is stored, which is the result of the Okada (1985) algorithm.
2. The script *Library_1d.py*, where all the functions needed are written
3. The script *1d_Transect.py*, to perform the superposition of unit contributions.

To launch the script:

> python3 1d_Transect.py --event_loc EVENT_LOC --component COMPONENT






