#!/bin/bash
echo "====================================================================="
echo "FOLLOW THESE STEPS: '$ cd external/dolfin/src/' and then run "
echo "Eg1: For 2D-onelevel: '$ python poisson2D_stochasticDDM_twolevel.py' "
echo "Eg2: For 3D-twolevel: '$ python poisson3D_stochasticDDM_twolevel.py' "
echo "Do the same for other python executables available in same folder "
echo "============================== PARALLE =============================="
echo "Eg1: mpiexec -n 4 python poisson2D_stochasticDDM_parallel_twolevel.py"
echo "Eg2: mpiexec -n 8 python poisson3D_stochasticDDM_parallel_twolevel.py"
echo "====================================================================="



docker start fenics_17_1 && docker exec -ti -u fenics fenics_17_1 /bin/bash -l

# docker exec -it wavebash bash

# docker exec -it wavebash cd ssfem_wave/external/dolfin/src/wave


# NP=$(<../../../data/meshData/num_partition.dat)
# echo "Using $NP processors"
# # echo "Enter 2 for 2D / 3 for 3D"
# # read input1

# # if [ $input1 == 2 ]
# # then




# mpiexec -n $NP python3 acoustic2D_detddm_parallel2L.py


# exit

# cd -


###*** Poisson2D/StochasticDDM/Twolevel
#python poisson3D_stochasticDDM_twolevel.py

## KLE/PCE data generator and importer
# !fenicsproject start dolfin
