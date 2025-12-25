#!/bin/bash
#SBATCH --account=def-asarkar
#SBATCH --nodes=1
##SBATCH --ntasks=32
#SBATCH --time=0-00:15
##SBATCH --mem-per-cpu=7700M
#SBATCH --tasks-per-node=32
#SBATCH --output=%j-3DPfenics.out
#SBATCH --mail-user=sudhipv@cmail.carleton.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## Import the required libraries to run fenics
source fenics_activate.sh


dijitso clean
export DIJITSO_CACHE_DIR='$SCRATCH/ssfem_wave/external/dolfin/src/static'
# In order to remove the cache error for NIAGARA machine only


mpiexec -n 2 python gmshData_fenicsXML_parallel_3D.py

mpiexec -n 2 python elasticity3D_stochasticDDM_parallel_twolevel.py


echo "###############################################"
echo "Trial Problem completed"
echo "###############################################"

module list

echo "###############################################"
echo "Mesh Conversion Starting"
echo "###############################################"

## Compile and execute the preprocessor code
## This will convert meshDataFile (*.dat) to fenicsDataFiles (*.xml) in ../data/
mpiexec -n 32 python gmshData_fenicsXML_parallel_3D.py


echo "###############################################"
echo "Mesh Conversion Finished"
echo "###############################################"
##############################################



module list

################################################

echo "###############################################"
echo "Matrix Assembly Starting"
echo "###############################################"
## Compile and execute the FEniCS/dolfin assembly routines
## This will create assemlby matrices-vectors in ../data/Amats/
# mpiexec -n 32 python -m mpi4py elasticity3D_stochasticDDM_parallel_twolevel.py

mpiexec -n 32 python elasticity3D_stochasticDDM_parallel_twolevel.py


echo "###############################################"
echo "Matrix Assembly Finsihed"
echo "###############################################"

exit
