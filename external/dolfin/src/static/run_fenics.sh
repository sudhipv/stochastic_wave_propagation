#!/bin/bash
#SBATCH --account=def-asarkar
#SBATCH --nodes=2
##SBATCH --ntasks=32
#SBATCH --time=0-00:40
#SBATCH --mem-per-cpu=3700M
#SBATCH --tasks-per-node=32
#SBATCH --output=/scratch/sudhipv/sudhipv/ssfem_wave/data/slurm/%j.out
#SBATCH --mail-user=sudhipv@cmail.carleton.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## Import the required libraries to run fenics
source fenics_activate.sh
#module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6
#source ~/software/dolfin_17_1/share/dolfin/dolfin.conf
#source ~/fenics_17_1/bin/activate

#PYTHONPATH=/home/sudhipv/fenics_17_r/lib/python3.5/site-packages:$PYTHONPATH


#echo $PATH

#echo $PYTHONPATH	



echo "Trial Probelm Starting"

mpiexec -n 2 python gmshData_fenicsXML_parallel_3D.py

mpiexec -n 2 python elasticity3D_stochasticDDM_parallel_twolevel.py


echo "Trial Problem Ended"

#PATH=$HOME/fenics_17_1/lib/python3.5/site-packages:$PATH
## Compile and execute the preprocessor code
## This will convert meshDataFile (*.dat) to fenicsDataFiles (*.xml) in ../data/
mpiexec -n 64 python gmshData_fenicsXML_parallel_3D.py

## Compile and execute the FEniCS/dolfin assembly routines
## This will create assemlby matrices-vectors in ../data/Amats/
mpiexec -n 64 python elasticity3D_stochasticDDM_parallel_twolevel.py

exit
