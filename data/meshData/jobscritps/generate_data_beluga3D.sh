#!/bin/bash
#SBATCH --time=0-00:30
#SBATCH --mem-per-cpu=3700M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --output=%j-3DP_Mesh.out
#SBATCH --mail-user=sudhipv@cmail.carleton.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

### first lc=current lc and second lc=desiredlc
### for cluster use lc=0.8 for 10K, lc=0.5 for 30K
#sed -i -e 's/lc=0.5;/lc=0.08;/g' fooServer3D.geo

## select the number of partitions
NP=32

module load StdEnv/2018.3
export PATH=${PATH}:/home/sudhipv/packages/gmsh/gmsh-3.0.4-source/install/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/sudhipv/packages/gmsh/gmsh-3.0.4-source/install/lib

#module --force purge
#module load arch/avx2 gcc/5.4.0
## create mesh file::
/home/sudhipv/packages/gmsh/gmsh-3.0.4-source/install/bin/./gmsh -3 foo.geo -part $NP -o foo3D.msh

# by loading gmsh installed in beluga due to error of local installation
#module load StdEnv/2020 intel/2020.1.217 gmsh/3.0.4
#gmsh -3 foo.geo -part $NP -o foo3D.msh

## create mesh data::
gfortran preprocmesh3D1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh3D2_AD_test.F90 -O2 -o ./a.out;./a.out

#module unload arch/avx2 gcc/5.4.0
#module load StdEnv/2018.3
## mesh visualization::
# gmsh gmsh.msh


