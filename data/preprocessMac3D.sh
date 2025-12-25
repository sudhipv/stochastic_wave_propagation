## KLE/PCE data generator and importer
## This take one command line input after number of dimensions (RVs)
# sh clean_all.sh

cd klePceData/
gfortran KLE_PCE_Data_commandLineArg.F90
./a.out
cd -

# sh clean_all.sh
## GMSH data generator
## Change to respective folder
cd meshData/

##### For elasticsity
# cp geofiles/slenderBeam.geo foo3D.geo


#### For layered Soil
cp geofiles/soil_layer_4.geo foo3D.geo

sed -i -e 's/lc1=1/lc1=1/g' foo3D.geo # Top part of soil
sed -i -e 's/lc2=1/lc2=2/g' foo3D.geo # middle
sed -i -e 's/lc3=1/lc3=25/g' foo3D.geo # Bottom

#sed -i -e 's/lc=0.1;/lc=0.08;/g' foo3D.geo

# For elastic wave
# cp geofiles/soil.geo foo3D.geo
# sed -i -e 's/lc1=1;/lc1=20;/g' foo3D.geo
# sed -i -e 's/lc2=1;/lc2=20;/g' foo3D.geo
## select the number of partitions

NP=8
## create mesh file::
gmsh -3 foo3D.geo -part $NP -o foo3D.msh

## create mesh data::
gfortran preprocmesh3D1_AD.F90 -O2 -o ./a.out;./a.out
gfortran preprocmesh3D2_AD.F90 -O2 -o ./a.out;./a.out

## remove unnecessory files created by "sed"
rm *.geo-e
cd -

## FEniCS data generator
# cd ../external/dolfin/src/wave/

# echo "Enter 1 for Parallel / 2 for Serial"
# read input1
# if [ $input1 == 1 ]
# then
#     export EVENT_NOKQUEUE=1
#     ##var=`cat ../../../meshData/num_partition.dat`
#     NP=$(<../../../../data/meshData/num_partition.dat)
#     echo "Using $NP processors"
#     /Users/sudhipv/documents/PETSc/petsc-3.7.5/arch-darwin-c-debug/bin/mpiexec -n $NP python3 gmshData_fenicsXML_parallel_3D.py
# elif [ $input1 == 2 ]
# then
#     python3 gmshData_fenicsXML_3D.py
# else
#     echo "STOP: Wrong Input"
# fi

# echo "NEXT: sh processFEniCS_Mac.sh"
# echo "=========================================================="
