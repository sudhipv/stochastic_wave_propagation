#!/bin/bash
echo "Enter 2 for 2D / 3 for 3D"
read input1

if [ $input1 == 2 ]
then
    echo "Entered in 2D"
    echo "Select one of the following"
    echo "1-Acoustic - Dterministic / 2 -Acoustic Stochastic"

    read input2

    if [ $input2 == 1 ]
    then
        echo "Selected: 2D Acoustic - Deterministic"
        cd ../src/2Dacoustic_detDDM/clusterMakefiles/
        cp makefile_mac ../makefile
        pwd
        cd ../
    elif [ $input2 == 2 ]
    then
        echo "Selected: 2D Acoustic - Stochastic"
        cd ../src/2Dacoustic_stoDDM/clusterMakefiles/
        # cd ../src/2Dpoisson_twolevel/clusterMakefiles/
        cp makefile_mac ../makefile
        pwd
        cd ../
    else
    echo "STOP: Wrong Input"
    fi

elif [ $input1 == 3 ]
then
    echo "Entered in 3D"
    echo "Select one of the following"
    echo "1- 3D elastic wave - Dterministic / 2 - 3D elastic wave Stochastic"

    read input3

    if [ $input3 == 1 ]
    then
        echo "Selected: 3D elastic wave - Deterministic"
        cd ../src/3Delasticwave_detDDM/clusterMakefiles/
        cp makefile_mac ../makefile
        pwd
        cd ../
    elif [ $input3 == 2 ]
    then
        echo "Selected: 3D elastic wave - Stochastic"
        cd ../src/3Delasticwave_stoDDM/clusterMakefiles/
        # cd ../src/2Dpoisson_twolevel/clusterMakefiles/
        cp makefile_mac ../makefile
        pwd
        cd ../
    else
    echo "STOP: Wrong Input"
    fi


else
    echo "STOP: Wrong Input"
fi

echo "Enter 0 for Nothing / 1 for VTK/Dat-outputs "
read input4

if [ $input4 == 1 ]
then
    echo "Selected VTK/Dat-outputs"
    pwd
    sed -i -e 's/outputFlag=0/outputFlag=1/g' main.F90
elif [ $input4 == 0 ]
then
    echo "No output selected"
    sed -i -e 's/outputFlag=1/outputFlag=0/g' main.F90
else
    echo "STOP: Wrong input"
    exit
fi


####*** Submit-Job method: same for all case
NP=$(<../../data/meshData/num_partition.dat)
echo "Using $NP processors"
rm ./a.out
make all
/Users/sudhipv/documents/PETSc/petsc-3.7.5/arch-darwin-c-debug/bin/mpiexec -np $NP ./a.out -log_view


# -info
# -mat_view ::ascii_info  - Prints info on matrix at conclusion of MatAssemblyEnd()
# -mat_view ::ascii_info_detail   - Prints more detailed info
# -mat_view   - Prints matrix in ASCII format
# -mat_view ::ascii_matlab    - Prints matrix in Matlab format
# -mat_view draw  - PetscDraws nonzero structure of matrix, using MatView() and PetscDrawOpenX().
# -display <name> - Sets display name (default is host)
# -draw_pause <sec>   - Sets number of seconds to pause after display

echo "====================================================================="
echo "FOR 3D Poisson or 3D elasticity"
echo "$ sh postprocessMac3D.sh"
echo "FOR 2D check *.vtk files in vtkOutputs"
echo "====================================================================="


cd ../../data/
