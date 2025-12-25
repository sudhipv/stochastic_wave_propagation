#!/bin/bash
cd /home/a/asarkar/sudhipv/packages/fenics_intel
rm -rf .cache
cd -

cd ../../external/dolfin/src/jobscripts/wave

echo "FEniCS: Enter 1 for Niagara"
read input1

if [ $input1 == 1 ]
then
    echo "FEniCS: Selected Niagara"
    cp -f fenics_activate_niagara_intel2017.sh ../../wave/fenics_activate.sh
    cp -f fenics_activate_niagara_intel2017.sh ../../static/fenics_activate.sh
else
    echo "STOP: Wrong input"
    exit
fi

echo "FEniCS: Select one of the following options: "
echo " 1 for 2D-Acoustic Deterministic two-level parallel "
echo " 2 for 2D-Acoustic Stochastic two-level parallel  "
echo " 3 for 3D-Elasticity Stochastic two-level parallel "
echo " 4 for 3D-Elastic wave Deterministic two-level parallel "
echo " 5 for 3D-Elastic wave Stochastic two-level parallel "
echo " 6 for 3D-Poisson Stochastic two-level parallel "

read input2
if [ $input2 == 1 ]
then
    echo "** FEniCS: 2D-Poisson Two-level Parallel Selected **"
    cp -f run_fenics_acoustic2D_deterministicDDM_parallel_twolevel_NIAGARA.sh ../../wave/run_fenics.sh
    cd ../../wave
elif [ $input2 == 2 ]
then
    echo "** FEniCS: 2D-Acoustic Stochastic Two-level Parallel Selected **"
    cp -f run_fenics_acoustic2D_stochasticDDM_parallel_twolevel_NIAGARA.sh ../../wave/run_fenics.sh
    cd ../../wave
elif [ $input2 == 3 ]
then
    echo "** FEniCS: 3D-Elasticity Stochastic Two-level Parallel Selected **"
    cp -f ./../static/run_fenics_elasticity3D_stochasticDDM_parallel_twolevel_NIAGARA.sh ../../static/run_fenics.sh
    cd ../../static
elif [ $input2 == 4 ]
then
    echo "** FEniCS: 3D-Elastic wave Deterministic Two-level Parallel Selected **"
    # cp -f run_fenics_elasticwave_deterministicDDM_parallel_twolevel_NIAGARA.sh ../../wave/run_fenics.sh

    ##### To run the layered version uncomment the line below.
    cp run_fenics_elasticwave_deterministicDDM_parallel_twolevel_layered_NIAGARA.sh ../../wave/run_fenics.sh

    cd ../../wave
elif [ $input2 == 5 ]
then
    echo "** FEniCS: 3D-Elastic wave Stochastic Two-level Parallel Selected **"
    cp -f run_fenics_elasticwave_stochasticDDM_parallel_twolevel_NIAGARA.sh ../../wave/run_fenics.sh
    cd ../../wave
elif [ $input2 == 6 ]
then
    echo "** FEniCS: 3D-Poisson Stochastic Two-level Parallel Selected **"
    cp -f ./../static/run_fenics_poisson3D_stochasticDDM_parallel_twolevel_NIAGARA.sh ../../static/run_fenics.sh
    cd ../../static
else
    echo "STOP: Wrong input"
    exit
fi

### Use next line to increase run-time for FEniCS code
sed -i -e 's/time=0-00:15/time=00:15:00/g' run_fenics.sh
sed -i -e 's/nodes=1/nodes=4/g' run_fenics.sh
sed -i -e 's/--tasks-per-node=32/--ntasks-per-node=40/g' run_fenics.sh
sed -i -e 's/-n 32/-n 160/g' run_fenics.sh
sed -i -e 's#output=%j-3DPfenics.out#output=/scratch/a/asarkar/sudhipv/ssfem_wave/data/slurm/fenics.out#g' run_fenics.sh

### To avoid recreating FEniCS-xml files uncomment next line gmshData_fenicsXML_3D.py
#sed -i -e 's/python gmshData_fenicsXML_3D.py/#python gmshData_fenicsXML_3D.py/g' run_fenics3D.sh

###*** Submit-Job method: same for all cases


################################TRIAL PROBLEM TO REMOVE ERROR FOR NIAGARA ############################
##############################################

# In order to remove the cache error for NIAGARA machine only

if [ $input2 == 4 ] || [ $input2 == 5 ] || [ $input2 == 6 ]
then
    source fenics_activate.sh

    dijitso clean
    export DIJITSO_CACHE_DIR='$SCRATCH/ssfem_wave/external/dolfin/src/wave'

    cd /home/a/asarkar/sudhipv/packages/fenics_new
    rm -rf .cache
    cd -

    # instant-clean
    echo "Trial Problem Started"
    mpiexec -n 2 python gmshData_fenicsXML_parallel_3D.py

    if [ $input2 == 4 ]
    then
        mpiexec -n 2 python elasticwave3D_detddm_parallel2L_layer.py
    elif [ $input2 == 5 ]
    then
        mpiexec -n 2 python elasticwave3D_stoddm_parallel2L.py
        #mpiexec -n 2 python elasticbeam3D_stoddm_parallel2L.py
    elif [ $input2 == 6 ]
    then
        mpiexec -n 2 python poisson3D_stochasticDDM_parallel_twolevel.py
    fi



    echo "Trial Problem Ended"


fi

################################################


sbatch run_fenics.sh
