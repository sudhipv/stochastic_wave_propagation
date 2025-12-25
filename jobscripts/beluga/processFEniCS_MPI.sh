#!/bin/bash
cd ../../external/dolfin/src/jobscripts/wave

echo "FEniCS: Enter 1 for Beluga"
read input1

if [ $input1 == 1 ]
then
    echo "FEniCS: Selected Beluga"
    # cp fenics_activate_cedarcopy.sh ../../wave/fenics_activate.sh
    # cp fenics_activate_cedarcopy.sh ../../static/fenics_activate.sh

### new installation
    cp fenics_activate_beluga_2023.sh ../../wave/fenics_activate.sh
    cp fenics_activate_beluga_2023.sh ../../static/fenics_activate.sh
else
    echo "STOP: Wrong input"
    exit
fi

echo "FEniCS: Select one of the following options: "
echo " 1 for 2D-Acoustic Deterministic two-level parallel "
echo " 2 for 2D-Acoustic Stochastic two-level parallel "
echo " 3 for 3D-Elasticity two-level parallel "
echo " 4 for 3D-Elastic wave Deterministic two-level parallel "
echo " 5 for 3D-Elastic wave Stochastic two-level parallel "
echo " 6 for 3D-Poisson Stochastic two-level parallel "
read input2
if [ $input2 == 1 ]
then
    echo "** FEniCS: 2D-Poisson Two-level Parallel Selected **"
    cp run_fenics_acoustic2D_deterministicDDM_parallel_twolevel.sh ../../wave/run_fenics.sh
    cd ../../wave
elif [ $input2 == 2 ]
then
    echo "** FEniCS: 2D-Acoustic Stochastic two-level parallel **"
    cp run_fenics_acoustic2D_stochasticDDM_parallel_twolevel.sh ../../wave/run_fenics.sh
    cd ../../wave
elif [ $input2 == 3 ]
then
    echo "** FEniCS: 3D-Elasticity Two-level Parallel Selected **"
    cp ./../static/run_fenics_elasticity3D_stochasticDDM_parallel_twolevel.sh ../../static/run_fenics.sh
    cd ../../static
elif [ $input2 == 4 ]
then
    echo "** FEniCS: 3D-Elastic wave Deterministic Two-level Parallel Selected **"
    cp run_fenics_elasticwave_deterministicDDM_parallel_twolevel.sh ../../wave/run_fenics.sh

##### To run the layered version uncomment the line below.
    # cp run_fenics_elasticwave_deterministicDDM_parallel_twolevel_layered.sh ../../wave/run_fenics.sh

    cd ../../wave
elif [ $input2 == 5 ]
then
    echo "** FEniCS: 3D-Elastic wave Stochastic Two-level Parallel Selected **"
    # cp -f run_fenics_elasticwave_stochasticDDM_parallel_twolevel.sh ../../wave/run_fenics.sh

    cp -f run_fenics_elasticwave_stochasticDDM_parallel_twolevel_layered.sh ../../wave/run_fenics.sh

    cd ../../wave
elif [ $input2 == 6 ]
then
    echo "** FEniCS: 3D-Poisson Stochastic Two-level Parallel Selected **"
    cp -f ./../static/run_fenics_poisson3D_stochasticDDM_parallel_twolevel.sh ../../static/run_fenics.sh
    cd ../../static
else
    echo "STOP: Wrong input"
    exit
fi

### Use next line to increase run-time for FEniCS code
sed -i -e 's/time=0-00:15/time=00:30:00/g' run_fenics.sh
sed -i -e 's/mem-per-cpu=7700M/mem-per-cpu=3700M/g' run_fenics.sh
sed -i -e 's/nodes=1/nodes=20/g' run_fenics.sh
sed -i -e 's/tasks-per-node=32/ntasks-per-node=40/g' run_fenics.sh
sed -i -e 's/-n 32/-n 800/g' run_fenics.sh
sed -i -e 's#output=%j-3DPfenics.out#output=/scratch/sudhipv/ssfem_wave/data/slurm/%j_stochasticlayer_fenics.out#g' run_fenics.sh
### To avoid recreating FEniCS-xml files uncomment next line gmshData_fenicsXML_3D.py
#sed -i -e 's/python gmshData_fenicsXML_3D.py/#python gmshData_fenicsXML_3D.py/g' run_fenics3D.sh

###*** Submit-Job method: same for all cases
sbatch run_fenics.sh
