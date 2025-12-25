cd ../../src

echo " MPI-Job: Select one of the following options: "
echo " 1 for 2D-Acoustic Deterministic"
echo " 2 for 2D-Acoustic Stochastic"
echo " 3 for 3D-Elasticity two-level"
echo " 4 for 3D-Elastic wave Deterministic"
echo " 5 for 3D-Elastic wave stochastic"
echo " 6 for 3D-Poisson stochastic"

read input2
if [ $input2 == 1 ]
then
    echo "MPI-Job: 2D-Acoustic Deterministic Two-level Selected"
    cd 2Dacoustic_detDDM/clusterMakefiles/
elif [ $input2 == 2 ]
then
    echo "MPI-Job: 2D-Acoustic Stochastic"
    cd 2Dacoustic_stoDDM/clusterMakefiles/
elif [ $input2 == 3 ]
then
    echo "MPI-Job: 3D-Elasticity Two-level Selected"
    cd 3Delasticity_twolevel/clusterMakefiles/
elif [ $input2 == 4 ]
then
    echo "MPI-Job: 3D-Elastic wave Deterministic Two-level Selected"
    cd 3Delasticwave_detDDM/clusterMakefiles/
elif [ $input2 == 5 ]
then
    echo "MPI-Job: 3D-Elastic wave stochastic Two-level Selected"
    cd 3Delasticwave_stoDDM/clusterMakefiles/
elif [ $input2 == 6 ]
then
    echo "MPI-Job: 3D-Poisson stochastic Two-level Selected"
    cd 3Dpoisson_twolevel/clusterMakefiles/
else
    echo "STOP: Wrong input"
    exit
fi


echo "MPI-Job: Enter 1 for Beluga"
read input1

if [ $input1 == 1 ]
then
    echo "Selected Beluga"
    cp makefile_beluga ../makefile
    cp run_ddm_beluga.sh ../run_ddm.sh
else
    echo "STOP: Wrong input"
    exit
fi


###*** Personalize each job based on requirements
cd ..

sed -i -e 's/time=0-00:10/time=11:30:30/g' run_ddm.sh
sed -i -e 's/mem-per-cpu=7700M/mem-per-cpu=3700M/g' run_ddm.sh
sed -i -e 's/tasks-per-node=32/tasks-per-node=40/g' run_ddm.sh
sed -i -e 's/nodes=1/nodes=20/g' run_ddm.sh
sed -i -e 's/-np 32/-np 800/g' run_ddm.sh
## Below d->dims, lc-> mesh density parameter, p->partitions
sed -i -e 's/job-name=DP_nRV_nNodes_nParts/job-name=3Dwavestochasticlayer_processmodel_90k/g' run_ddm.sh
sed -i -e 's#output=%x-%j.out#output=/scratch/sudhipv/ssfem_wave/data/slurm/%j-%x.out#g' run_ddm.sh
#sed -i -e 's/a.out/a.out -ksp_converged_reason/g' run_ddm.sh

echo "MPI-Job: Enter 0 for Nothing / 1 for VTK/Dat-outputs "
read input3

if [ $input3 == 1 ]
then
    echo "** MPI-Job: VTK/Dat-outputs selected **"
    sed -i -e 's/outputFlag=0/outputFlag=1/g' main.F90
elif [ $input3 == 0 ]
then
    echo "** MPI-Job: No VTK/Dat outputs selected **"
    sed -i -e 's/outputFlag=1/outputFlag=0/g' main.F90
else
    echo "STOP: Wrong input"
    exit
exit
fi

###*** Submit-Job method: same for all cases
sbatch run_ddm.sh
