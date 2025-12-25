#### Script file for installing petsc in remote machine

#### tested in niagara

module --force purge --all


module load CCEnv
module load StdEnv/2016.4
### loading StdEnv automatically loads other modules, but ensure these modules are loaded.

module load nixpkgs/16.09 intel/2016.4 openmpi/2.1.1 imkl/11.3.4.258

## This is also saved inside ./doc/petsc/

wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.7.5.tar.gz


tar -xvzf petsc-3.7.5.tar.gz

cd petsc-3.7.5

./configure --prefix=$HOME/software/petsctest/install/ --ignoreWarnings=1 --with-shared-libraries --with-scalar-type=real --with-blas-lapack-dir=${MKLROOT}/lib/intel64 --with-mpi-shared-libraries=1 --with-x=0 --with-x11=0 --with-mpi-dir=$EBROOTOPENMPI --with-debugging=no --with-64-bit-indices=0 --with-cuda=0 --with-c2html=0


# Follow guidelines after each step on screen by PETSc

