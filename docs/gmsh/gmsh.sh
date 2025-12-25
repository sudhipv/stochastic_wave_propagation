

####### Step by step instruction to install gmsh


cd $HOME

#  screen - command preferred if you know how to use it

mkdir packages

cd packages

mkdir gmsh

cd gmsh

#  Binary install
cp $SCRATCH/sudhipv/ssfem_wave/docs/gmsh/gmsh-3.0.4-Linux64.tgz ./

tar -xvzf gmsh-3.0.4-Linux64.tgz

# Executable is inside gmsh-3.0.4-Linux/bin/gmsh



# For source Install
### Needs BLAS and LAPCK - you can load module load StdEnv/2020.

# cp $SCRATCH/sudhipv/ssfem_wave/docs/gmsh/gmsh-3.0.4-source.tgz ./

# tar -xvzf gmsh-3.0.4-source.tgz

# cd gmsh-3.0.4-source

# cmake -DCMAKE_INSTALL_PREFIX=$HOME/packages/gmsh/install/ -DENABLE_FLTK=0

# make -j 32

# make install



