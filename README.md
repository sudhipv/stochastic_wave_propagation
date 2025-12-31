# Scalable Solvers for Spectral Stochastic FEM : Acoustic Wave Propagation and 3D Elastic Waves 

This solver have capabilities to solve SSFEM using Two level Domain decomposition based preconditioners for 2D acoustic wave and 3D elastic wave propagation. 

<img width="1000" height="1000" alt="mean_30" src="https://github.com/user-attachments/assets/c411d668-98cd-4b96-8f6c-129a4349f54b" />
<img width="1000" height="1000" alt="mean_60" src="https://github.com/user-attachments/assets/65d07d31-89d0-49e7-8e4a-d0033d6fdbbc" />


### Components
1. FEniCS: Employed for deterministic FE matrix-vector assembly
2. PETSc: Employed for (local) sparse storage of stochastic FE matrix,vector and associated parallel linear algebraic calculations
3. MPI: Employed for parallel processing of PETSc local routines

### Folder structure :

* jobscripts   :: USER EDITABLE FILES FOR BATCH SUBMISSION : NODES/TASKS/TIME
     :: niagara/cedar/graham/beluga - script files for respective machines in these folders  

    - preprocess.sh           ::  Generates KLE/PCE data and GMSH/DDM data
    - processFEniCS.sh    :: generates FEniCS/DDM assembly serially
    - processFEniCS_MPI.sh :: generates FEniCS/DDM assembly simultaneously in each core for each subdomain by parallel code
    - processMPI.sh          :: MPI/PETSc compile and execute 
    - loadmodules.sh         :: helper file to "module load" the exact PETSc/FEniCS/MPI stacks required by each cluster
    - salloc.sh              :: example SLURM salloc command for interactive debugging on ComputeCanada clusters
       
* src    :: executable : No need to edit anything here unless you have to change solver properties
    - 3Delasticity_twolevel 
    - 3Dpoisson_twolevel
    - 2Dpoisson_twolevel
    - 2Dacoustic_detDDM
    - 2Dacoustic_stoDDM
    - 3Delasticwave_detDDM
    - 3Delasticwave_stoDDM
    - README.md :: short summary of the PETSc/FEniCS driver directories
    
* external :: No need to edit unless you need to change the input characteristics like gaussian mean or sd  / fenics compilation / assembly 
    - puffin     :: external matlab FEM package
    - dolfin     :: external python/C++ FEM package
    - external/dolfin/src/static     :: Static Poisson/elasticity
    - external/dolfin/src/wave       :: wave problems
    - external/dolfin/src/jobscripts/static :: job scripts for static Poisson/elasticity
    - external/dolfin/src/jobscripts/wave  :: job scripts for wave problems
    
* data   :: Folder for mesh/kle/output data storage
    - klePceData :: KLE and PCE data
    - mallocData :: malloc allocation data
    - meshData   :: preprocessed mesh data
    - vtkOutputs :: outputs of simulation
    - slurm      :: All log files from runs
    - solution   :: interface/interior solution dumps + iteration counts from PETSc runs
    - preprocessMac2D.sh / preprocessMac3D.sh :: local workflow scripts that chain the KLE/PCE generator, gmsh preprocessing and mesh partitioning
    - processFEniCS_Mac.sh / processMPI_Mac.sh / postprocessMac3D.sh :: laptop automation helpers for launching FEniCS (Docker) and PETSc solvers
     
* docs :: PETSc/Fenics/gmsh subfolders
    - contains project related documents
    - Installation files for specific components
    - Error file
    - docker/fenics_docker.txt :: step-by-step guide for recreating the tested Fenics 2017.1.0 Docker container
    - laptop_instructions.txt :: quick reference for running the entire workflow on a personal workstation
    - Submitting a SLURM Job Script.pdf :: reminder of scheduler policies across ComputeCanada clusters


### Typical Workflow :: Preprocess -> Assembly -> Solve -> Postprocess

1. **Preprocess random field + mesh data**
    * Remote clusters : `$ sh jobscripts/<cluster>/preprocess.sh` automatically builds KLE/PCE data via `data/klePceData/KLE_PCE_Data_commandLineArg.F90`, updates `data/meshData/num_partition.dat`, and partitions the gmsh meshes with the requested `lc` and partition counts.
    * Local testing : `$ sh data/preprocessMac2D.sh` or `$ sh data/preprocessMac3D.sh` performs the same sequence (copy geo template, tweak `lc`, run gmsh, call `preprocmesh*.F90`). Edit `data/meshData/geofiles/*.geo` and `data/meshData/jobscritps` to point to the gmsh executable you compiled.
    * Mesh artefacts live in `data/meshData/*.msh`, the stochastic inputs live under `data/klePceData/`, and the number of partitions requested later by PETSc/FEniCS is stored in `data/meshData/num_partition.dat`.
2. **FEniCS / DDM assembly**
    * Remote clusters : `$ sh jobscripts/<cluster>/processFEniCS.sh` orchestrates serial assembly; `$ sh jobscripts/<cluster>/processFEniCS_MPI.sh` runs the same Python drivers with MPI so that each subdomain is assembled in parallel.
    * Local testing : `$ sh data/processFEniCS_Mac.sh` attaches to the validated Docker image (`fenics_17_1`) described in `docs/docker/fenics_docker.txt` and prints example commands such as `python poisson3D_stochasticDDM_twolevel.py` or `mpiexec -n $NP python poisson2D_stochasticDDM_parallel_twolevel.py`. Those scripts live in `external/dolfin/src/static` and `external/dolfin/src/wave`.
3. **MPI/PETSc solve**
    * Remote clusters : `$ sh jobscripts/<cluster>/processMPI.sh` loads the desired solver (`src/2Dacoustic_detDDM`, `src/3Delasticwave_stoDDM`, etc.), copies the matching makefile from `clusterMakefiles/`, and runs `mpiexec -np $(cat data/meshData/num_partition.dat)` on the generated `a.out`.
    * Local testing : `$ sh data/processMPI_Mac.sh` mirrors the cluster script, lets you toggle deterministic vs stochastic solvers, switches `outputFlag` between `.dat/.vtk` outputs, and finally executes PETSc using the location configured inside the script.
    * Solver output is deposited under `data/vtkOutputs/` (ParaView), `data/solution/` (Matlab post-processing), and `data/slurm/` (scheduler logs).
4. **Postprocess**
    * `data/postprocessMac3D.sh` copies the vtk/solution files to a workstation and launches Matlab/ParaView commands listed at the bottom of that script.
    * `data/preprocessMac2D.sh` already prints reminders for the next command (`sh processFEniCS_Mac.sh`), and you can review `docs/laptop_instructions.txt` for a concise checklist. Use `data/clean_all.sh` and `data/mesh_clean.sh` whenever you need to reset intermediary artefacts.


### Installation of Required packages
* GMSH       :: 3.0.4 (Use only this version) : Binary installation preferred
    * Binary : http://gmsh.info/bin/   : download the 3.0.4 archive that matches your OS and place it under `docs/gmsh/` (archives are no longer kept inside the repo to keep git pushes under 100 MB)
    * Source : https://gitlab.onelab.info/gmsh/gmsh/blob/master/README.txt     : download the matching source `.tgz` manually and adapt `docs/gmsh/gmsh.sh` if you prefer to build from sources
    * Install steps : `./docs/gmsh/gmsh.sh` (expects the downloaded archive to be next to the script or gmsh already available on `PATH`)
* PETSc      :: 3.7.5 (Use only this version) : Source Install Recommended
    * Local Install  : https://www.mcs.anl.gov/petsc/documentation/installation.html
    * Remote Install : Refer to 'PETSC_Install.pdf' under docs folder
    * Install Steps : petsc_install.sh
* FEniCS     :: 2017.1.0 (Use only this version) : Source Install Recommended
    * Local installation : https://fenics.readthedocs.io/en/latest/installation.html#
    * Remote Install: Follow instructions from "fenics_install.pdf" available in docs/fenics folder. One could also check corresponding script files for installation and modify accordingly as below.
    * Install Steps : Niagara : fenics_install_niagara_intel.sh, CEDAR/GRAHAM/BELUGA : fenics_install.sh
    *  For FEniCS installation issues and updates refer to [computeCanada wiki page for FEniCS](https://docs.computecanada.ca/wiki/FEniCS)
    *  Docker recipe : `docs/docker/fenics_docker.txt` rebuilds the tested `fenics_17_1` container for laptop use
    

### Other Packages which can be loaded from remote machines 
* MPI
* PYTHON : 3
## For visualization
* PARAVIEW   :: 5.4.0 : 
* MATLAB     :: 2017


## Running Code LOCALLY (Laptop / Docker testing)

All helper scripts mentioned below live under `data/` and mirror the cluster workflow so you can validate meshes, assemblies and PETSc solves on a workstation before submitting to SLURM. Detailed checkpoints are listed in `docs/laptop_instructions.txt`.

* $ sh data/preprocessMac2D.sh   :: builds KLE/PCE data, regenerates 2D meshes from `data/meshData/geofiles/*.geo`, lets you adjust the `lc` parameter, partitions with gmsh, and compiles the `preprocmesh*.F90` stages. Use `preprocessMac3D.sh` for 3D.
* $ sh data/processFEniCS_Mac.sh :: displays example serial/parallel Python invocations, attaches to the Fenics Docker container (`docker start fenics_17_1 && docker exec ...`), and leaves you inside `external/dolfin/src/` so you can run the desired deterministic/stochastic assembly scripts manually.
* $ sh data/processMPI_Mac.sh    :: selects the deterministic vs stochastic solver, copies the right `clusterMakefiles/makefile_mac`, flips `outputFlag` if you want `.vtk/.dat` dumps, recompiles the Fortran driver, and runs PETSc with `mpiexec -np $(cat data/meshData/num_partition.dat) ./a.out -log_view`.
* $ sh data/postprocessMac3D.sh  :: post-process 3D cases on laptops (copy vtk outputs, trigger Matlab scripts, thin out data folders). For 2D runs simply open the `vtkOutputs/out_*.vtk` files with ParaView.
* Remember to update the gmsh path within `data/meshData/jobscritps/*` and the PETSc path hard-coded in `data/processMPI_Mac.sh` so they match your local installation. Use `data/fenics_clean.sh`, `data/mesh_clean.sh`, and `data/clean_all.sh` to reset intermediate files between tests.


## Running Code in REMOTE CLUSTERS (CEDAR, GRAHAM, NIAGARA, BELUGA)

No special modules need to be loaded by user before running the code since this is automatically done in seperate script files as needed.
Code uses GFORTRAN for compiling fortran codes. Generally most HPC clusters contain gfortran installled on root. If not load the modules appropriately.

### Step-1 :: Get the repository from bitbucket (Note: private repo, needs permission to access)
* $ git clone https://sudhipv@bitbucket.org/sudhipv/ssfem_wave.git
* $ git fetch && git checkout "master"
* $ cd ssfem_wave/
* NOTE: Replace the URL with your GitHub (or other) remote if you mirrored the project elsewhere; every script assumes the root directory is still called `ssfem_wave`.


### Step-2 :: Change executable paths

*   Change the path to your GMSH installation inside : ssfem_wave/data/meshData/jobscritps (all files)
    Replace here : /home/sudhipv/project/sudhipv/packages/gmsh/gmsh-3.0.4-Linux/bin/./gmsh /path/to/your/executable
*   Change the path to your fenics installtion inside : "ssfem_wave/external/dolfin/src/jobscripts/static" and "ssfem_wave/external/dolfin/src/jobscripts/wave"
    Replace these lines by your path inside all files starting with "fenics_...": source ~/software/dolfin_intel/share/dolfin/dolfin.conf,  source ~/packages/fenics_intel/bin/activate.
    Replace mail info inside all files starting with "run_fenics_...."
*   Change your mail id and required log info inside all "run_ddm_..." files : "src/3Delasticity_twolevel/clusterMakefiles/" - This can be done for all folders inside ./src/
*   If you are not running the code first time, you might want to delete all the previous data file inside data using the command "sh ./data/clean_all.sh"
    

### Step-3 :: preprocess & process:: GMSH/FEniCS/MPI/PETSc
* $ cd jobscripts/name_of_cluster (niagara/cedar/graham/beluga)        :: from top level project directory
* $ source preprocess.sh      :: This script file takes input from the user for number of random variable to use for stochastic expansion as well as the physical dimension for the problem.     
                             One could edit out the parameters inside this script file as in steps 2, 3,4.
      
    1. It promts for one command line input for nDim (number of RVs): enter the desired integer (Eg: 2 for 2RVs)
     
     generates KLE/PCE data for selected nDim (RVs) case. Data generated inside ssfem_fenics/data/klePceData/
       
    2. Control mesh density using lc parameter: Eg: lc=0.1/lc=0.09  , here lc=0.1 changes to lc=0.09 (NOTE: keep first "lc" always fixed to "0.1")
    3. Control type of coarse grid using "wb" parameter: Eg: wb=1/wb=0 , changes coarse grid from wire-basket(wb=1) to vertex grid(wb=0)
    4. Control number of partitions using "NP" parameter: Eg: NP=32/NP=64 , changes numper of partitions from 32 to 64 (default=32)

    generates GMSH/DDM data for selected LC and NP
        
* $ source processFEniCS_MPI.sh    ::  Generates FEniCS mesh data and DDM assembly Mat-Vec data
    1. Control request time for job using "time" parameter: Eg: time=0-00:15/time=0-00:30 , changes job time from 15mins to 30mins (default=15mins)
    2. Control number of nodes and tasks-per-node with respective parameter values.
    
    
    Note : By default mail comes only when the process ends or fails. Changes to mail id and default value have to be changed inside other script files inside : ssfem_fenics/external/dolfin/src/jobscripts.
    Depending upon the problem you have selected this script file changes. The name of this file can be found from commands inside processFEniCS_MPI.sh. ex : run_fenics_poisson2D_stochasticDDM_parallel_twolevel.sh, 
    run_fenics_elasticity3D_stochasticDDM_parallel_twolevel.sh

    In many cases the process fails because of an issue with the way MPI Job is initialized. Please ensure your job is successfull by looking at the output folder you selected while running the job
    folder -default : /scratch/sudhipv/sudhipv/ssfem_fenics/data/slurm/.
    In general a successful job prints out details of random variable, applying BCs, and finally success message.

    A sample Output looks like :
 
    * Running FEniCS for 3D stochastic Poisson/twoLevel
    * number of partitions: 320
    * number of dimensions: 3
    * number of pceInputss: 10
    * number of pceOutputs: 20
    * ==========================Success============================



* $ sh processMPI.sh          :: compiles FORTRAN executables and run MPI/PETSc simulation
    1. Control request time for job using "time" parameter: Eg: time=0-00:10/time=0-00:15 , changes job time from 10mins to 15mins (default=10mins)
    2. Control request tasks per node using "tasks-per-node" parameter: Eg: tasks-per-node=32/tasks-per-node=16 , changes from 32 to 16 (default=32)
    3. Control request number of nodes using "nodes" parameter: Eg: nodes=1/nodes=2 , changes requested node from 1 to 2 (default=1)
    4. Contorl number of MPI processor to compile and run using "np" parameter: Eg: -np 32/-np 16 , changes from 32 to 16 (Note: np=nodes*tasks-per-node)
    5. Option to control job-name and output-file name eg: job-name=test/job-name=d2_lc05_p32 ,  (here d->dims, lc-> mesh density parameter, p->partitions)
    
    
### Step-3 :: postprocess3D :: Copy data to HomePC and use Matlab/ParaView for postprocessing
* from top level project directory @HomePC
* $ cd data/vtkOutputs
* $ copy "posprocess3D.sh" to your local pc and edit the file locations    
* $ Then run 'sh posprocess3D.sh'  :: should do all of the steps outlined below (or follow step by step guide)


#### Use ParaView for visualization
* $ paraview out_*.vtk      :: open *.vtk files with paraview to check outputs



## NOTE:: Additional Useful Commands
* $ sacct  :: to check status for all of your recent jobs submitted using slrum
* $ more slrum-JobID.out  :: for outputs info of perticular JobID
* $ git fetch && git checkout graham
* $ git checkout -f   :: this will discard any local changes (useful to pull latest changes and discard old one)



### Submit the job ###
* $ sbatch run_ddm.sh     !! submit the job
* $ squeue -u userID     !! check status
* $ scancel jobID          !! kill the job

### Once status is complete then check ###
* $ more outputfile*
* $ more errorfile*

### Output visualization : ###
* $ paraview out_deterministic.vtk               

### Clean compilation/output data ###
* $ sh clean.sh             :: cleans everything

### Who do I talk to? ###
* Sudhi Sharma P V : sudhi.pv@cmail.carleton.ca


### Reference ##


> **[Scalable Domain Decomposition Methods for Nonlinear and Time-Dependent Stochastic Systems](https://doi.org/10.22215/etd/2023-15817)**

**Authors:** Vasudevan, Padillath and Sharma, Sudhi  
**Institution:** Carleton University (2023)  
**DOI:** [10.22215/etd/2023-15817](https://doi.org/10.22215/etd/2023-15817)

<details>
<summary><b>Click to expand BibTeX citation</b></summary>

```bibtex
@phdthesis{vasudevan2023scalable,
  title={Scalable Domain Decomposition Methods for Nonlinear and Time-Dependent Stochastic Systems},
  author={Vasudevan, Padillath and Sharma, Sudhi},
  year={2023},
  school={Carleton University},
  doi={10.22215/etd/2023-15817}
}
\```


