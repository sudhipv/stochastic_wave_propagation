!! Main Code: By Ajit Desai, Augest/2015, PETSc Version : March/2016
!! purpose: is to define inputs, load FEM/DDM data and call DDM solvers
!




!!! EDITED BY SUDHI P V - ACOUSTIC WAVE PROPAGATION - STOCHASTIC DDM - TWO LEVEL
!!! 2020/06/19



!input:
!  PCE:KLE-Data:
!      ndim  : random dimension of KLE/PCE expansion used for input represenations
!      nord  : order of basis polynomials used for PCE expansion
!      npcein: number of PCE terms : npceout
!      Apce  : stochastic-stiffness matrix
!      Fpce  : stochastic-source vector
!      Cijk  : multidimentional moments : tripple product
!      sigma : σ for underlying Guassian process : Strenght of uncertainty
!      omegas: w's of KLE expansion terms : Rogers Book : Scale of uncertainty
!      multipliers  : √λ / √(a1 = (sin(2*w_i a1)/2*w_i a_1)): KLE expansion terms : Rogers Book
!
!  MESH-Data:
!      p e t : points: elements: triangles: for local or each subdomain
!      np ne nt     : number of points, elements & triangles for local or each subdomain
!      pg eg tg     : points: elements: triangles: for global or whole domain
!      npg neg ntg  : number of points, elements & triangles for global or whole domain
!      nb nci nri   : number of nodes on boundary, corner & remaining for local or each subdomain
!      nbg nbgc nbgr: number of nodes on boundary, corner & remaining for global or whole domain
!      dbounds : physical domain limits : x=[0, 1], y=[0,1] : for unit square
!
!  MPI-Data:
!      pid   : Processor ID
!      nproc : number of processors
!      ierr  : MPI errors
!
!  FEM-Data:
!      Amat  : FEM assembly matrix
!      Aadv  : advection matrix
!      Adif  : diffusion matrix
!      fvec  : FEM assembly vector
!
!  DDM-Data:
!      Aii   : interior-interior sub-block of Amat
!      Aig   : interior-boundary sub-block of Amat
!      Air   : interior-remainig sub-block of Amat
!      Aic   : interior-corner sub-block of Amat
!      Fi Fg : interior & boundary sub-block of Fvec
!
! output:
!      Ui Ur Uc : local interior, remaining & corner nodes solution vector
!      Ub Ub_g  : local & global boundary nodes solution vector
!!
!!---------------------------------------------------------------------------------------
PROGRAM main

    use common
    use myCommon
    use assembly
    use PETScSolvers

    implicit none

include 'mpif.h'

!!------------------------------------------------------
!! Two-Level-NNC-PCGM with PETSc Sparse Iterative Solver
!!------------------------------------------------------

!!MPI-Data:
    integer :: ierr, pid

!!Exit-Criteria
    integer, parameter :: maxiter = 100
    double precision, parameter :: tol = 1.0d-5
    integer, parameter :: outputFlag=1

    integer            :: nord,ndim,npceout,npcein

!!Mallocs-Data:
    !! To precalculate "mallocsCals = 1" : First time run
    !! To use calculated "mallocsCals = 0" : Repeated runs
    integer, parameter :: mallocsCals = 1

!!KLE-Data-Inputs
    integer          :: i, ncijk  !,indexi,j, k,
    integer, allocatable, dimension(:,:)          :: ijk
    double precision, allocatable, dimension(:)   :: cijk

!!MESH-Data-Inputs:
    integer    :: npg, neg, ntg, nbg, nbgc, nbgr, nParts !nhg,
    integer    :: np, ne, nt, nb, nci, nri  !nh,
    character(len=255) :: label, extension, filename
    !! integer, allocatable, dimension(:,:)          :: hg
    double precision, allocatable, dimension(:,:) :: pg
    integer, allocatable, dimension(:,:)          :: eg
    integer, allocatable, dimension(:,:)          :: tg

!!FEM & DDM-Data-Outputs:
    double precision, allocatable, dimension(:)   :: Ui, Ub, Ub_g, U_gi, U_g


!! Newmark Beta/Time discretization Parameters

    double precision :: T, deltaT, beta_NB, gamma_NB

!!PCKF-Data-Inputs/Outputs                                                       !!** PCKF ***
    integer         :: nmeasg, One  !nmeas
    integer         :: tcount

    double precision, allocatable, dimension(:,:) :: Af, Afb
    double precision, allocatable, dimension(:)   :: var_g

!!Time Calculation Data:
    integer :: c1,c2,cr
    double precision :: time1
    double precision :: time2

!!-----------------------------------------------------------------------------------------------
!! MPI initiation
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,pid,ierr)

!!-----------------------------------------------------------------------------------------------
!!!!!!!!!!> 1: Pre-Processing: Memory-Allocation / Assembly / Distribution of Amat
!!-----------------------------------------------------------------------------------------------
!!*** Pre-Processing Time T1 ***!!!
    call system_clock(count_rate=cr)
    call cpu_time(time1)
    call system_clock(c1)

!!-----------------------------------------------------------------------------------------------
    One = 1
    if (pid .eq. 0) then
        call readGlobMeshDim(nParts)

        !! Read mesh data : local & global
        extension = ''
        nmeasg = 0
        call readmeshdim(npg,neg,ntg,nbg,extension)
        call readgcrdim(nbgc,nbgr)
        allocate(pg(2,npg),eg(3,neg),tg(3,ntg))
        call readmeshdata(pg,eg,tg,npg,neg,ntg,extension)
        call readNB_Para(T,deltaT,beta_NB,gamma_NB,extension)

        call readPceData(nord, ndim, npceout, npcein)
        call readncijk(ncijk)

        print*, '---------- MESH INFO -----------------'
        print*, 'nParts             =',nParts
        print*, 'nNodes             =',npg
        print*, 'nBoundary Nodes    =',nbg
        print*, 'nCorner Nodes      =',nbgc
        print*, 'nRemaing Nodes     =',nbgr

        print*, '---------- PCE INFO ------------------'

        print*, 'Order     =', nord
        print*, 'Dimen     =', ndim
        print*, 'nPCEin    =', npcein
        print*, 'nPCE      =', npceout
        print*, 'nCijk     =',ncijk


        print*, '---------- TIME INFO ------------------'

        print*, 'Total Time, T      =',T
        print*, 'deltaT             =',deltaT

    end if


    CALL MPI_BCAST(nbg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nbgr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nbgc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nmeasg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(nord,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ndim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npceout,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npcein,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ncijk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    CALL MPI_BCAST(T,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(deltaT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(beta_NB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(gamma_NB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


!!-----------------------------------------------------------------------------------------------

    call int2str(extension,pid+1,4)
    extension = trim(extension)
    !!call readmeshdim3D(np,ne,nt,nh,nb,extension)          !!**
    call readmeshdim(np,ne,nt,nb,extension)
    call readcrdim(nci,nri,extension)
    !call readnmeas(nmeas,extension)                                               !!** PCKF ***


    tcount = NINT(T/deltaT)    ! Total time steps
!!-----------------------------------------------------------------------------------------------
!! Dynamic array allocation
    allocate(Ub(nb*npceout),ijk(ncijk,3),cijk(ncijk)) !!p(2,np),e(3,ne),t(3,nt),
    allocate(Ui((np-nb)*npceout),Ub_g(nbg*npceout))

!!_______________________________________________________________________________________________

    if (pid .eq. 0) then
        !!!!!!  ijk are the indices in same order and cijk is the value - Only nonzero cijk values are stored
        call readcijk(ncijk,ijk,cijk)
    end if
    CALL MPI_BCAST(ijk,ncijk*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(cijk,ncijk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


!!-----------------------------------------------------------------------------------------------
!!*** Pre-Processing Time T2 ***!!!
    if (pid .eq. 0) then
        call cpu_time(time2)
        call system_clock(c2)
        if (pid .eq. 0) print*, '--------------------------------------------------------------'
        print*, 'Time taken for pre-processing', dble(c2-c1)/dble(cr),'seconds'
        print*, 'CPU-Time for pre-processing', time2-time1,'seconds'
        if (pid .eq. 0) print*, '--------------------------------------------------------------'
    end if

!!!-----------------------------------------------------------------------------------------------
!!*** Processing Time T1 ***!!! This includes assembly
    call cpu_time(time1)
    call system_clock(c1)

!!-----------------------------------------------------------------------------------------------
!!!!!!!!> 3: Main Process: Calling Solver
!!-----------------------------------------------------------------------------------------------

!    call PETSc_stonncpcgm(pid,p,e,t,np,ne,nt,nb,nbg,ndim,npcein,npceout,nri,nci,nbgc,nomga,ncijk,&
!                    ijk,cijk,mIndex,sIndex,casep,dbounds,const_diff,omegas,multipliers,sigma,amp,&
!                    Ui,Ub,Ub_g,mallocsCals,maxiter,tol)

    ! print*, 'process',pid,'has ni',np-nb,'and size of Ui', size(Ui)


    call PETSc_stonncpcgm(pid,pg,tg,np,nb,npg,ntg,nbg,npcein,npceout,nri,nci,nbgc,ncijk,ijk,cijk,Ui,Ub,Ub_g, &
                         mallocsCals,maxiter,tol,T, tcount, deltaT, beta_NB, gamma_NB,outputFlag)


!!-----------------------------------------------------------------------------------------------
!!!  The difference is whether the time is spent in user space or kernel space.
!    User CPU time is time spent on the processor running your program's code (or code in libraries);
!!   system CPU time is the time spent running code in the operating system kernel on behalf of your program.
!!
!!!-----------------------------------------------------------------------------------------------
!!*** Processing Time T2 ***!!!
    if (pid .eq. 0) then
        call cpu_time(time2)
        call system_clock(c2)
        if (pid .eq. 0) print*, '---------------------------------------------------------'
        print*, 'Time taken by solver', dble(c2-c1)/dble(cr),'seconds'
        print*, 'CPU-Time taken by olver',time2-time1,'seconds'
        if (pid .eq. 0) print*, '---------------------------------------------------------'
    end if



!!-----------------------------------------------------------------------------------------------
!! Dynamic array de-allocation
   if (pid .eq. 0) then
      deallocate(pg,eg,tg)
   end if

    deallocate(Ub,Ui,Ub_g)
    !!deallocate(omegas,multipliers,xi,mIndex,sIndex)
    ! deallocate(var_g) !,uveci,uvecb,,measnodes,d)                          !!** PCKF ***

   CALL MPI_FINALIZE(ierr)


END PROGRAM main
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










