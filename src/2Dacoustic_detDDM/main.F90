

!!! EDITED BY SUDHI P V - ACOUSTIC WAVE PROPAGATION - DETERMINISTIC DDM - TWO LEVEL
!!! 2020/05/24



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
    integer, parameter :: maxiter = 1000
    double precision, parameter :: tol = 1.0d-5
    integer, parameter :: outputFlag=1

!!PCE-Data:
    !! Select casep
    !! 2: for automatic
    !! 3: for 3dim, 3ord
    !! 4: for 4dim, 3ord
    !! 5: for 5dim, 3ord
    !integer, parameter :: casep = 2
    !integer, parameter :: nomga = 9
    ! integer            :: nord,ndim,npceout,npcein

!!Mallocs-Data:
    !! To precalculate "mallocsCals = 1" : First time run
    !! To use calculated "mallocsCals = 0" : Repeated runs
    ! integer, parameter :: mallocsCals = 1

!!KLE-Data-Inputs
    ! integer          :: i, ncijk  !,indexi,j, k,
    ! integer, allocatable, dimension(:,:)          :: ijk
    ! double precision, allocatable, dimension(:)   :: cijk
    !double precision :: amp,const_diff,sigma,gamma
    !double precision, dimension(2,2)              :: dbounds
    !double precision, dimension(2)                   :: centre
    !double precision, allocatable, dimension(:)   :: omegas, multipliers, xi
    !integer, allocatable, dimension(:,:)          :: mIndex, sIndex

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
    integer         :: tcount, i

    double precision, allocatable, dimension(:,:) :: Af, Afb

    ! double precision, allocatable, dimension(:,:) :: U_wave

    ! double precision, allocatable, dimension(:)   :: var_g !uveci, uvecb
    !integer, allocatable, dimension(:) :: measnodes  !!,measnodesg
    !double precision, allocatable, dimension(:)   :: d

!!Time Calculation Data:
    integer :: c1,c2,cr
    double precision :: time1
    double precision :: time2
    !integer :: seed_size
    !integer, dimension (:), allocatable :: seed

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
        !!call readmeshdim3D(npg,neg,ntg,nhg,nbg,extension)
        call readgcrdim(nbgc,nbgr)
        allocate(pg(2,npg),eg(3,neg),tg(3,ntg)) !!,measnodesg(nmeasg))
        call readmeshdata(pg,eg,tg,npg,neg,ntg,extension)
        call readNB_Para(T,deltaT,beta_NB,gamma_NB,extension)

        print*, '--------------------------------------'
        print*, 'nParts             =',nParts
        print*, 'nNodes             =',npg
        print*, 'nBoundary Nodes    =',nbg
        print*, 'nCorner Nodes      =',nbgc
        print*, 'nRemaing Nodes     =',nbgr


        print*,'***************************************************'
        print*,'******** Time discretization Parameters ***********'
        print*,'***************************************************'

        print*, 'Total Time, T      =',T
        print*, 'deltaT             =',deltaT

    end if


    CALL MPI_BCAST(nbg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(npg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nbgr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nbgc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(nmeasg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)               !!** PCKF ***

    CALL MPI_BCAST(T,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(deltaT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(beta_NB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(gamma_NB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


    call int2str(extension,pid+1,4)
    extension = trim(extension)
    !!call readmeshdim3D(np,ne,nt,nh,nb,extension)          !!**
    call readmeshdim(np,ne,nt,nb,extension)
    call readcrdim(nci,nri,extension)
    !call readnmeas(nmeas,extension)                                               !!** PCKF ***


    tcount = NINT(T/deltaT)    ! Total time steps
!!-----------------------------------------------------------------------------------------------
!! Dynamic array allocation
    allocate(Ub(nb)) !!p(2,np),e(3,ne),t(3,nt),
    allocate(Ui(np-nb),Ub_g(nbg),U_g(npg),U_gi(npg))
    allocate(Af(np,1),Afb(nbg,1))
    !!allocate(uveci(np-nb),uvecb(nbg),d(nmeas),measnodes(nmeas))                   !!** PCKF ***

    !! allocate(Ur(nri*npceout),Uc(nci*npceout),UiUr((np-nci)*npceout),&
    !! FETIuveci(np-nci),FETIuvecb(nbgc),FETIAfb(nbgc,npceout+nmeasg),Uc_g(nbgc*npceout))

!!-----------------------------------------------------------------------------------------------
    ! if (pid .eq. 0) then
    !     call readcijk(ncijk,ijk,cijk)
    ! end if
    ! CALL MPI_BCAST(ijk,ncijk*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    ! CALL MPI_BCAST(cijk,ncijk,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    !call readmeshdata(p,e,t,np,ne,nt,extension)
    !call readmeasnodes(nmeas,measnodes,extension)                                 !!** PCKF ***

!!!-----------------------------------------------------------------------------------------------
!!!!!!!!!> 2: Calling Solver Only for Mallocs Calculations
!!!-----------------------------------------------------------------------------------------------
!! If you want to precalculate mallocs "mallocsCals = 1'      : First time run : Necessory
!! If have calculated mallocs then just use "mallocsCals = 0" : Repeated runs
!! mallocsCals = 1

!    if (mallocsCals .eq. 1) then
!
!    call PETSc_stonncpcgm(pid,p,e,t,np,ne,nt,nb,nbg,ndim,npcein,npceout,nri,nci,nbgc,nomga,ncijk,&
!                    ijk,cijk,mIndex,sIndex,casep,dbounds,const_diff,omegas,multipliers,sigma,amp,&
!                    Ui,Ub,Ub_g,mallocsCals,maxiter,tol)
!    end if

!!!-----------------------------------------------------------------------------------------------
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


    call PETSc_detnncpcgm(pid,np,pg,eg,tg,npg,ntg,nb,nbg,nri,nci,nbgc,Ui,Ub,Ub_g,maxiter,tol,T, tcount, deltaT, beta_NB, gamma_NB,outputFlag)


    ! call PETSc_stonncpcgm(pid,np,nb,nbg,npcein,npceout,nri,nci,nbgc,ncijk,ijk,cijk,Ui,Ub,Ub_g, &
    !                     mallocsCals,maxiter,tol)


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

!-----------------------------------------------------------------------------------------------
!! PCGM-Post-Processing
!!-----------------------------------------------------------------------------------------------

!!!!!!!!! For taking just one time step solution !!!!!!!!!!!!!
! if (outputFlag .eq. 1) then

!         call create_Af_pckfddm(nmeasg,One,np,nb,Ui,Ub,Af)
!         if (pid .eq. 0) then
!           call create_Afb_pckfddm(nmeasg,One,nbg,Ub_g,Afb)
!         end if

! !        xi(:) = 0.0d0
! !        xi(1) = -1.0d0

! !!!***** Mean *****!!!
!        call construct_coln(pid,Af(1:(np-nb),1),U_gi)
!        CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
!        if (pid .eq. 0) then
!           call construct_coln_boundary(Afb(:,1),U_g)
!           label='../../data/vtkOutputs/out_mean_prior_pce.vtk'
!           call writevtk(label,pg,tg,U_g,npg,ntg)
!           print*, 'Max-mean =', maxval(U_g)
!        end if


! end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! if (pid .eq. 0) then
!     print*, U_wave
!     stop 123
! end if

! if (outputFlag .eq. 1) then

!        if (pid .eq. 0) then

!               filename = '../../data/vtkOutputs/U_wave.dat'  !!**
!               open(unit=2,file=filename,status='replace')

!             do i= 1,5   !tcount
!               call int2str(extension,i,2)

!               ! U_wtemp = U_wave(:,i)

!               label='../../data/vtkOutputs/out_wave_' // trim(extension) // '.vtk'
!               call writevtk(label,pg,tg,U_wave(:,i),npg,ntg)

!               write(unit=2,fmt='(1ES17.8E3)') U_wave(:,i)

!             end do

!             close(2)

!        end if


! end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! if (pid .eq. 0) print*, 'U_g is', U_g



!!-----------------------------------------------------------------------------------------------
!! Dynamic array de-allocation
   if (pid .eq. 0) then
      deallocate(pg,eg,tg)
   end if

    deallocate(Ub,Ui,Ub_g,U_g,U_gi, Af, Afb)
    !!deallocate(omegas,multipliers,xi,mIndex,sIndex)
    ! deallocate(var_g) !,uveci,uvecb,,measnodes,d)                          !!** PCKF ***

   CALL MPI_FINALIZE(ierr)


END PROGRAM main
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
