!!Solvers Module: by Ajit Desai, 7/2015: PETSc Version : March/2016
!!contains: various solver for stochastic-DDM based PCGM with different Preconditioners
!
!input :
!    DDM-Data  : Descretized blocks of matrices & vectors
!    Mesh-Data : Global & local moes data
!    MPI-Data  : For MPI routines
!    PCE-Data  : For stochastic formulation
!
!output:
!    Ui, Ub, Urs, Ucs : Local solution vectors
!    Ub_g, Uc         : Global solution vectors
!
!The following include statements are required for KSP Fortran programs:
!    petscsys.h       : base PETSc routines
!    petscvec.h       : vectors
!    petscmat.h       : matrices
!    petscksp.h       : Krylov subspace methods
!    petscpc.h        : preconditioners
!
!Record of revisions:
!     Date        Programmer     Description of change
!  ==========     ==========     =====================
!  23/05/2014        AD          solve_lpcgm: Stochastic-Lumped-PCGM solver(L)
!  12/06/2014        AD          solve_wlpcgm: Stochastic-Weighted-Lumped-PCGM Solver (WL)
!  16/06/2014        AD          solve_nnpcgm: Stochastic-OneLevel-PCGM Solver added(NN)
!  23/07/2014        AD          solve_nncpcgm: Stochastic-TwoLevel-PCGM Solver added(NNC)
!  07/08/2015        AD          solve_coarsepcgm: Coarse-PCGM subroutine added (D-CG)
!  12/08/2015        AD          solve_fetidp: FETI-DP with Dirichlet-PCGM Solver (FETIDP)
!  15/03/2016        AD          PETSc_stonncpcgm: Sparse PETSc-Two-Level-PCGM Solver (NNC)
!  17/03/2016        AD          solve_coarsepcgm: Sparse PETSc-Coarse-PCGM Solver (S-CG)




!!! EDITED BY SUDHI P V - ACOUSTIC WAVE PROPAGATION - STOCHASTIC DDM - TWO LEVEL
!!! 2020/05/24

!!!! P infront of a quantity represent PETSc vec or matrices
!!!! i, b at the end of a quantity represent interior and interface nodes
!!
!!----------------------------------------------------------------------------------------------
MODULE PETScSolvers

use PETScommon
use PETScAssembly
!!use common
!!use myCommon

implicit none

contains


!!!----------------------------------------------------------------------------------------------
!!!!!!!!!! Two-Level Neumann-Neumann-Coarse (NNC) PCGM Solver : Algorithm 5 from report !!!!!!!!!
!!!----------------------------------------------------------------------------------------------
!!!! DEVELOPER  NOTES  :

!!!! ALL THE QUANTITIES ARE STOCHASTIC UNLESS ITS LOADED FROM FENICS....
!!!! THERE MIGHT BE DIFFERENCE IN NAMING CONVENTION USED FROM DETERMINISTIC CODE SINCE THE VARIABLE NAMES BECOMES TOO BIG
!!!!  IN ORDER TO EASILY CONVERT THE CODE TO TIME DEPENDANT CASE AND USE THE LEGACY CODE, MANY VARIABLES ARE KEPT WITH THE SAME NAME

!!!! VARIABLE WITH P PREFIX MEANS ITS PETSC OBJECT, BUT SOME OTHERS LIKE Asgg, Asgi etc are also PETSC MATRICES BUT DON'T HAVE P PREFIX TO HAVE
!!!!! CONSISTENCY WITH THE PREVIOUSLY WRITTEN CODE




SUBROUTINE PETSc_stonncpcgm(pid,pg,tg,np,nb,npg,ntg,nbg,npcein,npceout,nr,nc,nbgc,ncijk,ijk,cijk,Ui,Ub,Ub_g, &
                         mallocsCals,maxiter,tol,T, tcount, deltaT, beta_NB, gamma_NB,outputFlag)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

!!!-------------------------------------------------------------------------------------------
!! PETSc-Variables Declaration
    Vec              PetVeci,PetVecir,PetVecc,PetVecb,PetVecr   !! PETSc Vectors
    Vec              SolVeci,SolVecir,SolVecc,SolVecb,SolVecr
    Vec              PetVecbg,PetVeccg,PetVecc2,PetVecc3, PPb

    Vec              Fsi,Fsg          !!! Force vector
    Vec              Pu0i, Pv0i, Pa0i    !!! Initial displacement, velocity and acceleration
    Vec              Pu0b, Pv0b, Pa0b    !!! Initial displacement, velocity and acceleration
    Vec              PMMuvi, PMMuvb   !! Mass Multiplier
    Vec              PCMuvi, PCMuvb   !! Damping Multiplier

    Vec              Mni, Mnb !!! Force vector contribution from Mass times previous soultions
    Vec              Cni, Cnb !!! Force vector contribution from Damping times previous soultions

    Mat              RsMat,DsMat,BcMat,RrMat,RcMat              !! PETSc Matrix
    Mat              Asmat,Asii,Asgg,Asrr,Ascc
    Mat              Asgi,Asri,Asci,Ascr

    Mat              Msmat,Msii,Msgg,Msrr,Mscc             !!! PETSc Decomposed Mass and Damping Matrices
    Mat              Msgi,Msri,Msci,Mscr

    Mat              Csmat,Csii,Csgg,Csrr,Cscc
    Mat              Csgi,Csri,Csci,Cscr

    KSP              kspAm,kspAi,kspAc                          !! PETSc solver context
    PC               pcAm,pcAi,pcAc                             !! PETSc preconditioner

    PetscErrorCode   ierr

    PetscLogStage    stages(5)
!!!-------------------------------------------------------------------------------------------
!! Fortran-Variables Declaration
    integer :: pid, i, maxiter, mallocsCals
    integer :: np, ni, nb, nbg, npg, ntg
    integer :: nip, nbp, nrp, ncp, nirp, nbgp, nbgcp
    integer :: nc, nr, nbgc, nir, npceout, ncijk, npcein!casep,
    integer :: One, Zero, outputFlag

    integer, dimension((np-nb)*npceout) :: tempi                !! Temp arrays
    integer, dimension((np-nc)*npceout) :: tempir               !! require for VecGetVales
    integer, dimension(nc*npceout)      :: tempc
    integer, dimension(nr*npceout)      :: tempr
    integer, dimension(nb*npceout)      :: tempb
    integer, dimension(nbg*npceout)     :: tempbg
    integer, dimension(nbgc*npceout)    :: tempcg

    integer, dimension(ncijk,3)         :: ijk

    ! integer, allocatable, dimension(:) :: tb            !! Array of indices for boundary nodes below interior nodes
    !integer, dimension(ncijk,3)         :: ijk
    ! integer, dimension(3,ne)            :: e
    ! integer, dimension(3,nt)            :: t

    !!double precision, dimension(nomga)           :: omegas, multipliers
    double precision, dimension(nbg*npceout)     :: Ub_g, rb_g, Qb_g, Utemp_g
    double precision, dimension(nbg*npceout)     :: Pb_g, Zb_g, RZb
    double precision, dimension(nbgc*npceout)    :: Dc, Zc, DcBc
    double precision, dimension(nr*npceout)      :: Vri, FrS
    !!double precision, dimension(2,2)             :: dbounds
    double precision, dimension(ncijk)           :: cijk
    double precision, dimension((np-nc)*npceout)  :: ybc
    double precision, dimension(nb*npceout)      :: Ub
    double precision, dimension((np-nb)*npceout)  :: Ui

!!!! Vectors used for assembling global solition
    double precision, allocatable, dimension(:)    :: U_g
    double precision, allocatable, dimension(:)    :: U_gi
    ! double precision, dimension(2,np)            :: p

    double precision :: NegOne, PosOne   !!const_diff, sigma, amp, N
    double precision :: rho_next, rho_curr, alpha, beta, err, tol !! residual

    !!integer, dimension(ndim,npceout) :: mIndex
    !!integer, dimension(2,ndim)       :: sIndex


    !! Newmark Beta/Time discretization Parameters

    double precision :: T, deltaT, beta_NB, gamma_NB

    double precision, allocatable, dimension(:) :: u0i, u0b
    double precision, allocatable, dimension(:) :: v0i, v0b
    double precision, allocatable, dimension(:) :: a0i, a0b

    double precision, dimension(:,:) :: pg
    ! integer, dimension(:,:)          :: eg
    integer, dimension(:,:)          :: tg

    integer         :: tcount, nbcount

    ! double precision, allocatable, dimension(npg,tcount) :: U_wave                                 !!! Wave solution for all time steps

    double precision :: tolc = 1.0d-6                           !! Coarse-solver tolerence
    double precision, allocatable, dimension(:,:) :: Af, Afb
    double precision, allocatable, dimension(:)   :: var_g

    character(len=255) :: label, extension, filename, extension_pc

    ! stop 123
!!!-------------------------------------------------------------------------------------------
!! PETSc-Initialize
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)


    ni   = np-nb           !! det-interior nodes
    nir  = np-nc           !! det-interior+remaining nodes
    nip  = (ni*npceout)    !! sto-interior nodes
    nbp  = (nb*npceout)    !! sto-local-boundary nodes
    ncp  = (nc*npceout)    !! sto-local-corner nodes
    nrp  = (nr*npceout)    !! sto-local-remaining nodes
    nirp = (nir*npceout)   !! sto-interior+remaining nodes
    nbgp = (nbg*npceout)   !! sto-global-boundary nodes
    nbgcp= (nbgc*npceout)  !! sto-global-corner nodes

    allocate(u0i(ni*npceout), u0b(nb*npceout), v0i(ni*npceout), v0b(nb*npceout), a0i(ni*npceout), a0b(nb*npceout))
    allocate(U_g(npg),U_gi(npg))
    allocate(Af(np,npceout),Afb(nbg,npceout),var_g(npg))
    ! allocate(U_wave(npg, tcount))

!! Initiate Array's
    rho_next = 0.0d0
    rho_curr = 0.0d0
    NegOne   = -1.0
    PosOne   = 1.0
    One      = 1
    Zero     = 0

!!!-------------------------------------------------------------------------------------------
!! Assemble only Subdomain-Level Deterministic matrices to calculate Mallocs

    IF (mallocsCals .eq. 1) then

        if (pid .eq. 0) print*, '------------------------------------------------------'
        if (pid .eq. 0) print*, 'Initializing Two-Level-NNC-PCGM Mallocs Calculation...'
        if (pid .eq. 0) print*, '------------------------------------------------------'

        !call GetMallocs(pid,p,e,t,np,ne,nt,nb,nc,nr,ndim,npcein,npceout,nomga,ncijk,ijk,&
        !                cijk,ni,nip,nbp,ncp,nrp,nirp,dbounds, const_diff, omegas, &
        !                mIndex, sIndex, multipliers, sigma, casep, ierr)

        call GetMallocsShort(pid,np,nb,nc,nr,npceout,nip,nbp,ncp,nrp,nirp,ierr)


    Ui(:)  = 0.0d0
    Ub(:)  = 0.0d0
    Ub_g(:)= 0.0d0

    END IF

!!!---------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, '------------------------------------------------'
    if (pid .eq. 0) print*, 'Initializing STOCHASTIC-Two-Level-NNC-PCGM...'
    if (pid .eq. 0) print*, '------------------------------------------------'



!!!! Only for Coarse Level Time calculation
    call PetscLogStageRegister("Coarse Solver",stages(2),ierr)
    call PetscLogStageRegister("Post Processing",stages(3),ierr)
    call PetscLogStageRegister("MPI-Main",stages(4),ierr)
    call PetscLogStageRegister("MPI-Coarse",stages(5),ierr)
!!!!!!!!!!!!!!!!!!!!!!!!

    ! call GetMatSeqTwoLevel(pid,PAmat,PAii,PAgg,PAgi,PAri,PArr,PAcc,PAci,PAcr,PMii,PMgg,PMgi,PMri,PMrr,PMcc,PMci,PMcr,&
                           ! PCii,PCgg,PCgi,PCri,PCrr,PCcc,PCci,PCcr,np,nb,nc,nr,nbgc,ierr)

      call GetStoMatSeq( pid, Asmat, Asii, Asgg, Asgi, Asri, Asrr, Ascc, Asci, Ascr, Msmat, Msii,Msgg,Msgi,Msri,Msrr,Mscc,Msci,Mscr, &
                             Csmat,Csii,Csgg,Csgi,Csri,Csrr,Cscc,Csci,Cscr,np,nb,nc,nr,npcein,npceout,ncijk,ijk,cijk,ni,nip,nbp,ncp,nrp,nirp,&
                             beta_NB, gamma_NB, deltaT,ierr)



!!!!!! ************ TEST CASE FOR STATIC CASE********************!!!!!!!!!!!!!!!!!!!!
      ! call StoMatSeqOneLevelShort( pid, Asmat, Asii, Asgg, Asgi, Asri, Asrr, Ascc, Asci, Ascr, &
      !                      np,nb,nc,nr,npcein,npceout,ncijk,ijk,cijk,ni,nip,nbp,ncp,nrp,nirp,ierr)




      ! call GetStoForceVecSeq(pid, nbcount, Fsi, Fsg, ni, nb, nip, nbp, ierr)


!!!!!! **************                          ********************!!!!!!!!!!!!!!!!!!!!


!!!!!!!                                                !!!!!!!!!!!!!!!!!

      ! call GetKTstochastic(pid, Asmat, Asii, Asgg, Asgi, Asri, Asrr, Ascc, Asci, Ascr, Msmat,Msii,Msgg,Msgi,Msri,Msrr,Mscc,Msci,Mscr, &
      !                        Csmat,Csii,Csgg,Csgi,Csri,Csrr,Cscc,Csci,Cscr,np,nb,nc,nr,ni,nip,nbp,ncp,nrp,nirp,beta_NB, deltaT,gamma_NB, ierr)

     call GetStoInitialCond(pid, Pu0i, Pv0i, Pa0i,Pu0b, Pv0b, Pa0b, np, nb, npceout, ierr)


    call GetVecSeqTemp2(nirp,PetVecir,Solvecir,tempir,ierr)
    call GetVecSeqTemp2(nip,PetVeci,Solveci,tempi,ierr)
    call GetVecSeqTemp2(nbp,PetVecb,Solvecb,tempb,ierr)
    call GetVecSeqTemp2(ncp,PetVecc,Solvecc,tempc,ierr)
    call GetVecSeqTemp2(nrp,PetVecr,Solvecr,tempr,ierr)


    call GetVecSeqDummy(nbgcp,PetVeccg,tempcg,ierr)
    call GetVecSeqDummy(nbgp,PetVecbg,tempbg,ierr)


    call GetVecSeqTemp1(ncp,PetVecc2,ierr)
    call GetVecSeqTemp1(ncp,PetVecc3,ierr)
    call GetVecSeqTemp1(nbp,PetVecb,ierr)
    call GetVecSeqTemp1(nip,PetVeci,ierr)
    call GetVecSeqTemp1(nbp,PPb,ierr)


    call GetStoInitialCondSubdomain(pid,ni,nb,npceout,tempi,tempb,Pu0i,Pv0i,Pa0i,Pu0b,Pv0b,Pa0b,u0i,u0b,v0i,v0b,a0i,a0b, ierr)


    if (pid .eq. 0) print*, 'total time steps is ',tcount

    !!!      Newmark Beta Loop Starts        !!!!!

DO nbcount = 1,tcount

    if (pid .eq. 0) print*, '---------------------------------------------------------------'
    if (pid .eq. 0) print*, '       Newmark Beta Iteration number :',nbcount
    if (pid .eq. 0) print*, '---------------------------------------------------------------'

    ! print*, nbg, pid, T, deltaT, beta_NB, gamma_NB

    call GetStoMassDampMultiplier(pid,nbcount,nip,nbp,deltaT,beta_NB,gamma_NB, u0i,u0b,v0i,v0b,a0i,a0b,tempi,tempb, PMMuvi, PMMuvb, PCMuvi, PCMuvb,ierr)

    ! print*, 'process id', pid
    ! print*, u0b
    ! stop 123
    ! call VecView(PMMuvi,PETSC_VIEWER_STDOUT_SELF, ierr);


    call GetStoForceVecSeq(pid, nbcount, Fsi, Fsg, ni, nb, nip, nbp, ierr)


    call GetStoTransientForce(pid,nbcount,nip,nbp,Msii,Msgg,Msgi,Csii,Csgg,Csgi, PMMuvi, PMMuvb, PCMuvi, PCMuvb, Mni, Mnb, Cni, Cnb, Fsi, Fsg,ierr)


    !!!! MILESTONE VERIFIED Fsi -SAME AS DETERMINISTIC FOR FIRST ITERATION

    !!!! MILESTONE VERIFIED Fsg -SAME AS DETERMINISTIC FOR FIRST ITERATION

    ! if (pid .eq. 0) then
    !     ! if(nbcount .eq. 2) then
    !         print*, 'procees is', pid
    !         print*, 'ii is', nip
    !         print*, 'Asii is'
    !         call MatView(Asii,PETSC_VIEWER_STDOUT_SELF, ierr);
    !         stop 123
    !     ! endif
    ! end if

!!!--------------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, '---------------------------------------------------------------'
    if (pid .eq. 0) print*, 'Using Sparse-Iterative-PETSc-KSP Solver For Interior Problem...'
    if (pid .eq. 0) print*, '---------------------------------------------------------------'

!!---------------------------------------------------------------------------------------------
!!  Uncomment this section for Execution Time Calculations !!!!!
!!---------------------------------------------------------------------------------------------
    !integer :: c1,c2,cr
    !double precision :: time1
    !double precision :: time2
    !call system_clock(count_rate=cr)

!!!---------------------------------------------------------------------------------------------
!! Pre-Iteration : 0th Iteration : Step 1 to 7. - Finding the Initial Residual r_gamma = g_gamma - S u_gamma0
!!!---------------------------------------------------------------------------------------------
    !! Step 3 : Solve
    call PETSc_KSP(kspAi,pcAi,Asii,Fsi,SolVeci,ierr)

    !! Step 4 : Compute
    call MatMult(Asgi,SolVeci,SolVecb,ierr)       !! A(1) x(2) = b(3)
    call VecAYPX(SolVecb,NegOne,Fsg,ierr);        !! (y,a,x):y = x + a*y

    call GetRs(pid,nb,nbg,npceout,RsMat)                   !! Sparse Rs Operator

    ! call MatView(RsMat, PETSC_VIEWER_STDOUT_SELF, ierr)
    ! stop 123

    call MatMultTranspose(RsMat,SolVecb,PetVecbg,ierr)     !! A(1)'x(2) = b(3)
    call VecGetValues(PetVecbg, nbgp, tempbg, RZb, ierr)   !! PETSc to Fortran Vec


    call PetscLogStagePush(stages(4),ierr)

    !! STEP 5 : Gather                                     !! Use PETSc_Reduce in future
    CALL MPI_ALLREDUCE(RZb,rb_g,nbgp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)


    call PetscLogStagePop(ierr)
    !!!! MILESTONE Rbg VERIFIED for first iteration - same as deterministic for mean term

    ! if(pid .eq. 0) then
    !     print*,'rbg is', rb_g
    !     stop 123
    ! end if



    ! stop 123
! !!-----------------------------------------------------------------------------------------------
! !!  Uncomment this section for Residual Calculation !!!!!
! !!-----------------------------------------------------------------------------------------------
!     !if (pid .eq. 0) then
!     !residual = SQRT(dot_product(rb_g,rb_g))
!     !open (unit=7,file="residual.dat",action="write",status="replace")
!     !print*,residual
!     !write (7,*), residual
!     !close (7)
!     !end if


! !!-----------------------------------------------------------------------------------------------
! !!  STEP 6 : Two-Level NNC : Algorithm 6 from the report,
!     ALGORITHM 4 ajit thesis, step 2 - 5
!!!!! In order to do Two-Level PC we need to solve course problem which is Fcc uc = dc  Eq. 3.26
!!!!   dc = Eq. 3.28
!!!!! Finding out the residual force vector for remaining F_r and corner nodes F_c found below.
! !!-----------------------------------------------------------------------------------------------
    ! SUB-STEP 3 : Scatter
    call SetVecSeq(nbgp,rb_g,PetVecbg,ierr)                !! Fortran to PETSc Vec
    call MatMult(RsMat,PetVecbg,SolVecb,ierr)              !! Sparse PETSc Mat*Vec

    !! SUB-STEP 4 : Weighting
    call GetDs(pid,nb,npceout,DsMat)                       !! Sparse Ds Operator
    call MatMult(DsMat,SolVecb,PetVecb,ierr)

    !! SUB-STEP 5 & 6 : Compute
    call getRr(pid, nb, nr, npceout, RrMat)                !! Sparse RrS Operator
    call MatMult(RrMat,PetVecb,PetVecr,ierr)               !!


    ! print*, 'process id', pid
    ! print*, nr
    ! call VecView(PetVecr,PETSC_VIEWER_STDOUT_SELF, ierr);
    ! print*,'end vector'
    ! stop 123

    call VecGetValues(PetVecr, nrp, tempr, FrS, ierr)      !! Convert FrS to PETScVec-
                                                           !!-For future modification
    call getRc(pid, nb, nc, npceout, RcMat)                !! Sparse RcS Operator
    call MatMult(RcMat,PetVecb,PetVecc2,ierr)            !!!! F_c Subdomain-Level

! !!-------------------------This is the solving stage v1 = s_rr \f_r -------------------------------------------
    !! SUB-STEP 7 : Solve
    ybc(:) = 0.0d0                                         !! Check if necessory or not
    ybc((nip+1):nirp) = FrS

    ! stop 123

    ! print*, 'process id', pid
    ! print*, nir
    ! call MatView(PAmat,PETSC_VIEWER_STDOUT_SELF, ierr);
    ! print*,'end matrix'
    ! stop 123

    call SetVecSeq(nirp,ybc,PetVecir,ierr)
    call PETSc_KSP(kspAm,pcAm,Asmat,PetVecir,SolVecir,ierr)
    call VecGetValues(SolVecir, nirp, tempir, ybc, ierr)   !! Re-using ybc : Check for error

    Vri = ybc((nip+1):nirp)                                !! Re-using Vr1 : Check for error

!!!------------------- PMV:Sec1: Scr*v1 = [Acr - (Aci *inv(Aii)*Air)]*v1- -------------------!!!
 !!!! Part of Eq. 3.28 Ajits thesis !!!!!!!!!!!!!

    ! SUB-STEP 8 : Compute
    call SetVecSeq(nrp,Vri,PetVecr,ierr)
    call MatMultTranspose(Asri,PetVecr,PetVeci,ierr)       !! A(1)'x(2) = b(3)
    call KSPSolve(kspAi,PetVeci,SolVeci,ierr)
    call MatMult(Asci,SolVeci,SolVecc,ierr)
    call MatMult(Ascr,PetVecr,PetVecc,ierr)
    call VecAXPY(PetVecc,NegOne,SolVecc,ierr)              !! (y,a,x): y = y-a*x

    !! SUB-STEP 9 : Update.  fc - S_cr * v1
    call VecAYPX(PetVecc,NegOne,PetVecc2,ierr)             !! (y,a,x): y = x-a*y

    !! SUB-STEP 10 : Gather
    call GetBc(pid, nc, nbgc, npceout, BcMat)              !! Sparse RcS Operator
    call MatMultTranspose(BcMat,PetVecc,PetVeccg,ierr)
    call VecGetValues(PetVeccg, nbgcp, tempcg, DcBc, ierr)


    call PetscLogStagePush(stages(4),ierr)
    CALL MPI_ALLREDUCE(DcBc,Dc,nbgcp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call PetscLogStagePop(ierr)
    ! stop 123
! !!-----------------------------------------------------------------------------------------------
! !! SUB-STEP 12 : Coarse Problem: LUMPED-PCGM : Algorithm 7 from the report.
! !!                          Fcc uc = dc
! !!
! !!-----------------------------------------------------------------------------------------------

    call SetPETScKSP(kspAc,pcAc,Ascc,ierr)

    call PetscLogStagePush(stages(2),ierr)

    call solve_coarsepcgm(pid,nirp,ncp,nrp,nip,nbgcp,BcMat,kspAm,kspAi,kspAc,pcAm,pcAi,pcAc, &
                        tempir, tempr, tempcg, PetVecc2, PetVecc3, PetVeccg, ybc, Vri,DcBc, &
                        PetVeci,SolVeci,PetVecir,SolVecir,PetVecc,Solvecc,PetVecr,SolVecr, &
                        Asci,Asri,Ascr,Ascc,Dc,Zc,tolc, stages(5))

    call PetscLogStagePop(ierr)

! !!-----------------------------------------------------------------------------------------------


        !! SUB-STEP 14 : Scatter.  Algorithm 4 - Step 12-18
    call PetscLogStagePush(stages(4),ierr)
    CALL MPI_BCAST(Zc,nbgcp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call PetscLogStagePop(ierr)
    ! if (nbcount .eq. 1) then
    !     if(pid .eq. 0) then
    !         print*, 'process id', pid
    !         print*, 'nbgcp is', nbgcp
    !         print*, 'Zc is', Zc
    !     end if
    ! ! end if
    ! stop 123

    call SetVecSeq(nbgcp,Zc,PetVeccg,ierr)
    call MatMult(BcMat,PetVeccg,PetVecc,ierr)

!!!------------------- PMV:Sec1: Src*v1 = [Arc - (Ari *inv(Aii)*Aic)]*v1- -------------------!!!
    !! SUB-STEP 15 : Compute
    call MatMultTranspose(Asci,PetVecc,PetVeci,ierr)
    call KSPSolve(kspAi,PetVeci,SolVeci,ierr)
    call MatMult(Asri,SolVeci,SolVecr,ierr)
    call MatMultTranspose(Ascr,PetVecc,PetVecr,ierr)
    call VecAXPY(PetVecr,NegOne,SolVecr,ierr)
    call VecGetValues(PetVecr, nrp, tempr, Vri, ierr)

    ! stop 123

    !! SUB-STEP 16 : Update
    Vri = FrS - Vri

!!-----------------------------------------------------------------------------------------------
!! Solving Naumann-Naumann Problem: solve(nrp,1,Srr,ZrS,v2)
!!-----------------------------------------------------------------------------------------------
    !! SUB-STEP 17 : Solve
    ybc(:) = 0.0d0
    ybc((nip+1):nirp) = Vri

    call SetVecSeq(nirp,ybc,PetVecir,ierr)
    call KSPSolve(kspAm,PetVecir,SolVecir,ierr)
    call VecGetValues(SolVecir, nirp, tempir, ybc, ierr)

    Vri = ybc((nip+1):nirp)

!!-----------------------------------------------------------------------------------------------
    !!!! Step 17 Algorithm 4.
!!-----------------------------------------------------------------------------------------------
    !! SUB-STEP 18 & 19 : Compute
    call SetVecSeq(nrp,Vri,PetVecr,ierr)
    call MatMultTranspose(RrMat,PetVecr,PetVecb,ierr)
    call MatMultTranspose(RcMat,PetVecc,SolVecb,ierr)

    !! SUB-STEP 20 : Update
    call VecAXPY(SolVecb,PosOne,PetVecb,ierr)

    !! SUB-STEP 21 : Weighting
    call MatMult(DsMat,SolVecb,PetVecb,ierr)

    !! SUB-STEP 22 : Gather
    call MatMultTranspose(RsMat,PetVecb,PetVecbg,ierr)
    call VecGetValues(PetVecbg, nbgp, tempbg, RZb, ierr)

    call PetscLogStagePush(stages(4),ierr)
    CALL MPI_REDUCE(RZb,Zb_g,nbgp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call PetscLogStagePop(ierr)
    ! stop 123
!!-----------------------------------------------------------------------------------------------

!!!!! Zb_g is the Two level Neuamann Neumann Preconditioned residual
!!-----------------------------------------------------------------------------------------------
   !! STEP 8 & 9 : Compute
    if (pid .eq. 0) Pb_g = Zb_g
    if (pid .eq. 0) rho_next = dot_product(rb_g,Zb_g)  ! rho_next is actually delta in Algorithm 1 Ajits Thesis

!!!---------------------------------------------------------------------------------------------
!! PCGM Iteration : Main For Loop for each iteration: Step 10 to 26
!!!---------------------------------------------------------------------------------------------
    Ub_g(:)  = 0.0d0
    DO i = 1,maxiter

        !! STEP 12 : Scatter


        call PetscLogStagePush(stages(4),ierr)
        CALL MPI_BCAST(Pb_g,nbgp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)

        call SetVecSeq(nbgp,Pb_g,PetVecbg,ierr)
        call MatMult(RsMat,PetVecbg,PPb,ierr)

!!!! Algorithm 3 - step 4 - 7
        !! STEP 13 : Solve
        call MatMultTranspose(Asgi,PPb,PetVeci,ierr)
        call KSPSolve(kspAi,PetVeci,SolVeci,ierr)
        call MatMult(Asgi,SolVeci,SolVecb,ierr)
        call MatMult(Asgg,PPb,PetVecb,ierr)
        call VecAXPY(PetVecb,NegOne,SolVecb,ierr)

        !! STEP 15 : Gather step 9 Algorithm 3
        call MatMultTranspose(RsMat,PetVecb,PetVecbg,ierr)
        call VecGetValues(PetVecbg, nbgp, tempbg, RZb, ierr)

        call PetscLogStagePush(stages(4),ierr)
        CALL MPI_REDUCE(RZb,Qb_g,nbgp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)
!!!!!!!!        Qb_g is the Parallel Matrix Vector product output in the Outer PCGM Main iteration

    !!-----------------------------------------------------------------------------------------------
        !! STEP 16, 17, 18 & 19 : Compute & Update
        if (pid .eq. 0) then
            rho_curr = rho_next
            alpha = rho_curr/dot_product(Qb_g,Pb_g)
            err = alpha*alpha*dot_product(Pb_g,Pb_g)/dot_product(Ub_g,Ub_g)
            Ub_g = Ub_g + alpha*Pb_g
            rb_g = rb_g - alpha*Qb_g

!!-----------------------------------------------------------------------------------------------
!!  Uncomment this section for Residual Calculation !!!!!
!!-----------------------------------------------------------------------------------------------
            !residual = SQRT(dot_product(rb_g,rb_g))
            !print*,residual
            !open(unit=7,file="residual.dat", position="append", status='old')
            !write (7,*), residual
            !close (7)
!!-----------------------------------------------------------------------------------------------

            !! STEP 20 : Exit
            print*, '----------------'
            print*, 'Main Iteration #',i, ', relative error of ',err
            print*, '----------------'
            !!print*, 'and sum of Ub of ', sum(Ub_g)
        end if

        call PetscLogStagePush(stages(4),ierr)
        CALL MPI_BCAST(err,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)

        if (err .lt. tol) exit

!!-----------------------------------------------------------------------------------------------
!!  STEP 21 : Two-Level NNC : Algorithm 4 from the Ajits Thesis
!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 3 : Scatter. steps 3 to 5
        call PetscLogStagePush(stages(4),ierr)
        CALL MPI_BCAST(rb_g,nbgp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)

        call SetVecSeq(nbgp,rb_g,PetVecbg,ierr)
        call MatMult(RsMat,PetVecbg,SolVecb,ierr)

        !! SUB-STEP 4 : Weighting
        call MatMult(DsMat,SolVecb,PetVecb,ierr)

        !! SUB-STEP 5 & 6 : Compute
        call MatMult(RrMat,PetVecb,PetVecr,ierr)
        call VecGetValues(PetVecr, nrp, tempr, FrS, ierr)
        call MatMult(RcMat,PetVecb,PetVecc2,ierr)

!!-----------------------------------------------------------------------------------------------

!!!

!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 7 : Solve
        ybc(:) = 0.0d0
        ybc((nip+1):nirp) = FrS   ! Force vector which needs to be used after the coarse solve

!!!!    Local Neumann Solve -- Algorithm 5
        call SetVecSeq(nirp,ybc,PetVecir,ierr)
        call KSPSolve(kspAm,PetVecir,SolVecir,ierr)
        call VecGetValues(SolVecir, nirp, tempir, ybc, ierr)

        Vri = ybc((nip+1):nirp)

!!------------------- PMV:Sec1: Scr*v1 = [Acr - (Aci *inv(Aii)*Air)]*v1- -------------------!!!
 !!!! Part of Eq. 3.28 Ajits thesis !!!!!!!!!!!!! Algorithm 4

        !! SUB-STEP 8 : Compute
        call SetVecSeq(nrp,Vri,PetVecr,ierr)
        call MatMultTranspose(Asri,PetVecr,PetVeci,ierr)
        call KSPSolve(kspAi,PetVeci,SolVeci,ierr)
        call MatMult(Asci,SolVeci,SolVecc,ierr)
        call MatMult(Ascr,PetVecr,PetVecc,ierr)
        call VecAXPY(PetVecc,NegOne,SolVecc,ierr)

        !! SUB-STEP 10 : Update. step 8
        call VecAYPX(PetVecc,NegOne,PetVecc2,ierr)     !!! PetVecc2 is actually F_c Residual force in corner nodes

        !! SUB-STEP 10 : Gather
        call MatMultTranspose(BcMat,PetVecc,PetVeccg,ierr)
        call VecGetValues(PetVeccg, nbgcp, tempcg, DcBc, ierr)

        call PetscLogStagePush(stages(4),ierr)
        CALL MPI_ALLREDUCE(DcBc,Dc,nbgcp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)

        ! stop 123
!!-----------------------------------------------------------------------------------------------

!!-----------------------------------------------------------------------------------------------
!!  SUB-STEP 12 : Coarse Problem: LUMPED-PCGM : Algorithm 9 from the report
!!-----------------------------------------------------------------------------------------------
        call PetscLogStagePush(stages(2),ierr)

        call solve_coarsepcgm(pid,nirp,ncp,nrp,nip,nbgcp,BcMat,kspAm,kspAi,kspAc,pcAm,pcAi,pcAc,&
                            tempir, tempr, tempcg, PetVecc2, PetVecc3, PetVeccg, ybc, Vri,DcBc,&
                            PetVeci,SolVeci,PetVecir,SolVecir,PetVecc,Solvecc,PetVecr,SolVecr,&
                            Asci,Asri,Ascr,Ascc,Dc,Zc,tolc,stages(5))

        call PetscLogStagePop(ierr)

        ! stop 123
!!-----------------------------------------------------------------------------------------------

        !!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 14 : Scatter
        call PetscLogStagePush(stages(4),ierr)
        CALL MPI_BCAST(Zc,nbgcp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)


        call SetVecSeq(nbgcp,Zc,PetVeccg,ierr)
        call MatMult(BcMat,PetVeccg,PetVecc,ierr)

!!!------------------- PMV:Sec1: Src*v1 = [Arc - (Ari *inv(Aii)*Aic)]*v1- -------------------!!!
        !! SUB-STEP 15 : Compute
        call MatMultTranspose(Asci,PetVecc,PetVeci,ierr)
        call KSPSolve(kspAi,PetVeci,SolVeci,ierr)
        call MatMult(Asri,SolVeci,SolVecr,ierr)
        call MatMultTranspose(Ascr,PetVecc,PetVecr,ierr)
        call VecAXPY(PetVecr,NegOne,SolVecr,ierr)
        call VecGetValues(PetVecr, nrp, tempr, Vri, ierr)

        !! SUB-STEP 16 : Update.  step 15 Algorithm 4
        Vri = FrS - Vri

!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 17 : Solve
        ybc(:) = 0.0d0
        ybc((nip+1):nirp) = Vri

        call SetVecSeq(nirp,ybc,PetVecir,ierr)
        call KSPSolve(kspAm,PetVecir,SolVecir,ierr)
        call VecGetValues(SolVecir, nirp, tempir, ybc, ierr)

        Vri = ybc((nip+1):nirp)

!!-----------------------------------------------------------------------------------------------
        !! SUB-STEP 18 & 19 : Compute
        call SetVecSeq(nrp,Vri,PetVecr,ierr)
        call MatMultTranspose(RrMat,PetVecr,PetVecb,ierr)
        call MatMultTranspose(RcMat,PetVecc,SolVecb,ierr)

        !! SUB-STEP 20 : Update
        call VecAXPY(SolVecb,PosOne,PetVecb,ierr)

        !! SUB-STEP 21 : Weighting
        call MatMult(DsMat,SolVecb,PetVecb,ierr)

        !! SUB-STEP 22 : Gather
        call MatMultTranspose(RsMat,PetVecb,PetVecbg,ierr)
        call VecGetValues(PetVecbg, nbgp, tempbg, RZb, ierr)

        call PetscLogStagePush(stages(4),ierr)
        CALL MPI_REDUCE(RZb,Zb_g,nbgp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)
!!-----------------------------------------------------------------------------------------------
        !! STEP 22, 23 & 24 : Compute & Update
        if (pid .eq. 0) then
            rho_next = dot_product(rb_g,Zb_g)
            beta = rho_next/rho_curr
            Pb_g = Zb_g + beta*Pb_g
        end if


    END DO


    ! stop 123
!!!---------------------------------------------------------------------------------------------
!! Post Iteration : To calculate local interior/interface solutions : Step 27 to 30
!!!---------------------------------------------------------------------------------------------
    !! STEP 28 : Scatter. Eq. 3.17 Ajits Thesis ,
    call PetscLogStagePush(stages(4),ierr)
    CALL MPI_BCAST(Ub_g,nbgp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call PetscLogStagePop(ierr)

    call SetVecSeq(nbgp,Ub_g,PetVecbg,ierr)
    call MatMult(RsMat,PetVecbg,PetVecb,ierr)

    !! STEP 29 : Compute & Solve
    call MatMultTranspose(Asgi,PetVecb,PetVeci,ierr)
    call VecAYPX(PetVeci,NegOne,Fsi,ierr)          !! VERIFY THIS STEP
    call KSPSolve(kspAi,PetVeci,SolVeci,ierr)

    call VecGetValues(SolVeci, nip, tempi, Ui, ierr)
    call VecGetValues(PetVecb, nbp, tempb, Ub, ierr)    !!!! In order to do the update at subdomain level


    ! if (nbcount .eq. 1) then
    !      if(pid .eq. 0) then
    !          print*, 'process id', pid
    !          print*, Ub_g
    !          stop 123
    !      end if
    !     ! call VecView(PetVecbg,PETSC_VIEWER_STDOUT_SELF, ierr);

    ! end if

    !
    ! print*,'end vector'
    ! stop 123

    !!!!!!!! Forming the Solution Vector for display !!!!!!!!!!!!!!


     call PetscLogStagePush(stages(3),ierr)

if (outputFlag .eq. 1) then

    ! if(mod(nbcount,5) == 0) then

        Af(:,1)  = 0.0d0
        Afb(:,1) = 0.0d0

        call create_Af_pckfddm(Zero,npceout,np,nb,Ui,Ub,Af)  !!!! Function to create full solution but Ui is only for a particular subdomain
        if (pid .eq. 0) then
          call create_Afb_pckfddm(Zero,npceout,nbg,Ub_g,Afb)  !!!! Temporary storing of boundary solution needed for final creation after reduce
        end if


!!!***** Mean *****!!!
       call construct_coln(pid,Af(1:(np-nb),1),U_gi)                   !!!! Putting all the interior nodes in order for display
       CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)        !!!! Summing up all the nodes from all subdoamins which
!         ! results in sum of boundary nodes too which has to be replaced . Interior nodes when summed doesn't effetc each other
       if (pid .eq. 0) then
        ! print*, 'inside mean creation'
          call int2str(extension,nbcount,2)
          call construct_coln_boundary(Afb(:,1),U_g)                !!! Replacing the summed up boundary nodes to correct one
          label='../../data/vtkOutputs/out_wave_mean_' // trim(extension) //'.vtk'
          call writevtk(label,pg,tg,U_g,npg,ntg)
          print*, 'Max-mean =', maxval(U_g)
       end if


       ! print*, 'finished mean '

!!!***** SD *****!!!
       var_g(:) = 0.0d0
       do i = 2,npceout
          call construct_coln(pid,Ui((i-1)*(np-nb)+1:i*(np-nb)),U_gi)
          call MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if (pid .eq. 0) then
             call construct_coln_boundary(Ub_g((i-1)*nbg+1:i*nbg),U_g)
             var_g = var_g + U_g**2
          end if
       end do
       if (pid .eq. 0) then

          call int2str(extension,nbcount,2)
          var_g = sqrt(var_g)
          label='../../data/vtkOutputs/out_wave_sd_' // trim(extension) //'.vtk'
          call writevtk(label,pg,tg,var_g,npg,ntg)
          print*, 'Max-SD   =', maxval(var_g)
       end if


!!***** Higher Order PCE Coefficients *****!!!
    ! do i = 2,npceout

    !     if (mod(i,2) == 0) then

    !       call construct_coln(pid,Af(1:(np-nb),i),U_gi)
    !       CALL MPI_REDUCE(U_gi,U_g,npg,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    !       if (pid .eq. 0) then
    !          call construct_coln_boundary(Afb(:,i),U_g)
    !          call int2str(extension_pc,nbcount,2)
    !          call int2str(extension,i,2)
    !          label='../../data/vtkOutputs/out_pce_' // trim(extension) // '_' // trim(extension_pc) // '.vtk'
    !          call writevtk(label,pg,tg,U_g,npg,ntg)
    !       end if

    !     endif

    ! end do


   ! endif

end if



!     end if
     call PetscLogStagePop(ierr)

    ! print*, 'process id', pid
    ! print*, 'nbcount is', nbcount
    ! print*, v0b
    ! stop 123

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call NewmarkBetaUpdate(pid, nip, nbp, deltaT, beta_NB, gamma_NB, Ui, Ub, u0i, u0b, v0i, v0b, a0i, a0b, ierr)




END DO


    deallocate(u0i, u0b, v0i, v0b, a0i, a0b)
    deallocate(U_g,U_gi)
    deallocate(Af,Afb,var_g)
! !!-----------------------------------------------------------------------------------------------
!!! PETSc-Destroy
    call VecDestroy(PetVeci,ierr)
    call VecDestroy(PetVecir,ierr)
    call VecDestroy(PetVecc,ierr)
    call VecDestroy(PetVeccg,ierr)
    call VecDestroy(PetVecb,ierr)
    call VecDestroy(PetVecbg,ierr)
    call VecDestroy(PetVecr,ierr)
    call VecDestroy(SolVeci,ierr)
    call VecDestroy(SolVecir,ierr)
    call VecDestroy(SolVecc,ierr)
    call VecDestroy(SolVecb,ierr)
    call VecDestroy(SolVecr,ierr)
    call VecDestroy(PetVecc2,ierr)
    call VecDestroy(PetVecc3,ierr)
    call VecDestroy(Fsi,ierr)
    call VecDestroy(Fsg,ierr)
    call VecDestroy(PPb,ierr)
    Call MatDestroy(RsMat,ierr)
    Call MatDestroy(DsMat,ierr)
    Call MatDestroy(BcMat,ierr)
    Call MatDestroy(RrMat,ierr)
    Call MatDestroy(RcMat,ierr)
    Call MatDestroy(Asii,ierr)
    Call MatDestroy(Asgi,ierr)
    Call MatDestroy(Asgg,ierr)
    Call MatDestroy(Ascc,ierr)
    Call MatDestroy(Asrr,ierr)
    Call MatDestroy(Asci,ierr)
    Call MatDestroy(Ascr,ierr)
    Call MatDestroy(Asri,ierr)
    Call MatDestroy(Asmat,ierr)
    Call KspDestroy(kspAi,ierr)
    Call KspDestroy(kspAm,ierr)
    Call KspDestroy(kspAc,ierr)

!!-----------------------------------------------------------------------------------------------
!!! PETSc-Finalize
call PetscFinalize(ierr)



END SUBROUTINE PETSc_stonncpcgm


! !!!!---------------------------------------------------------------------------------------------
! !!!!!!!!!!!!!!! Coarse Problem: Lumped-PCGM Solver : Algorithm 7 from the report !!!!!!!!!!!!!!
! !!!!---------------------------------------------------------------------------------------------
 SUBROUTINE solve_coarsepcgm(pid,nirp,ncp,nrp,nip,nbgcp,BcMat,kspAm,kspAi,kspAc,pcAm,pcAi,pcAc, &
                        tempir, tempr, tempcg, PetVecc2, PetVecc3, PetVeccg, ybc, Vri,DcBc, &
                        PetVeci,SolVeci,PetVecir,SolVecir,PetVecc,Solvecc,PetVecr,SolVecr, &
                        Asci,Asri,Ascr,Ascc,Dc,Zc,tolc,stageC)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

!!----------------------------------------------------------------------------------------
!! PETSc-Variables Declaration
    Vec              PetVeci,PetVecir,PetVecc,PetVecr
    Vec              SolVeci,SolVecir,SolVecc,SolVecr
    Vec              PetVecc2,PetVecc3,PetVeccg

    Mat              Asri,Asci,Ascr,Ascc,BcMat

    KSP              kspAm, kspAi, kspAc
    PC               pcAm, pcAi, pcAc

    PetscScalar      NegOne

    PetscLogStage    stageC

!!----------------------------------------------------------------------------------------
!! Fortran-Variables Declaration
    integer :: maxiter = 100
    integer :: pid,ierr,j
    integer :: nirp,ncp,nrp,nip,nbgcp

    integer, dimension(nirp)  :: tempir
    integer, dimension(nrp)   :: tempr
    integer, dimension(nbgcp) :: tempcg

    double precision, dimension(nbgcp) :: Pb_gc,Zb_gc,Qb_gc
    double precision, dimension(nbgcp) :: Dc,Zc,DcBc
    double precision, dimension(nrp)   :: Vri
    double precision, dimension(nirp)  :: ybc

    double precision :: rho_nextc, rho_currc, alphac, betac, errc, tolc

!!----------------------------------------------------------------------------------------
!! Initiate Array's
    rho_nextc = 0.0d0
    rho_currc = 0.0d0
    NegOne    = -1.0


!!!---------------------------------------------------------------------------------------
!! Pre-Iteration : 0th Iteration :step 2 algorithm 7
!   Algorithm 9 : L-PP
!!!---------------------------------------------------------------------------------------
!!!---------------- PPE: zi = inv(M)*ri -----------------------------------------------!!!
    !!STEP 3 : Scatter
    call SetVecSeq(nbgcp,Dc,PetVeccg,ierr)
    call MatMult(BcMat,PetVeccg,PetVecc,ierr)


    !STEP 4 : Solve      Acc\rc , where rc = Bc * dc
    call KSPSolve(kspAc,PetVecc,SolVecc,ierr)

    !!STEP 5 : Gather
    call MatMultTranspose(BcMat,SolVecc,PetVeccg,ierr)
    call VecGetValues(PetVeccg, nbgcp, tempcg, DcBc, ierr)

    call PetscLogStagePush(stageC,ierr)
    CALL MPI_REDUCE(DcBc,Zb_gc,nbgcp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call PetscLogStagePop(ierr)

!!!                 Gathered the initially preconditioned residual to all subdomains.

    !STEP 7 & 8 : Compute. !!!!! rho_nextc is actually delta in Ajits thesis - Algorithm 7 step 4
    if (pid .eq. 0) Pb_gc = Zb_gc
    if (pid .eq. 0) rho_nextc = dot_product(Dc,Zb_gc)

!!!---------------------------------------------------------------------------------------
!! PCGM-Iteration : Main For Loop for each iteration: Step 6 to 15 - Ajit's thesis Algorithm 7
!!!---------------------------------------------------------------------------------------
    Zc(:) = 0.0d0
    do j = 1,maxiter

        !!STEP 11 : Scatter.

        !!!! Eq. 3.64 : Parallel Schur Mat-Vec Product   Algorithm 8 Ajits Thesis
        call PetscLogStagePush(stageC,ierr)
        CALL MPI_BCAST(Pb_gc,nbgcp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)

        call SetVecSeq(nbgcp,Pb_gc,PetVeccg,ierr)
        call MatMult(BcMat,PetVeccg,PetVecc,ierr)

!!!---------------- PMV:Sec1: v1 = Src*p1 = [Arc - (Ari *inv(Aii)*Aic)]*p1 -----------------!!!
        !! STEP 12 : Compute
        call MatMultTranspose(Asci,PetVecc,PetVeci,ierr)
        call KSPSolve(kspAi,PetVeci,SolVeci,ierr)
        call MatMult(Asri,SolVeci,SolVecr,ierr)

        call MatMultTranspose(Ascr,PetVecc,PetVecr,ierr)
        call VecAXPY(PetVecr,NegOne,SolVecr,ierr);
        call VecGetValues(PetVecr, nrp, tempr, Vri, ierr)

! !!!---------------- PMV:Sec2: v2 = inv(Srr)*v1 ---------------------------------------------!!!
        !!STEP 13 : Solve
        ybc(:) = 0.0d0
        ybc((nip+1):nirp)=Vri

        call SetVecSeq(nirp,ybc,PetVecir,ierr)
        call KSPSolve(kspAm,PetVecir,SolVecir,ierr)
        call VecGetValues(SolVecir, nirp, tempir, ybc, ierr)

        Vri = ybc((nip+1):nirp)

!!!---------------- PMV:Sec3: v3 = Scr*v2 = [Acr - (Aci *inv(Aii)*Air)]*v2 ---------------!!!
        !!STEP 14 : Compute
        call SetVecSeq(nrp,Vri,PetVecr,ierr)
        call MatMultTranspose(Asri,PetVecr,PetVeci,ierr)
        call KSPSolve(kspAi,PetVeci,SolVeci,ierr)
        call MatMult(Asci,SolVeci,SolVecc,ierr)
        call MatMult(Ascr,PetVecr,PetVecc2,ierr)
        call VecAXPY(PetVecc2,NegOne,SolVecc,ierr)

!!!!---------------- PMV:Sec4: v4 = Scc*p1 = [Acc - (Aci *inv(Aii)*Aic)]*v3 ----------------!!!
        !!STEP 15 : Compute
        call MatMultTranspose(Asci,PetVecc,PetVeci,ierr)

        call KSPSolve(kspAi,PetVeci,SolVeci,ierr)

        call MatMult(Asci,SolVeci,SolVecc,ierr)

        ! print*, 'process id', pid
        ! call VecView(PetVecc, PETSC_VIEWER_STDOUT_SELF, ierr)

        ! stop 123
        call MatMult(Ascc,PetVecc,PetVecc3,ierr)

        call VecAXPY(PetVecc3,NegOne,SolVecc,ierr);

        ! stop 123
! !!----------------------------------------------------------------------------------------
!         !!STEP 16 : Update.   step 8 in algorithm 8.  q1 = v4 - v3
        call VecAYPX(PetVecc2,NegOne,PetVecc3,ierr)
!         !!STEP 17 : Gather
        call MatMultTranspose(BcMat,PetVecc2,PetVeccg,ierr)
        call VecGetValues(PetVeccg, nbgcp, tempcg, DcBc, ierr)

        call PetscLogStagePush(stageC,ierr)
        CALL MPI_REDUCE(DcBc,Qb_gc,nbgcp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)
!!!                     Found out Qb_gc - Matrix Vector Output for PCGM Coarse Level
        !!STEP 19, 20, 21 & 22 : Compute & Update
        if (pid .eq. 0) then
            rho_currc = rho_nextc
            alphac = rho_currc/dot_product(Qb_gc,Pb_gc)
            errc = alphac*alphac*dot_product(Pb_gc,Pb_gc)/dot_product(Zc,Zc)
            Zc = Zc + alphac*Pb_gc  !!!!! Its the solution vector
            Dc = Dc - alphac*Qb_gc
            ! print*, 'iteration # ',j,'relative error of', errc
        end if

        !!STEP 23 : Exit
        call PetscLogStagePush(stageC,ierr)
        CALL MPI_BCAST(errc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)

        if (errc .lt. tolc) exit

! !!----------------------------------------------------------------------------------------
! !!!---------------- PPE: zi = inv(M)*ri -----------------------------------------------!!!
        !!STEP 25 : Scatter
        call PetscLogStagePush(stageC,ierr)
        CALL MPI_BCAST(Dc,nbgcp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)

        call SetVecSeq(nbgcp,Dc,PetVeccg,ierr)
        call MatMult(BcMat,PetVeccg,PetVecc,ierr)

        !!STEP 26 : Solve
        call KSPSolve(kspAc,PetVecc,SolVecc,ierr)

        !!STEP 27 : Gather
        call MatMultTranspose(BcMat,SolVecc,PetVeccg,ierr)
        call VecGetValues(PetVeccg, nbgcp, tempcg, DcBc, ierr)

        call PetscLogStagePush(stageC,ierr)
        CALL MPI_REDUCE(DcBc,Zb_gc,nbgcp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call PetscLogStagePop(ierr)

        !!STEP 29, 30 to 31 : Compute & Update
        if (pid .eq. 0) then
            rho_nextc = dot_product(Dc,Zb_gc)
            betac = rho_nextc/rho_currc
            Pb_gc = Zb_gc + betac*Pb_gc
        end if
     end do


 END SUBROUTINE solve_coarsepcgm


END MODULE PETScSolvers
