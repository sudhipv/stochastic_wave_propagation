!! This MODULE is a part of HPSFEM package : Ajit Desai : Dec 2015
!! This is the common MODULE contains various subrountines for PETSc
!! This is mainly written to make use of PETSc for solving FEM or DDM problems
!! If you want to include any other PETSc subroutines, then this is the place
!! Contains:
!! 1. PETScKSP     : Construct PETSc KSP setup and solve using iterative solver
!! 2. PETScMUMPS   : Construct PETSc MUMPS setup and slover using direct solver
!! 3. SetPETScKSP  : Only to construct PETSc KSP set-up, no solve
!! 4. GetVecSeq    : Construct PETSc Vector : Sequential  : Serial
!! 5. GetMatSeq    : Construct PETSc Matrix : Sequential  : Serial
!! 6. GetVecSeqDummy:To setup dummy PETSc-Vec with temp   : Serial
!! 7. GetVecSeqTemp1:To setup 1-dymmy PETSc-Vec (no temp) : Serial
!! 8. GetVecSeqTemp2:To setup 2-dymmy PETSc-Vec with temp : Serial
!! 9. StoVecSeqOneLevel: To call Sparse Stochastic Vector Assembly
!!10. StoMatSeqOneLevel: To call Sparse Stochastic Matrix Assemly
!!11. GetRs,GetBc  : To construct Rs/Bc restriction operator : Remove in future
!!12. GetDs        : To construct bolock-diagonal weighted matrix
!!13. getRc,getRr  : To consturct RcS/RrS restriction operators
!!14. GetVecMPI    : Construct PETSc Vector : MPI format  : Parallel
!!15. PETScSolver  : Call PETSc assembled Mat, Vec and solver Ax=b using ksp
!!
!!-----------------------------------------------------------------------------------

MODULE PETScommon

 use PETScAssembly
 use assembly
 use common


IMPLICIT NONE


CONTAINS


!!!*********************************************************

!!!*********************************************************
!!! Subroutine: To Call PETSc KSP solver in general way

!!!!Error: Cannot change attributes of USE-associated symbol petscksp at (1)

!!! Name being changed to PETSc_KSP for PETSc - 3.9.2
!!! Subroutine: To Call PETSc KSP solver in general way
SUBROUTINE PETSc_KSP(ksp,pc,PetMat,PetVec,Solvec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

    Vec              PetVec,SolVec
    Mat              PetMat
    KSP              ksp
    PC               pc

    PetscInt         maxits
    PetscReal        tol
    PetscErrorCode   ierr


    call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
    call KSPSetOperators(ksp,PetMat,PetMat,ierr)
    call KSPSetFromOptions(ksp,ierr)   !! position changed : moved down

    call KSPSetType(ksp,KSPCG,ierr)
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCBJACOBI,ierr)
    !call PCSetType(pc,PCGAMG,ierr)
    !!call PCSetType(pc,PCKSP,ierr)

    tol = .0000001
    maxits = 200
    call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,maxits,ierr)

    call KSPSolve(ksp,PetVec,SolVec,ierr)

END SUBROUTINE PETSc_KSP



!!!!*********************************************************
!!!! Subroutine : To call PETSc MUMPS solver in general way
!SUBROUTINE PETScMUMPS(ksp,pc,PetMat,PetVec,Solvec,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscpc.h>
!
!    Vec              PetVec,SolVec
!    Mat              PetMat
!    KSP              ksp
!    PC               pc
!
!    Mat              F
!    PetscInt         ival,icntl
!
!    PetscErrorCode   ierr
!
!    call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
!    call KSPSetOperators(ksp,PetMat,PetMat,ierr)
!
!    call KSPSetType(ksp,KSPPREONLY,ierr)
!    call KSPGetPC(ksp,pc,ierr)
!    !call PCSetType(pc,PCLU,ierr)           !! LU Factorization
!    call PCSetType(pc,PCCHOLESKY,ierr)    !! Cholesky Factorization
!    call PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS,ierr)
!    call PCFactorSetUpMatSolverPackage(pc,ierr)
!    call PCFactorGetMatrix(pc,F,ierr)
!
!    !! sequential ordering
!    icntl = 7
!    ival  = 2
!    call MatMumpsSetIcntl(F,icntl,ival,ierr)
!
!    call KSPSetFromOptions(ksp,ierr)
!    call KSPGetPC(ksp,pc,ierr)
!
!    call KSPSolve(ksp,PetVec,SolVec,ierr)
!
!
!END SUBROUTINE PETScMUMPS

!!!*********************************************************
!!! Subroutine: Just To Setup PETSc-KSP solver in general way
SUBROUTINE SetPETScKSP(ksp,pc,PetMat,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

    Mat              PetMat
    KSP              ksp
    PC               pc

    PetscErrorCode   ierr

    call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
    call KSPSetOperators(ksp,PetMat,PetMat,ierr)
    call KSPSetFromOptions(ksp,ierr)   !! position changed : moved down


    call KSPSetType(ksp,KSPCG,ierr)
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCBJACOBI,ierr)
    !call PCSetType(pc,PCGAMG,ierr)
    !! call PCSetType(pc,PCKSP,ierr)

    !! Just to set KSP: no need to call KSPSolve
    !! call KSPSolve(ksp,PetVec,SolVec,ierr)

END SUBROUTINE SetPETScKSP



!!!!!!!!!!!! Sudhi P V !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE GetMatSeqTwoLevel(pid,PAmat,PAii,PAgg,PAgi,PAri,PArr,PAcc,PAci,PAcr,PMii,PMgg,PMgi,PMri,PMrr,PMcc,PMci,PMcr,&
                             PCii,PCgg,PCgi,PCri,PCrr,PCcc,PCci,PCcr,np,nb,nc,nr,nbgc,nVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    Mat              PAmat,PAii,PAgg,PArr,PAcc
    Mat              PAgi,PAri,PAci,PAcr

    Mat              PMii,PMgg,PMrr,PMcc
    Mat              PMgi,PMri,PMci,PMcr

    Mat              PCii,PCgg,PCrr,PCcc
    Mat              PCgi,PCri,PCci,PCcr

    PetscInt         nb, nir, nc, nr, nbV,niV,nrV,nirV,ncV
    PetscInt         ni,m,nz,one
    PetscInt         ii, jj
    PetscInt         id, jd
    PetscScalar      Amat, Mmat, Cmat
    PetscErrorCode   ierr

    integer :: np, nbgc, nVec
    double precision :: const_diff

    integer              :: pid
    character(len=255)   :: str1, str2, str3

    double precision, dimension(2,2) :: dbounds

    double precision, allocatable, dimension(:,:)  :: Adii,Adgg,Adrr,Adcc
    double precision, allocatable, dimension(:,:)  :: Adgi,Adri,Adci,Adcr

    double precision, allocatable, dimension(:,:)  :: Mdii,Mdgg,Mdrr,Mdcc
    double precision, allocatable, dimension(:,:)  :: Mdgi,Mdri,Mdci,Mdcr

    double precision, allocatable, dimension(:,:)  :: Cdii,Cdgg,Cdrr,Cdcc
    double precision, allocatable, dimension(:,:)  :: Cdgi,Cdri,Cdci,Cdcr

    REAL*8, PARAMETER :: eps = 1d-16

    ni = np-nb
    nir = np-nc !! Interior plus remaining nodes

    niV = ni *nVec
    nirV = nir * nVec
    nbV = nb * nVec
    ncV = nc *nVec
    nrV = nr * nVec


    !!-----------------------------------------------------------------------------------------------
    allocate(Adii(niV,niV), Adgg(nbV,nbV), Adrr(nrV,nrV), Adcc(ncV,ncV),&
    Adgi(nbV,niV), Adri(nrV,niV), Adci(ncV,niV), Adcr(ncV,nrV))

    allocate(Mdii(niV,niV), Mdgg(nbV,nbV), Mdrr(nrV,nrV), Mdcc(ncV,ncV),&
    Mdgi(nbV,niV), Mdri(nrV,niV), Mdci(ncV,niV), Mdcr(ncV,nrV))

    allocate(Cdii(niV,niV), Cdgg(nbV,nbV), Cdrr(nrV,nrV), Cdcc(ncV,ncV),&
    Cdgi(nbV,niV), Cdri(nrV,niV), Cdci(ncV,niV), Cdcr(ncV,nrV))


!!!! will be error from PETSc if memory allocated for numbe of non zeros is less
    nz = niV*3  !! number of non-zeros per row for matrix pre-allocation

    call MatCreateSeqAIJ(PETSC_COMM_SELF, niV, niV, nz, PETSC_NULL_INTEGER, PAii, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbV, nbV, nz, PETSC_NULL_INTEGER, PAgg, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbV, niV, nz, PETSC_NULL_INTEGER, PAgi, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nirV, nirV, nz, PETSC_NULL_INTEGER, PAmat, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrV, niV, nz, PETSC_NULL_INTEGER, PAri, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, niV, nz, PETSC_NULL_INTEGER, PAci, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrV, nrV, nz, PETSC_NULL_INTEGER, PArr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, ncV, nz, PETSC_NULL_INTEGER, PAcc, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, nrV, nz, PETSC_NULL_INTEGER, PAcr, ierr)


    call MatCreateSeqAIJ(PETSC_COMM_SELF, niV, niV, nz, PETSC_NULL_INTEGER, PMii, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbV, nbV, nz, PETSC_NULL_INTEGER, PMgg, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbV, niV, nz, PETSC_NULL_INTEGER, PMgi, ierr)
    ! call MatCreateSeqAIJ(PETSC_COMM_SELF, nir, nir, nz, PETSC_NULL_INTEGER, PMmat, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrV, niV, nz, PETSC_NULL_INTEGER, PMri, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, niV, nz, PETSC_NULL_INTEGER, PMci, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrV, nrV, nz, PETSC_NULL_INTEGER, PMrr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, ncV, nz, PETSC_NULL_INTEGER, PMcc, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, nrV, nz, PETSC_NULL_INTEGER, PMcr, ierr)


    call MatCreateSeqAIJ(PETSC_COMM_SELF, niV, niV, nz, PETSC_NULL_INTEGER, PCii, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbV, nbV, nz, PETSC_NULL_INTEGER, PCgg, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbV, niV, nz, PETSC_NULL_INTEGER, PCgi, ierr)
    ! call MatCreateSeqAIJ(PETSC_COMM_SELF, nir, nir, nz, PETSC_NULL_INTEGER, PCmat, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrV, niV, nz, PETSC_NULL_INTEGER, PCri, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, niV, nz, PETSC_NULL_INTEGER, PCci, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrV, nrV, nz, PETSC_NULL_INTEGER, PCrr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, ncV, nz, PETSC_NULL_INTEGER, PCcc, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncV, nrV, nz, PETSC_NULL_INTEGER, PCcr, ierr)


!!!!!!!!!! Loading the FEniCS Assembled Matrices !!!!!!!!!!!!!!!!!!!!!!!


!!!!! Loaded A matrix is actually K_T - KEPT SAME NAME FOR EASINESS

    call int2str(str1,pid+1,1)
    ! call int2str(str2,1,1)
    ! print*,'str1 is', str1

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1)// '/ADii' // trim(str1)// '.dat'
    ! print*,'str3 is', str3
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adii
    close(1)

    !print*, str3, SHAPE(Adii)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgg' // trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adgg
    close(2)

    !print*, str3, SHAPE(Adgg)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADig' // trim(str1) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adgi
    close(3)

    !print*, str3, SHAPE(Adgi)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADic' // trim(str1) // '.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adci
    close(1)

    !print*, str3, SHAPE(Adci)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADir' // trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adri
    close(2)

    !print*, str3, SHAPE(Adri)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADcc' // trim(str1) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adcc
    close(3)

    !print*, str3, SHAPE(Adcc)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADrr' // trim(str1) // '.dat'
    open(unit=4,file=str3,status='old')
    read(unit=4,fmt=*) Adrr
    close(4)

    !print*, str3, SHAPE(Adrr)


    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADrc' // trim(str1) // '.dat'
    open(unit=5,file=str3,status='old')
    read(unit=5,fmt=*) Adcr
    close(5)

    !print*, str3, SHAPE(Adcr)


!!!!!!!!!!! Mass Matrix !!!!!!!!!!!!!!!!


    call int2str(str1,pid+1,1)
    ! call int2str(str2,1,1)
    ! print*,'str1 is', str1

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1)// '/MDii' // trim(str1)// '.dat'
    ! print*,'str3 is', str3
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Mdii
    close(1)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/MDgg' // trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Mdgg
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/MDig' // trim(str1) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Mdgi
    close(3)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/MDic' // trim(str1) // '.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Mdci
    close(1)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/MDir' // trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Mdri
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/MDcc' // trim(str1) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Mdcc
    close(3)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/MDrr' // trim(str1) // '.dat'
    open(unit=4,file=str3,status='old')
    read(unit=4,fmt=*) Mdrr
    close(4)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/MDrc' // trim(str1) // '.dat'
    open(unit=5,file=str3,status='old')
    read(unit=5,fmt=*) Mdcr
    close(5)

    ! print*, str3, SHAPE(Mdcr)



!!!!!!!!! Damping Matrix !!!!!!!!!!!!!!!


    call int2str(str1,pid+1,1)
    ! call int2str(str2,1,1)
    ! print*,'str1 is', str1

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1)// '/CDii' // trim(str1)// '.dat'
    ! print*,'str3 is', str3
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Cdii
    close(1)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/CDgg' // trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Cdgg
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/CDig' // trim(str1) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Cdgi
    close(3)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/CDic' // trim(str1) // '.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Cdci
    close(1)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/CDir' // trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Cdri
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/CDcc' // trim(str1) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Cdcc
    close(3)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/CDrr' // trim(str1) // '.dat'
    open(unit=4,file=str3,status='old')
    read(unit=4,fmt=*) Cdrr
    close(4)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/CDrc' // trim(str1) // '.dat'
    open(unit=5,file=str3,status='old')
    read(unit=5,fmt=*) Cdcr
    close(5)

    ! print*, str3, SHAPE(Cdcr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    one    = 1

    !print*, ni

    !print*,'Before assembly'

    ! call MatView(Adii, PETSC_VIEWER_STDOUT_SELF)

!!!!!! PAii && PAmat

!!!!! Fortran and Fenics uses different type of ordering for matrices, row order and column order, so all the fenics matrices taken here will get inverted.
!!!!!! Here Because of the way fenics assembled the matrices ii, cc, rr, cc are also not symmetric,,so have to invert the way they are taken
!!!!!! assemble it correctly. For other matrices like cr, ig, ic etc...inversion has been done while loading itself.

    do id = 1,niV
        do jd = 1,niV

            if ( ABS(Adii(jd,id)) > eps ) then

                Amat = Adii(jd,id)

                call MatSetValues(PAii,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)
                call MatSetValues(PAmat,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)

            end if

            if (ABS(Mdii(jd,id)) > eps ) then

                Mmat = Mdii(jd,id)
                call MatSetValues(PMii,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)

            end if

            if (ABS(Cdii(jd,id)) > eps ) then

                Cmat = Cdii(jd,id)
                call MatSetValues(PCii,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)

            end if


        end do
    end do


    ! if (pid .eq. 0) then
    !     print*, 'procees is', pid
    !     print*, 'Mdii  is', Mdii(18,17)
    !     ! call MatView(PMii, PETSC_VIEWER_STDOUT_SELF, ierr)
    !     stop 123
    ! end if

!!!!!!!!!!!!!!!!!!!! DOING three SEPERATE LOOPS TO ENSURE !!!!!!!!!!!!!!!

    ! do id = 1,ni
    !     do jd = 1,ni
    !         if (Mdii(id,jd) .ne. 0) then

    !             Mmat = Mdii(id,jd)

    !             call MatSetValues(PMii,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)
    !             call MatSetValues(PMmat,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)

    !         end if
    !     end do
    ! end do



    ! do id = 1,ni
    !     do jd = 1,ni
    !         if (Cdii(id,jd) .ne. 0) then

    !             Cmat = Cdii(id,jd)

    !             call MatSetValues(PCii,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)
    !             call MatSetValues(PCmat,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)

    !         end if
    !     end do
    ! end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!! PAgg

    do id = 1,nbV
        do jd = 1,nbV

            if (ABS(Adgg(jd,id)) > eps ) then

                Amat = Adgg(jd,id)
                call MatSetValues(PAgg,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)

            end if

            if (ABS(Mdgg(jd,id)) > eps) then

                Mmat = Mdgg(jd,id)
                call MatSetValues(PMgg,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)

            end if

            if (ABS(Cdgg(jd,id)) > eps) then

                Cmat = Cdgg(jd,id)
                call MatSetValues(PCgg,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)

            end if

        end do
    end do

!!!!!! PAgi

    do id = 1,nbV
        do jd = 1,niV

            if (ABS(Adgi(id,jd)) > eps ) then

                Amat = Adgi(id,jd)
                call MatSetValues(PAgi,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)

            end if

            if (ABS(Mdgi(id,jd)) > eps ) then

                Mmat = Mdgi(id,jd)
                call MatSetValues(PMgi,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)

            end if

            if (ABS(Cdgi(id,jd)) > eps ) then

                Cmat = Cdgi(id,jd)
                call MatSetValues(PCgi,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)

            end if

        end do
    end do

!!!!!!! PAri

    ! print*, 'assembling PAmat for ri', nir
    ! print*, 'ni plus nr is', ni+nr
    !!!! Have to confirm whether its done correctly or not ...potential error point

    do id = 1,nrV
        do jd = 1,niV

            if (ABS(Adri(id,jd)) > eps) then

                Amat = Adri(id,jd)
                call MatSetValues(PAri,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)
                call MatSetValues(PAmat,one,niV+id-1,one,jd-1,Amat,INSERT_VALUES,ierr)
                call MatSetValues(PAmat,one,jd-1,one,niV+id-1,Amat,INSERT_VALUES,ierr)

            end if

            if (ABS(Mdri(id,jd)) > eps) then

                Mmat = Mdri(id,jd)
                call MatSetValues(PMri,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)

            end if

            if (ABS(Cdri(id,jd)) > eps) then

                Cmat = Cdri(id,jd)
                call MatSetValues(PCri,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)

            end if

        end do
    end do

!!!!! PAci

    do id = 1,ncV
        do jd = 1,niV
            if (ABS(Adci(id,jd)) > eps) then
                Amat = Adci(id,jd)
                call MatSetValues(PAci,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)
            end if

            if (ABS(Mdci(id,jd)) > eps) then
                Mmat = Mdci(id,jd)
                call MatSetValues(PMci,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)
            end if

            if (ABS(Cdci(id,jd)) > eps) then
                Cmat = Cdci(id,jd)
                call MatSetValues(PCci,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)
            end if

        end do
    end do

!!!!! PArr
    ! print*, 'assembling PAmat for rr', nir
    ! print*, 'ni plus nr is', ni+nr

    do id = 1,nrV
        do jd = 1,nrV

            if (ABS(Adrr(jd,id)) > eps) then
                Amat = Adrr(jd,id)
                call MatSetValues(PArr,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)
                call MatSetValues(PAmat,one,niV+id-1,one,niV+jd-1,Amat,INSERT_VALUES,ierr)
            end if

            if (ABS(Mdrr(jd,id)) > eps) then
                Mmat = Mdrr(jd,id)
                call MatSetValues(PMrr,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)
            end if

            if (ABS(Cdrr(jd,id)) > eps) then
                Cmat = Cdrr(jd,id)
                call MatSetValues(PCrr,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)
            end if

        end do
    end do

!!!!! PAcc

    do id = 1,ncV
        do jd = 1,ncV

            if (ABS(Adcc(jd,id)) > eps) then
                Amat = Adcc(jd,id)
                call MatSetValues(PAcc,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)
            end if

            if (ABS(Mdcc(jd,id))> eps) then
                Mmat = Mdcc(jd,id)
                call MatSetValues(PMcc,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)
            end if

            if (ABS(Cdcc(jd,id)) > eps) then
                Cmat = Cdcc(jd,id)
                call MatSetValues(PCcc,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)
            end if

        end do
    end do

!!!!! PAcr

    do id = 1,ncV
        do jd = 1,nrV

            if (ABS(Adcr(id,jd)) > eps) then
                Amat = Adcr(id,jd)
                call MatSetValues(PAcr,one,id-1,one,jd-1,Amat,INSERT_VALUES,ierr)
            end if

            if (ABS(Mdcr(id,jd)) > eps) then
                Mmat = Mdcr(id,jd)
                call MatSetValues(PMcr,one,id-1,one,jd-1,Mmat,INSERT_VALUES,ierr)
            end if

            if (ABS(Cdcr(id,jd)) > eps) then
                Cmat = Cdcr(id,jd)
                call MatSetValues(PCcr,one,id-1,one,jd-1,Cmat,INSERT_VALUES,ierr)
            end if

        end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    deallocate(Adii,Adgg,Adrr,Adcc,Adgi,Adri,Adci,Adcr)
    deallocate(Mdii,Mdgg,Mdrr,Mdcc,Mdgi,Mdri,Mdci,Mdcr)
    deallocate(Cdii,Cdgg,Cdrr,Cdcc,Cdgi,Cdri,Cdci,Cdcr)


    call MatAssemblyBegin(PAii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PAgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PAgi,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PAmat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PAri,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PAci,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PArr,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PAcc,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PAcr,MAT_FINAL_ASSEMBLY,ierr)


    call MatAssemblyEnd(PAii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PAgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PAgi,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PAmat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PAri,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PAci,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PArr,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PAcc,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PAcr,MAT_FINAL_ASSEMBLY,ierr)

!!!!!!!!!!!!!!! Mass matrices !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call MatAssemblyBegin(PMii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PMgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PMgi,MAT_FINAL_ASSEMBLY,ierr)
    ! call MatAssemblyBegin(PMmat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PMri,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PMci,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PMrr,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PMcc,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PMcr,MAT_FINAL_ASSEMBLY,ierr)


    call MatAssemblyEnd(PMii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PMgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PMgi,MAT_FINAL_ASSEMBLY,ierr)
    ! call MatAssemblyEnd(PMmat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PMri,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PMci,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PMrr,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PMcc,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PMcr,MAT_FINAL_ASSEMBLY,ierr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call MatAssemblyBegin(PCii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PCgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PCgi,MAT_FINAL_ASSEMBLY,ierr)
    ! call MatAssemblyBegin(PCmat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PCri,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PCci,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PCrr,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PCcc,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PCcr,MAT_FINAL_ASSEMBLY,ierr)


    call MatAssemblyEnd(PCii,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PCgg,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PCgi,MAT_FINAL_ASSEMBLY,ierr)
    ! call MatAssemblyEnd(PCmat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PCri,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PCci,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PCrr,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PCcc,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PCcr,MAT_FINAL_ASSEMBLY,ierr)


!!!! PMii, PAgg, PAgi        verified

    ! print*, 'process id', pid
    ! call MatView(PAgi, PETSC_VIEWER_STDOUT_SELF, ierr)
    ! print*, 'end matrix'
    ! stop 123



    ! if (pid .eq. 0) then
    !     print*, 'procees is', pid
    !     print*, 'PMii  is'
    !     call MatView(PMii, PETSC_VIEWER_STDOUT_SELF, ierr)
    !     stop 123
    ! end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! print*,'Finished PETSc Mat assembly'

    ! call MatView(PAii, PETSC_VIEWER_STDOUT_SELF)
    if (pid .eq. 0) print*, '....................................................................'
    if (pid .eq. 0) print*, '...............Successfully assembled PETSc matrices................'
    if (pid .eq. 0) print*, '....................................................................'

END SUBROUTINE GetMatSeqTwoLevel




SUBROUTINE GetInitialCond(pid, ni,np,u0i,u0b,v0i,v0b,a0i,a0b,nVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    ! Vec              Pu0, Pv0, Pa0
    PetscInt         ni, nb, np,npV
    PetscScalar      uvec,vvec,avec
    PetscErrorCode   ierr

    ! double precision, allocatable, dimension(:) :: u0, v0, a0

    double precision, allocatable, dimension(:) :: u0i, u0b
    double precision, allocatable, dimension(:) :: v0i, v0b
    double precision, allocatable, dimension(:) :: a0i, a0b

    ! double precision, dimension(:,:)            :: U_wave

    character(len=255)      :: str1, str3
    integer :: pid, nid, id, one, ii, nVec,j


    ni= np-nb
    npV = np * nVec
    !!-----------------------------------------------------------------------------------------------
    !! Allocate memory to read pre-assembled vectors
    ! allocate(u0i(ni*nVec), u0b(nb*nVec), v0i(ni*nVec), v0b(nb*nVec), a0i(ni*nVec), a0b(nb*nVec))

    ! call VecCreateSeq(PETSC_COMM_SELF, npV, Pu0, ierr)
    ! call VecCreateSeq(PETSC_COMM_SELF, npV, Pv0, ierr)
    ! call VecCreateSeq(PETSC_COMM_SELF, npV, Pa0, ierr)

    !!-----------------------------------------------------------------------------------------------
    !! FEniCS based Vec-Assembly procedure (using preassembled vecs)
    one= 1

    ! bi = 0.0d0
    ! bg = 0.0d0

    ! print*, 'Assembling Force Vectors'
    ! print*, 'process id',pid
    ! print*, 'ni is', ni
    ! print*, 'nb is', nb
    call int2str(str1,pid+1,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/u0i_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) u0i
    close(2)


    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/v0i_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) v0i
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/a0i_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) a0i
    close(2)


    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/u0b_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) u0b
    close(2)


    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/v0b_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) v0b
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/a0b_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) a0b
    close(2)

    !---------------------------------------->fsi



    ! print*,'Shape of bi',SHAPE(bi)




    ! print*,'Shape of bg',SHAPE(bg)
    ! call VecView(PFg,PETSC_VIEWER_STDOUT_SELF);



    ! print*, 'Assembling Initial Conditions'
    ! print*, 'process id',pid
    ! print*, 'np is',np
    ! print*,'Shape of u0',SHAPE(u0)


    ! do j = 1,nVec
    !     do id = 1,ni
    !             ii = (j-1)*ni+id-1
    !             uvec = u0((j-1)*np+id)
    !             vvec = v0((j-1)*np+id)
    !             avec = a0((j-1)*np+id)
    !             call VecSetValues(Pu0,one,ii,uvec,INSERT_VALUES,ierr)
    !             call VecSetValues(Pv0,one,ii,vvec,INSERT_VALUES,ierr)
    !             call VecSetValues(Pa0,one,ii,avec,INSERT_VALUES,ierr)
    !     end do
    ! end do

    ! if (pid .eq. 1) then
    !     print*,'process ID is', pid
    !     call VecView(Pu0,PETSC_VIEWER_STDOUT_SELF,ierr)
    ! end if


    ! do j = 1,nVec
    !     do id = 1,nb
    !             ii = nVec*ni+id-1
    !             uvec = u0((j-1)*ni+id)
    !             vvec = v0((j-1)*ni+id)
    !             avec = a0((j-1)*ni+id)
    !             call VecSetValues(Pu0,one,ii,uvec,INSERT_VALUES,ierr)
    !             call VecSetValues(Pv0,one,ii,vvec,INSERT_VALUES,ierr)
    !             call VecSetValues(Pa0,one,ii,avec,INSERT_VALUES,ierr)
    !     end do
    ! end do


    ! if (pid .eq. 1) then
    !     print*,'process ID is', pid
    !     call VecView(Pu0,PETSC_VIEWER_STDOUT_SELF,ierr)
    ! end if


    ! do id = 1,npV
    !         ii = id-1
    !         ivec = v0(id)
    !         call VecSetValues(Pv0,one,ii,ivec,INSERT_VALUES,ierr)
    ! end do

    ! do id = 1,npV
    !         ii = id-1
    !         ivec = a0(id)
    !         call VecSetValues(Pa0,one,ii,ivec,INSERT_VALUES,ierr)
    ! end do


    ! call VecAssemblyBegin(Pu0,ierr)
    ! call VecAssemblyBegin(Pv0,ierr)
    ! call VecAssemblyBegin(Pa0,ierr)

    ! call VecAssemblyEnd(Pu0,ierr)
    ! call VecAssemblyEnd(Pv0,ierr)
    ! call VecAssemblyEnd(Pa0,ierr)

    ! if (pid .eq. 1) then
    !     print*,'process ID is', pid
    !     call VecView(Pu0,PETSC_VIEWER_STDOUT_SELF,ierr)
    ! end if

    ! print*,"process id is", pid
    ! print*,"u0i is", u0i

    ! print*,"process id is", pid
    ! print*,"u0b is", u0b

! DEALLOCATE(u0i,v0i,a0i,u0b,v0b,a0b)

    if (pid .eq. 1) print*, '....................................................................'
    if (pid .eq. 1) print*, '...............Successfully assembled PETSc vectors................'
    if (pid .eq. 1) print*, '....................................................................'


END SUBROUTINE GetInitialCond


SUBROUTINE GetInitialCondSubdomain(pid,ni,nb,u0i,u0b,v0i,v0b,a0i,a0b,nVec, ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
    !!!!! Better implementation might be possible. Used PetscScalar


    ! Vec              Pu0, Pv0, Pa0
    PetscErrorCode   ierr

    integer :: ni, nb, i, nVec,pid
    character(len=255)      :: str1, str3

    ! integer, dimension((ni*nVec)) :: tempi
    ! integer, allocatable, dimension(:) :: tb

    double precision, dimension(:) :: u0i, u0b
    double precision, dimension(:) :: v0i, v0b
    double precision, dimension(:) :: a0i, a0b


    call int2str(str1,pid+1,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/u0i_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) u0i
    close(2)


    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/v0i_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) v0i
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/a0i_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) a0i
    close(2)


    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/u0b_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) u0b
    close(2)


    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/v0b_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) v0b
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/a0b_'// trim(str1) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) a0b
    close(2)



    ! print*,"process id is", pid
    ! print*,"u0i is", u0i

    ! print*,"process id is", pid
    ! print*,"u0i is", u0b

END SUBROUTINE GetInitialCondSubdomain


SUBROUTINE GetMassDampMultiplier(pid,nbcount,ni,nb,deltaT,beta_NB,gamma_NB, u0i,u0b,v0i,v0b,a0i,a0b,tempi,tempb, PMMuvi, PMMuvb, PCMuvi, PCMuvb,nVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    ! Vec              Pu0, Pv0, Pa0
    Vec              PMMuvi,PMMuvb
    Vec              PCMuvi,PCMuvb

    PetscInt         s1, s2, s3
    PetscScalar      ivec
    PetscInt         ni, nb
    ! PetscScalar, pointer :: xx_v(:)

    double precision, dimension(:) :: u0i, u0b
    double precision, dimension(:) :: v0i, v0b
    double precision, dimension(:) :: a0i, a0b

    PetscScalar, allocatable, dimension(:) :: MMuvi, MMuvb
    PetscScalar, allocatable, dimension(:) :: CMuvi, CMuvb

    PetscErrorCode   ierr

    integer :: one, nbcount, nVec
    integer :: i, ii,id
    double precision :: deltaT, beta_NB, gamma_NB

    integer, dimension(ni*nVec) :: tempi
    integer, dimension(nb*nVec) :: tempb

    ! integer, allocatable, dimension(:) :: tb

    integer              :: pid
    character(len=255)   :: str1, str2, str3

    ! allocate(tb(nb), u0i(ni), u0b(nb), v0i(ni), v0b(nb), a0i(ni), a0b(nb))
    allocate(MMuvi(ni*nVec),MMuvb(nb*nVec), CMuvi(ni*nVec), CMuvb(nb*nVec))


    one = 1

    ! print*, 'Inside Mass Multiplier'
    ! print*,'process ID is', pid
    ! print*, 'u0i is', u0i
    ! print*, 'u0b is', u0b



    !!! VERIFIED UOI AND UOB for 2 subdomains giving correct values

    ! if (nbcount .eq. 2) print*,'process ID is', pid

    ! if (nbcount .eq. 2) print*, u0b

    ! call VecView(Pu0,PETSC_VIEWER_STDOUT_SELF);


    !!!! MMuvb verified with MATLAB for process 8
    !!!! MMuvi verified with MATLAB for process 2
    !!!! CMuvi verified with MATLAB for process 1
    !!!! CMuvb verified with MATLAB for process 2

    MMuvi = (u0i + deltaT * v0i)/(beta_NB*deltaT**2) + ((1-2*beta_NB)*a0i/(2*beta_NB))

    MMuvb = (u0b + deltaT * v0b)/(beta_NB*deltaT**2) + ((1-2*beta_NB)*a0b/(2*beta_NB))


    CMuvi = ( (gamma_NB*deltaT*MMuvi) - v0i - ((1-gamma_NB)*a0i*deltaT) )
    CMuvb = ( (gamma_NB*deltaT*MMuvb) - v0b - ((1-gamma_NB)*a0b*deltaT) )


!!!! MILESTONE MMuvi,MMuvb,CMuvi,CMuvb checked for nbcount 2...
    ! if (pid .eq. 0) then
    !     if(nbcount .eq. 2) then
    !         print*, 'procees is', pid
    !         print*, 'CMuvb  is', CMuvb
    !         stop 123
    !     endif
    ! end if

    ! print*, tempb

!!!!!!!!! Converting Back to PETSc vectors

    call VecCreateSeq(PETSC_COMM_SELF, ni*nVec, PMMuvi, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, nb*nVec, PMMuvb, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, ni*nVec, PCMuvi, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, nb*nVec, PCMuvb, ierr)


    call VecSetValues(PMMuvi,ni*nVec,tempi,MMuvi,INSERT_VALUES,ierr)
    call VecSetValues(PMMuvb,nb*nVec,tempb,MMuvb,INSERT_VALUES,ierr)
    call VecSetValues(PCMuvi,ni*nVec,tempi,CMuvi,INSERT_VALUES,ierr)
    call VecSetValues(PCMuvb,nb*nVec,tempb,CMuvb,INSERT_VALUES,ierr)


    call VecAssemblyBegin(PMMuvi,ierr)
    call VecAssemblyBegin(PMMuvb,ierr)
    call VecAssemblyBegin(PCMuvi,ierr)
    call VecAssemblyBegin(PCMuvb,ierr)

    call VecAssemblyEnd(PMMuvi,ierr)
    call VecAssemblyEnd(PMMuvb,ierr)
    call VecAssemblyEnd(PCMuvi,ierr)
    call VecAssemblyEnd(PCMuvb,ierr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DEALLOCATE(MMuvi,MMuvb, CMuvi, CMuvb)




if (pid .eq. 0) print*, '....................................................................'
if (pid .eq. 0) print*, '...............Mass and Damping Multipliers Calculated..............'
if (pid .eq. 0) print*, '....................................................................'


END SUBROUTINE GetMassDampMultiplier


SUBROUTINE GetForceVecSeq(pid, nbcount,tempi, tempb, PFi, PFg, ni, nb, nVec, ierr)


#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>


    Vec              PFi, PFg
    PetscScalar      ivec
    PetscErrorCode   ierr


    integer :: ni, nb, ii,id, pid, nbcount, one, nVec

    character(len=255)      :: str1, str3

    double precision, allocatable, dimension(:) :: bi, bg

    integer, dimension(ni*nVec) :: tempi
    integer, dimension(nb*nVec) :: tempb

    double precision :: zero

    allocate(bi(ni*nVec), bg(nb*nVec))

    one = 1
    zero = 0.0d0

    if (nbcount .eq. 1) then

        call VecCreateSeq(PETSC_COMM_SELF, ni*nVec, PFi, ierr)
        call VecCreateSeq(PETSC_COMM_SELF, nb*nVec, PFg, ierr)

        call int2str(str1,pid+1,1)
        str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/bi'// trim(str1) // '.dat'
        open(unit=1,file=str3,status='old')
        read(unit=1,fmt=*) bi
        close(1)

        str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/bg'// trim(str1) // '.dat'
        open(unit=2,file=str3,status='old')
        read(unit=2,fmt=*) bg
        close(2)



        do id = 1,ni*nVec
                ii = id-1
                ivec = bi(id)
                call VecSetValues(PFi,one,ii,ivec,INSERT_VALUES,ierr)
        end do

        ! print*, 'process id', pid
        ! print*, nb

        do id = 1,nb*nVec
                ii = id-1
                ivec = bg(id)
                call VecSetValues(PFg,one,ii,ivec,INSERT_VALUES,ierr)
        end do


        call VecAssemblyBegin(PFi,ierr)
        call VecAssemblyBegin(PFg,ierr)


        call VecAssemblyEnd(PFi,ierr)
        call VecAssemblyEnd(PFg,ierr)

    else



        call VecZeroEntries(PFi,ierr)
        call VecZeroEntries(PFg,ierr)
        ! do id = 1,ni
        !     ii = id-1
        !     ivec = zero
        !     call VecSetValues(PFi,one,ii,ivec,INSERT_VALUES,ierr)
        ! end do

        ! ! print*, 'process id', pid
        ! ! print*, nb

        ! do id = 1,nb
        !     ii = id-1
        !     ivec = 0.0d0
        !     call VecSetValues(PFg,one,ii,ivec,INSERT_VALUES,ierr)
        ! end do

        ! call VecSetValues(PFi, ni, tempi, PETSC_NULL_SCALAR, INSERT_VALUES,ierr)
        ! call VecSetValues(PFg, nb, tempb, PETSC_NULL_SCALAR, INSERT_VALUES,ierr)

    endif




    ! if (pid .eq. 0) then
    !     if(nbcount .eq. 2)then
    !         print*, 'procees is', pid
    !         print*, 'PFi  is'
    !         call VecView(PFi, PETSC_VIEWER_STDOUT_SELF, ierr)
    !         stop 123
    !     endif
    ! end if


DEALLOCATE(bi,bg)

END SUBROUTINE GetForceVecSeq

SUBROUTINE GetTransientForce(pid,nbcount,ni,nb,PMii,PMgg,PMgi,PCii,PCgg,PCgi, PMMuvi, PMMuvb, PCMuvi, PCMuvb, Mni, Mnb, Cni, Cnb,PFi, PFg,nVec, ierr)


#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
!!! Function to calculate the transient force vector which is the sum of force vector at n+1 th step and contributions from previous steps

    !!!! F_T = b_n+1 + M_n *Mass Multiplier + C_n * Damping Multiplier

    Mat              PMii,PMgg,PMgi
    Mat              PCii,PCgg,PCgi

    Vec              PMMuvi,PMMuvb
    Vec              PCMuvi,PCMuvb

    Vec              PFi,PFg

    PetscScalar      PosOne
    PetscErrorCode   ierr

    Vec              out1,out2, out3, out4
    Vec              Mni, Mnb !!! Force vector contribution from Mass times previous soultions
    Vec              Cni, Cnb !!! Force vector contribution from Damping times previous soultions

    integer             :: ni, nb, pid, nbcount, nVec

    ! PetscViewer       PETSC_VIEWER_STDOUT_SELF

    PosOne = 1

    call GetVecSeqTemp1(ni*nVec,out1,ierr)
    call GetVecSeqTemp1(nb*nVec,out2,ierr)
    call GetVecSeqTemp1(ni*nVec,out3,ierr)
    call GetVecSeqTemp1(nb*nVec,out4,ierr)

    call GetVecSeqTemp1(ni*nVec,Mni,ierr)
    call GetVecSeqTemp1(nb*nVec,Mnb,ierr)
    call GetVecSeqTemp1(ni*nVec,Cni,ierr)
    call GetVecSeqTemp1(nb*nVec,Cnb,ierr)



    ! print*,'process ID is', pid

    ! print*, 'PMMuvi is'

    ! call VecView(PMMuvi,PETSC_VIEWER_STDOUT_SELF,ierr);

    ! ! call MatView(PMii,PETSC_VIEWER_STDOUT_SELF);
    ! call VecView(Mni,PETSC_VIEWER_STDOUT_SELF,ierr);

    ! print*, 'out1 is', out1


    ! ........Mass Component of Force ................

    call MatMult(PMii,PMMuvi,out1, ierr)                  !!!!! out1 = PMii * PMMuvi              !!! M_ii * mmuv_i = out1

!!!!! MILESTONE checked out1 for nbcount = 2
    ! if (pid .eq. 0) then

    !     if(nbcount .eq. 2)then
    !         print*, 'procees is', pid
    !         print*, 'out1  is'
    !         call VecView(out1, PETSC_VIEWER_STDOUT_SELF, ierr)
    !         stop 123
    !     endif
    ! end if

    ! print*, 'process id', pid
    ! call VecView(PMMuvi, PETSC_VIEWER_STDOUT_SELF, ierr)
    ! print*, 'end vector'
    ! stop 123

    call MatMultTranspose(PMgi,PMMuvb,Mni, ierr)          !!!!! PMgi^T * PMMuvb = Mni              !!!! MatMultTransposeAdd doesn't work in Fortran

    call VecAXPY(Mni,PosOne,out1, ierr)                   !!!!! Mni = Mni + out1


!!!!!! Mni verified but values with very high e -16, 27 etc are not true ????



    call MatMult(PMgg,PMMuvb,out2,ierr)                  !!!!! out2 = PMgg * PMMuvb              !!! M_gg * mmuv_g = out2
    call MatMultAdd(PMgi,PMMuvi,out2,Mnb, ierr)           !!!!! PMgi* PMMuvi + out2 = Mnb         !!! M_gi * mmuv_i + out2 = Mn_b



    ! if (pid .eq. 0) then
    !     if(nbcount .eq. 2)then
    !         print*, 'procees is', pid
    !         print*, 'Mnb  is'
    !         call VecView(Mnb, PETSC_VIEWER_STDOUT_SELF, ierr)
    !         stop 123
    !     endif
    ! end if


!!!! MILESTONE checked Mni, Mnb for nbcount = 2



!...........Damping Component of Force  ..........................


    call MatMult(PCii,PCMuvi,out3, ierr)                  !!!!! out1 = PCii * PCMuvi              !!! C_ii * cmuv_i = out3
    call MatMultTranspose(PCgi,PCMuvb,Cni, ierr)          !!!!! PCgi^T * PCMuvb = Cni              !!!! MatMultTransposeAdd doesn't work in Fortran
    call VecAXPY(Cni,PosOne,out3, ierr)                   !!!!! Cni = Cni + out3

    ! call VecView(PMMuvb,PETSC_VIEWER_STDOUT_SELF, ierr);

    call MatMult(PCgg,PCMuvb,out4,ierr)                  !!!!! out4 = PCgg * PCMuvb              !!! C_gg * cmuv_g = out4
    call MatMultAdd(PCgi,PCMuvi,out4,Cnb, ierr)           !!!!! PCgi* PCMuvi + out4 = Cnb         !!! C_gi * cmuv_i + out4 = Cn_b


!!!!                 Adding the Components to Force Vector of Current Time step                                      !!!


    call VecAXPY(Cni,PosOne,Mni, ierr)               ! Cni = Cni + 1 * Mni


    ! if (pid .eq. 0) then
    !     if(nbcount .eq. 2)then
    !         print*, 'procees is', pid
    !         print*, 'PFi  is'
    !         call VecView(PFi, PETSC_VIEWER_STDOUT_SELF, ierr)
    !         stop 123
    !     endif
    ! end if


    call VecAXPY(PFi,PosOne,Cni, ierr)               ! PFi = PFi + 1 * Cni

    call VecAXPY(Cnb,PosOne,Mnb, ierr)               ! Cnb = Cnb + 1 * Mnb
    call VecAXPY(PFg,PosOne,Cnb, ierr)               ! PFg = PFg + 1 * Cnb


    ! if (pid .eq. 0) then
    !     if(nbcount .eq. 2) then
    !         print*, 'procees is', pid
    !         call VecView(PFg, PETSC_VIEWER_STDOUT_SELF, ierr)
    !         stop 123
    !     end if
    ! end if

    ! if (pid .eq. 0) then
    !     if(nbcount .eq. 2)then
    !         print*, 'procees is', pid
    !         print*, 'PFi  is'
    !         call VecView(PFi, PETSC_VIEWER_STDOUT_SELF, ierr)
    !         stop 123
    !     endif
    ! end if





END SUBROUTINE GetTransientForce


































!!!*********************************************************
!! Subroutine : for stochastic-sparse-vector assembly
SUBROUTINE StoVecSeqOneLevel(Fsi,Fsg,p,e,t,np,ne,nt,nb,ndim,npceout,nomga,nip,nbp,amp,dbounds, &
                             mIndex, sIndex, ierr) !mIndex

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              Fsi, Fsg
    PetscInt         nip, nbp
    PetscErrorCode   ierr

    integer :: np, ne, nt, nb, ndim, npceout,nomga
    double precision                  :: amp
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(2,2)  :: dbounds

    integer, dimension(ndim,npceout) :: mIndex
    integer, dimension(2,ndim) :: sIndex

    !!integer, dimension(nip) :: tempi
    !!integer, dimension(nbp) :: tempb
    !!nip = ((np-nb)*npceout)
    !!nbp = nb*npceout

    call VecCreateSeq(PETSC_COMM_SELF, nip, Fsi, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, nbp, Fsg, ierr)
    !!call VecDuplicate(Fsi,SolVeci2,ierr)
    !!call VecDuplicate(Fsg,SolVecb2,ierr)

    call PETscSubVecAssembly(Fsi, Fsg, p, e, t, np, ne, nt, nb,ndim,npceout,nomga, &
                             amp, dbounds, mIndex, sIndex) !mIndex

    call VecAssemblyBegin(Fsi,ierr)
    call VecAssemblyBegin(Fsg,ierr)
    call VecAssemblyEnd(Fsi,ierr)
    call VecAssemblyEnd(Fsg,ierr)


END SUBROUTINE StoVecSeqOneLevel


!!!*********************************************************
!!! Subroutine: To construc PETSc Matrix in most general way
!! Subroutine : for stochastic-sparse-matrix assembly
SUBROUTINE StoMatSeqOneLevel( pid, Asmat, Asii, Asgg, Asgi, Asri, Asrr, Ascc, Asci, Ascr, &
                           p,e,t,np,ne,nt,nb,nci,nri,ndim,npcein,npceout,nomga,ncijk,ijk, &
                           cijk,nid,nip,nbp, ncp, nrp, nirp, dbounds, const_diff, omegas, &
                           mIndex, sIndex, multipliers,sigma,casep,ierr) !mIndex

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              Asmat,Asii,Asgg,Asrr,Ascc
    Mat              Asgi,Asri,Asci,Ascr

    PetscInt         nip,nbp,one
    PetscInt         nirp, nrp, ncp
    PetscInt         ii, jj, ii2, jj2
    PetscInt         id, jd, nni, nnj
    PetscScalar      imat
    PetscErrorCode   ierr

    integer :: i, j, k, indexi, ncijk, ndim, npceout,npcein
    integer :: np, ne, nt,nb,nid,nci,nri,casep,nomga
    double precision :: const_diff, sigma

    integer, dimension(ncijk,3)       :: ijk
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(ncijk):: cijk
    double precision, dimension(2,2)  :: dbounds
    double precision, dimension(nomga):: omegas, multipliers

    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    integer              :: pid
    character(len=255)   :: str1
    PetscInt,ALLOCATABLE :: nnzi(:),nnzb(:),nnzbi(:),nnzri(:),nnzci(:)
    PetscInt,ALLOCATABLE :: nnzcr(:), nnzcc(:), nnzir(:), nnzrr(:)

    double precision, allocatable, dimension(:,:)  :: Adii,Adgg,Adrr,Adcc
    double precision, allocatable, dimension(:,:)  :: Adgi,Adri,Adci,Adcr

!!-----------------------------------------------------------------------------------------------
    allocate(Adii(nid,nid), Adgg(nb,nb), Adrr(nri,nri), Adcc(nci,nci),&
             Adgi(nb,nid), Adri(nri,nid), Adci(nci,nid), Adcr(nci,nri))

!!-----------------------------------------------------------------------------------------------
    allocate(nnzi(nip),nnzb(nbp),nnzbi(nbp),nnzri(nrp),nnzrr(nrp))
    allocate(nnzci(ncp), nnzcr(ncp), nnzcc(ncp), nnzir(nirp))

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzir' // trim(str1) // '.dat'
    open(unit=9,file=str1,status='old')
    read(unit=9,fmt=*) nnzir
    close(9)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzi' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    read(unit=1,fmt=*) nnzi
    close(1)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzb' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    read(unit=2,fmt=*) nnzb
    close(2)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzbi' // trim(str1) // '.dat'
    open(unit=3,file=str1,status='old')
    read(unit=3,fmt=*) nnzbi
    close(3)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzri' // trim(str1) // '.dat'
    open(unit=4,file=str1,status='old')
    read(unit=4,fmt=*) nnzri
    close(4)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzci' // trim(str1) // '.dat'
    open(unit=5,file=str1,status='old')
    read(unit=5,fmt=*) nnzci
    close(5)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzrr' // trim(str1) // '.dat'
    open(unit=66,file=str1,status='old')
    read(unit=66,fmt=*) nnzrr
    close(66)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzcr' // trim(str1) // '.dat'
    open(unit=7,file=str1,status='old')
    read(unit=7,fmt=*) nnzcr
    close(7)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzcc' // trim(str1) // '.dat'
    open(unit=8,file=str1,status='old')
    read(unit=8,fmt=*) nnzcc
    close(8)

    call MatCreateSeqAIJ(PETSC_COMM_SELF, nirp, nirp, 0, nnzir, Asmat,ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nip, 0, nnzi,  Asii, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nip, 0, nnzbi, Asgi, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, 0, nnzb,  Asgg, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nip, 0, nnzri, Asri, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nip, 0, nnzci, Asci, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nrp, 0, nnzrr, Asrr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nrp, 0, nnzcr, Ascr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, ncp, 0, nnzcc, Ascc, ierr)

!!!-----------------------------------------------------------------------------------------------
!!! Need to optimize pre-memory allocations
!!! For nOrd = 1; nDim = 50; nPCE=51;   nMesh=1.5K; nzi = 600,  nzb = 300,  nzc = 200
!!! For nOrd = 2; nDim = 15; nPCE=136;  nMesh=1.5K; nzi = 1600, nzb = 1100, nzc = 700
!!! For nOrd = 2; nDim = 20; nPCE=231;  nMesh=1.5K; nzi = 2500, nzb = 1500, nzc = 1000
!!! For nOrd = 2; nDim = 25; nPCE=351;  nMesh=1.5K; nzi = 3500, nzb = 2500, nzc = 1500
!!! For nOrd = 2; nDim = 50; nPCE=1326; nMesh=1.5K; nzi = 12000, nzb = 9000, nzc = 6000
!
!!! HPCluster:
!!! For nOrd = 2; nDim = 25; nPCE=351;  nMesh=13K; nzi = 5000, nzb = 4000, nzc = 3000
!!! For nOrd = 2; nDim = 50; nPCE=351;  nMesh=5K;  nzi =12000, nzb = 8000, nzc = 5000
!!! For nOrd = 2; nDim = 50; nPCE=351;  nMesh=5K;  nzi =12000, nzb = 8000, nzc = 6000 (Acr,Acc)
!
!
!    nzi = 500 !! number of non-zeros per row for matrix mallocs
!    nzb = 300 !! for ndim = 3/4: nzi = 500, nzb = 300, nzc = 200
!    nzc = 200 !! for ndim = 5;
!
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nirp, nirp, nzi, PETSC_NULL_INTEGER, Asmat,ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nip, nzi, PETSC_NULL_INTEGER, Asii, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, nzb, PETSC_NULL_INTEGER, Asgg, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nip, nzb, PETSC_NULL_INTEGER, Asgi, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nip, nzb, PETSC_NULL_INTEGER, Asri, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nrp, nzb, PETSC_NULL_INTEGER, Asrr, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nip, nzb, PETSC_NULL_INTEGER, Asci, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nrp, nzb, PETSC_NULL_INTEGER, Ascr, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, ncp, nzc, PETSC_NULL_INTEGER, Ascc, ierr)
!
!    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nbp, nzb, PETSC_NULL_INTEGER, Asig, ierr)
!    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nrp, nzb, PETSC_NULL_INTEGER, Asir, ierr)
!    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, ncp, nzb, PETSC_NULL_INTEGER, Asic, ierr)
!    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, ncp, nzb, PETSC_NULL_INTEGER, Asrc, ierr)
!
!
!    !!call MatSetOption(Asii,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
!    !!call MatSetOption(Asgg,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
!    !!call MatSetOption(Asii,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
!    !!call MatSetOption(Asgg,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr)
!
!!-----------------------------------------------------------------------------------------------
!! Method : 3
one    = 1
indexi = 1

do k = 1,npcein

    !! Deterministic Matices
    Adii(:,:) = 0.0d0
    Adgg(:,:) = 0.0d0
    Adgi(:,:) = 0.0d0
    Adri(:,:) = 0.0d0
    Adci(:,:) = 0.0d0
    Adcr(:,:) = 0.0d0
    Adcc(:,:) = 0.0d0
    Adrr(:,:) = 0.0d0
    !!Adig(:,:) = 0.0d0
    !!Adir(:,:) = 0.0d0
    !!Adic(:,:) = 0.0d0
    !!Adrc(:,:) = 0.0d0

    if (k .eq. 1) then
        !! Det-Advection Matrix
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,1,&
                        0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
                        dbounds,mIndex,sIndex,0,0)
        !! Det-Diffusion Matrix
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,2,&
                const_diff,omegas,multipliers,sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,1)
    else
        !! Sto-Diffusion Matrices
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,2,&
                const_diff,omegas,multipliers,sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,k)
    end if


!!-----------------------------------------------------------------------------------------------
    do i = 1,npceout
    do j = 1,npceout

        if ((k .eq. ijk(indexi,1)) .and. (i .eq. ijk(indexi,2)) .and. (j .eq. ijk(indexi,3))) then

    !!---------------------------------------->Asii/Asmat
        nni = ((i-1)*nid)
        nnj = ((j-1)*nid)

        do id = 1,nid
        do jd = 1,nid
            if (Adii(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adii(id,jd)
                call MatSetValues(Asii,one,ii,one,jj,imat,ADD_VALUES,ierr)
                call MatSetValues(Asmat,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!---------------------------------------->Asgg
        nni = ((i-1)*nb)
        nnj = ((j-1)*nb)

        do id = 1,nb
        do jd = 1,nb
            if (Adgg(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adgg(id,jd)
                call MatSetValues(Asgg,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        !!!---------------------------------------->Asig
        !    nni = ((i-1)*nid)
        !    nnj = ((j-1)*nb)
        !
        !    do id = 1,nid
        !    do jd = 1,nb
        !        if (Adig(id,jd) .ne. 0) then
        !            ii = (nni+(id-1))
        !            jj = (nnj+(jd-1))
        !            imat = cijk(indexi)*Adig(id,jd)
        !            call MatSetValues(Asig,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !        end if
        !    end do
        !    end do

    !!---------------------------------------->Asgi
        nni = ((i-1)*nb)
        nnj = ((j-1)*nid)

        do id = 1,nb
        do jd = 1,nid
            if (Adgi(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adgi(id,jd)
                call MatSetValues(Asgi,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!!---------------------------------------->Asir/Asmat
        !    nni = ((i-1)*nid)
        !    nnj = ((j-1)*nri)
        !
        !    do id = 1,nid
        !    do jd = 1,nri
        !        if (Adir(id,jd) .ne. 0) then
        !            ii = (nni+(id-1))
        !            jj = (nnj+(jd-1))
        !            imat = cijk(indexi)*Adir(id,jd)
        !            call MatSetValues(Asir,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !
        !            jj2 = nip + jj
        !            call MatSetValues(Asmat,one,ii,one,jj2,imat,ADD_VALUES,ierr)
        !
        !        end if
        !    end do
        !    end do

    !!---------------------------------------->Asri/Asmat
        nni = ((i-1)*nri)
        nnj = ((j-1)*nid)

        do id = 1,nri
        do jd = 1,nid
            if (Adri(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adri(id,jd)
                call MatSetValues(Asri,one,ii,one,jj,imat,ADD_VALUES,ierr)

                ii2 = nip + ii
                call MatSetValues(Asmat,one,ii2,one,jj,imat,ADD_VALUES,ierr)
                call MatSetValues(Asmat,one,jj,one,ii2,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!---------------------------------------->Asrr/Asmat
        nni = ((i-1)*nri)
        nnj = ((j-1)*nri)

        do id = 1,nri
        do jd = 1,nri
            if (Adrr(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adrr(id,jd)
                call MatSetValues(Asrr,one,ii,one,jj,imat,ADD_VALUES,ierr)

                ii2 = nip + ii
                jj2 = nip + jj
                call MatSetValues(Asmat,one,ii2,one,jj2,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!!---------------------------------------->Asic
        !    nni = ((i-1)*nid)
        !    nnj = ((j-1)*nci)
        !
        !    do id = 1,nid
        !    do jd = 1,nci
        !        if (Adic(id,jd) .ne. 0) then
        !            ii = (nni+(id-1))
        !            jj = (nnj+(jd-1))
        !            imat = cijk(indexi)*Adic(id,jd)
        !            call MatSetValues(Asic,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !        end if
        !    end do
        !    end do

    !!---------------------------------------->Asci
        nni = ((i-1)*nci)
        nnj = ((j-1)*nid)

        do id = 1,nci
        do jd = 1,nid
            if (Adci(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adci(id,jd)
                call MatSetValues(Asci,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!!---------------------------------------->Asrc
        !    nni = ((i-1)*nri)
        !    nnj = ((j-1)*nci)
        !
        !    do id = 1,nri
        !    do jd = 1,nci
        !        if (Adrc(id,jd) .ne. 0) then
        !            ii = (nni+(id-1))
        !            jj = (nnj+(jd-1))
        !            imat = cijk(indexi)*Adrc(id,jd)
        !            call MatSetValues(Asrc,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !        end if
        !    end do
        !    end do

    !!---------------------------------------->Ascr
        nni = ((i-1)*nci)
        nnj = ((j-1)*nri)

        do id = 1,nci
        do jd = 1,nri
            if (Adcr(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adcr(id,jd)
                call MatSetValues(Ascr,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do


    !!---------------------------------------->Ascc
        nni = ((i-1)*nci)
        nnj = ((j-1)*nci)

        do id = 1,nci
        do jd = 1,nci
            if (Adcc(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adcc(id,jd)
                call MatSetValues(Ascc,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        indexi = indexi+1
        end if

    end do
    end do
end do

deallocate(Adii,Adgg,Adrr,Adcc,Adgi,Adri,Adci,Adcr)
deallocate(nnzb,nnzi,nnzbi,nnzri,nnzci,nnzrr,nnzcr,nnzcc,nnzir)

call MatAssemblyBegin(Asmat,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asii,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asgg,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asgi,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asrr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asri,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Ascc,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asci,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Ascr,MAT_FINAL_ASSEMBLY,ierr)

call MatAssemblyEnd(Asii,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asgg,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asgi,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asrr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asri,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Ascc,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asci,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Ascr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asmat,MAT_FINAL_ASSEMBLY,ierr)

!!call MatAssemblyBegin(Asig,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asir,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asic,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asrc,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asig,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asir,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asic,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asrc,MAT_FINAL_ASSEMBLY,ierr)


END SUBROUTINE StoMatSeqOneLevel



!!!*********************************************************
!! PETSc-Vec assembly using FEniCS assembled ddm Vecs # NEW
SUBROUTINE StoVecSeqOneLevelShort(pid, Fsi, Fsg, np, nb, npceout, nip, nbp, nVec, ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              Fsi, Fsg
    PetscInt         nip, nbp
    PetscScalar      ivec
    PetscErrorCode   ierr

    integer :: np, nb, npceout

    double precision, allocatable, dimension(:) :: bi, bg
    character(len=255)      :: str2, str3
    integer :: pid, id, ii, one, nVec, nid, nbd

    REAL*8, PARAMETER :: EPSILON = 1d-14

    !!-----------------------------------------------------------------------------------------------
    !! Allocate memory to read pre-assembled vectors
    nid = (np-nb)*nVec
    nbd = nb*nVec
    one= 1

    allocate(bi(nid), bg(nbd))

    !!-----------------------------------------------------------------------------------------------
    !! FEniCS based Vec-Assembly procedure (using preassembled vecs)
    call VecCreateSeq(PETSC_COMM_SELF, nip, Fsi, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, nbp, Fsg, ierr)

    bi = 0.0d0
    bg = 0.0d0

    call int2str(str2,pid+1,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str2) // '/bi1.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) bi
    close(1)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str2) // '/bg1.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) bg
    close(2)

    !---------------------------------------->fsi
    do id = 1,nid
        if (ABS( bi(id) ) > EPSILON) then
        ii = (id-1)
        ivec = bi(id)
        call VecSetValues(Fsi,one,ii,ivec,ADD_VALUES,ierr)
        end if
    end do

    !!---------------------------------------->fsg
    do id = 1,nbd
        if (ABS( bg(id) ) > EPSILON) then
        ii = (id-1)
        ivec = bg(id)
        call VecSetValues(Fsg,one,ii,ivec,ADD_VALUES,ierr)
        end if
    end do

    call VecAssemblyBegin(Fsi,ierr)
    call VecAssemblyBegin(Fsg,ierr)
    call VecAssemblyEnd(Fsi,ierr)
    call VecAssemblyEnd(Fsg,ierr)


DEALLOCATE(bi,bg)

END SUBROUTINE StoVecSeqOneLevelShort



!!!*********************************************************
!!! Subroutine: To construc PETSc Matrix in most general way
!! Subroutine : for stochastic-sparse-matrix assembly
!SUBROUTINE StoMatSeqOneLevelShort( pid, Asmat, Asii, Asgg, Asgi, Asri, Asrr, Ascc, Asci, Ascr, &
!                            p,e,t,np,ne,nt,nb,nci,nri,ndim,npcein,npceout,nomga,ncijk,ijk, &
!                            cijk,nid,nip,nbp, ncp, nrp, nirp, dbounds, const_diff, omegas, &
!                            mIndex, sIndex, multipliers,sigma,casep,ierr) !mIndex

SUBROUTINE StoMatSeqOneLevelShort( pid, Asmat, Asii, Asgg, Asgi, Asri, Asrr, Ascc, Asci, Ascr, &
                 np,nb,nc,nr,npcein,npceout,ncijk,ijk,cijk,ni,nip,nbp,ncp,nrp,nirp,nVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              Asmat,Asii,Asgg,Asrr,Ascc
    Mat              Asgi,Asri,Asci,Ascr

    PetscInt         nip,nbp,one
    PetscInt         nirp, nrp, ncp
    PetscInt         ii, jj, ii2, jj2
    PetscInt         id, jd, nni, nnj
    PetscScalar      imat
    PetscErrorCode   ierr

    integer :: i, j, k, pid, indexi, ncijk, npceout, npcein  !!ndim
    integer :: np, nb, ni, nc, nr, nid, nbd, ncd, nrd, nVec  !ne, nt,casep,nomga
    !!double precision :: const_diff, sigma

    integer, dimension(ncijk,3)       :: ijk
    double precision, dimension(ncijk):: cijk
    !integer, dimension(3,ne)          :: e
    !integer, dimension(3,nt)          :: t
    !double precision, dimension(2,np) :: p
    !double precision, dimension(2,2)  :: dbounds
    !double precision, dimension(nomga):: omegas, multipliers
    !integer, dimension(ndim,npceout)  :: mIndex
    !integer, dimension(2,ndim)        :: sIndex

    character(len=255)   :: str1, str2, str3
    PetscInt,ALLOCATABLE :: nnzi(:),nnzb(:),nnzbi(:),nnzri(:),nnzci(:)
    PetscInt,ALLOCATABLE :: nnzcr(:), nnzcc(:), nnzir(:), nnzrr(:)

    double precision, allocatable, dimension(:,:)  :: Adii,Adgg,Adrr,Adcc
    double precision, allocatable, dimension(:,:)  :: Adgi,Adri,Adci,Adcr

    REAL*8, PARAMETER :: EPSILON = 1d-14

    ncd = nc*nVec
    nrd = nr*nVec
    nid = ni*nVec
    nbd = nb*nVec

    !!-----------------------------------------------------------------------------------------------
    allocate(Adii(nid,nid), Adgg(nbd,nbd), Adrr(nrd,nrd), Adcc(ncd,ncd),&
    Adgi(nbd,nid), Adri(nrd,nid), Adci(ncd,nid), Adcr(ncd,nrd))

    !!-----------------------------------------------------------------------------------------------
    allocate(nnzi(nip),nnzb(nbp),nnzbi(nbp),nnzri(nrp),nnzrr(nrp))
    allocate(nnzci(ncp), nnzcr(ncp), nnzcc(ncp), nnzir(nirp))

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzir' // trim(str1) // '.dat'
    open(unit=9,file=str1,status='old')
    read(unit=9,fmt=*) nnzir
    close(9)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzi' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    read(unit=1,fmt=*) nnzi
    close(1)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzb' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    read(unit=2,fmt=*) nnzb
    close(2)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzbi' // trim(str1) // '.dat'
    open(unit=3,file=str1,status='old')
    read(unit=3,fmt=*) nnzbi
    close(3)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzri' // trim(str1) // '.dat'
    open(unit=4,file=str1,status='old')
    read(unit=4,fmt=*) nnzri
    close(4)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzci' // trim(str1) // '.dat'
    open(unit=5,file=str1,status='old')
    read(unit=5,fmt=*) nnzci
    close(5)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzrr' // trim(str1) // '.dat'
    open(unit=66,file=str1,status='old')
    read(unit=66,fmt=*) nnzrr
    close(66)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzcr' // trim(str1) // '.dat'
    open(unit=7,file=str1,status='old')
    read(unit=7,fmt=*) nnzcr
    close(7)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzcc' // trim(str1) // '.dat'
    open(unit=8,file=str1,status='old')
    read(unit=8,fmt=*) nnzcc
    close(8)

!!-----------------------------------------------------------------------------------------------

    call MatCreateSeqAIJ(PETSC_COMM_SELF, nirp, nirp, 0, nnzir, Asmat,ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nip, 0, nnzi,  Asii, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nip, 0, nnzbi, Asgi, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, 0, nnzb,  Asgg, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nip, 0, nnzri, Asri, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nip, 0, nnzci, Asci, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nrp, 0, nnzrr, Asrr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nrp, 0, nnzcr, Ascr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, ncp, 0, nnzcc, Ascc, ierr)

!! Method : 3
one    = 1
indexi = 1

do k = 1,npcein

    !! Deterministic Matices
    Adii(:,:) = 0.0d0
    Adgg(:,:) = 0.0d0
    Adgi(:,:) = 0.0d0
    Adri(:,:) = 0.0d0
    Adci(:,:) = 0.0d0
    Adcr(:,:) = 0.0d0
    Adcc(:,:) = 0.0d0
    Adrr(:,:) = 0.0d0
    !!Adig(:,:) = 0.0d0
    !!Adir(:,:) = 0.0d0
    !!Adic(:,:) = 0.0d0
    !!Adrc(:,:) = 0.0d0

    !!!-----------------------------------------------------------------------------------------------
    !!! Python-dolfin
    !!! The pressembled FEniCS DD Matrices for each subdomain for each PCE mode
    call int2str(str1,pid+1,1)
    call int2str(str2,k,1)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADii' // trim(str2) // '.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adii
    close(1)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgg' // trim(str2) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adgg
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgi' // trim(str2) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adgi
    close(3)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADci' // trim(str2) // '.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adci
    close(1)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADri' // trim(str2) // '.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adri
    close(2)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADcc' // trim(str2) // '.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adcc
    close(3)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADrr' // trim(str2) // '.dat'
    open(unit=4,file=str3,status='old')
    read(unit=4,fmt=*) Adrr
    close(4)

    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADcr' // trim(str2) // '.dat'
    open(unit=5,file=str3,status='old')
    read(unit=5,fmt=*) Adcr
    close(5)


    !!-----------------------------------------------------------------------------------------------
    do i = 1,npceout
    do j = 1,npceout

        if ((k .eq. ijk(indexi,1)) .and. (i .eq. ijk(indexi,2)) .and. (j .eq. ijk(indexi,3))) then

        !!---------------------------------------->Asii/Asmat
        nni = ((i-1)*nid)
        nnj = ((j-1)*nid)

        do id = 1,nid
        do jd = 1,nid
            if (ABS( Adii(id,jd) ) > EPSILON) then
            ii = (nni+(id-1))
            jj = (nnj+(jd-1))
            imat = cijk(indexi)*Adii(id,jd)
            call MatSetValues(Asii,one,ii,one,jj,imat,ADD_VALUES,ierr)
            call MatSetValues(Asmat,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        !!---------------------------------------->Asgg
        nni = ((i-1)*nbd)
        nnj = ((j-1)*nbd)

        do id = 1,nbd
        do jd = 1,nbd
            if (ABS( Adgg(id,jd) ) > EPSILON) then
            ii = (nni+(id-1))
            jj = (nnj+(jd-1))
            imat = cijk(indexi)*Adgg(id,jd)
            call MatSetValues(Asgg,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        !!!---------------------------------------->Asig
        !    nni = ((i-1)*nid)
        !    nnj = ((j-1)*nbd)
        !
        !    do id = 1,nid
        !    do jd = 1,nbd
        !        if (ABS( Adig(id,jd) ) > EPSILON)
        !            ii = (nni+(id-1))
        !            jj = (nnj+(jd-1))
        !            imat = cijk(indexi)*Adig(id,jd)
        !            call MatSetValues(Asig,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !        end if
        !    end do
        !    end do

        !!---------------------------------------->Asgi
        nni = ((i-1)*nbd)
        nnj = ((j-1)*nid)

        do id = 1,nbd
        do jd = 1,nid
            if (ABS( Adgi(id,jd) ) > EPSILON) then
            ii = (nni+(id-1))
            jj = (nnj+(jd-1))
            imat = cijk(indexi)*Adgi(id,jd)
            call MatSetValues(Asgi,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        !!!---------------------------------------->Asir/Asmat
        !    nni = ((i-1)*nid)
        !    nnj = ((j-1)*nrd)
        !
        !    do id = 1,nid
        !    do jd = 1,nrd
        !        if (ABS( Adir(id,jd) ) > EPSILON) then
        !            ii = (nni+(id-1))
        !            jj = (nnj+(jd-1))
        !            imat = cijk(indexi)*Adir(id,jd)
        !            call MatSetValues(Asir,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !
        !            jj2 = nip + jj
        !            call MatSetValues(Asmat,one,ii,one,jj2,imat,ADD_VALUES,ierr)
        !
        !        end if
        !    end do
        !    end do

        !!---------------------------------------->Asri/Asmat
        nni = ((i-1)*nrd)
        nnj = ((j-1)*nid)

        do id = 1,nrd
        do jd = 1,nid
            if (ABS( Adri(id,jd) ) > EPSILON) then
            ii = (nni+(id-1))
            jj = (nnj+(jd-1))
            imat = cijk(indexi)*Adri(id,jd)
            call MatSetValues(Asri,one,ii,one,jj,imat,ADD_VALUES,ierr)

            ii2 = nip + ii
            call MatSetValues(Asmat,one,ii2,one,jj,imat,ADD_VALUES,ierr)
            call MatSetValues(Asmat,one,jj,one,ii2,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        !!---------------------------------------->Asrr/Asmat
        nni = ((i-1)*nrd)
        nnj = ((j-1)*nrd)

        do id = 1,nrd
        do jd = 1,nrd
            if (ABS( Adrr(id,jd) ) > EPSILON) then
            ii = (nni+(id-1))
            jj = (nnj+(jd-1))
            imat = cijk(indexi)*Adrr(id,jd)
            call MatSetValues(Asrr,one,ii,one,jj,imat,ADD_VALUES,ierr)

            ii2 = nip + ii
            jj2 = nip + jj
            call MatSetValues(Asmat,one,ii2,one,jj2,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        !!!---------------------------------------->Asic
        !    nni = ((i-1)*nid)
        !    nnj = ((j-1)*ncd)
        !
        !    do id = 1,nid
        !    do jd = 1,ncd
        !        if (Adic(id,jd) .ne. 0) then
        !            ii = (nni+(id-1))
        !            jj = (nnj+(jd-1))
        !            imat = cijk(indexi)*Adic(id,jd)
        !            call MatSetValues(Asic,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !        end if
        !    end do
        !    end do

        !!---------------------------------------->Asci
        nni = ((i-1)*ncd)
        nnj = ((j-1)*nid)

        do id = 1,ncd
        do jd = 1,nid
            if (ABS( Adci(id,jd) ) > EPSILON) then
            ii = (nni+(id-1))
            jj = (nnj+(jd-1))
            imat = cijk(indexi)*Adci(id,jd)
            call MatSetValues(Asci,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        !!!---------------------------------------->Asrc
        !    nni = ((i-1)*nrd)
        !    nnj = ((j-1)*ncd)
        !
        !    do id = 1,nrd
        !    do jd = 1,ncd
        !        if (Adrc(id,jd) .ne. 0) then
        !            ii = (nni+(id-1))
        !            jj = (nnj+(jd-1))
        !            imat = cijk(indexi)*Adrc(id,jd)
        !            call MatSetValues(Asrc,one,ii,one,jj,imat,ADD_VALUES,ierr)
        !        end if
        !    end do
        !    end do

        !!---------------------------------------->Ascr
        nni = ((i-1)*ncd)
        nnj = ((j-1)*nrd)

        do id = 1,ncd
        do jd = 1,nrd
            if (ABS( Adcr(id,jd) ) > EPSILON) then
            ii = (nni+(id-1))
            jj = (nnj+(jd-1))
            imat = cijk(indexi)*Adcr(id,jd)
            call MatSetValues(Ascr,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do


        !!---------------------------------------->Ascc
        nni = ((i-1)*ncd)
        nnj = ((j-1)*ncd)

        do id = 1,ncd
        do jd = 1,ncd
            if (ABS( Adcc(id,jd) ) > EPSILON) then
            ii = (nni+(id-1))
            jj = (nnj+(jd-1))
            imat = cijk(indexi)*Adcc(id,jd)
            call MatSetValues(Ascc,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        indexi = indexi+1
        end if

    end do
    end do
end do

deallocate(Adii,Adgg,Adrr,Adcc,Adgi,Adri,Adci,Adcr)
deallocate(nnzb,nnzi,nnzbi,nnzri,nnzci,nnzrr,nnzcr,nnzcc,nnzir)

call MatAssemblyBegin(Asmat,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asii,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asgg,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asgi,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asrr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asri,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Ascc,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asci,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Ascr,MAT_FINAL_ASSEMBLY,ierr)

call MatAssemblyEnd(Asii,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asgg,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asgi,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asrr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asri,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Ascc,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asci,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Ascr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asmat,MAT_FINAL_ASSEMBLY,ierr)

!!call MatAssemblyBegin(Asig,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asir,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asic,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asrc,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asig,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asir,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asic,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asrc,MAT_FINAL_ASSEMBLY,ierr)


END SUBROUTINE StoMatSeqOneLevelShort






!!!*********************************************************
!!! Subroutine: To construc PETSc Matrix in most general way
!! Subroutine : for stochastic-sparse-matrix assembly
SUBROUTINE StoMatSeqOneLevelFullMallocs(Asmat, Asii, Asgg, Asgi, Asri, Asrr, Ascc, Asci, Ascr, &
                             p,e,t,np,ne,nt,nb,nci,nri,ndim,npceout,nomga,ncijk,ijk,cijk, &
                             nid,nip,nbp,ncp,nrp,nirp,dbounds,const_diff,omegas, &
                             mIndex, sIndex, multipliers,sigma,casep,ierr) !mIndex

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              Asmat,Asii,Asgg,Asrr,Ascc
    Mat              Asgi,Asri,Asci,Ascr

    PetscInt         nip,nbp,nzi,nzb,nzc,one
    PetscInt         nirp, nrp, ncp
    PetscInt         ii, jj, ii2, jj2
    PetscInt         id, jd, nni, nnj
    PetscScalar      imat
    PetscErrorCode   ierr

    integer :: i, j, k, indexi, ncijk, ndim, npceout
    integer :: np, ne, nt, nb, nid, nci, nri, casep, nomga
    double precision :: const_diff, sigma

    integer, dimension(ncijk,3)       :: ijk
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(ncijk):: cijk
    double precision, dimension(2,2)  :: dbounds
    double precision, dimension(nomga):: omegas, multipliers

    double precision, allocatable, dimension(:,:)  :: Adii,Adgg,Adrr,Adcc
    double precision, allocatable, dimension(:,:)  :: Adgi,Adri,Adci,Adcr

    integer, dimension(ndim,npceout) :: mIndex
    integer, dimension(2,ndim) :: sIndex

!!-----------------------------------------------------------------------------------------------
    allocate(Adii(nid,nid), Adgg(nb,nb), Adrr(nri,nri), Adcc(nci,nci),&
             Adgi(nb,nid), Adri(nri,nid), Adci(nci,nid), Adcr(nci,nri))

!!-----------------------------------------------------------------------------------------------
!! Need to optimize pre-memory allocations
!! For nOrd = 1; nDim = 50; nPCE=51;   nMesh=1.5K; nzi = 600,  nzb = 300,  nzc = 200
!! For nOrd = 2; nDim = 15; nPCE=136;  nMesh=1.5K; nzi = 1600, nzb = 1100, nzc = 700
!! For nOrd = 2; nDim = 20; nPCE=231;  nMesh=1.5K; nzi = 2500, nzb = 1500, nzc = 1000
!! For nOrd = 2; nDim = 25; nPCE=351;  nMesh=1.5K; nzi = 3500, nzb = 2500, nzc = 1500
!! For nOrd = 2; nDim = 50; nPCE=1326; nMesh=1.5K; nzi = 12000, nzb = 9000, nzc = 6000

!! HPCluster:
!! For nOrd = 2; nDim = 25; nPCE=351;  nMesh=13K; nzi = 5000, nzb = 4000, nzc = 3000
!! For nOrd = 2; nDim = 50; nPCE=351;  nMesh=5K;  nzi =12000, nzb = 8000, nzc = 5000
!! For nOrd = 2; nDim = 50; nPCE=351;  nMesh=5K;  nzi =12000, nzb = 8000, nzc = 6000 (Acr,Acc)

    nzi = 3500 !! number of non-zeros per row for matrix mallocs
    nzb = 2500 !! for ndim = 3/4: nzi = 500, nzb = 300, nzc = 200
    nzc = 1500 !! for ndim = 5;

    call MatCreateSeqAIJ(PETSC_COMM_SELF, nirp, nirp, nzi, PETSC_NULL_INTEGER, Asmat,ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nip, nzi, PETSC_NULL_INTEGER, Asii, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, nzb, PETSC_NULL_INTEGER, Asgg, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nip, nzb, PETSC_NULL_INTEGER, Asgi, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nip, nzb, PETSC_NULL_INTEGER, Asri, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nrp, nzb, PETSC_NULL_INTEGER, Asrr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nip, nzb, PETSC_NULL_INTEGER, Asci, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nrp, nzb, PETSC_NULL_INTEGER, Ascr, ierr)
    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, ncp, nzc, PETSC_NULL_INTEGER, Ascc, ierr)

    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nbp, nzb, PETSC_NULL_INTEGER, Asig, ierr)
    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, nrp, nzb, PETSC_NULL_INTEGER, Asir, ierr)
    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nip, ncp, nzb, PETSC_NULL_INTEGER, Asic, ierr)
    !!call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, ncp, nzb, PETSC_NULL_INTEGER, Asrc, ierr)

!!-----------------------------------------------------------------------------------------------
!! Method : 3
one    = 1
indexi = 1

do k = 1,npceout

    !! Deterministic Matices
    Adii(:,:) = 0.0d0
    Adgg(:,:) = 0.0d0
    Adgi(:,:) = 0.0d0
    Adri(:,:) = 0.0d0
    Adci(:,:) = 0.0d0
    Adcr(:,:) = 0.0d0
    Adcc(:,:) = 0.0d0
    Adrr(:,:) = 0.0d0
    !!Adig(:,:) = 0.0d0
    !!Adir(:,:) = 0.0d0
    !!Adic(:,:) = 0.0d0
    !!Adrc(:,:) = 0.0d0

!! UnderStudy: July 7, 2016

    if (k .eq. 1) then
        !! Det-Advection Matrix
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,1,&
                        0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
                        dbounds,mIndex,sIndex,0,0)
        !! Det-Diffusion Matrix
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,2,&
                const_diff,omegas,multipliers,sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,1)
    else
        !! Sto-Diffusion Matrices
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,2,&
                const_diff,omegas,multipliers,sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,k)
    end if


!!-----------------------------------------------------------------------------------------------
    do i = 1,npceout
    do j = 1,npceout

        if ((k .eq. ijk(indexi,1)) .and. (i .eq. ijk(indexi,2)) .and. (j .eq. ijk(indexi,3))) then

    !!---------------------------------------->Asii/Asmat
        nni = ((i-1)*nid)
        nnj = ((j-1)*nid)

        do id = 1,nid
        do jd = 1,nid
            if (Adii(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adii(id,jd)
                call MatSetValues(Asii,one,ii,one,jj,imat,ADD_VALUES,ierr)
                call MatSetValues(Asmat,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!---------------------------------------->Asgg
        nni = ((i-1)*nb)
        nnj = ((j-1)*nb)

        do id = 1,nb
        do jd = 1,nb
            if (Adgg(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adgg(id,jd)
                call MatSetValues(Asgg,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

!        !!---------------------------------------->Asig
!            nni = ((i-1)*nid)
!            nnj = ((j-1)*nb)
!
!            do id = 1,nid
!            do jd = 1,nb
!                if (Adig(id,jd) .ne. 0) then
!                    ii = (nni+(id-1))
!                    jj = (nnj+(jd-1))
!                    imat = cijk(indexi)*Adig(id,jd)
!                    call MatSetValues(Asig,one,ii,one,jj,imat,ADD_VALUES,ierr)
!                end if
!            end do
!            end do

    !!---------------------------------------->Asgi
        nni = ((i-1)*nb)
        nnj = ((j-1)*nid)

        do id = 1,nb
        do jd = 1,nid
            if (Adgi(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adgi(id,jd)
                call MatSetValues(Asgi,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

!        !!---------------------------------------->Asir/Asmat
!            nni = ((i-1)*nid)
!            nnj = ((j-1)*nri)
!
!            do id = 1,nid
!            do jd = 1,nri
!                if (Adir(id,jd) .ne. 0) then
!                    ii = (nni+(id-1))
!                    jj = (nnj+(jd-1))
!                    imat = cijk(indexi)*Adir(id,jd)
!                    call MatSetValues(Asir,one,ii,one,jj,imat,ADD_VALUES,ierr)
!
!                    jj2 = nip + jj
!                    call MatSetValues(Asmat,one,ii,one,jj2,imat,ADD_VALUES,ierr)
!
!                end if
!            end do
!            end do

    !!---------------------------------------->Asri/Asmat
        nni = ((i-1)*nri)
        nnj = ((j-1)*nid)

        do id = 1,nri
        do jd = 1,nid
            if (Adri(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adri(id,jd)
                call MatSetValues(Asri,one,ii,one,jj,imat,ADD_VALUES,ierr)

                ii2 = nip + ii
                call MatSetValues(Asmat,one,ii2,one,jj,imat,ADD_VALUES,ierr)
                call MatSetValues(Asmat,one,jj,one,ii2,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

    !!---------------------------------------->Asrr/Asmat
        nni = ((i-1)*nri)
        nnj = ((j-1)*nri)

        do id = 1,nri
        do jd = 1,nri
            if (Adrr(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adrr(id,jd)
                call MatSetValues(Asrr,one,ii,one,jj,imat,ADD_VALUES,ierr)

                ii2 = nip + ii
                jj2 = nip + jj
                call MatSetValues(Asmat,one,ii2,one,jj2,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

!        !!---------------------------------------->Asic
!            nni = ((i-1)*nid)
!            nnj = ((j-1)*nci)
!
!            do id = 1,nid
!            do jd = 1,nci
!                if (Adic(id,jd) .ne. 0) then
!                    ii = (nni+(id-1))
!                    jj = (nnj+(jd-1))
!                    imat = cijk(indexi)*Adic(id,jd)
!                    call MatSetValues(Asic,one,ii,one,jj,imat,ADD_VALUES,ierr)
!                end if
!            end do
!            end do

    !!---------------------------------------->Asci
        nni = ((i-1)*nci)
        nnj = ((j-1)*nid)

        do id = 1,nci
        do jd = 1,nid
            if (Adci(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adci(id,jd)
                call MatSetValues(Asci,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

!        !!---------------------------------------->Asrc
!            nni = ((i-1)*nri)
!            nnj = ((j-1)*nci)
!
!            do id = 1,nri
!            do jd = 1,nci
!                if (Adrc(id,jd) .ne. 0) then
!                    ii = (nni+(id-1))
!                    jj = (nnj+(jd-1))
!                    imat = cijk(indexi)*Adrc(id,jd)
!                    call MatSetValues(Asrc,one,ii,one,jj,imat,ADD_VALUES,ierr)
!                end if
!            end do
!            end do

    !!---------------------------------------->Ascr
        nni = ((i-1)*nci)
        nnj = ((j-1)*nri)

        do id = 1,nci
        do jd = 1,nri
            if (Adcr(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adcr(id,jd)
                call MatSetValues(Ascr,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do


    !!---------------------------------------->Ascc
        nni = ((i-1)*nci)
        nnj = ((j-1)*nci)

        do id = 1,nci
        do jd = 1,nci
            if (Adcc(id,jd) .ne. 0) then
                ii = (nni+(id-1))
                jj = (nnj+(jd-1))
                imat = cijk(indexi)*Adcc(id,jd)
                call MatSetValues(Ascc,one,ii,one,jj,imat,ADD_VALUES,ierr)
            end if
        end do
        end do

        indexi = indexi+1
        end if

    end do
    end do
end do


deallocate(Adii,Adgg,Adrr,Adcc,Adgi,Adri,Adci,Adcr)

call MatAssemblyBegin(Asmat,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asii,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asgg,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asgi,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asrr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asri,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Ascc,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Asci,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyBegin(Ascr,MAT_FINAL_ASSEMBLY,ierr)

call MatAssemblyEnd(Asii,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asgg,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asgi,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asrr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asri,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Ascc,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asci,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Ascr,MAT_FINAL_ASSEMBLY,ierr)
call MatAssemblyEnd(Asmat,MAT_FINAL_ASSEMBLY,ierr)

!!call MatAssemblyBegin(Asig,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asir,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asic,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyBegin(Asrc,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asig,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asir,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asic,MAT_FINAL_ASSEMBLY,ierr)
!!call MatAssemblyEnd(Asrc,MAT_FINAL_ASSEMBLY,ierr)


END SUBROUTINE StoMatSeqOneLevelFullMallocs


!!!*********************************************************
!!! Subroutine: To extract non-zero elements from PETSc Matrix Assembly
SUBROUTINE GetMallocsShort(pid,np,nb,nc,nr,npceout,nip,nbp,ncp,nrp,nirp,nVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


PetscInt         nip, nbp, one, iNZ  !!nzi,nzb,
PetscInt         nirp, nrp, ncp, nVec
PetscInt         id, jd, kd !!ii, jj, nni, nnj
PetscErrorCode   ierr


    PetscInt,ALLOCATABLE      :: idxmi(:),idxni(:),idxmb(:),idxnb(:),idxnir(:)
    PetscInt,ALLOCATABLE      :: idxmbi(:), idxnbi(:), idxmrr(:), idxnrr(:)
    PetscInt,ALLOCATABLE      :: idxmri(:), idxnri(:), idxmci(:), idxnci(:)
    PetscInt,ALLOCATABLE      :: idxmcr(:), idxncr(:), idxmcc(:), idxncc(:)
    PetscInt,ALLOCATABLE      :: idxmi2(:), idxni2(:), idxmi3(:), idxni3(:)

    integer :: k, kk, nid, pid
    integer :: np, nb, nbd, npceout
    integer :: nc, nr, ncd, nrd

    integer, dimension(npceout)       :: nZcijk
    character(len=255)                :: str1, str3

    double precision, allocatable, dimension(:,:)  :: Adii,Adgg,Adrr,Adcc
    double precision, allocatable, dimension(:,:)  :: Adgi,Adri,Adci,Adcr

    REAL*8, PARAMETER :: EPSILON = 1d-14

    !!-----------------------------------------------------------------------------------------------
    nid = (np-nb)*nVec
    ncd = nc*nVec
    nrd = nr*nVec
    nbd = nb*nVec
    one = 1

    !!-----------------------------------------------------------------------------------------------
    ALLOCATE (Adii(nid,nid), Adgg(nbd,nbd), Adrr(nrd,nrd), Adcc(ncd,ncd))
    ALLOCATE (Adgi(nbd,nid), Adri(nrd,nid), Adci(ncd,nid), Adcr(ncd,nrd))
    ALLOCATE (idxmi(nid),idxni(nip),idxmb(nbd),idxnb(nbp),idxmbi(nbd),idxnbi(nbp))
    ALLOCATE (idxmri(nrd),idxnri(nrp),idxmci(ncd),idxnci(ncp),idxmrr(nrd),idxnrr(nrp))
    ALLOCATE (idxmcr(ncd),idxncr(ncp),idxmcc(ncd),idxncc(ncp))
    ALLOCATE (idxmi2(nid), idxni2(nip), idxmi3(nrd), idxni3(nrp),idxnir(nirp))

    !!-----------------------------------------------------------------------------------------------
    !! Precalculated nZijk for each PCE order and Dimension
    open(unit=2,file='../../data/klePceData/nZijk.dat',status='old')
    read(unit=2,fmt='(I8)') nZcijk
    close(2)

    !!-----------------------------------------------------------------------------------------------
    !! non-zeros(ADii)
    kk = 1   !! kk=1, because the non-zero structure is assumed to be same for each PCE mode
    call int2str(str1,pid+1,1) !!call int2str(str2,kk,1)

    !! For Dolfin we use ADii***2.dat: because 1-matrix has different non-zero structure
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADii2.dat'
    open(unit=1,file=str3,status='old')
    read(unit=1,fmt=*) Adii
    close(1)

    !!-----------------------------------------------------------------------------------------------
    !! non-zeros(ADgg)
    call int2str(str1,pid+1,1)   !!call int2str(str2,kk,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgg2.dat'
    open(unit=2,file=str3,status='old')
    read(unit=2,fmt=*) Adgg
    close(2)

    !!-----------------------------------------------------------------------------------------------
    !! non-zeros(ADgi)
    call int2str(str1,pid+1,1) !!call int2str(str2,kk,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADgi2.dat'
    open(unit=3,file=str3,status='old')
    read(unit=3,fmt=*) Adgi
    close(3)

    !!-----------------------------------------------------------------------------------------------
    !! non-zeros(ADri)
    call int2str(str1,pid+1,1) !!call int2str(str2,kk,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADri2.dat'
    open(unit=4,file=str3,status='old')
    read(unit=4,fmt=*) Adri
    close(4)

    !!-----------------------------------------------------------------------------------------------
    !! non-zeros(ADci)
    call int2str(str1,pid+1,1) !!call int2str(str2,kk,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADci2.dat'
    open(unit=5,file=str3,status='old')
    read(unit=5,fmt=*) Adci
    close(5)

    !!-----------------------------------------------------------------------------------------------
    !! non-zeros(ADri)
    call int2str(str1,pid+1,1) !!call int2str(str2,kk,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADrr2.dat'
    open(unit=5,file=str3,status='old')
    read(unit=5,fmt=*) Adrr
    close(5)

    !!-----------------------------------------------------------------------------------------------
    !! non-zeros(ADci)
    call int2str(str1,pid+1,1) !!call int2str(str2,kk,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADcc2.dat'
    open(unit=7,file=str3,status='old')
    read(unit=7,fmt=*) Adcc
    close(7)

    !!-----------------------------------------------------------------------------------------------
    !! non-zeros(ADci)
    call int2str(str1,pid+1,1) !!call int2str(str2,kk,1)
    str3 = '../../external/dolfin/data/Amats/subdom' // trim(str1) // '/ADcr2.dat'
    open(unit=8,file=str3,status='old')
    read(unit=8,fmt=*) Adcr
    close(8)

    !!-----------------------------------------------------------------------------------------------
    !! Here, using first level non-Zero & cijk structure we extrac second level non-Zero structure
    DO k = 1,npceout

        !!-----------------------------------------------------------------------------------------------
        do id = 1,nid
        iNZ = 0
            do jd = 1,nid
                if (ABS (Adii(id,jd) ) > EPSILON) then
                !if (Adii(id,jd) .ne. 0) then
                !if (abs(Adii(id,jd)) .gt. 0.000000001) then
                iNZ = iNZ+1
                end if
            end do
            idxmi(id) = (iNZ)
            !print*, iNZ

            do kd = 1,nrd
                if (ABS (Adri(kd,id) ) > EPSILON) then
                !if (Adri(kd,id) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
            idxmi2(id) = (iNZ)
        !Print*, iNZ
        !print*,'--------**---------'
        end do
        !print*,'-----------------'
        !print*,nid
        !print*,'--------$$-------'
        !print*,idxmi
        idxni(((k-1)*nid+1):(k*nid)) = nZcijk(k)*idxmi
        idxni2(((k-1)*nid+1):(k*nid)) = nZcijk(k)*idxmi2


        do id = 1,nbd
        iNZ = 0
            do jd = 1,nbd
                if (ABS (Adgg(id,jd) ) > EPSILON) then
                !if (Adgg(id,jd) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
            idxmb(id) = (iNZ)
        end do
        idxnb(((k-1)*nbd+1):(k*nbd)) = nZcijk(k)*idxmb


        do id = 1,nbd
        iNZ = 0
            do jd = 1,nid
                if (ABS (Adgi(id,jd) ) > EPSILON) then
                !if (Adgi(id,jd) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
            idxmbi(id) = (iNZ)
        end do
        idxnbi(((k-1)*nbd+1):(k*nbd)) = nZcijk(k)*idxmbi


        do id = 1,nrd
        iNZ = 0
            do jd = 1,nid
                if (ABS (Adri(id,jd) ) > EPSILON) then
                !if (Adri(id,jd) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
            idxmri(id) = (iNZ)
            do kd = 1,nrd
                if (ABS (Adrr(id,kd) ) > EPSILON) then
                !if (Adrr(id,kd) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
        idxmi3(id) = (iNZ)
        end do
        idxnri(((k-1)*nrd+1):(k*nrd)) = nZcijk(k)*idxmri
        idxni3(((k-1)*nrd+1):(k*nrd)) = nZcijk(k)*idxmi3


        do id = 1,ncd
        iNZ = 0
            do jd = 1,nid
                if (ABS (Adci(id,jd) ) > EPSILON) then
                !if (Adci(id,jd) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
        idxmci(id) = (iNZ)
        end do
        idxnci(((k-1)*ncd+1):(k*ncd)) = nZcijk(k)*idxmci


        do id = 1,nrd
        iNZ = 0
            do jd = 1,nrd
                if (ABS (Adrr(id,jd) ) > EPSILON) then
                !if (Adrr(id,jd) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
        idxmrr(id) = (iNZ)
        end do
        idxnrr(((k-1)*nrd+1):(k*nrd)) = nZcijk(k)*idxmrr


        do id = 1,ncd
        iNZ = 0
            do jd = 1,nrd
                if (ABS (Adcr(id,jd) ) > EPSILON) then
                !if (Adcr(id,jd) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
        idxmcr(id) = (iNZ)
        end do
        idxncr(((k-1)*ncd+1):(k*ncd)) = nZcijk(k)*idxmcr


        do id = 1,ncd
        iNZ = 0
            do jd = 1,ncd
                if (ABS (Adcc(id,jd) ) > EPSILON) then
                !if (Adcc(id,jd) .ne. 0) then
                iNZ = iNZ+1
                end if
            end do
        idxmcc(id) = (iNZ)
        end do
        idxncc(((k-1)*ncd+1):(k*ncd)) = nZcijk(k)*idxmcc

    END DO

    !!-----------------------------------------------------------------------------------------------
    !! Write mallocs files for respective matrices  !! Symmetry structure exploited
    idxnir(1:nip) = idxni2
    idxnir(nip+1:nirp) = idxni3

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzir' // trim(str1) // '.dat'
    open(unit=9,file=str1,status='replace')
    write(9,*) idxnir
    close(9)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzi' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='replace')
    write(1,*) idxni
    close(1)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzb' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='replace')
    write(2,*) idxnb
    close(2)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzbi' // trim(str1) // '.dat'
    open(unit=3,file=str1,status='replace')
    write(3,*) idxnbi
    close(3)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzri' // trim(str1) // '.dat'
    open(unit=4,file=str1,status='replace')
    write(4,*) idxnri
    close(4)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzci' // trim(str1) // '.dat'
    open(unit=5,file=str1,status='replace')
    write(5,*) idxnci
    close(5)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzrr' // trim(str1) // '.dat'
    open(unit=66,file=str1,status='replace')
    write(66,*) idxnrr
    close(66)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzcr' // trim(str1) // '.dat'
    open(unit=7,file=str1,status='replace')
    write(7,*) idxncr
    close(7)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzcc' // trim(str1) // '.dat'
    open(unit=8,file=str1,status='replace')
    write(8,*) idxncc
    close(8)

    DEALLOCATE(Adii,Adgg,Adrr,Adcc,Adgi,Adri,Adci,Adcr)
    DEALLOCATE(idxnb,idxmb,idxni,idxmi,idxnbi,idxmbi,idxmci,idxnci)
    DEALLOCATE(idxmri,idxnri,idxmrr,idxnrr,idxmcr,idxncr,idxmcc,idxncc)
    DEALLOCATE(idxmi2,idxni2,idxmi3,idxni3,idxnir)

END SUBROUTINE GetMallocsShort



!!!*********************************************************
!!! Subroutine: To non-zero elements from PETSc Matrix Assembly
SUBROUTINE GetMallocs(pid,p,e,t,np,ne,nt,nb,nci,nri,ndim,npcein,npceout,nomga,ncijk,ijk,&
                      cijk,nid,nip,nbp,ncp,nrp,nirp, dbounds, const_diff, omegas, &
                      mIndex,sIndex,multipliers,sigma,casep,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    PetscInt         nip,nbp,one
    PetscInt         nirp, nrp, ncp
    PetscInt         id, jd, kd
    PetscErrorCode   ierr

    integer                   :: iNZ,pid
    character(len=255)        :: str1
    integer,dimension(npceout):: nZcijk

    PetscInt,ALLOCATABLE      :: idxmi(:),idxni(:),idxmb(:),idxnb(:),idxnir(:)
    PetscInt,ALLOCATABLE      :: idxmbi(:), idxnbi(:), idxmrr(:), idxnrr(:)
    PetscInt,ALLOCATABLE      :: idxmri(:), idxnri(:), idxmci(:), idxnci(:)
    PetscInt,ALLOCATABLE      :: idxmcr(:), idxncr(:), idxmcc(:), idxncc(:)
    PetscInt,ALLOCATABLE      :: idxmi2(:), idxni2(:), idxmi3(:), idxni3(:)

    integer :: k, indexi, ncijk, ndim, npceout, npcein
    integer :: np, ne, nt, nb, nid, nci, nri, casep, nomga
    double precision :: const_diff, sigma

    integer, dimension(ncijk,3)       :: ijk
    integer, dimension(3,ne)          :: e
    integer, dimension(3,nt)          :: t
    double precision, dimension(2,np) :: p
    double precision, dimension(ncijk):: cijk
    double precision, dimension(2,2)  :: dbounds
    double precision, dimension(nomga):: omegas, multipliers
    integer, dimension(ndim,npceout)  :: mIndex
    integer, dimension(2,ndim)        :: sIndex

    double precision, allocatable, dimension(:,:)  :: Adii,Adgg,Adrr,Adcc
    double precision, allocatable, dimension(:,:)  :: Adgi,Adri,Adci,Adcr

!!-----------------------------------------------------------------------------------------------
    ALLOCATE (Adii(nid,nid), Adgg(nb,nb), Adrr(nri,nri), Adcc(nci,nci),&
              Adgi(nb,nid), Adri(nri,nid), Adci(nci,nid), Adcr(nci,nri))

    ALLOCATE (idxmi(nid),idxni(nip),idxmb(nb),idxnb(nbp),idxmbi(nb),idxnbi(nbp))
    ALLOCATE (idxmri(nri),idxnri(nrp),idxmci(nci),idxnci(ncp),idxmrr(nri),idxnrr(nrp))
    ALLOCATE (idxmcr(nci),idxncr(ncp),idxmcc(nci),idxncc(ncp))
    ALLOCATE (idxmi2(nid), idxni2(nip), idxmi3(nri), idxni3(nrp),idxnir(nirp))

    open(unit=2,file='../../data/klePceData/nZijk.dat',status='old')
    read(unit=2,fmt='(I8)') nZcijk
    close(2)
    !!print*,nZcijk


!!-----------------------------------------------------------------------------------------------
!! Method : 3
one    = 1
indexi = 1

do k = 1,npceout

    !! Deterministic Matices
    Adii(:,:) = 0.0d0
    Adgg(:,:) = 0.0d0
    Adgi(:,:) = 0.0d0
    Adri(:,:) = 0.0d0
    Adci(:,:) = 0.0d0
    Adcr(:,:) = 0.0d0
    Adcc(:,:) = 0.0d0
    Adrr(:,:) = 0.0d0
    !!Adig(:,:) = 0.0d0
    !!Adir(:,:) = 0.0d0
    !!Adic(:,:) = 0.0d0
    !!Adrc(:,:) = 0.0d0

    !! UnderStudy: July 7, 2016

    if (k .eq. 1) then
        !! Det-Advection Matrix
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,1,&
        0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
        dbounds,mIndex,sIndex,0,0)
        !! Det-Diffusion Matrix
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,2,&
        const_diff,omegas,multipliers,sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,1)
    else
        !! Sto-Diffusion Matrices
        call SubAssembleMatrix(Adii,Adgg,Adgi,Adri,Adci,Adcr,Adcc,Adrr,p,e,t,np,ne,nt,nb,nci,ndim,npceout,nomga,2,&
        const_diff,omegas,multipliers,sigma,[0.0d0,0.0d0,0.0d0,0.0d0],dbounds,mIndex,sIndex,casep,k)
    end if

    !!-----------------------------------------------------------------------------------------------
    do id = 1,nid
        iNZ = 0
        do jd = 1,nid
        if (Adii(id,jd) .ne. 0) then
        !if (abs(Adii(id,jd)) .gt. 0.000000001) then
        iNZ = iNZ+1
        end if
        end do
        idxmi(id) = (iNZ)
        !print*, iNZ

        do kd = 1,nri
        if (Adri(kd,id) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmi2(id) = (iNZ)
        !Print*, iNZ
        !print*,'--------**---------'
    end do
    !print*,'-----------------'
    !print*,nid
    !print*,'--------$$-------'
    !print*,idxmi
    idxni(((k-1)*nid+1):(k*nid)) = nZcijk(k)*idxmi
    idxni2(((k-1)*nid+1):(k*nid)) = nZcijk(k)*idxmi2


    do id = 1,nb
        iNZ = 0
        do jd = 1,nb
        if (Adgg(id,jd) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmb(id) = (iNZ)
    end do
    idxnb(((k-1)*nb+1):(k*nb)) = nZcijk(k)*idxmb


    do id = 1,nb
        iNZ = 0
        do jd = 1,nid
        if (Adgi(id,jd) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmbi(id) = (iNZ)
    end do
    idxnbi(((k-1)*nb+1):(k*nb)) = nZcijk(k)*idxmbi


    do id = 1,nri
    iNZ = 0
    do jd = 1,nid
        if (Adri(id,jd) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmri(id) = (iNZ)
        do kd = 1,nri
        if (Adrr(id,kd) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmi3(id) = (iNZ)
    end do
    idxnri(((k-1)*nri+1):(k*nri)) = nZcijk(k)*idxmri
    idxni3(((k-1)*nri+1):(k*nri)) = nZcijk(k)*idxmi3


    do id = 1,nci
        iNZ = 0
        do jd = 1,nid
        if (Adci(id,jd) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmci(id) = (iNZ)
    end do
    idxnci(((k-1)*nci+1):(k*nci)) = nZcijk(k)*idxmci


    do id = 1,nri
        iNZ = 0
        do jd = 1,nri
        if (Adrr(id,jd) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmrr(id) = (iNZ)
    end do
    idxnrr(((k-1)*nri+1):(k*nri)) = nZcijk(k)*idxmrr


    do id = 1,nci
        iNZ = 0
        do jd = 1,nri
        if (Adcr(id,jd) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmcr(id) = (iNZ)
    end do
    idxncr(((k-1)*nci+1):(k*nci)) = nZcijk(k)*idxmcr


    do id = 1,nci
        iNZ = 0
        do jd = 1,nci
        if (Adcc(id,jd) .ne. 0) then
        iNZ = iNZ+1
        end if
        end do
        idxmcc(id) = (iNZ)
    end do
    idxncc(((k-1)*nci+1):(k*nci)) = nZcijk(k)*idxmcc

end do
!!-----------------------------------------------------------------------------------------------

    idxnir(1:nip) = idxni2
    idxnir(nip+1:nirp) = idxni3

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzir' // trim(str1) // '.dat'
    open(unit=9,file=str1,status='replace')
    write(9,*) idxnir
    close(9)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzi' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='replace')
    write(1,*) idxni
    close(1)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzb' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='replace')
    write(2,*) idxnb
    close(2)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzbi' // trim(str1) // '.dat'
    open(unit=3,file=str1,status='replace')
    write(3,*) idxnbi
    close(3)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzri' // trim(str1) // '.dat'
    open(unit=4,file=str1,status='replace')
    write(4,*) idxnri
    close(4)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzci' // trim(str1) // '.dat'
    open(unit=5,file=str1,status='replace')
    write(5,*) idxnci
    close(5)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzrr' // trim(str1) // '.dat'
    open(unit=66,file=str1,status='replace')
    write(66,*) idxnrr
    close(66)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzcr' // trim(str1) // '.dat'
    open(unit=7,file=str1,status='replace')
    write(7,*) idxncr
    close(7)

    call int2str(str1,pid+1,4)
    str1 = '../../data/mallocData/nnzcc' // trim(str1) // '.dat'
    open(unit=8,file=str1,status='replace')
    write(8,*) idxncc
    close(8)


    DEALLOCATE(Adii,Adgg,Adrr,Adcc,Adgi,Adri,Adci,Adcr)
    DEALLOCATE(idxnb,idxmb,idxni,idxmi,idxnbi,idxmbi,idxmci,idxnci)
    DEALLOCATE(idxmri,idxnri,idxmrr,idxnrr,idxmcr,idxncr,idxmcc,idxncc)
    DEALLOCATE(idxmi2,idxni2,idxmi3,idxni3,idxnir)


END SUBROUTINE GetMallocs
!!!!*********************************************************


!!!*********************************************************
!!! Subroutine: To construc New-PETSc Vector in Sequential Form
SUBROUTINE GetVecSeq(n,fvec,PetVec,Solvec,temp1,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec,SolVec
    PetscInt         i,ii
    PetscInt         n,one
    PetscScalar      ivec
    PetscErrorCode   ierr

    double precision, dimension(n) :: fvec
    integer, dimension(n)          :: temp1

    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)

    call VecDuplicate(PetVec,SolVec,ierr)

    one = 1
    do i = 1,n
        ii = i-1
        temp1(i) = i-1
        ivec = fvec(i)
        call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
    end do

    call VecAssemblyBegin(PetVec,ierr)
    call VecAssemblyEnd(PetVec,ierr)


END SUBROUTINE GetVecSeq



!!!*********************************************************
!!! Subroutine: To construc allready created Sequential PETSc Vector
SUBROUTINE SetVecSeq(n,fvec,PetVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec
    PetscInt         i,ii
    PetscInt         n,one
    PetscScalar      ivec
    PetscErrorCode   ierr

    double precision, dimension(n) :: fvec
    !!integer, dimension(n)          :: temp1
    !!call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)
    !!call VecDuplicate(PetVec,SolVec,ierr)

    one = 1
    do i = 1,n
        ii = i-1
        !!temp1(i) = i-1
        ivec = fvec(i)
        call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
    end do

    call VecAssemblyBegin(PetVec,ierr)
    call VecAssemblyEnd(PetVec,ierr)

END SUBROUTINE SetVecSeq


!!!*********************************************************
!!! Subroutine: To construc New PETSc Matrix in Sequential Form
SUBROUTINE GetMatSeq(n,m,Amat,PetMat,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              PetMat
    PetscInt         i,j,ii,jj
    PetscInt         n,m,nz,one
    PetscScalar      imat
    PetscErrorCode   ierr

    double precision, dimension(n,m) :: Amat

    nz = 12  !! number of non-zeros per row for matrix pre-allocation
    call MatCreateSeqAIJ(PETSC_COMM_SELF,n,m,nz,PETSC_NULL_INTEGER,PetMat,ierr)

    one = 1
    do i = 1,n
    do j = 1,m
        if (Amat(i,j) .ne. 0) then
            ii = i-1
            jj = j-1
            imat = Amat(i,j)
            call MatSetValues(PetMat,one,ii,one,jj,imat,INSERT_VALUES,ierr)
        end if
    end do
    end do

    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)


END SUBROUTINE GetMatSeq

!!!*********************************************************
!!! Subroutine: To construc empty sequential PETSc Vector
SUBROUTINE GetVecSeqDummy(n,PetVec,temp1,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec
    PetscInt         i, n
    PetscErrorCode   ierr

    integer, dimension(n) :: temp1

    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)

    do i = 1,n
       temp1(i) = i-1
    end do


END SUBROUTINE GetVecSeqDummy


!!!*********************************************************
!!! Subroutine: To construc one-empty sequential PETSc Vector
SUBROUTINE GetVecSeqTemp1(n,PetVec,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec
    PetscInt         n
    PetscErrorCode   ierr

    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)


END SUBROUTINE GetVecSeqTemp1


!!!*********************************************************
!!! Subroutine: To construc two-empty sequential PETSc Vector
!! with it's position indieces defined in array temp
SUBROUTINE GetVecSeqTemp2(n,PetVec,Solvec,temp1,ierr)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

    Vec              PetVec,SolVec
    PetscInt         i, n
    PetscErrorCode   ierr

    integer, dimension(n) :: temp1

    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)
    call VecDuplicate(PetVec,SolVec,ierr)

    do i = 1,n
        temp1(i) = i-1
    end do

END SUBROUTINE GetVecSeqTemp2


!!!!*********************************************************
subroutine GetRs(pid,nb,nbg,npceout,nVec,RsMat)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>


    Mat              RsMat
    PetscScalar      imat
    PetscErrorCode   ierr
    PetscInt         one, nbp, nbgp !, ii, jj
    PetscInt         ii1, jj1, ii2, jj2, ii3, jj3


    character(len=255) :: str1
    integer :: nb, nbg, npceout, i, j, pid, temp1, nVec, nbd, nbgd
    !!integer, dimension(nb) :: temp2

    one  = 1
    nbd  = nb*nVec
    nbgd = nbg*nVec
    nbp  = nbd*npceout
    nbgp = nbgd*npceout

    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbgp, one, PETSC_NULL_INTEGER, RsMat,ierr)

    call int2str(str1,pid+1,4)
    str1 = '../../data/meshData/bnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')

    imat = 1.0
    do i = 1,nb
        !!temp1 = temp2(i)
        read(unit=1,fmt=*) temp1

        do j = 1,npceout

            ii1 = (i+(j-1)*nbd)-1
            jj1 = (temp1+(j-1)*nbgd)-1

            ii2 = nb + (i+(j-1)*nbd)-1
            jj2 = nbg + (temp1+(j-1)*nbgd)-1

            ii3 = 2*nb + (i+(j-1)*nbd)-1
            jj3 = 2*nbg + (temp1+(j-1)*nbgd)-1

            call MatSetValues(RsMat,one,ii1,one,jj1,imat,INSERT_VALUES,ierr)
            call MatSetValues(RsMat,one,ii2,one,jj2,imat,INSERT_VALUES,ierr)
            call MatSetValues(RsMat,one,ii3,one,jj3,imat,INSERT_VALUES,ierr)

        end do
    end do

    close(1)

    call MatAssemblyBegin(RsMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(RsMat,MAT_FINAL_ASSEMBLY,ierr)

end subroutine GetRs


!!**************************************************************************
!! Subroutine to constrct Block-diagonal-weighting matrix 'D_s'
!! The diagonal entries of the Det-Weighting-Mat are equal to the reciprocal
!! of the number of subdomains which contain the same interface node
SUBROUTINE GetDs(pid,nb,npceout,nVec,DsMat)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    Mat              DsMat
    PetscScalar      imat
    PetscInt         one, nbp, ii1, ii2, ii3
    PetscErrorCode   ierr


    CHARACTER(len=255) :: extension, str1

    INTEGER :: npg,neg,ntg,nbg,ndom,np,ne,nt,nb,npceout
    INTEGER :: i, j, k, Ri, Rj, Ri2, Rj2, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: GlobalBN, bnodes12, bnodes13

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CountR
    INTEGER :: nh, nhg, nbd, nVec

    one = 1
    nbd = nb*nVec
    nbp = nbd*npceout

    open(unit=1,file='../../data/meshData/meshdim.dat',status='old')
    read(unit=1,fmt='(6I8)') npg, neg, ntg, nhg, nbg, ndom
    close(1)


    ALLOCATE(GlobalBN(nbg),CountR(nbg))
    open(unit=1,file='../../data/meshData/boundary_nodes.dat',status='old')
    READ(1,*) GlobalBN
    close(1)

    CountR(:) = 0.0d0
    DO k = 1,ndom
        CALL int2str(extension,k,4)               ! there is scope for optimization here
        extension = trim(extension)
        !CALL readmeshdim(np,ne,nt,nb,extension)
        CALL readmeshdim3D(np,ne,nt,nh,nb,extension)
        ALLOCATE(bnodes12(nb))

        CALL int2str(str1,k,4)
        str1 = '../../data/meshData/nbnodes' // trim(str1) // '.dat'  ! changing from nbnodes to nodes wont work
        OPEN(unit=1,file=str1,status='old')       ! Need to confirm this change
        READ(1,*) bnodes12                        ! Date Feb 1, 2014
        CLOSE(1)

        DO i = 1,nbg
           Ri = GlobalBN(i)
           DO j = 1,nb
              Rj = bnodes12(j)
              IF (Ri == Rj) THEN
                 CountR(i) = CountR(i) + 1
              END IF
           END DO
        END DO

        DEALLOCATE(bnodes12)
    END DO

    CALL int2str(extension,pid+1,4)
    extension = trim(extension)
    !CALL readmeshdim(np,ne,nt,nb,extension)
    CALL readmeshdim3D(np,ne,nt,nh,nb,extension)
    ALLOCATE(bnodes13(nb))

    CALL int2str(str1,pid+1,4)
    str1 = '../../data/meshData/nbnodes' // trim(str1) // '.dat'
    OPEN(unit=1,file=str1,status='old')
    READ(1,*) bnodes13
    CLOSE(1)

    call MatCreateSeqAIJ(PETSC_COMM_SELF, nbp, nbp, one, PETSC_NULL_INTEGER, DsMat,ierr)

    Do i = 1,nb
       Ri2 = bnodes13(i)
       DO j = 1,nbg
          Rj2 = GlobalBN(j)
          IF (Ri2 == Rj2) THEN

             !!Dmat(i,i) = 1/CountR(j)
             imat = 1/CountR(j)
             do k = 1,npceout

                ii1 = (i+(k-1)*nbd)-1
                ii2 = nb + (i+(k-1)*nbd)-1
                ii3 = 2*nb + (i+(k-1)*nbd)-1

                call MatSetValues(DsMat,one,ii1,one,ii1,imat,INSERT_VALUES,ierr)
                call MatSetValues(DsMat,one,ii2,one,ii2,imat,INSERT_VALUES,ierr)
                call MatSetValues(DsMat,one,ii3,one,ii3,imat,INSERT_VALUES,ierr)

             end do

          END IF
       END DO
    END DO

    call MatAssemblyBegin(DsMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(DsMat,MAT_FINAL_ASSEMBLY,ierr)

    DEALLOCATE(bnodes13, GlobalBN, CountR)

END SUBROUTINE GetDs


!!****************************************************************
!! Subroutine to constrct Block-diagonal-restriction matrix 'Bc_s'
!! restriction operator that maps the global corner unknown 'Uc'
!! into its local 'Uc_s' components
SUBROUTINE GetBc(pid, nc, nbgc, npceout, nVec, BcMat)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    Mat              BcMat
    PetscScalar      imat
    PetscInt         one, ncp, nbgcp, ncd, nbgcd !, ii, jj
    PetscInt         ii1, jj1, ii2, jj2, ii3, jj3
    PetscErrorCode   ierr


    CHARACTER(len=255) :: str1 !!,extension, filename

    INTEGER :: nc,nbgc,npceout, nVec
    INTEGER :: i, j, k, Ri, Rj, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: cnode12, GlobalCN

    one   = 1
    ncd   = nc*nVec
    nbgcd = nbgc*nVec
    ncp   = ncd*npceout
    nbgcp = nbgcd*npceout

    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nbgcp, one, PETSC_NULL_INTEGER, BcMat,ierr)

    ALLOCATE(GlobalCN(nbgc))
    open(unit=2,file='../../data/meshData/corner_nodes.dat',status='old')
    read(2,*) GlobalCN
    close(2)

    ALLOCATE(cnode12(nc))
    call int2str(str1,pid+1,4)
    str1 = '../../data/meshData/cnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    READ(2,*) cnode12
    close(2)

    imat = 1.0
    DO i = 1,nc
       Ri = cnode12(i)
       DO j = 1,nbgc
          Rj = GlobalCN(j)
          IF (Ri == Rj) THEN

            !!BcSmat(i,j) = 1
            do k = 1,npceout

                ii1 = (i+(k-1)*ncd)-1
                jj1 = (j+(k-1)*nbgcd)-1

                ii2 = nc + (i+(k-1)*ncd)-1
                jj2 = nbgc + (j+(k-1)*nbgcd)-1

                ii3 = 2*nc + (i+(k-1)*ncd)-1
                jj3 = 2*nbgc + (j+(k-1)*nbgcd)-1

                call MatSetValues(BcMat,one,ii1,one,jj3,imat,INSERT_VALUES,ierr)
                call MatSetValues(BcMat,one,ii2,one,jj2,imat,INSERT_VALUES,ierr)
                call MatSetValues(BcMat,one,ii3,one,jj1,imat,INSERT_VALUES,ierr)
            end do

          END IF
       END DO
    END DO

    call MatAssemblyBegin(BcMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(BcMat,MAT_FINAL_ASSEMBLY,ierr)

    DEALLOCATE(cnode12, GlobalCN)

END SUBROUTINE GetBc


!!*****************************************************************
!! Subroutine to constrct Boolean operators 'Rr_s' that extract the
!! local interface unknowns 'Ub_s' into corner 'Ur_s'
SUBROUTINE getRr(pid, nb, nr, npceout, nVec, RrMat)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    Mat              RrMat
    PetscScalar      imat
    PetscInt         one, nrp, nbp, nrd, nbd  !, ii, jj
    PetscInt         ii1,jj1,ii2,jj2,ii3,jj3
    PetscErrorCode   ierr


    CHARACTER(len=255) :: str1 !!,extension, filename,

    INTEGER :: nb,nr,npceout, nVec
    INTEGER :: i, j, k, Ri, Rj, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: rnode12, bnodes12

    ALLOCATE(rnode12(nr), bnodes12(nb))

    one  = 1
    nrd = nr*nVec
    nbd = nb*nVec
    nbp = nbd*npceout
    nrp = nrd*npceout

    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrp, nbp, one, PETSC_NULL_INTEGER, RrMat,ierr)

    call int2str(str1,pid+1,4)
    str1 = '../../data/meshData/rnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    READ(1,*) rnode12
    close(1)

    CALL int2str(str1,pid+1,4)
    str1 = '../../data/meshData/nbnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    READ(2,*) bnodes12
    close(2)

    imat = 1.0
    Do i = 1,nr
       Ri = rnode12(i)
       DO j = 1,nb
          Rj = bnodes12(j)
          IF (Ri == Rj) THEN

            !!RrSmat(i,j) = 1
            do k = 1,npceout

                ii1 = (i+(k-1)*nrd)-1
                jj1 = (j+(k-1)*nbd)-1

                ii2 = nr + (i+(k-1)*nrd)-1
                jj2 = nb + (j+(k-1)*nbd)-1

                ii3 = 2*nr + (i+(k-1)*nrd)-1
                jj3 = 2*nb + (j+(k-1)*nbd)-1

                call MatSetValues(RrMat,one,ii1,one,jj1,imat,INSERT_VALUES,ierr)
                call MatSetValues(RrMat,one,ii2,one,jj2,imat,INSERT_VALUES,ierr)
                call MatSetValues(RrMat,one,ii3,one,jj3,imat,INSERT_VALUES,ierr)

            end do

          END IF
       END DO
    END DO

    call MatAssemblyBegin(RrMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(RrMat,MAT_FINAL_ASSEMBLY,ierr)

    DEALLOCATE(rnode12, bnodes12)

END SUBROUTINE getRr

!!*****************************************************************
!! Subroutine to constrct Boolean operators 'Rc_S' that extract the
!! local corner 'Uc_s' from local interface unknowns 'Ub_s'
SUBROUTINE getRc(pid, nb, nc, npceout, nVec, RcMat)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>

    Mat              RcMat
    PetscScalar      imat
    PetscInt         one, ncp, nbp, ncd, nbd  !!, ii, jj
    PetscInt         ii1,jj1,ii2,jj2,ii3,jj3
    PetscErrorCode   ierr

    CHARACTER(len=255) :: str1 !,extension, filename,

    INTEGER :: nb,nc,npceout, nVec
    INTEGER :: i, j, k, Ri, Rj, pid
    INTEGER, DIMENSION(:), ALLOCATABLE :: cnode12, bnodes12

    ALLOCATE(cnode12(nc), bnodes12(nb))

    one  = 1
    nbd = nb*nVec
    ncd = nc*nVec
    nbp = nbd*npceout
    ncp = ncd*npceout

    call MatCreateSeqAIJ(PETSC_COMM_SELF, ncp, nbp, one, PETSC_NULL_INTEGER, RcMat,ierr)

    call int2str(str1,pid+1,4)
    str1 = '../../data/meshData/cnodes' // trim(str1) // '.dat'
    open(unit=1,file=str1,status='old')
    READ(1,*) cnode12
    close(1)

    CALL int2str(str1,pid+1,4)
    str1 = '../../data/meshData/nbnodes' // trim(str1) // '.dat'
    open(unit=2,file=str1,status='old')
    READ(2,*) bnodes12
    close(2)

    imat = 1.0
    Do i = 1,nc
        Ri = cnode12(i)
        DO j = 1,nb
        Rj = bnodes12(j)
            IF (Ri == Rj) THEN

            !! RcSmat(i,j) = 1
            do k = 1,npceout

                ii1 = (i+(k-1)*ncd)-1
                jj1 = (j+(k-1)*nbd)-1

                ii2 = nc + (i+(k-1)*ncd)-1
                jj2 = nb + (j+(k-1)*nbd)-1

                ii3 = 2*nc + (i+(k-1)*ncd)-1
                jj3 = 2*nb + (j+(k-1)*nbd)-1

                call MatSetValues(RcMat,one,ii1,one,jj1,imat,INSERT_VALUES,ierr)
                call MatSetValues(RcMat,one,ii2,one,jj2,imat,INSERT_VALUES,ierr)
                call MatSetValues(RcMat,one,ii3,one,jj3,imat,INSERT_VALUES,ierr)
            end do

            END IF
        END DO
    END DO

    call MatAssemblyBegin(RcMat,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(RcMat,MAT_FINAL_ASSEMBLY,ierr)


    DEALLOCATE(cnode12, bnodes12)

END SUBROUTINE getRc



SUBROUTINE NewmarkBetaUpdate(pid, ni, nb, deltaT, beta_NB, gamma_NB, Ui, Ub, u0i, u0b, v0i, v0b, a0i, a0b, ierr)


#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>


    PetscScalar      PosOne
    PetscErrorCode   ierr


    integer :: ni, nb, pid
    double precision :: deltaT, beta_NB, gamma_NB, ccc


    !!!! Previous step displacement, velocity and acceleration vectors decomposed

    double precision, dimension(ni) :: u0i, v0i, a0i
    double precision, dimension(nb) :: u0b, v0b , a0b

!!!! New displacement vector

    double precision, dimension(nb) :: Ub, Ab, Vb, aaa, bbb
    double precision, dimension(ni) :: Ui, Ai, Vi



! anew_interface = ((unew_sub_g - un_sub_g - deltaT*pn_sub_g)/(0.5*deltaT^2) - ((1-(2*beta))*an_sub_g))/(2*beta);

! pnew_interface = pn_sub_g + deltaT * ((1-gamma)*an_sub_g + gamma*anew_interface);



! anew_interior = ((unew_interior - un_interior - deltaT*pn_interior)/(0.5*deltaT^2) - ((1-(2*beta))*an_interior))/(2*beta);

! pnew_interior = pn_interior + deltaT * ((1-gamma)*an_interior + gamma*anew_interior);



    Ab = ((Ub - u0b - deltaT*v0b)/(0.5*deltaT**2) -  ((1-2*beta_NB)*a0b))/(2*beta_NB)

    ! bbb = Ub - u0b

    ! ccc = (0.5*deltaT**2)

    ! aaa = (Ub - u0b - deltaT*v0b)/(0.5*deltaT**2)

    Ai = ((Ui - u0i - deltaT*v0i)/(0.5*deltaT**2) -  ((1-2*beta_NB)*a0i))/(2*beta_NB)

    Vb = v0b + deltaT * ((1-gamma_NB)*a0b + gamma_NB * Ab)

    Vi = v0i + deltaT * ((1-gamma_NB)*a0i + gamma_NB * Ai)


! if (pid .eq. 0) then

    ! print*, 'process id', pid
    ! print*, 'Ui is', Ui

    ! print*, 'process id', pid
    ! print*, 'u0i is', u0i

    ! print*,'bbb is', bbb

    ! print*, 'ccc is', ccc

    ! print*, 'process id', pid
    ! print*, 'aaa is', aaa

!     print*, 'process id', pid
!     print*, 'Ai is', Ai

!     stop 123
! end if

    u0i = Ui
    u0b = Ub

    v0i = Vi
    v0b = Vb

    a0i = Ai
    a0b = Ab


END SUBROUTINE NewmarkBetaUpdate


!!!!*********************************************************
!!!! Subroutine: To call PETSc Assembly Vector in most general way
!SUBROUTINE GetPetVec(n,PetVec,SolVec,pg,eg,tg,npg,neg,ntg,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec,SolVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: npg, neg, ntg
!    double precision :: amp
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!    call VecCreate(PETSC_COMM_WORLD,PetVec,ierr)
!    call VecSetFromOptions(PetVec,ierr)
!    call VecSetSizes(PetVec,PETSC_DECIDE,n,ierr)
!    call VecDuplicate(PetVec,SolVec,ierr)
!
!    call PETScAssemblevector(PetVec, SolVec, pg, eg, tg, npg, neg, ntg, amp, dbounds)
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetPetVec
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Vector in Sequential Form
!SUBROUTINE GetPetVecSeq(n,fvec,PetVec,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    double precision, dimension(n) :: fvec
!    !integer, dimension(n)          :: temp1
!
!    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)
!    !call VecDuplicate(PetVec,SolVec,ierr)
!
!    one = 1
!    do i = 1,n
!        ii = i-1
!!        temp1(i) = i-1
!        ivec = fvec(i)
!        call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
!    end do
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetPetVecSeq
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!SUBROUTINE GetPetMat(n,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!    Mat              PetMat
!    PetscInt         i,j,ii,jj
!    PetscInt         n,one
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    integer :: npg, neg, ntg
!    double precision :: amp, meanc
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!    call MatCreate(PETSC_COMM_WORLD,PetMat,ierr)
!    call MatSetFromOptions(PetMat,ierr)
!    call MatSetSizes(PetMat,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!    call MatSetUp(PetMat,ierr)
!
!    call PETScAssembleADVmatrix(PetMat,pg,eg,tg,npg,neg,ntg,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!                               [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call PETScAssembleDIFmatrix(PetMat,pg,eg,tg,npg,neg,ntg,2,meanc,[0.0d0,0.0d0,0.0d0,0.0d0],&
!                               [0.0d0,0.0d0,0.0d0,0.0d0], 0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetPetMat
!
!
!!!!*********************************************************
!!!! Subroutine: To call PETSc Assembly Vector in most general way
!SUBROUTINE GetPetVecSelf(n,PetVec,SolVec,pg,eg,tg,npg,neg,ntg,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec,SolVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: npg, neg, ntg
!    double precision :: amp
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!    call VecCreateSeq(PETSC_COMM_SELF, n, PetVec, ierr)
!    call VecDuplicate(PetVec,SolVec,ierr)
!
!    call PETScAssemblevector(PetVec, SolVec, pg, eg, tg, npg, neg, ntg, amp, dbounds)
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetPetVecSelf
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!SUBROUTINE GetPetMatSelf(n,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!    Mat              PetMat
!    PetscInt         i,j,ii,jj
!    PetscInt         n,one
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    integer :: npg, neg, ntg
!    double precision :: amp, meanc
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, n,PETSC_NULL_INTEGER, PetMat, ierr)
!
!    call PETScAssembleADVmatrix(PetMat,pg,eg,tg,npg,neg,ntg,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!    [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call PETScAssembleDIFmatrix(PetMat,pg,eg,tg,npg,neg,ntg,2,meanc,[0.0d0,0.0d0,0.0d0,0.0d0],&
!    [0.0d0,0.0d0,0.0d0,0.0d0], 0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetPetMatSelf
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Vector in most general way
!SUBROUTINE GetVec(n,fvec,PetVec,Solvec,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec,SolVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    double precision, dimension(n) :: fvec
!
!    call VecCreate(PETSC_COMM_WORLD,PetVec,ierr)
!    call VecSetFromOptions(PetVec,ierr)
!    call VecSetSizes(PetVec,PETSC_DECIDE,n,ierr)
!    call VecDuplicate(PetVec,SolVec,ierr)
!
!    one = 1
!    do i = 1,n
!    ii = i-1
!    ivec = fvec(i)
!    call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
!    end do
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetVec
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!SUBROUTINE GetMat(n,AMAT,PetMat,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!    Mat              PetMat
!    PetscInt         i,j,ii,jj
!    PetscInt         n,one
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    double precision, dimension(n,n) :: Amat
!
!    call MatCreate(PETSC_COMM_WORLD,PetMat,ierr)
!    call MatSetFromOptions(PetMat,ierr)
!    call MatSetSizes(PetMat,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!    call MatSetUp(PetMat,ierr)
!
!    one = 1
!    do i = 1,n
!    do j = 1,n
!    if (Amat(i,j) .ne. 0) then
!    ii = i-1
!    jj = j-1
!    imat = Amat(i,j)
!    call MatSetValues(PetMat,one,ii,one,jj,imat,INSERT_VALUES,ierr)
!    end if
!    end do
!    end do
!
!    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetMat



!!!!*********************************************************
!SUBROUTINE StoVecSeqOneLevel(Fsi,Fsg,SolVeci,SolVecb,tempi,tempb,p,e,t,np,ne,nt,nb,npceout,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              Fsi,Fsg,SolVeci,SolVecb
!    PetscInt         i,ii
!    PetscInt         one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb, npceout, nip, nbp, id
!    double precision                  :: amp
!    integer, dimension(3,ne)          :: e
!    integer, dimension(3,nt)          :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2)  :: dbounds
!
!    integer, dimension(np-nb*npceout) :: tempi
!    integer, dimension(nb*npceout)    :: tempb
!
!    double precision, dimension((np-nb)*npceout)  :: Fdi
!    double precision, dimension(nb*npceout)       :: Fdg
!
!    nip = ((np-nb)*npceout)
!    nbp = nb*npceout
!
!    call VecCreateSeq(PETSC_COMM_SELF, nip, Fsi, ierr)
!    call VecCreateSeq(PETSC_COMM_SELF, nbp, Fsg, ierr)
!
!    call VecDuplicate(Fsi,SolVeci,ierr)
!    call VecDuplicate(Fsg,SolVecb,ierr)
!
!    !call PETScOneLevelAssemVec(Fdi, Fdg, p, e, t, np, ne, nt, nb, amp, dbounds,tempi,tempb)
!    call SubAssembleVector(Fdi, Fdg, p, e, t, np, ne, nt, nb, npceout, amp, dbounds)
!
!    do id = 1,nip
!        tempi(id) = (id-1)
!        if (Fdi(id) .ne. 0) then
!            ii = (id-1)
!            ivec = Fdi(id)
!            call VecSetValues(Fsi,one,ii,ivec,ADD_VALUES,ierr)
!        end if
!    end do
!
!    do id = 1,nbp
!        tempb(id) = id-1
!        if (Fdi(id) .ne. 0) then
!            ii = (id-1)
!            ivec = Fdg(id)
!            call VecSetValues(Fsg,one,ii,ivec,ADD_VALUES,ierr)
!        end if
!    end do
!
!    call VecAssemblyBegin(Fsi,ierr)
!    call VecAssemblyBegin(Fsg,ierr)
!    call VecAssemblyEnd(Fsi,ierr)
!    call VecAssemblyEnd(Fsg,ierr)
!
!
!END SUBROUTINE StoVecSeqOneLevel

!!!!*********************************************************
!SUBROUTINE GetVecSeqOneLevel(PFi,PFg,SolVeci,SolVecb,tempi,tempb,p,e,t,np,ne,nt,nb,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PFi,PFg,SolVeci,SolVecb
!    PetscInt         i,ii
!    PetscInt         ni,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb
!    double precision :: amp
!    integer, dimension(3,ne) :: e
!    integer, dimension(3,nt) :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2) :: dbounds
!
!    integer, dimension(np-nb) :: tempi
!    integer, dimension(nb) :: tempb
!
!    ni = np-nb
!    call VecCreateSeq(PETSC_COMM_SELF, ni, PFi, ierr)
!    call VecCreateSeq(PETSC_COMM_SELF, nb, PFg, ierr)
!
!    call VecDuplicate(PFi,SolVeci,ierr)
!    call VecDuplicate(PFg,SolVecb,ierr)
!
!    call PETScOneLevelAssemVec(PFi, PFg, p, e, t, np, ne, nt, nb, amp, dbounds,tempi,tempb)
!
!    call VecAssemblyBegin(PFi,ierr)
!    call VecAssemblyBegin(PFg,ierr)
!    call VecAssemblyEnd(PFi,ierr)
!    call VecAssemblyEnd(PFg,ierr)
!
!
!END SUBROUTINE GetVecSeqOneLevel
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!!SUBROUTINE GetMatSeqOneLevel(n,m,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!SUBROUTINE GetMatSeqOneLevel(PAii,PAgg,PAig,PAgi, p, e, t, np, ne, nt, nb, dbounds, const_diff, ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!    Mat              PAii,PAgg,PAig,PAgi
!    PetscInt         i,j,ii,jj
!    PetscInt         ni,m,nz,one
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb
!    double precision :: const_diff
!    integer, dimension(3,ne) :: e
!    integer, dimension(3,nt) :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2) :: dbounds
!
!    ni = np-nb
!    nz = 12  !! number of non-zeros per row for matrix pre-allocation
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ni, ni, nz, PETSC_NULL_INTEGER, PAii, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nb, nb, nz, PETSC_NULL_INTEGER, PAgg, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ni, nb, nz, PETSC_NULL_INTEGER, PAig, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nb, ni, nz, PETSC_NULL_INTEGER, PAgi, ierr)
!
!    call PETScOneLevelAssemADVmatrix(PAii,PAgg,PAig,PAgi,p,e,t,np,ne,nt,nb,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
!    [0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call PETScOneLevelAssemDIFmatrix(PAii,PAgg,PAig,PAgi,p,e,t,np,ne,nt,nb,2,const_diff,[0.0d0,0.0d0,0.0d0,0.0d0],&
!    [0.0d0,0.0d0,0.0d0,0.0d0], 0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!    call MatAssemblyBegin(PAii,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAgg,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAig,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAgi,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAii,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAgg,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAig,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAgi,MAT_FINAL_ASSEMBLY,ierr)
!
!
!END SUBROUTINE GetMatSeqOneLevel
!
!
!!!!*********************************************************
!SUBROUTINE GetVecSeqTwoLevel(PFi,PFg,PFc,SolVeci,SolVecb,SolVecc,tempi,tempb,tempc,p,e,t,np,ne,nt,nb,nci,amp,dbounds,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PFi,PFg,PFc
!    Vec              SolVeci,SolVecb,SolVecc
!    PetscInt         i,ii
!    PetscInt         ni,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb, nci
!    double precision :: amp
!    integer, dimension(3,ne) :: e
!    integer, dimension(3,nt) :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2) :: dbounds
!
!    integer, dimension(np-nb) :: tempi
!    integer, dimension(nb)    :: tempb
!    integer, dimension(nci)   :: tempc
!
!    ni = np-nb
!
!    call VecCreateSeq(PETSC_COMM_SELF, ni, PFi, ierr)
!    call VecCreateSeq(PETSC_COMM_SELF, nb, PFg, ierr)
!    call VecCreateSeq(PETSC_COMM_SELF, nci, PFc, ierr)
!
!    call VecDuplicate(PFi,SolVeci,ierr)
!    call VecDuplicate(PFg,SolVecb,ierr)
!    call VecDuplicate(PFc,SolVecc,ierr)
!
!    call PETScTwoLevelAssemVec(PFi, PFg, PFc, p, e, t, np, ne, nt, nb, nci,amp,dbounds,tempi,tempb,tempc)
!
!    call VecAssemblyBegin(PFi,ierr)
!    call VecAssemblyBegin(PFg,ierr)
!    call VecAssemblyBegin(PFc,ierr)
!    call VecAssemblyEnd(PFi,ierr)
!    call VecAssemblyEnd(PFg,ierr)
!    call VecAssemblyEnd(PFc,ierr)
!
!
!END SUBROUTINE GetVecSeqTwoLevel
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Matrix in most general way
!!SUBROUTINE GetMatSeqOneLevel(n,m,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!SUBROUTINE GetMatSeqTwoLevel(PAmat,PAcc,PAic,PAci,PAir,PAri,PArc,PAcr,p, e, t, &
!                             np, ne, nt, nb, nci, dbounds, const_diff, ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!    Mat              PAmat,PAcc,PAic,PAci
!    Mat              PAir,PAri,PArc,PAcr
!    PetscInt         i,j,ii,jj
!    PetscInt         ni,m,nz,one
!    PetscErrorCode   ierr
!
!    integer :: np, ne, nt, nb, nci, nc, nr
!    double precision :: const_diff
!    integer, dimension(3,ne) :: e
!    integer, dimension(3,nt) :: t
!    double precision, dimension(2,np) :: p
!    double precision, dimension(2,2) :: dbounds
!
!    ni = np-nb
!    nc = np-nci
!    nr = nb-nci
!
!    nz = 12  !! number of non-zeros per row for matrix pre-allocation
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nc, nc, nz, PETSC_NULL_INTEGER, PAmat, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nci, nci, nz, PETSC_NULL_INTEGER, PAcc, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ni, nci, nz, PETSC_NULL_INTEGER, PAic, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nci, ni, nz, PETSC_NULL_INTEGER, PAci, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, ni, nr, nz, PETSC_NULL_INTEGER, PAir, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nr, ni, nz, PETSC_NULL_INTEGER, PAri, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nr, nci, nz, PETSC_NULL_INTEGER, PArc, ierr)
!    call MatCreateSeqAIJ(PETSC_COMM_SELF, nci, nr, nz, PETSC_NULL_INTEGER, PAcr, ierr)
!
!    call PETScTwoLevelAssemADVmatrix(PAmat,PAcc,PAic,PAci,PAir,PAri,PArc,PAcr,p,e,t,np,ne,nt,nb,nci,1,0.0d0, &
!    [0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!    call PETScTwoLevelAssemDIFmatrix(PAmat,PAcc,PAic,PAci,PAir,PAri,PArc,PAcr,p,e,t,np,ne,nt,nb,nci,2,const_diff, &
!    [0.0d0,0.0d0,0.0d0,0.0d0],[0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!    call MatAssemblyBegin(PAmat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAcc,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAic,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAci,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAir,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAri,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PArc,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyBegin(PAcr,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAmat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAcc,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAic,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAci,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAir,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAri,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PArc,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PAcr,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetMatSeqTwoLevel
!
!
!!!!*********************************************************
!!!! Subroutine: To construc PETSc Vector in Parallel MPI Form
!SUBROUTINE GetVecMPI(n,fvec,PetVec,Solvec,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!
!    Vec              PetVec,SolVec
!    PetscInt         i,ii
!    PetscInt         n,one
!    PetscScalar      ivec
!    PetscErrorCode   ierr
!
!    double precision, dimension(n) :: fvec
!
!    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,PetVec,ierr);
!
!    call VecDuplicate(PetVec,SolVec,ierr)
!
!    one = 1
!    do i = 1,n
!    ii = i-1
!    ivec = fvec(i)
!    call VecSetValues(PetVec,one,ii,ivec,INSERT_VALUES,ierr)
!    end do
!
!    call VecAssemblyBegin(PetVec,ierr)
!    call VecAssemblyEnd(PetVec,ierr)
!
!END SUBROUTINE GetVecMPI
!
!
!!!*********************************************************
!!! Subroutine: To construc PETSc Matrix in Sequential Form
!!! Not working : Dec 10,2015
!SUBROUTINE GetMatMPI(n,Amat,PetMat,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!
!    Mat              PetMat
!    PetscInt         i,j,ii,jj
!    PetscInt         n,one
!    PetscInt         Istart,Iend
!    PetscScalar      imat
!    PetscErrorCode   ierr
!
!    double precision, dimension(n,n) :: Amat
!
!    call MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,&
!     & n,n,0,PETSC_NULL_INTEGER,0,PETSC_NULL_INTEGER,PetMat,ierr);
!
!    call MatGetOwnershipRange(PetMat,Istart,Iend,ierr)
!
!    one = 1
!    do i = Istart,Iend-1
!    do j = Istart,Iend-1
!        if (Amat(i,j) .ne. 0) then
!            ii = i
!            jj = j
!            imat = Amat(i,j)
!            call MatSetValues(PetMat,one,ii,one,jj,imat,INSERT_VALUES,ierr)
!        end if
!    end do
!    end do
!
!    call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!    call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE GetMatMPI
!
!!*********************************************************
!! Subroutine : Call PETSc Matrix, Vector and Ksp
!! and solve the linear problem of the form Ax = b
!SUBROUTINE PETScSolver(n, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, temp1, temp2)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscpc.h>
!
!    Vec              PetVec,SolVec
!    Mat              PetMat
!    KSP              ksp
!    PC               pc
!
!    PetscInt         rank,size
!    PetscBool        flg
!    PetscReal        tol
!    PetscInt         n
!    PetscErrorCode   ierr
!
!    double precision, dimension(n)   :: temp2
!    integer, dimension(n)            :: temp1
!
!    integer :: npg, neg, ntg
!    double precision :: amp, meanc
!    integer, dimension(3,neg) :: eg
!    integer, dimension(3,ntg) :: tg
!    double precision, dimension(2,npg) :: pg
!    double precision, dimension(2,2) :: dbounds
!
!!!----------------------------------------------------------------
!!! PETSc-Initialize
!    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
!    call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
!    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!
!!! PETSc-Hello World
!    WRITE(*,*) 'Hello from PETSc World, rank:',rank, ' of total',size,' processes.'
!
!!!----------------------------------------------------------------
!!!PETSc-Assembled-Vector
!    call GetPetVec(n,PetVec,SolVec,pg,eg,tg,npg,neg,ntg,amp,dbounds,ierr)
!
!!!PETSc-Assembled-Matrix
!    call GetPetMat(n,PetMat, pg, eg, tg, npg, neg, ntg, amp, dbounds, meanc, ierr)
!
!!!PETSc-KSP
!    call PETScKSP(ksp,pc,PetMat,PetVec,Solvec,ierr)
!
!!!----------------------------------------------------------------
!!!PETSc-Output to Fortran Output
!    call VecGetValues(SolVec, n, temp1, temp2, ierr)
!
!!!----------------------------------------------------------------
!!! PETSc-Finalize
!    call VecDestroy(PetVec,ierr)
!    call VecDestroy(SolVec,ierr)
!    Call MatDestroy(PetMat,ierr)
!    Call KspDestroy(ksp,ierr)
!
!    call PetscFinalize(ierr)
!
!END SUBROUTINE PETScSolver
!
!!!*********************************************************
!!! Subroutine : Call PETSc Matrix, Vector and Ksp
!!! and solve the linear problem of the form Ax = b
!SUBROUTINE PETScSolverPreAssembled(n,fvec,Amat,temp1,temp2)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscpc.h>
!
!Vec              PetVec,SolVec
!Mat              PetMat
!KSP              ksp
!PC               pc
!
!PetscInt         rank,size
!PetscBool        flg
!PetscReal        tol
!PetscInt         n
!PetscErrorCode   ierr
!
!double precision, dimension(n,n) :: Amat
!double precision, dimension(n)   :: fvec
!double precision, dimension(n)   :: U
!double precision, dimension(n)   :: temp2
!integer, dimension(n)            :: temp1
!
!    !!----------------------------------------------------------------
!    !! PETSc-Initialize
!    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
!    call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr)
!    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!
!    !! PETSc-Hello World
!    WRITE(*,*) 'Hello from PETSc World, rank:',rank, ' of total',size,' processes.'
!
!    !!----------------------------------------------------------------
!    !!PETSc-Vector
!    call GetVec(n,fvec,PetVec,Solvec,ierr)      !! Generalize Form
!    !call GetVecSeq(n,fvec,PetVec,Solvec,ierr)   !! Sequential Form
!    !call GetVecMPI(n,fvec,PetVec,Solvec,ierr)   !! ParallelMPI Form
!
!    !!PETSc-Matrix
!    call GetMat(n,Amat,PetMat,ierr)             !! Generalize Form
!    !call GetMatSeq(n,Amat,PetMat,ierr)          !! Sequential Form
!    !call GetMatMPI(n,Amat,PetMat,ierr)          !! ParallelMPI Form : ERROR
!
!    !!PETSc-KSP
!    call PETScKSP(ksp,pc,PetMat,PetVec,Solvec,ierr)
!
!    !!----------------------------------------------------------------
!    !!PETSc-Output to Fortran Output
!    call VecGetValues(SolVec, n, temp1, temp2, ierr)
!
!    !!----------------------------------------------------------------
!    !! PETSc-Finalize
!    call VecDestroy(PetVec,ierr)
!    call VecDestroy(SolVec,ierr)
!    Call MatDestroy(PetMat,ierr)
!    Call KspDestroy(ksp,ierr)
!
!call PetscFinalize(ierr)
!
!END SUBROUTINE PETScSolverPreAssembled




END MODULE PETScommon
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%** END **%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!!*********************************************************
!!! Subroutine :
!SUBROUTINE PETScMatrix(PetMat,pg,eg,tg,npg,neg,ntg,dbounds,const_diff,ierr)
!
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!
!Mat              PetMat
!PetscInt         n
!PetscErrorCode   ierr
!
!double precision :: const_diff
!double precision, dimension(2,2) :: dbounds
!double precision, dimension(2,npg) :: pg
!integer, dimension(3,ntg) :: tg
!integer, dimension(3,neg) :: eg
!integer :: npg, neg, ntg
!
!n=npg
!call MatCreate(PETSC_COMM_WORLD,PetMat,ierr)
!call MatSetFromOptions(PetMat,ierr)
!call MatSetSizes(PetMat,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!call MatSetUp(PetMat,ierr)
!
!call PETScAssembleADVmatrix(PetMat,pg,eg,tg,npg,neg,ntg,1,0.0d0,[0.0d0,0.0d0,0.0d0,0.0d0],&
![0.0d0,0.0d0,0.0d0,0.0d0],0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,0,0)
!call PETScAssembleDIFmatrix(PetMat,pg,eg,tg,npg,neg,ntg,2,const_diff,[0.0d0,0.0d0,0.0d0,0.0d0],&
![0.0d0,0.0d0,0.0d0,0.0d0], 0.0d0,[0.0d0,0.0d0,0.0d0],dbounds,1,1)
!
!call MatAssemblyBegin(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!call MatAssemblyEnd(PetMat,MAT_FINAL_ASSEMBLY,ierr)
!
!END SUBROUTINE PETScMatrix


