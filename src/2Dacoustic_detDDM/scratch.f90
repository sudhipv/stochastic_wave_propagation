
SUBROUTINE PETSc_detnncpcgm(pid,p,e,t,np,ne,nt,nb,nbg,nr,nc,nbgc,Ui,Ub,Ub_g,maxiter,tol)

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>

!!!-------------------------------------------------------------------------------------------
!! PETSc-Variables Declaration
    Vec              PetVeci,PetVecir,PetVecc,PetVecb,PetVecr   !! PETSc Vectors
    Vec              SolVeci,SolVecir,SolVecc,SolVecb,SolVecr
    Vec              PetVecbg,PetVeccg,PetVecc2,PetVecc3
    Vec              Fsi,Fsg,PPb

    Mat              RsMat,DsMat,BcMat,RrMat,RcMat              !! PETSc Matrix
    Mat              PAmat,PAii,PAgg,PArr,PAcc
    Mat              PAgi,PAri,PAci,PAcr

    KSP              kspAm,kspAi,kspAc                          !! PETSc solver context
    PC               pcAm,pcAi,pcAc                             !! PETSc preconditioner

!!!-------------------------------------------------------------------------------------------
!! Fortran-Variables Declaration
    integer :: pid, ierr, i, maxiter
    integer :: np, ni, nb, nbg, ne, nt
    ! integer :: nip, nbp, nrp, ncp, nirp, nbgp, nbgcp
    integer :: nc, nr, nbgc, nir !casep,

    integer, dimension((np-nb)) :: tempi                !! Temp arrays
    integer, dimension((np-nc)) :: tempir               !! require for VecGetVales
    integer, dimension(nc)      :: tempc
    integer, dimension(nr)      :: tempr
    integer, dimension(nb)      :: tempb
    integer, dimension(nbg)     :: tempbg
    integer, dimension(nbgc)    :: tempcg
    !integer, dimension(ncijk,3)         :: ijk
    integer, dimension(3,ne)            :: e
    integer, dimension(3,nt)            :: t

    !!double precision, dimension(nomga)           :: omegas, multipliers
    double precision, dimension(nbg)     :: Ub_g, rb_g, Qb_g
    double precision, dimension(nbg)     :: Pb_g, Zb_g, RZb
    double precision, dimension(nbgc)    :: Dc, Zc, DcBc
    double precision, dimension(nr)      :: Vri, FrS
    !!double precision, dimension(2,2)             :: dbounds
    !double precision, dimension(ncijk)           :: cijk
    double precision, dimension(np-nc) :: ybc
    double precision, dimension(nb)      :: Ub
    double precision, dimension(np-nb) :: Ui
    double precision, dimension(2,np)            :: p

    double precision :: NegOne, PosOne    !!const_diff, sigma, amp, N
    double precision :: rho_next, rho_curr, alpha, beta, err, tol !! residual

    !!integer, dimension(ndim,npceout) :: mIndex
    !!integer, dimension(2,ndim)       :: sIndex

    double precision :: tolc = 1.0d-6                           !! Coarse-solver tolerence

!!!-------------------------------------------------------------------------------------------
!! PETSc-Initialize
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    ni   = np-nb           !! det-interior nodes
    nir  = np-nc           !! det-interior+remaining nodes
    ! nip  = (ni*npceout)    !! sto-interior nodes
    ! nbp  = (nb*npceout)    !! sto-local-boundary nodes
    ! ncp  = (nc*npceout)    !! sto-local-corner nodes
    ! nrp  = (nr*npceout)    !! sto-local-remaining nodes
    ! nirp = (nir*npceout)   !! sto-interior+remaining nodes
    ! nbgp = (nbg*npceout)   !! sto-global-boundary nodes
    ! nbgcp= (nbgc*npceout)  !! sto-global-corner nodes

!! Initiate Array's
    rho_next = 0.0d0
    rho_curr = 0.0d0
    NegOne   = -1.0
    PosOne   = 1.0

!!!-------------------------------------------------------------------------------------------
!! Assemble only Subdomain-Level Deterministic matrices to calculate Mallocs

    ! IF (mallocsCals .eq. 1) then

    !     if (pid .eq. 0) print*, '------------------------------------------------------'
    !     if (pid .eq. 0) print*, 'Initializing Two-Level-NNC-PCGM Mallocs Calculation...'
    !     if (pid .eq. 0) print*, '------------------------------------------------------'

    !     !call GetMallocs(pid,p,e,t,np,ne,nt,nb,nc,nr,ndim,npcein,npceout,nomga,ncijk,ijk,&
    !     !                cijk,ni,nip,nbp,ncp,nrp,nirp,dbounds, const_diff, omegas, &
    !     !                mIndex, sIndex, multipliers, sigma, casep, ierr)

    !     call GetMallocsShort(pid,np,nb,nc,nr,npceout,nip,nbp,ncp,nrp,nirp,ierr)

    Ui(:)  = 0.0d0
    Ub(:)  = 0.0d0
    Ub_g(:)= 0.0d0

    ! END IF

!!!---------------------------------------------------------------------------------------
    if (pid .eq. 0) print*, '------------------------------------------------'
    if (pid .eq. 0) print*, 'Initializing Stochastic-Two-Level-NNC-PCGM...'
    if (pid .eq. 0) print*, '------------------------------------------------'

!!!-------------------------------------------------------------------------------------------
!! Assemble Subdomain-Level Matrices directly as PETSc Matrices
!
!    call StoMatSeqOneLevel( pid, PAmat, PAii, PAgg, PAgi, PAri, PArr, PAcc, PAci, PAcr, &
!                           p,e,t,np,ne,nt,nb,nc,nr,ndim,npcein,npceout,nomga,ncijk,ijk, &
!                           cijk,ni,nip,nbp,ncp, nrp, nirp, dbounds, const_diff, omegas, &
!                           mIndex, sIndex, multipliers, sigma, casep, ierr)


    call GetMatSeqTwoLevel(PAmat,PAii,PAgg,PAgi,PAri,PArr,PAcc,PAci,PAcr,np,nb,nc,nr,nbgc,ierr)

    stop 123
