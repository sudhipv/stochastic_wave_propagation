# Copyright (C) 2022 Sudhi P V Ph.D. Candidate sudhipv@cmail.carleton.ca

### Elastic wave propagation in a 3D soil medium
#### Application of Boundary condition may be wrong
################################################################
## Sudhi P V: Modified May-2022
################################################################

from dolfin import *
import numpy as np
import time
import os as os
from scipy import linalg

from mpi4py import MPI

## Initialize MPI for python
comm = MPI.COMM_WORLD
ip = comm.Get_rank()
## print('My rank is ',ip)

if ip==0:
    print("==========================================================")
    print("Running FEniCS in Parallel for 3D deterministic Elastic wave Assembly...")

##### VERY IMPORTANT - TO AVOID FENICS FROM RE ORDERING THE NODES
parameters['reorder_dofs_serial'] = False


## Global deterministic assembly for each KLE mode (here it's only mean term)
def detAssembly(a,m,l):

    if (ip == 0):
        print("Inside Deterministic Assembly")


    ## Dummy problem to define LAYOUT of global problem (the easy way)
    A_g = assemble(Constant(0.)* inner(sigma(u), epsilon(v))*dx)
    M_g = assemble(Constant(0.)*(rho * dot(u,v)*dx))
    L_g = assemble(Constant(0.)*dot(f, v)*dx )

    ## Get dofmap to construct cell-to-dof connectivity
    dofmap = V.dofmap()

    ## Perform assembling
    for cell in cells(mesh):
        dof_idx = dofmap.cell_dofs(cell.index())

        ## Assemble local rhs and lhs
        a_local  = assemble_local(a, cell)
        m_local  = assemble_local(m, cell)
        l_local  = assemble_local(l, cell)

        ## Assemble into global system
        A_g.add_local(a_local,dof_idx, dof_idx)
        M_g.add_local(m_local,dof_idx, dof_idx)
        L_g.add_local(l_local,dof_idx)

    ## Finalize assembling
    A_g.apply("add"), M_g.apply("add"), L_g.apply("add")

#################################################
    # if (ip == 0):
    #     print("Setting boundary conditions")

#### X = 200, Y = 114, Z = 200 ##
#### x[0] = 200, x[1] = 114, x[2] = 200 ##

    tol = 1E-14
    def left(x, on_boundary):
        return on_boundary and x[0] < tol

    def right(x, on_boundary):
        return on_boundary and near(x[0],200.0)

    def bottom(x, on_boundary):
        return on_boundary and x[1] < tol

    def front(x, on_boundary):
        return on_boundary and x[2] < tol

    def back(x, on_boundary):
        return on_boundary and near(x[2],200.0)

    bc_L = DirichletBC(V, Constant((0, 0, 0)), left)
    bc_R = DirichletBC(V, Constant((0, 0, 0)), right)
    bc_BT = DirichletBC(V, Constant((0, 0, 0)), bottom)
    bc_F = DirichletBC(V, Constant((0, 0, 0)), front)
    bc_BK = DirichletBC(V, Constant((0, 0, 0)), back)


    bcs = [bc_L,bc_R,bc_BT,bc_F,bc_BK]

    for bc in bcs:
        bc.apply(A_g)
        #bc.apply(M_g)
        bc.apply(L_g)

    A = A_g.array()
    M = M_g.array()
    b = L_g.get_local()

    return A,M,b


##---------------------------------------------------------------------------
## One-level deterministic assembly matrix decomposer
### Rearranging the degrees of freedom to have all interior and interface nodes together for all displacement components


def OneLevelDDMat(A, M, C, b,u0, v0, a0, dt,beta,gamma, nP, nB,nComp):


    K_T = (M/(beta*dt*dt)) + (C *gamma/(beta*dt))  + A

    # for bc in bcs:
    #     bc.apply(K_T)
    #     #bc.apply(M_g)
    #     bc.apply(b)



    nI = nP - nB
    bi = np.zeros([nI*nComp])
    bg = np.zeros([nB*nComp])

    u0 = u0.vector().array()
    v0 = v0.vector().array()
    a0 = a0.vector().array()

    u0i = np.zeros([nI*nComp])
    u0b = np.zeros([nB*nComp])

    v0i = np.zeros([nI*nComp])
    v0b = np.zeros([nB*nComp])

    a0i = np.zeros([nI*nComp])
    a0b = np.zeros([nB*nComp])

    ADii = np.zeros([nI*nComp,nI*nComp])
    ADgg = np.zeros([nB*nComp,nB*nComp])
    ADig = np.zeros([nI*nComp,nB*nComp])

    MDii = np.zeros([nI*nComp,nI*nComp])
    MDgg = np.zeros([nB*nComp,nB*nComp])
    MDig = np.zeros([nI*nComp,nB*nComp])

    CDii = np.zeros([nI*nComp,nI*nComp])
    CDgg = np.zeros([nB*nComp,nB*nComp])
    CDig = np.zeros([nI*nComp,nB*nComp])

    for j in range(nComp):
        bi[j*nI:(j+1)*nI] = b[nP*j:(nP*j)+nI]
        bg[j*nB:(j+1)*nB] = b[(nP*j)+nI:(nP*j)+nP]

        u0i[j*nI:(j+1)*nI] = u0[nP*j:(nP*j)+nI]
        u0b[j*nB:(j+1)*nB] = u0[(nP*j)+nI:(nP*j)+nP]

        v0i[j*nI:(j+1)*nI] = v0[nP*j:(nP*j)+nI]
        v0b[j*nB:(j+1)*nB] = v0[(nP*j)+nI:(nP*j)+nP]

        a0i[j*nI:(j+1)*nI] = a0[nP*j:(nP*j)+nI]
        a0b[j*nB:(j+1)*nB] = a0[(nP*j)+nI:(nP*j)+nP]



        for k in range(nComp):
            ADii[j*nI:(j+1)*nI, k*nI:(k+1)*nI] = K_T[nP*j:(nP*j)+nI, nP*k:(nP*k)+nI]
            ADgg[j*nB:(j+1)*nB, k*nB:(k+1)*nB] = K_T[(nP*j)+nI:(nP*j)+nP, (nP*k)+nI:(nP*k)+nP];
            ADig[j*nI:(j+1)*nI, k*nB:(k+1)*nB] = K_T[nP*j:(nP*j)+nI, (nP*k)+nI:(nP*k)+nP];

            MDii[j*nI:(j+1)*nI, k*nI:(k+1)*nI] = M[nP*j:(nP*j)+nI, nP*k:(nP*k)+nI]
            MDgg[j*nB:(j+1)*nB, k*nB:(k+1)*nB] = M[(nP*j)+nI:(nP*j)+nP, (nP*k)+nI:(nP*k)+nP];
            MDig[j*nI:(j+1)*nI, k*nB:(k+1)*nB] = M[nP*j:(nP*j)+nI, (nP*k)+nI:(nP*k)+nP];

            CDii[j*nI:(j+1)*nI, k*nI:(k+1)*nI] = C[nP*j:(nP*j)+nI, nP*k:(nP*k)+nI]
            CDgg[j*nB:(j+1)*nB, k*nB:(k+1)*nB] = C[(nP*j)+nI:(nP*j)+nP, (nP*k)+nI:(nP*k)+nP];
            CDig[j*nI:(j+1)*nI, k*nB:(k+1)*nB] = C[nP*j:(nP*j)+nI, (nP*k)+nI:(nP*k)+nP];


    return ADii, ADig, ADgg, MDii, MDig, MDgg, CDii, CDig, CDgg, bi, bg, u0i, v0i, a0i,u0b, v0b, a0b

##---------------------------------------------------------------------------
##---------------------------------------------------------------------------
## Two-level deterministic assembly matrix decomposer
### Rearranging the degrees of freedom to have all remaining and corner nodes together for all displacement components

def TwoLevelDDMat(A, M, C, b, nP, nB, nC, nComp):

    K_T = (M/(beta*dt*dt)) + (C *gamma/(beta*dt)) + A


    nI = nP - nB
    nR = nB - nC
    nIR = nP - nC

    ADic = np.zeros([nI*nComp,nC*nComp])
    ADir = np.zeros([nI*nComp,nR*nComp])
    ADcc = np.zeros([nC*nComp,nC*nComp])
    ADrr = np.zeros([nR*nComp,nR*nComp])
    ADrc = np.zeros([nR*nComp,nC*nComp])


    MDic = np.zeros([nI*nComp,nC*nComp])
    MDir = np.zeros([nI*nComp,nR*nComp])
    MDcc = np.zeros([nC*nComp,nC*nComp])
    MDrr = np.zeros([nR*nComp,nR*nComp])
    MDrc = np.zeros([nR*nComp,nC*nComp])


    CDic = np.zeros([nI*nComp,nC*nComp])
    CDir = np.zeros([nI*nComp,nR*nComp])
    CDcc = np.zeros([nC*nComp,nC*nComp])
    CDrr = np.zeros([nR*nComp,nR*nComp])
    CDrc = np.zeros([nR*nComp,nC*nComp])

    for j in range(nComp):
        for k in range(nComp):
            ADir[j*nI:(j+1)*nI, k*nR:(k+1)*nR] = K_T[nP*j:(nP*j)+nI, (nP*k)+nI:(nP*k)+nIR]
            ADic[j*nI:(j+1)*nI, k*nC:(k+1)*nC] = K_T[nP*j:(nP*j)+nI, (nP*k)+nIR:(nP*k)+nP]
            ADrr[j*nR:(j+1)*nR, k*nR:(k+1)*nR] = K_T[(nP*j)+nI:(nP*j)+nIR, (nP*k)+nI:(nP*k)+nIR]
            ADcc[j*nC:(j+1)*nC, k*nC:(k+1)*nC] = K_T[(nP*j)+nIR:(nP*j)+nP, (nP*k)+nIR:(nP*k)+nP]
            ADrc[j*nR:(j+1)*nR, k*nC:(k+1)*nC] = K_T[(nP*j)+nI:(nP*j)+nIR, (nP*k)+nIR:(nP*k)+nP]

            MDir[j*nI:(j+1)*nI, k*nR:(k+1)*nR] = M[nP*j:(nP*j)+nI, (nP*k)+nI:(nP*k)+nIR]
            MDic[j*nI:(j+1)*nI, k*nC:(k+1)*nC] = M[nP*j:(nP*j)+nI, (nP*k)+nIR:(nP*k)+nP]
            MDrr[j*nR:(j+1)*nR, k*nR:(k+1)*nR] = M[(nP*j)+nI:(nP*j)+nIR, (nP*k)+nI:(nP*k)+nIR]
            MDcc[j*nC:(j+1)*nC, k*nC:(k+1)*nC] = M[(nP*j)+nIR:(nP*j)+nP, (nP*k)+nIR:(nP*k)+nP]
            MDrc[j*nR:(j+1)*nR, k*nC:(k+1)*nC] = M[(nP*j)+nI:(nP*j)+nIR, (nP*k)+nIR:(nP*k)+nP]

            CDir[j*nI:(j+1)*nI, k*nR:(k+1)*nR] = C[nP*j:(nP*j)+nI, (nP*k)+nI:(nP*k)+nIR]
            CDic[j*nI:(j+1)*nI, k*nC:(k+1)*nC] = C[nP*j:(nP*j)+nI, (nP*k)+nIR:(nP*k)+nP]
            CDrr[j*nR:(j+1)*nR, k*nR:(k+1)*nR] = C[(nP*j)+nI:(nP*j)+nIR, (nP*k)+nI:(nP*k)+nIR]
            CDcc[j*nC:(j+1)*nC, k*nC:(k+1)*nC] = C[(nP*j)+nIR:(nP*j)+nP, (nP*k)+nIR:(nP*k)+nP]
            CDrc[j*nR:(j+1)*nR, k*nC:(k+1)*nC] = C[(nP*j)+nI:(nP*j)+nIR, (nP*k)+nIR:(nP*k)+nP]

    # return ADci, ADri, ADcc, ADrr, ADcr
    return ADic, ADir, ADcc, ADrr, ADrc, MDic, MDir, MDcc, MDrr, MDrc, CDic, CDir, CDcc, CDrr, CDrc

##---------------------------------------------------------------------------


################### MAIN CODE STARTS ##################################
## Path to the global meshdims.dat
datapath = './../../../../data/'
meshdimpath = datapath + 'meshData/meshdim.dat'
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[5].astype(int)   ## For 3D


if (ip == 0):
    print("number of partitions:", nParts)

# Nsteps  = 20
# dt = Constant(T/Nsteps)
dt = 0.005;
T = 50*dt;
beta = 0.25;
gamma = 0.5;

if (ip == 0):
    apath = "./../NBparam.dat"
    np.savetxt(apath, (T,dt, beta, gamma),newline=" ",fmt='%1.4e')


nComp = 3    ## 3 for 3D % For assembly


### Deterministic test case
E = 64E+6
nu = 0.3
mu = Constant(E/(2.0*(1.0+nu)))
lambda_ = Constant((E*nu)/((1+nu)*(1-(2.0*nu)))) ##(2.556*0.2778)/((1+0.2778)*(1-(2.0*0.2778)))=1.25
rho = 1700


#velocity = sqrt((lambda_+2*mu)/rho)
#print('Velocity of wave',velocity)
#c_deltat = velocity*dt
#print('Distance in one timestep',c_deltat)


## Get DDM mesh dimension data  ## NOTE: Zfill=4 for both Matlab & Fortran DDM data
mdpath = datapath + "meshData/meshdim"+str(ip+1).zfill(4)+".dat"
meshdim = np.genfromtxt(mdpath)
nNodes = meshdim[0].astype(int)     ## For 3D
nBnodes = meshdim[4].astype(int)    ## For 3D

mcrdpath = datapath + "meshData/dimcrn"+str(ip+1).zfill(4)+".dat"
meshdim = np.genfromtxt(mcrdpath)
nCnodes = meshdim[0].astype(int)
nRnodes = meshdim[1].astype(int)

## Create mesh and define function space
mpath = "../../data/foo"+str(ip+1)+".xml"   ## path to mesh files
if ip == 0:
    print("Input XML-Mesh:", mpath)
mesh = Mesh(mpi_comm_self(),mpath)      ## mesh = Mesh(mpath)
## Initialize the connectivity between facets and cells
tdim = mesh.topology().dim()
mesh.init(tdim-1, tdim)


V = VectorFunctionSpace(mesh, 'P', 1)
gfun = Function(V)

## Define strain and stress
def epsilon(u):
    return sym(nabla_grad(u))

def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)


# Define variational problem
u = TrialFunction(V)
# Initialization of fields (displacement, velocity, acceleration)
u0, v0, a0 = Function(V), Function(V), Function(V)

d = gfun.ufl_domain().geometric_dimension()# space dimension #
v = TestFunction(V)
f = Constant((0, 0, 0))
Tr = Constant((0, 0, 0))


### Gaussian Pulse Initail Condition
ui  = Expression(("0.05*exp(-0.01*(pow(x[0]-100.0,2)+pow(x[1]-90.0,2)+pow(x[2]-100.0,2) ))","0.05*exp(-0.01*(pow(x[0]-100.0,2)+pow(x[1]-90.0,2)+pow(x[2]-100.0,2)))","0.05*exp(-0.01*(pow(x[0]-100.0,2)+pow(x[1]-90.0,2)+pow(x[2]-100.0,2)))"),degree=2)
# u0 = interpolate(ui, V)

# ui = Expression(("0.0","0.0","0.0"), degree=2)

vi = Expression(("0.0","0.0","0.0"), degree=0)
ai = Expression(("0.0","0.0","0.0"), degree=0)


#######

u0 = interpolate(ui,V)

v0 = interpolate(vi, V)

a0 = interpolate(ai, V)


if (ip == 0):
    print("Creating Variational Forms")


# Variation Form
a = Form(inner(sigma(u), epsilon(v))*dx)
m = Form((rho * dot(u,v)*dx))
l = Form(dot(f, v)*dx )


## Save deterministic assembly matrices
subpath = "../../data/Amats/subdom000".rstrip("0")+str(ip+1)

if not os.path.exists(subpath):
    os.makedirs(subpath)

## Invoke Deterministic Assembly procedure
A,M,b = detAssembly(a,m,l)


### Damped
C = 0.0403 * M + 0.0025 * A;

#C = 0.0 * M + 0.0 * A;

print("Saving subdomain level Matrices for part",ip+1)
print(np.shape(A))
print(np.shape(M))
print(np.shape(b))
Apath = subpath+"/A000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(Apath, A)

Mpath = subpath+"/M000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(Mpath, M)

bpath = subpath+"/b000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(bpath, b)


ADii, ADig, ADgg,MDii, MDig, MDgg, CDii, CDig, CDgg, bi, bg,  u0i, v0i, a0i,u0b, v0b, a0b = OneLevelDDMat(A, M, C, b,u0, v0, a0,dt,beta,gamma, nNodes, nBnodes, nComp)


ADic, ADir, ADcc, ADrr, ADrc, MDic, MDir, MDcc, MDrr, MDrc, CDic, CDir, CDcc, CDrr, CDrc = TwoLevelDDMat(A, M, C, b, nNodes, nBnodes, nCnodes,nComp)


ADiipath = subpath+"/ADii000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(ADiipath, ADii)
ADigpath = subpath+"/ADig000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(ADigpath, ADig)
ADggpath = subpath+"/ADgg000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(ADggpath, ADgg)

MDiipath = subpath+"/MDii000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(MDiipath, MDii)
MDigpath = subpath+"/MDig000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(MDigpath, MDig)
MDggpath = subpath+"/MDgg000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(MDggpath, MDgg)

CDiipath = subpath+"/CDii000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(CDiipath, CDii)
CDigpath = subpath+"/CDig000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(CDigpath, CDig)
CDggpath = subpath+"/CDgg000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(CDggpath, CDgg)


fipath = subpath+"/bi000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(fipath, bi)
fgpath = subpath+"/bg000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(fgpath, bg)

upath = subpath+"/u0i_"+str(ip+1)+".dat"
np.savetxt(upath, u0i)
vpath = subpath+"/v0i_"+str(ip+1)+".dat"
np.savetxt(vpath, v0i)
apath = subpath+"/a0i_"+str(ip+1)+".dat"
np.savetxt(apath, a0i)


upath = subpath+"/u0b_"+str(ip+1)+".dat"
np.savetxt(upath, u0b)
vpath = subpath+"/v0b_"+str(ip+1)+".dat"
np.savetxt(vpath, v0b)
apath = subpath+"/a0b_"+str(ip+1)+".dat"
np.savetxt(apath, a0b)

##################Two Level ####################################


ADicpath = subpath+"/ADic000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(ADicpath, ADic)
ADirpath = subpath+"/ADir000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(ADirpath, ADir)
ADccpath = subpath+"/ADcc000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(ADccpath, ADcc)
ADrrpath = subpath+"/ADrr000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(ADrrpath, ADrr)
ADrcpath = subpath+"/ADrc000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(ADrcpath, ADrc)


MDicpath = subpath+"/MDic000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(MDicpath, MDic)
MDirpath = subpath+"/MDir000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(MDirpath, MDir)
MDccpath = subpath+"/MDcc000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(MDccpath, MDcc)
MDrrpath = subpath+"/MDrr000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(MDrrpath, MDrr)
MDrcpath = subpath+"/MDrc000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(MDrcpath, MDrc)


CDicpath = subpath+"/CDic000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(CDicpath, CDic)
CDirpath = subpath+"/CDir000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(CDirpath, CDir)
CDccpath = subpath+"/CDcc000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(CDccpath, CDcc)
CDrrpath = subpath+"/CDrr000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(CDrrpath, CDrr)
CDrcpath = subpath+"/CDrc000".rstrip("0")+str(ip+1)+".dat"
np.savetxt(CDrcpath, CDrc)


# fcpath = subpath+"/bc000".rstrip("0")+str(ip+1)+".dat"
# np.savetxt(fcpath, bc)
# frpath = subpath+"/br000".rstrip("0")+str(ip+1)+".dat"
# np.savetxt(frpath, br)



##############################################################

if (ip == 0):
    print("==========================Success============================")

########################  ############################





















