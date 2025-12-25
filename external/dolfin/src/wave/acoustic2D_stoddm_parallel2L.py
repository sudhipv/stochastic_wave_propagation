

### Acoustic Wave Propagation Problem ######

#### Copyright (C) Sudhi P V #####

### June -17 - 2020


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
    print("Running FEniCS in Parallel for 2D stochastic Acoustic Assembly...")

##### VERY IMPORTANT - TO AVOID FENICS FROM RE ORDERING THE NODES
parameters['reorder_dofs_serial'] = False


def OneLevelDDMat(A, M, C, b,u0, v0, a0, dt,beta,gamma, nNodes, nBnodes):


    # K_T = (M/(beta*dt*dt)) + (C *gamma/(beta*dt))  + A

    # K_T can be calculated only after stochastic projection

    ADii = A[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    ADig = A[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    ADgg = A[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]

    MDii = M[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    MDig = M[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    MDgg = M[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]

    CDii = C[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    CDig = C[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    CDgg = C[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]


    bi   = b[0:(nNodes-nBnodes)]
    bg   = b[(nNodes-nBnodes):nNodes]


    u0_vec = u0.vector().array()
    u0_vecI = u0_vec[0:(nNodes-nBnodes)]
    u0_vecB = u0_vec[(nNodes-nBnodes):nNodes]

    v0_vec = v0.vector().array()
    v0_vecI = v0_vec[0:(nNodes-nBnodes)]
    v0_vecB = v0_vec[(nNodes-nBnodes):nNodes]

    a0_vec = a0.vector().array()
    a0_vecI = a0_vec[0:(nNodes-nBnodes)]
    a0_vecB = a0_vec[(nNodes-nBnodes):nNodes]


    if(k ==0):
        return ADii, ADig, ADgg, MDii, MDig, MDgg, CDii, CDig, CDgg, bi, bg, u0_vecI, v0_vecI, a0_vecI, u0_vecB, v0_vecB, a0_vecB
    else:
        return ADii, ADig, ADgg



def TwoLevelDDMat(A, M, C, b, nNodes, nBnodes, nCnodes):

    # K_T = (M/(beta*dt*dt)) + (C *gamma/(beta*dt))  + A


    ADic = A[0:(nNodes-nBnodes),(nNodes-nCnodes):nNodes]
    ADir = A[0:(nNodes-nBnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    ADcc = A[(nNodes-nCnodes):nNodes,(nNodes-nCnodes):nNodes]
    ADrr = A[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    ADrc = A[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nCnodes):nNodes]


    MDic = M[0:(nNodes-nBnodes),(nNodes-nCnodes):nNodes]
    MDir = M[0:(nNodes-nBnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    MDcc = M[(nNodes-nCnodes):nNodes,(nNodes-nCnodes):nNodes]
    MDrr = M[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    MDrc = M[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nCnodes):nNodes]


    CDic = C[0:(nNodes-nBnodes),(nNodes-nCnodes):nNodes]
    CDir = C[0:(nNodes-nBnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    CDcc = C[(nNodes-nCnodes):nNodes,(nNodes-nCnodes):nNodes]
    CDrr = C[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nBnodes):(nNodes-nCnodes)]
    CDrc = C[(nNodes-nBnodes):(nNodes-nCnodes),(nNodes-nCnodes):nNodes]


    bc   = b[(nNodes-nCnodes):nNodes]
    br   = b[(nNodes-nBnodes):(nNodes-nCnodes)]


    ##### Mass and Damping matrices are same for each PC coefficients
    ##### Initial conditions are deterministic so non zero for only mean term and zero for all other input PC terms.

    if(k ==0):
        return ADic, ADir, ADcc, ADrr, ADrc, MDic, MDir, MDcc, MDrr, MDrc, CDic, CDir, CDcc, CDrr, CDrc, br, bc
    else:
        return ADic, ADir, ADcc, ADrr, ADrc



## Global deterministic assembly for each KLE mode (here it's only mean term)
def detAssembly(a,m,l):

    # if (ip == 0):
    #     print("Inside Deterministic Assembly")
    ## Dummy problem to define LAYOUT of global problem (the easy way)
    A_g = assemble(Constant(0.)* c * inner(nabla_grad(u1),nabla_grad(w))*dx)
    M_g = assemble(Constant(0.)*u1*w*dx)
    L_g = assemble(Constant(0.)*f*w*dx)


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
    # if (k == 0):
    #     print("Setting boundary conditions")

    def boundary_S(x, on_boundary):
        tol = 1E-14
        return on_boundary


    # Mark boundary subdomians (2D domain (0,1) x (0,1)
    left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
    right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
    top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0)
    bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)

    bc_L = DirichletBC(V, Constant(0), left)
    bc_R = DirichletBC(V, Constant(0), right)
    bc_T = DirichletBC(V, Constant(0), top)
    bc_B = DirichletBC(V, Constant(0), bottom)

    ## Apply boundary conditions
    bcs = [bc_L,bc_R,bc_T,bc_B]


    if k==0:
        for bc in bcs:
            bc.apply(A_g)
            bc.apply(M_g)
            bc.apply(L_g)

    A = A_g.array()
    M = M_g.array()
    b = L_g.get_local()

    ##### Matrices created using assemble functions don't preserve symmetry and
    ##### thus are not completely symmetric. either use assemble_system or do manual code to make them symmetric
    ### Refer : Chapter 6 Fenics book 2011-10-27 section 6.3

    return A,M,b

################################################################################

####              Stochastic Part of Code                        #

################################################################################


### Ajit's Part Used
## Spatialy varying KLE/PCE modes
class MyExpression(Expression):
    def __init__(self, params, **kwargs):
        self.index  = params[0]
        self.ndim   = params[1]
        self.sIndex = params[2].astype(int)
        self.mIndex = params[3].astype(int)

    def eval(self, values, x):
        ## For bx=by=1 good for upto 15 RVs cases (unit square)
        multipliers = [0.92184, 0.49248, 0.29374, 0.20437, 0.15576]
        omegas = [1.30654, 3.67319, 6.58462, 9.63168, 12.72324]

        # Log normal Mean and Standard Deviation of Underlaying Gaussian
        meanl = 1.005     # similar to meanc in Fortran package
        sigma = 0.1      # similar to sigma in Fortran package

### Trunctaed PCE of lognormal process
## Automation to n-RVs : ----------------------------------------------------------------
## KLE: obtain KLE terms
        g = []
        for i in range(self.ndim):
            Xindex = self.sIndex[i,0]
            # print('Xindex', Xindex)
            ### First Dimesnion
            if (Xindex % 2) == 0:
                ##print("even")
                gg1 = multipliers[Xindex-1] * (sin(omegas[Xindex-1]*(x[0]-0.5)))
            else:
                ##print("odd")
                gg1 = multipliers[Xindex-1] * (cos(omegas[Xindex-1]*(x[0]-0.5)))

            ### Second Dimension
            Yindex = self.sIndex[i,1]
            if (Yindex % 2) == 0:
                ##print("even")
                gg2 = multipliers[Yindex-1] * (sin(omegas[Yindex-1]*(x[1]-0.5)))
            else:
                ##print("odd")
                gg2 = multipliers[Yindex-1] * (cos(omegas[Yindex-1]*(x[1]-0.5)))

            g.append(sigma*gg1*gg2)


## PCE: obtain PCE terms
        newY = 1.0
        i = self.index
        for j in range(self.ndim):
            idx = mIndex[i,j]
            if idx == 0:
                yy = 1.0
            elif idx == 1:
                yy = g[j]
            else:
                varFacto = nfactorial(idx)
                yy = (g[j]**(idx))/np.sqrt(varFacto)

            ## Keep all indentations properly (there was error here due to indentation)
            newY =  newY*yy

        values[0] = meanl * newY

##---------------------------------------------------------------------------
## Stochastic variational form
def stoVartional(params, u1, w, f):
    kleID = params[0]
    # print('KLE id', kleID)
    cd = MyExpression(params, degree=0)
    m = u1*w*dx
    a = cd* inner(nabla_grad(u1),nabla_grad(w))*dx
    l = f*w*dx
    return a,m,l

##---------------------------------------------------------------------------

################### MAIN CODE STARTS ##################################


nComp = 1 # 3 for 3D, 2 for 2D

## Path to the global meshdims.dat
meshdimpath = '../../../../data/meshData/meshdim.dat'
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[4].astype(int)     ## For 2D

if (ip == 0):
    print("number of partitions:", nParts)

T = 0.5;
dt = 6.5e-3;
beta = 0.25;
gamma = 0.5;

if (ip == 0):
    apath = "./../NBparam.dat"
    np.savetxt(apath, (T,dt, beta, gamma),newline=" ",fmt='%1.4e')


## Get DDM mesh dimension data  ## NOTE: Zfill=4 for both Matlab & Fortran DDM data
mdpath = "../../../../data/meshData/meshdim"+str(ip+1).zfill(4)+".dat"
meshdim = np.genfromtxt(mdpath)
nNodes = meshdim[0].astype(int)    ## For 2D
nBnodes = meshdim[3].astype(int)   ## For 2D

mcrdpath = "../../../../data/meshData/dimcrn"+str(ip+1).zfill(4)+".dat"
meshdim = np.genfromtxt(mcrdpath)
nCnodes = meshdim[0].astype(int)
nRnodes = meshdim[1].astype(int)


pcedatapath = "../../../../data/klePceData/pcedata.dat"
pcedata = np.genfromtxt(pcedatapath)
order_u = pcedata[0].astype(int)     ## order of output PCE (input==2)
dim_L   = pcedata[1].astype(int)     ## number of KLE modes (RVs)
PCE_u   = pcedata[2].astype(int)     ## number of PCE outputs
PCE_A   = pcedata[3].astype(int)    ## number of PCE inputs

## Load KLE/PCE Indices for Automation
sIndexPath = "../../../../data/klePceData/sortIndex.dat"
sIndex = np.genfromtxt(sIndexPath)

mIndexPath = "../../../../data/klePceData/multiIndex.dat"
mIndex = np.genfromtxt(mIndexPath)

# print('sort index is', sIndex)


def nfactorial(nf):
    if nf == 0:
        return 1
    else:
        return nf * nfactorial(nf-1)

## Create mesh and define function space
mpath = "../../data/foo"+str(ip+1)+".xml"   ## path to mesh files
if ip == 0:
    print("Input XML-Mesh:", mpath)
mesh = Mesh(mpi_comm_self(),mpath)

# Initialize the connectivity between facets and cells
tdim = mesh.topology().dim()
mesh.init(tdim-1, tdim)
# print("cells are",mesh.cells())
# print("coordinates are",mesh.coordinates())

# meshfile = File("mesh_"+str(ip+1)+".pvd")
# meshfile << mesh


# Function space definition over the mesh
V = FunctionSpace(mesh, "CG", 1)   #Continuous Galerkin for the displacement field

# Test and Trial function
u1, w = TrialFunction(V), TestFunction(V)
# Initialization of fields (displacement, velocity, acceleration)
u0, v0, a0 = Function(V), Function(V), Function(V)

##### Initial condition ######

### Gaussian Pulse Initail Condition
ui  = Expression(("1*exp(-100*(pow(x[0]-0.7,2)+pow(x[1]-0.7,2)))"),degree=1)

# ui  = Expression(("1*exp(-50*(pow(x[0]-0.5,2)+pow(x[1]-0.5,2)))"),degree=1)
#######

# ui  = Expression(("sin(m*pi*x[0])*sin(n*pi*x[1])"),degree=1,m = 2,n = 1)

u0 = interpolate(ui, V)

v0 = interpolate(Constant(0.0), V)

a0 = interpolate(Constant(0.0), V)

f  = Constant((0.0))
# f = Expression(("5*pow(pi,2)*sin(m*pi*x[0])*sin(n*pi*x[1])"),degree=2,m = 1,n = 2)

c = 1 # wave velocity


for k in range(int(float(PCE_A))):

# Governing equation in the weak form

    if (ip == 0):
        print("Creating Variational Forms")

    params = ([k, dim_L, sIndex, mIndex])

    ##--------------------------------------------------------------------
    ## FEniCS (genearalize) assembly procedure ## NEW
    ## Invoke variation fomrulaiton for each KLE/PCE mode
    a,m,l = stoVartional(params, u1, w, f)

    ## Invoke Deterministic Assembly procedure for each KLE/PCE mode
    A,M,b = detAssembly(a,m,l)     ## New Assembly

    ## Save deterministic assembly matrices in ".mat" or ".txt" format
    subpath = "../../data/Amats/subdom000".rstrip("0")+str(ip+1)

    if not os.path.exists(subpath):
        os.makedirs(subpath)

    #### Rayleigh Damping
    if(k ==0):
        ### Damped
        # C = 0.5 * M + 0.01 * A;

        C = 0.5445 * M + 0.0174 * A

    # else

    #  Since the PC coefficients of damping matrix is product of algebraic constant, they are not seperately calculated and stored.
    #     C = 0.0174 * A

    # print("Saving subdomain level Matrices for part",ip+1)
    # print(np.shape(A))
    # print(np.shape(M))
    # print(np.shape(b))
    # Apath = subpath+"/A000".rstrip("0")+str(ip+1)+".dat"
    # np.savetxt(Apath, A)

    # Mpath = subpath+"/M000".rstrip("0")+str(ip+1)+".dat"
    # np.savetxt(Mpath, M)

    # bpath = subpath+"/b000".rstrip("0")+str(ip+1)+".dat"
    # np.savetxt(bpath, b)

    ## Decompose assembly matrices intor DD-blocks
    #### According to the random damping matrix PC approximation, mean coefficient is C = alpha M + beta K_0. and other coefficients are just beta*K_i.
    #### Thus only K_i is needed for constructing all the PC coefficients of damping matrix. Thats y we are saving only some terms of PC coefficients.

    if(k == 0):
        ADii, ADig, ADgg,MDii, MDig, MDgg, CDii, CDig, CDgg, bi, bg, u0_vecI, v0_vecI, a0_vecI, u0_vecB, v0_vecB, a0_vecB \
         = OneLevelDDMat(A, M, C, b,u0, v0, a0,dt,beta,gamma, nNodes, nBnodes)
    else:
        ADii, ADig, ADgg = OneLevelDDMat(A, M, C, b,u0, v0, a0,dt,beta,gamma, nNodes, nBnodes)

    if(k == 0):
        ADic, ADir, ADcc, ADrr, ADrc, MDic, MDir, MDcc, MDrr, MDrc, CDic, CDir, CDcc, CDrr, CDrc, br, bc = TwoLevelDDMat(A, M, C, b, nNodes, nBnodes, nCnodes)
    else:
        ADic, ADir, ADcc, ADrr, ADrc = TwoLevelDDMat(A, M, C, b, nNodes, nBnodes, nCnodes)


    print("Saving Decomposed subdomain level Matrices for part",ip+1)

    ADiipath = subpath+"/ADii000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADiipath, ADii)
    ADigpath = subpath+"/ADig000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADigpath, ADig)
    ADggpath = subpath+"/ADgg000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADggpath, ADgg)


    ADicpath = subpath+"/ADic000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADicpath, ADic)
    ADirpath = subpath+"/ADir000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADirpath, ADir)
    ADccpath = subpath+"/ADcc000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADccpath, ADcc)
    ADrrpath = subpath+"/ADrr000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADrrpath, ADrr)
    ADrcpath = subpath+"/ADrc000".rstrip("0")+str(k+1)+".dat"
    np.savetxt(ADrcpath, ADrc)


    if(k ==0):

        fipath = subpath+"/bi000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(fipath, bi)
        fgpath = subpath+"/bg000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(fgpath, bg)


        fcpath = subpath+"/bc000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(fcpath, bc)
        frpath = subpath+"/br000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(frpath, br)

        MDiipath = subpath+"/MDii000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(MDiipath, MDii)
        MDigpath = subpath+"/MDig000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(MDigpath, MDig)
        MDggpath = subpath+"/MDgg000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(MDggpath, MDgg)

        CDiipath = subpath+"/CDii000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(CDiipath, CDii)
        CDigpath = subpath+"/CDig000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(CDigpath, CDig)
        CDggpath = subpath+"/CDgg000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(CDggpath, CDgg)

        upath = subpath+"/u0I_"+str(k+1)+".dat"
        np.savetxt(upath, u0_vecI)
        vpath = subpath+"/v0I_"+str(k+1)+".dat"
        np.savetxt(vpath, v0_vecI)
        apath = subpath+"/a0I_"+str(k+1)+".dat"
        np.savetxt(apath, a0_vecI)

        upath = subpath+"/u0B_"+str(k+1)+".dat"
        np.savetxt(upath, u0_vecB)
        vpath = subpath+"/v0B_"+str(k+1)+".dat"
        np.savetxt(vpath, v0_vecB)
        apath = subpath+"/a0B_"+str(k+1)+".dat"
        np.savetxt(apath, a0_vecB)


        MDicpath = subpath+"/MDic000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(MDicpath, MDic)
        MDirpath = subpath+"/MDir000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(MDirpath, MDir)
        MDccpath = subpath+"/MDcc000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(MDccpath, MDcc)
        MDrrpath = subpath+"/MDrr000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(MDrrpath, MDrr)
        MDrcpath = subpath+"/MDrc000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(MDrcpath, MDrc)


        CDicpath = subpath+"/CDic000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(CDicpath, CDic)
        CDirpath = subpath+"/CDir000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(CDirpath, CDir)
        CDccpath = subpath+"/CDcc000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(CDccpath, CDcc)
        CDrrpath = subpath+"/CDrr000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(CDrrpath, CDrr)
        CDrcpath = subpath+"/CDrc000".rstrip("0")+str(k+1)+".dat"
        np.savetxt(CDrcpath, CDrc)


    ##############################################################

    if (ip == 0):
        print("==========================Success============================")

########################  ############################











