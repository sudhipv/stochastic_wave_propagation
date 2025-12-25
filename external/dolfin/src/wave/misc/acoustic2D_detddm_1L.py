

### Acoustic Wave Propagation Problem ######

#### Copyright (C) Sudhi P V #####

### MAY -24 - 2020

from dolfin import *
import numpy as np
import time
import os as os
from scipy import linalg
import matplotlib.pyplot as plt

##### VERY IMPORTANT - TO AVOID FENICS FROM RE ORDERING THE NODES

parameters['reorder_dofs_serial'] = False


def OneLevelDDMat(A, M, C, b,u0,v0,a0, dt,beta,gamma, nNodes, nBnodes):


    K_T = (M/(beta*dt*dt)) + (C *gamma/(beta*dt))  + A

    ADii = K_T[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    ADig = K_T[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    ADgg = K_T[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]

    MDii = M[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    MDig = M[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    MDgg = M[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]

    CDii = C[0:(nNodes-nBnodes),0:(nNodes-nBnodes)]
    CDig = C[0:(nNodes-nBnodes),(nNodes-nBnodes):nNodes]
    CDgg = C[(nNodes-nBnodes):nNodes,(nNodes-nBnodes):nNodes]


    bi   = b[0:(nNodes-nBnodes)]
    bg   = b[(nNodes-nBnodes):nNodes]

    u0 = u0.vector().array()
    v0 = v0.vector().array()
    a0 = a0.vector().array()

    return ADii, ADig, ADgg, MDii, MDig, MDgg, CDii, CDig, CDgg, bi, bg, u0, v0, a0


## Global deterministic assembly for each KLE mode (here it's only mean term)
def detAssembly(a,m,l):

    print("Inside Deterministic Assembly")
    ## Dummy problem to define LAYOUT of global problem (the easy way)
    A_g = assemble(Constant(0.)* c * c * inner(nabla_grad(u1),nabla_grad(w))*dx)
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

    print("Setting boundary conditions")

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

    for bc in bcs:
        bc.apply(A_g)
        bc.apply(M_g)
        bc.apply(L_g)

    A = A_g.array()
    M = M_g.array()
    b = L_g.get_local()

    return A,M,b



################### MAIN CODE STARTS ##################################


nComp = 1 # 3 for 3D, 2 for 2D

## Path to the global meshdims.dat
meshdimpath = '../../../../data/meshData/meshdim.dat'
nParts = np.genfromtxt(meshdimpath)
nParts = nParts[4].astype(int)     ## For 2D
print("number of partitions:", nParts)

T = 1;
dt = 0.01;
beta = 0.25;
gamma = 0.5;


apath = "./../NBparam.dat"
np.savetxt(apath, (T,dt, beta, gamma))

for ip in range(nParts):
    ## Get DDM mesh dimension data  ## NOTE: Zfill=4 for both Matlab & Fortran DDM data
    mdpath = "../../../../data/meshData/meshdim"+str(ip+1).zfill(4)+".dat"
    meshdim = np.genfromtxt(mdpath)
    nNodes = meshdim[0].astype(int)    ## For 2D
    nBnodes = meshdim[3].astype(int)   ## For 2D

    mcrdpath = "../../../../data/meshData/dimcrn"+str(ip+1).zfill(4)+".dat"
    meshdim = np.genfromtxt(mcrdpath)
    nCnodes = meshdim[0].astype(int)
    nRnodes = meshdim[1].astype(int)

    ## Create mesh and define function space
    mpath = "../../data/foo"+str(ip+1)+".xml"   ## path to mesh files
    if ip == 0:
        print("Input XML-Mesh:", mpath)
    mesh = Mesh(mpath)

    # Initialize the connectivity between facets and cells
    tdim = mesh.topology().dim()
    mesh.init(tdim-1, tdim)

    # Function space definition over the mesh
    V = FunctionSpace(mesh, "CG", 1)   #Continuous Galerkin for the displacement field

    # Test and Trial function
    u1, w = TrialFunction(V), TestFunction(V)
    # Initialization of fields (displacement, velocity, acceleration)
    u0, v0, a0 = Function(V), Function(V), Function(V)

    ##### Initial condition ######

    ui  = Expression(("sin(m*pi*x[0])*sin(n*pi*x[1])"),degree=2,m = 2,n = 1)

    u0 = interpolate(ui, V)

    v0 = interpolate(Constant(0.0), V)

    a0 = interpolate(Constant(0.0), V)

    f  = Constant((0.0))
    # f = Expression(("5*pow(pi,2)*sin(m*pi*x[0])*sin(n*pi*x[1])"),degree=2,m = 1,n = 2)
    c = 1 # wave velocity
    # Governing equation in the weak form

    print("Creating Variational Forms")


    m = u1*w*dx
    a = c* c * inner(nabla_grad(u1),nabla_grad(w))*dx
    l = f*w*dx

    ## Save deterministic assembly matrices in ".mat" or ".txt" format
    subpath = "../../data/Amats/subdom000".rstrip("0")+str(ip+1)

    if not os.path.exists(subpath):
        os.makedirs(subpath)

   ## Invoke Deterministic Assembly procedure for each KLE/PCE mode
    A,M,b = detAssembly(a,m,l) ## New assembly procedure

    #### Rayleigh Damping

    C = 0.0 * M + 0.0 * A;


    ### To check original Matrices ###

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

    #### A matrices here are K_T kept the same name for easier conversion of static code

    ADii, ADig, ADgg,MDii, MDig, MDgg, CDii, CDig, CDgg, bi, bg, u0, v0, a0 = OneLevelDDMat(A, M, C,b,u0,v0,a0,dt,beta,gamma, nNodes, nBnodes)
    # ADci, ADri, ADcc, ADrr, ADcr = TwoLevelDDMat(A, b, nNodes, nBnodes, nCnodes, nComp)

    print("Saving Decomposed subdomain level Matrices for part",ip+1)

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

    upath = subpath+"/u0_"+str(ip+1)+".dat"
    np.savetxt(upath, u0)
    vpath = subpath+"/v0_"+str(ip+1)+".dat"
    np.savetxt(vpath, v0)
    apath = subpath+"/a0_"+str(ip+1)+".dat"
    np.savetxt(apath, a0)


print("==========================Success============================")

########################  ############################











