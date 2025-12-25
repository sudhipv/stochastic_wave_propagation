

### Acoustic Wave Propagation Problem ######

#### Sudhi P V #####


from dolfin import *
import numpy as np
import time
import os as os
from scipy import linalg
import matplotlib.pyplot as plt


def newmark_update(u, u0, v0, a0, beta, gamma, dt):

    u_vec = u
    print("u0 is",u0)
    u0_vec = u0
    v0_vec = v0
    a0_vec = a0

 # Update acceleration and velocity
    a_vec = (1.0/(2.0*beta))*( (u_vec - u0_vec -
    v0_vec*dt)/(0.5*dt*dt) - (1.0-2.0*beta)*a0_vec )

    v_vec = dt*((1.0-gamma)*a0_vec + gamma*a_vec) + v0_vec

    v0 = v_vec
    a0 = a_vec
    u0 = u



def GetTransientForce(u0,v0,a0,nNodes,nBnodes,A,M,b):

    u0_vec = u0.vector()
    v0_vec, a0_vec = v0.vector(), a0.vector()
    fun_T = (1.0/(2.0*beta))*(((u0_vec+dt*v0_vec)/0.5*dt*dt) + (1-2*beta)*a0_vec)

    print("shape of u vec",np.shape(u0_vec))
    print("Number of nodes in this subdomain",nNodes)
    print("fun_T is", np.shape(fun_T))
    M_n = M.dot(fun_T)
    F_T = b + M_n

    print("M_n is", np.shape(M_n))
    print("F_T is", np.shape(F_T))

    bi   = F_T[0:(nNodes-nBnodes)]
    bg   = F_T[(nNodes-nBnodes):nNodes]



    return bi,bg

def GetInterior(nParts,nBg,nodesBG,U_gamma,U):


    for ip in range(nParts):

        npart = ip+1

        mdpath = "../../../../data/meshData/meshdim"+str(npart).zfill(4)+".dat"
        meshdim = np.genfromtxt(mdpath)
        nP = meshdim[0].astype(int)     ## For 3D
        nB = meshdim[3].astype(int)    ## For 3D
        print("Number of nodes",nP)
        print("Number of boundary nodes",nB)

        nodespath = "../../../../data/meshData/nodes"+str(npart).zfill(4)+".dat"
        nodei = np.genfromtxt(nodespath) # Nodes in the subdomain interior at top
        nodei = nodei.astype(int)

        nI = nP - nB # Number of interior nodes in the sub domain
        print("Number of interior nodes",nI)

        one = 1
        subpath = "../../data/Amats/subdom"+str(npart)
        ADiipath = subpath+"/ADii"+str(npart)+".dat"
        ADii = np.genfromtxt(ADiipath)

        ADigpath = subpath+"/ADig"+str(npart)+".dat"
        ADig = np.genfromtxt(ADigpath)

        bipath = subpath+"/bi"+str(npart)+".dat"
        bi = np.genfromtxt(bipath)

        Rpath = "../../../../data/meshData/R_"+str(npart)+".dat"
        R = np.genfromtxt(Rpath)
        print("shape of R", np.shape(R))

        Ub = R.dot(U_gamma)
        print("shape of Ub",np.shape(Ub))
        out1 = bi - ADig.dot(Ub)
        Ui = linalg.solve(ADii,out1)
        print("Interior Solution",Ui)
        npi = len(Ui)

        ######### Re arranging the nodes to form Global Solution #########

        for k in range(npi):

            nni = nodei[k]
            print("nni",nni)
            U[nni-1,:] = Ui[k]


            # elif(k < 2*nI):

            #     nni = nodei[k-nI]
            #     # print("nni",nni)
            #     U[(nPG-1)+nni,:] = Ui[k]

            # else:
            #     nni = nodei[k-(2*nI)]
            #     # print("nni",nni)
            #     U[(2*nPG-1)+nni,:] = Ui[k]


    subpath = '../../../../data/solution/'
    Uinterior_path = subpath+'interior_solu.dat'
    np.savetxt(Uinterior_path,U)
    # print("U interior is",U)


    npB = len(U_gamma)
    print("U_gamma is",U_gamma)
    print("U before putting gamma nodes is",U)

    for j in range(npB):

        tempG = nodesBG[j]
        print("tempG",tempG)
        U[tempG-1,:] = U_gamma[j]

        # elif(j < 2*nBG):

        #     tempG = nPG + nodesBG[j-nBG]
        #     # print("tempG",tempG)
        #     U[tempG-1,:] = U_gamma[j]

        # else:
        #     tempG = 2*nPG + nodesBG[j-(2*nBG)]
        #     # print("tempG",tempG)
        #     U[tempG-1,:] = U_gamma[j]

    subpath = '../../../../data/solution/'
    U_path = subpath+'final_solu.dat'
    np.savetxt(U_path,U)

    return U
#####################################################


#################################################



################### MAIN CODE STARTS ##################################

# Parameters for the Newmark beta method
beta, gamma = 0.25, 0.5


# Time stepping parameters
dt = 0.1#0.0005# timestep
T = 1
t = 0.0
nt = np.int(T/dt)

C = 1 # wave velocity

m = 2
n = 1


def dot_u(u):
    return (gamma/(beta*dt))*(u - u0) - (gamma/beta - 1.0)*v0 - dt*(gamma/(2.0*beta) - 1.0)*a0

def ddot_u(u):
    return (1.0/(beta*dt**2))*(u - u0 - dt*v0) - (1.0/(2.0*beta) - 1.0)*a0


## Path to the global meshdims.dat
meshdimpath = '../../../../data/meshData/meshdim.dat'
meshdat = np.genfromtxt(meshdimpath)
nParts = meshdat[4].astype(int)
print("number of partitions:", nParts)   ## For 2D
nPG = meshdat[0].astype(int)
print("number of points:", nPG)   ## For 2D
nBg = meshdat[3].astype(int)
print("number of Boundary nodes:", nBg)

meshdimpath = '../../../../data/meshData/boundary_nodes.dat'
Gbfile = np.genfromtxt(meshdimpath)
nodesBG = Gbfile.astype(int)

nDim = 1;


S_global = np.zeros([nBg,nBg])
g_global = np.zeros([nBg])

count = 0

while abs(t - T) > 1e-14:
    t += dt

############ DDM Based Solver ######################

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
        # u1, w = TrialFunction(V), TestFunction(V)


        if count==0:    # Initialization of fields (displacement, velocity, acceleration)
            u0, v0, a0 = Function(V), Function(V), Function(V)

            ##### Initial condition ######

            ui  = Expression(("sin(m*pi*x[0])*sin(n*pi*x[1])"),degree=1,m = 2,n = 1)

            u0 = interpolate(ui, V)

        # print("u0 vector is ", np.shape(u0.vector().array())


        #####Exact Solution ###############
        u_e = Expression(("sin(m*pi*x[0])*sin(n*pi*x[1])*cos(c*pi*sqrt(pow(m,2)+pow(n,2))*t)"),degree=1,c=C,m=m,n=n,t=0)


        print("u_initial is",u0.vector().array())

    ####################################


        kk = ip+1
        subpath = "../../data/Amats/subdom"+str(kk)

        Apath = subpath+"/A"+str(kk)+".dat"
        A = np.genfromtxt(Apath)

        Mpath = subpath+"/M"+str(kk)+".dat"
        M = np.genfromtxt(Mpath)

        bpath = subpath+"/b"+str(kk)+".dat"
        b = np.genfromtxt(bpath)

        ADiipath = subpath+"/ADii"+str(kk)+".dat"
        ADii = np.genfromtxt(ADiipath)

        ADggpath = subpath+"/ADgg"+str(kk)+".dat"
        ADgg = np.genfromtxt(ADggpath)

        ADigpath = subpath+"/ADig"+str(kk)+".dat"
        ADig = np.genfromtxt(ADigpath)

        # MDiipath = subpath+"/MDii"+str(kk)+".dat"
        # MDii = np.genfromtxt(MDiipath)

        # MDggpath = subpath+"/MDgg"+str(kk)+".dat"
        # MDgg = np.genfromtxt(MDggpath)

        # MDigpath = subpath+"/MDig"+str(kk)+".dat"
        # MDig = np.genfromtxt(MDigpath)

        # fipath = subpath+"/bi"+str(kk)+".dat"
        # bi = np.genfromtxt(fipath)

        # fgpath = subpath+"/bg"+str(kk)+".dat"
        # bg = np.genfromtxt(fgpath)

        subpath = "../../data/Amats/subdom000".rstrip("0")+str(ip+1)

        if not os.path.exists(subpath):
            os.makedirs(subpath)

        Rpath = "../../../../data/meshData/R_"+str(kk)+".dat"
        R = np.genfromtxt(Rpath)
        print("shape of R", np.shape(R))

        ADgi = ADig.T



        bi,bg = GetTransientForce(u0,v0,a0,nNodes,nBnodes,A,M,b)

        fipath = subpath+"/bi000".rstrip("0")+str(kk)+".dat"
        np.savetxt(fipath, bi)
        fgpath = subpath+"/bg000".rstrip("0")+str(kk)+".dat"
        np.savetxt(fgpath, bg)


        ###### Direct Schur Complement Solution ################

        ss1 = linalg.solve(ADii,ADig)
        S_local = ADgg - ADgi.dot(ss1)
        ss2 = S_local.dot(R)
        S_global = S_global + R.T.dot(ss2)

        gg1 = linalg.solve(ADii,bi)
        print("shape of bi", np.shape(bi))

        g_local = bg - ADgi.dot(gg1)
        print("shape of g global", np.shape(g_local))
        g_global = g_global + R.T.dot(g_local)



    u_gamma = linalg.solve(S_global,g_global)

###################### Interface Solution Found ###################
    print("shape of u_gamma is", np.shape(u_gamma))
    print("u_gamma is", u_gamma)

    subpath = '../../../../data/solution/'
    Ugamma_path = subpath+'interface_solu.dat'
    np.savetxt(Ugamma_path,u_gamma)

    ############## Interior Nodes Solution ######################

    U = np.zeros([nPG*nDim,1])
    print("Shape of U",np.shape(U))

    U = GetInterior(nParts,nBg,nodesBG,u_gamma,U)
##################################################################

    u_new = U
    print("u_new is", u_new)


#################################################################

    newmark_update(u_new, u0_vec, v0_vec, a0_vec, beta, gamma, dt)
    u_e.t = t

# Compute error at vertices
    u_exact = interpolate(u_e, V)
    print("u exact",u_exact.vector().array())
    print("unew",np.shape(u_new))
    # error = np.abs(u_exact.vector().array() - u_new)
    # print("t = %.2f: error = %.3g" % (t, error))
    file << u_new
    tseries[0:ndof,count]= u_new.compute_vertex_values()
    t_exact[0:ndof,count]= u_exact.compute_vertex_values()

    count+=1
    print(count)



print("==========================Success============================")








