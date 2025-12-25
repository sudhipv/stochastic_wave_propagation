


### Acoustic Wave Propagation Problem - Serial ######

##### Monte Carlo Simulation using Fortran Decomposed Mesh and Fenics Assembly

#### Copyright (C) Sudhi P V #####

### June -30 - 2020





from dolfin import *
import numpy as np
import time
import os as os
from scipy import linalg
import numpy.random as npr


## Initialize MPI for python
# comm = MPI.COMM_WORLD
# ip = comm.Get_rank()
# print('My rank is ',ip)

if ip==0:
    print("==========================================================")
    print("Running FEniCS in Parallel for 2D stochastic Acoustic Assembly...")

##### VERY IMPORTANT - TO AVOID FENICS FROM RE ORDERING THE NODES
parameters['reorder_dofs_serial'] = False

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


def newmark_update(u, u0, v0, a0, beta, gamma, dt):

    u_vec, u0_vec = u.vector(), u0.vector()
    v0_vec, a0_vec = v0.vector(), a0.vector()

 # Update acceleration and velocity
    a_vec = (1.0/(2.0*beta))*( (u_vec - u0_vec -
    v0_vec*dt)/(0.5*dt*dt) - (1.0-2.0*beta)*a0_vec )

    v_vec = dt*((1.0-gamma)*a0_vec + gamma*a_vec) + v0_vec
    v0.vector()[:], a0.vector()[:] = v_vec, a_vec
    u0.vector()[:] = u.vector()



###############################    Main Code Starts ##################################################################


T = 0.5;
dt = 0.01;
beta = 0.25;
gamma = 0.5;


# Reading Mesh information
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "file.h5", "r")
hdf.read(mesh, "/mesh", False)
cd = CellFunction("size_t", mesh)
hdf.read(cd, "/cd")
fd = FacetFunction("size_t", mesh)
hdf.read(fd, "/fd")

ndof = mesh.num_vertices()
print(ndof)

# Function space definition over the mesh
V = FunctionSpace(mesh, "CG", 1)   #Continuous Galerkin for the displacement field

# Test and Trial function
u1, w = TrialFunction(V), TestFunction(V)

# Initialization of fields (displacement, velocity, acceleration)
u0, v0, a0 = Function(V), Function(V), Function(V)


def dot_u(u):
  return (gamma/(beta*dt))*(u - u0) - (gamma/beta - 1.0)*v0 - dt*(gamma/(2.0*beta) - 1.0)*a0

def ddot_u(u):
  return (1.0/(beta*dt**2))*(u - u0 - dt*v0) - (1.0/(2.0*beta) - 1.0)*a0


  #Initial conditions
f  = Constant((0.0))

# ui  = Expression(("1*exp(-50*(pow(x[0]-0.5,2)+pow(x[1]-0.5,2)))"),degree=2)


##### Initial condition ######

ui  = Expression(("sin(m*pi*x[0])*sin(n*pi*x[1])"),degree=1,m = 2,n = 1)

u0 = interpolate(ui, V)

v0 = interpolate(Constant(0.0), V)

a0 = interpolate(Constant(0.0), V)

c = 1 # wave velocity

ndim = 2

nSamples = npr.randn(ndim)

print(nSamples)

# for k in range(nSamples):


#     params = ([k, dim_L, sIndex, mIndex])

#     ## Invoke variation fomrulaiton for each sample
#     a,m,l = stoVartional(params, u1, w, f)

#     ## Invoke Deterministic Assembly procedure for each sample
#     A,M,b = detAssembly(a,m,l)     ## New Assembly


#     ### Damped
#     C = 0.5 * M + 0.01 * A;


#     K_t = (M/(beta*dt*dt))  + A


#     u = Function(V)

#     tseries = np.zeros([ndof,nt+1])

#     while count < nt+1:

#        t += dt
#        F_t = b + M * ( ((u0.vector() + dt*v0.vector()) /(beta*dt*dt)) + ((1.0-2.0*beta)*a0.vector()/(2*beta)) )
#        solve(K_t,u.vector(),F_t)
#        newmark_update(u, u0, v0, a0, beta, gamma, dt)

#        file << u
#        tseries[0:ndof,count]= u.compute_vertex_values()

#        count+=1
#        print(count)

#     #######################################

#     # print(shape(tseries))
#     np.savetxt('./u_tseries.dat', tseries)
#     np.savetxt('./u_exact_tseries.dat', t_exact)



















