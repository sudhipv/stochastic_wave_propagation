#
# Copyright (C) 2017 Ajit Desai, PhD Candidate, ajit.ndesai@gmail.com
#
# Code to convert mesh data in points and cells format to dolfin-xml format
# Using meshio and lxml module https://github.com/nschloe/meshio

# First added:  2007-10-08
# meshio can be directly use to convert any mesh: meshio-convert input.msh output.vtu

# EDITED BY SUDHI - MAY 24 2020

import numpy as np
import meshio as ms                   ## pip3 install meshio
import lxml                           ## pip3 install lxml

print("========================================================")
print("Convert GMSH msh-mesh to FEniCS xml-mesh")


## Path to the global meshdims.dat       ## Fortran Decomposed
nParts = np.genfromtxt("../../../../data/meshData/meshdim.dat")
nParts = nParts[4].astype(int)
print("number of partitions:", nParts)

for i in range(nParts):

    ## Path to the global points*.dat    ## Fortran Decomposed
    ppath = "../../../../data/meshData/points"+str(i+1).zfill(4)+".dat"
    pts = np.genfromtxt(ppath)
    nPoints = len(pts)   ## Number of points (used for 2D)
    points = np.c_[pts, np.zeros(nPoints)]  ## adding zeros in 3rd column for 2D mesh

    if (i==0):
        print("Input GMSH: points:    ", ppath)

    ## Path to the global triangles*.dat ## Fortran Decomposed
    rpath = "../../../../data/meshData/triangles"+str(i+1).zfill(4)+".dat" ## Global
    elements = np.genfromtxt(rpath)
    cells = elements.astype(int)       ## converted all elements to integer

    cells = cells - 1                  ## index start with zero in python
    cells = {'triangle': cells}        ## meshio required input cells format

    if (i==0):
        print("Input GMSH: triangles: ", rpath)

    mpath = "../../data/foo"+str(i+1)+".xml" ## output file path
    mesh = ms.Mesh(points,cells)
    ms.write(mpath,mesh)   ## write mesh.xml for dolfin
    if (i==0):
        print("Output FEniCS-XML mesh:" , mpath)

print("==========================Success============================")
