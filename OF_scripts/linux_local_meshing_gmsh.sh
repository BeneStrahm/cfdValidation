#!/bin/bash
ncores=16 

# Create logMesh folder
mkdir logMesh

# Create logRun folder
mkdir logRun

# Converting the Mesh to OpenFOAM
gmshToFoam gmsh.msh > logMesh/gmshToFoam.log 

# Changing boundary entries to wall
changeDictionary -constant > logMesh/changeDictionary.log 

# Execute checkMesh
checkMesh > logMesh/checkMesh.log  