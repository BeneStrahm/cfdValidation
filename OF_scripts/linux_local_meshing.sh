#!/bin/bash
ncores=16 

# Create logMesh folder
mkdir logMesh

# Create logRun folder
mkdir logRun

# Execute blockMesh
blockMesh > logMesh/blockMesh.log 

# Execute checkMesh
checkMesh > logMesh/checkMesh.log  