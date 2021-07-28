#!/bin/bash
# Params
# ---------
ncores=16 

# # Cleaning
# # ---------

# rm -r logRun/*
# rm -r postProcessing/

# # Ramp-Up
# # ---------

# # Copy RampUp-Dictionaries
# cp system/controlDict.rampup system/controlDict
# cp constant/dynamicMeshDict.rampup constant/dynamicMeshDict

# # Execute decomposePar
# decomposePar > logRun/decomposePar_RampUp.log

# # Execute renumberMesh
# mpirun --allow-run-as-root -np $ncores renumberMesh -parallel -overwrite > logMesh/renumberMesh_RampUp.log

# # Run Simulation 
# mpirun --allow-run-as-root -np $ncores pimpleFoam -parallel > logRun/pimpleFoam_RampUp.log

# # Execute reconstructPar
# reconstructPar > logRun/reconstructPar_RampUp.log

# # Clean-up directory
# rm -r processor*

# # Deepsave latest timestep
# latestTime=$(foamListTimes -latestTime)
# cp -r $latestTime $latestTime.orig

# Run
# ---------
# Copy run-Dictionaries
cp system/controlDict.run system/controlDict
cp constant/dynamicMeshDict.run constant/dynamicMeshDict

# Execute decomposePar
decomposePar > logRun/decomposePar_Run.log

# Execute renumberMesh
mpirun --allow-run-as-root -np $ncores renumberMesh -parallel -overwrite > logMesh/renumberMesh_Run.log

# Run Simulation 
mpirun --allow-run-as-root -np $ncores pimpleFoam -parallel > logRun/pimpleFoam_Run.log

# FoamToVTK
foamToVTK -parallel > logRun/foamToVTK.log

# Execute reconstructPar
reconstructPar > logRun/reconstructPar_Run.log

# Clean-up directory
# rm -r processor*