Ich habe meine Cases immer aufgegliedert:
    - dynamicMeshDict verweist auf den Ordner amr, d.h. amr enthält die Inhalte die in das dynamicMeshDict gehören
    - controlDict verweist auf post, d.h. post enthält alle scripte für das postProcessing
    - der Ordner vars enthält diverse Variablen (bspw. T, dT ..) die in den verschiedenen dicts benötigt werden

Für das dynamicMeshDict verwende ich https://github.com/HenningScheufler/multiDimAMR/tree/v2012

Bisher lief der case 40s als LES mit einer CFL von 2 (.rampup), ab dann als LES mit einer CFL von 0.7 (.run) und amr

Dieser Ablauf ist in linux_local_run.sh hinterlegt.

Aktuell ist unter U/inlet timeVaryingMappedFixedValue hinterlegt, der case wird also so nicht laufen. 

Unter cfdValidation/inlet ist ein inlet hinterlegt, diesen nach constant/boundaryData/ linken, dann würde der Case probeweise laufen. Ist natürlich noch nicht der finale Inlet. 

