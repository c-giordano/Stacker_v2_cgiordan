#!/bin/bash

scan () {
    npoints=$1
    cd AsiScan2D
    combineTool.py -M MultiDimFit ../workspace.root --algo grid --points $npoints --cminDefaultMinimizerStrategy 0 \
        -P r_TTTT -P r_TTT -t -1 --floatOtherPOIs 1 -n .scan_asi --setParameters r_TTTT=1,r_TTT=1 --job-mode condor --task-name scan_2D_asi --split-points 50 \
        --setParameterRanges r_TTTT=0,4:r_TTT=0,50
    cd ..
    cd ObsScan2D
    combineTool.py -M MultiDimFit ../workspace.root --algo grid --points $npoints --cminDefaultMinimizerStrategy 0 \
        -P r_TTTT -P r_TTT --floatOtherPOIs 1 -n .scan_obs --setParameters r_TTTT=1,r_TTT=1 --job-mode condor --task-name scan_2D_obs --split-points 50 \
        --setParameterRanges r_TTTT=0,4:r_TTT=0,50
    cd ..

    combineTool.py -M MultiDimFit -d workspace.root --setParameterRanges r_TTTT=0,4:r_TTT=0,50 --cminDefaultMinimizerStrategy 0  -t -1 --setParameters r_TTTT=1,r_TTT=1 -n bf_asi -P r_TTTT -P r_TTT --floatOtherPOIs 1 --job-mode condor --task-name bf_asi
    combineTool.py -M FitDiagnostics -d workspace.root --setParameterRanges r_TTTT=0,4:r_TTT=0,50 --cminDefaultMinimizerStrategy 0  -t -1 --setParameters r_TTTT=1,r_TTT=1 --saveShapes -n asi_corr --redefineSignalPOIs r_TTTT,r_TTT --job-mode condor --task-name fitdiagnostics_asi

    combineTool.py -M MultiDimFit -d workspace.root --setParameterRanges r_TTTT=0,4:r_TTT=0,50 --cminDefaultMinimizerStrategy 0  --setParameters r_TTTT=1,r_TTT=1 -n bf_obs -P r_TTTT -P r_TTT --floatOtherPOIs 1 --job-mode condor --task-name bf_obs
    combineTool.py -M FitDiagnostics -d workspace.root --setParameterRanges r_TTTT=0,4:r_TTT=0,50 --cminDefaultMinimizerStrategy 0  --setParameters r_TTTT=1,r_TTT=1 --saveShapes -n obs_corr --redefineSignalPOIs r_TTTT,r_TTT --job-mode condor --task-name fitdiagnostics_obs
}

if [ "$1" == "init" ]; then
    combineCards.py y16=DC_SM_2016.txt y17=DC_SM_2017.txt y18=DC_SM_2018.txt > DC_Combined.txt

    text2workspace.py DC_Combined.txt -o workspace.root \
    -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO 'map=.*/tttt:r_TTTT[1,-5,5]' --PO 'map=.*/tttW:r_TTT[1,-50,50]' --PO 'map=.*/tttj:r_TTT[1,-50,50]' --PO 'map=.*/ttW:r_TTW[1,-5,5]' --PO 'map=.*/ttZ:r_TTZ[1,-5,5]' --channel-masks
    mkdir -p AsiScan2D ObsScan2D

    scan 10000
fi

if [ "$1" == "collect" ]; then
    mkdir -p pics

    hadd -f asiscan_2D.root AsiScan2D/*root
    hadd -f obsscan_2D.root ObsScan2D/*root
    python3 /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/plotCorrelations.py -i fitDiagnosticsobs_corr.root:fit_s -p r_TTT,r_TTTT
    python3 /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/get.py correlationMatrix.root r_TTTT r_TTT

    python3 /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/plot2DScan.py obsscan_2D.root -o Scan2D_obs -bf higgsCombinebf_obs.MultiDimFit.mH120.root  --x r_TTTT --y r_TTT
    python3 /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/plot2DScan.py asiscan_2D.root -o Scan2D_asi -bf higgsCombinebf_asi.MultiDimFit.mH120.root  --x r_TTTT --y r_TTT
fi
