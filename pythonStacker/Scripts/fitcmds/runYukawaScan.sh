#!/bin/bash

set -e

# define t2w command in a function:
t2w () {

    text2workspace.py $1 -o $2 -m 125 -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  \
        --PO verbose \
        --PO 'map=.*/NOTHING:at[1,-2.0,2.0]' \
        --PO 'map=.*/NOTHING:bt[0,-2.0,2.0]' \
        --PO 'map=.*/ttH:r_TTH=expr::r_TTH("at**2 + 0.46 * bt**2", at, bt)' \
        --PO 'map=.*/tttt:r_TTTT=expr::r_TTTT("1.0091 - 0.1249 * at**2 + 0.2931 * bt**2 + 0.2843 * at**2 * bt**2 + 0.1158 * at**4 + 0.1887 * bt**4", at, bt)' \
        --PO 'map=.*/tttj:r_TTTJ=expr::r_TTTJ("1.3720 - 1.1026 * at - 0.0016 * bt + 1.2721 * at**2 + 0.0480 * bt**2 + 0.0317 * at * bt - 0.7346 * at**3 + 0.0006 * bt**3 - 0.0332 * at**2 * bt - 0.2273 * at * bt**2 + 0.1931 * at**4 + 0.1659 * bt**4 + 0.4923 * at**2 * bt**2 + 0.0013 * at**3 * bt + 0.0013 * at * bt**3", at, bt)' \
        --PO 'map=.*/tttW:r_TTTW=expr::r_TTTW("1.1339 - 0.1499 * at - 0.0012 * bt - 0.0844 * at**2 + 0.0508 * bt**2 + 0.0004 * at * bt - 0.0217 * at**3 + 0.0005 * bt**3 + 0.0002 * at**2 * bt - 0.0153 * at * bt**2 + 0.1221 * at**4 + 0.1486 * bt**4 + 0.2962 * at**2 * bt**2 - 0.0010 * at**3 * bt - 0.0010 * at * bt**3", at, bt)' \
        --PO 'map=.*/ttW:r_TTW[1,-5,5]' \
        --PO 'map=.*/ttZ:r_TTZ[1,-5,5]'
}

fitting_asimov () {
    # mkdir AsiScan2D
    npoints_twodim=$2
    npoints_onedim=$3

    cd AsiScan2D
    combineTool.py -M MultiDimFit ../$1 --algo grid --points $npoints_twodim --cminDefaultMinimizerStrategy 0  -P at -P bt -t -1 --floatOtherPOIs 1 -n .scan_asi --setParameters at=1,bt=0 --job-mode condor --task-name scan_2D_asi --split-points 20
    cd ..

    combineTool.py -M MultiDimFit -d $1 --cminDefaultMinimizerStrategy 0  -t -1 --setParameters at=1,bt=0 -n bf_asi -P at -P bt --floatOtherPOIs 1 --job-mode condor --task-name bf_asi
    combineTool.py -M FitDiagnostics -d $1 --cminDefaultMinimizerStrategy 0  -t -1 --setParameters at=1,bt=0 --saveShapes -n asi_corr --job-mode condor --task-name fitdiagnostics_asi

    # # hadd results, then plot
    # python plotCorrelations.py -i fitDiagnosticsasi_corr.root:fit_s -p at,bt
    # python getCorrMatrix.py correlationMatrix.root at bt
    # python plot2DScan.py scan_2D_asi_limsyst.root 

    # mkdir AsiScan1D
    cd AsiScan1D
    combineTool.py -M MultiDimFit ../$1 --algo grid --points $npoints_onedim --cminDefaultMinimizerStrategy 0  -P at -t -1 --floatOtherPOIs 1 -n .at.scan_asi --setParameters at=1,bt=0 --job-mode condor --task-name at_scan_1D_asi --split-points 5
    combineTool.py -M MultiDimFit ../$1 --algo grid --points $npoints_onedim --cminDefaultMinimizerStrategy 0  -P bt -t -1 --floatOtherPOIs 1 -n .bt.scan_asi --setParameters at=1,bt=0 --job-mode condor --task-name bt_scan_1D_asi --split-points 5
    cd ..
}

fitting_observed() {
    # mkdir ObsScan2D

    npoints_twodim=$2
    npoints_onedim=$3

    cd ObsScan2D
    combineTool.py -M MultiDimFit ../$1 --algo grid --points $npoints_twodim --cminDefaultMinimizerStrategy 0  -P at -P bt --floatOtherPOIs 1 -n .scan_obs --setParameters at=1,bt=0 --job-mode condor --task-name scan_2D_obs --split-points 20
    cd ..

    combineTool.py -M MultiDimFit -d $1 --cminDefaultMinimizerStrategy 0 --setParameters at=1,bt=0 -n bf_obs -P at -P bt --floatOtherPOIs 1 --job-mode condor --task-name bf_obs
    combineTool.py -M FitDiagnostics -d $1 --cminDefaultMinimizerStrategy 0 --setParameters at=1,bt=0 --saveShapes -n obs_corr --job-mode condor --task-name fitdiagnostics_obs
    combineTool.py -M FitDiagnostics -d $1 --cminDefaultMinimizerStrategy 0 --redefineSignalPOIs at,bt --robustFit 1 --setParameters at=1,bt=0 --saveShapes -n obsPF --job-mode condor --task-name fitdiagnostics_obs

    # # hadd results, then plot
    # python plotCorrelations.py -i fitDiagnosticsasi_corr.root:fit_s -p at,bt
    # python getCorrMatrix.py correlationMatrix.root at bt
    # python plot2DScan.py scan_2D_asi_limsyst.root 

    # mkdir ObsScan1D
    cd ObsScan1D
    combineTool.py -M MultiDimFit ../$1 --algo grid --points $npoints_onedim --cminDefaultMinimizerStrategy 0  -P at --floatOtherPOIs 1 -n .at.scan_obs --setParameters at=1,bt=0 --job-mode condor --task-name at_scan_1D_obs --split-points 5
    combineTool.py -M MultiDimFit ../$1 --algo grid --points $npoints_onedim --cminDefaultMinimizerStrategy 0  -P bt --floatOtherPOIs 1 -n .bt.scan_obs --setParameters at=1,bt=0 --job-mode condor --task-name bt_scan_1D_obs --split-points 5
    cd ..
}

postprocess_scans () {
    mkdir -p pics

    hadd obs_at_scan_1D.root ObsScan1D/higgsCombine.at.scan_obs.POINTS.*
    hadd asi_at_scan_1D.root AsiScan1D/higgsCombine.at.scan_asi.POINTS.*
    hadd obs_bt_scan_1D.root ObsScan1D/higgsCombine.bt.scan_obs.POINTS.*
    hadd asi_bt_scan_1D.root AsiScan1D/higgsCombine.bt.scan_asi.POINTS.*
    hadd obs_scan_2D.root ObsScan2D/higgsCombine.scan_obs.POINTS.*
    hadd asi_scan_2D.root AsiScan2D/higgsCombine.scan_asi.POINTS.*

    python3 /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/plot1DScan.py obs_at_scan_1D.root -o pics/Scan1D_at --POI at --translate /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/POITranslate.json --others asi_at_scan_1D.root:Asimov:2
    python3 /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/plot1DScan.py obs_bt_scan_1D.root -o pics/Scan1D_bt --POI bt --translate /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/POITranslate.json --others asi_bt_scan_1D.root:Asimov:2
    python3 /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/plot2DScan.py asi_scan_2D.root -o Scan2D_asi -bf higgsCombinebf_asi.MultiDimFit.mH120.root -d pics/Scan2D_asi.root
    python3 /user/nivanden/plots/pythonStacker/Scripts/plotsandprocess/plot2DScan.py obs_scan_2D.root -o Scan2D_obs -bf higgsCombinebf_obs.MultiDimFit.mH120.root -d pics/Scan2D_obs.root
}

impacts_initial () {
    mkdir ImpactsAsi ImpactsObs
    cd ImpactsAsi
    cp ../workspace.root .
    combineTool.py -M Impacts -d workspace.root --setParameters r_TTW=1,r_TTZ=1,at=1,bt=0\
    --redefineSignalPOIs at,bt -t -1 --robustFit 1 --doInitialFit -m 125 --job-mode condor --task-name asi_initial_fit
    cd ../ImpactsObs
    cp ../workspace.root .
    combineTool.py -M Impacts -d workspace.root --setParameters r_TTW=1,r_TTZ=1,at=1,bt=0\
    --redefineSignalPOIs at,bt --robustFit 1 --doInitialFit -m 125 --job-mode condor --task-name obs_initial_fit
    cd ..
}

impacts_fits () {
    cd ImpactsAsi
    combineTool.py -M Impacts -d workspace.root --setParameters r_TTW=1,r_TTZ=1,at=1,bt=0\
    --redefineSignalPOIs at,bt -t -1 --robustFit 1 --doFits -m 125 --job-mode condor --task-name asi_initial_fit
    cd ../ImpactsObs
    combineTool.py -M Impacts -d workspace.root --setParameters r_TTW=1,r_TTZ=1\
    --redefineSignalPOIs at,bt --robustFit 1 --doFits -m 125 --job-mode condor --task-name obs_initial_fit
    cd ..
}

impacts_collect () {
    cd ImpactsAsi
    combineTool.py -M Impacts -d workspace.root --redefineSignalPOIs at,bt -t -1 -m 125 -o impacts.json
    plotImpacts.py -i impacts.json -o at_impacts --POI at
    plotImpacts.py -i impacts.json -o bt_impacts --POI bt

    cd ../ImpactsObs
    combineTool.py -M Impacts -d workspace.root --redefineSignalPOIs at,bt -m 125 -o impacts.json
    plotImpacts.py -i impacts.json -o at_impacts --POI at
    plotImpacts.py -i impacts.json -o bt_impacts --POI bt
    cd ..
}

gof () {
    cd gof
    # combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data \
    # --setParameters r_TTW=1,r_TTZ=1,at=1,bt=0 --redefineSignalPOIs at,bt
    # Create a loop to run 10 times, each time generate a random number as a seed for the GOF test
    mkdir highstat
    cd highstat
    for i in {1..100}
    do
        RANDOM_SEED=$(od -An -N2 -tu2 /dev/random | tr -d ' ' | awk '{print $1 % 100000}')
        combineTool.py -M GoodnessOfFit -d workspace.root -t 10 --algo saturated -m 125 --job-mode condor --task-name "gof_asi_$RANDOM_SEED" \
        --setParameters r_TTW=1,r_TTZ=1,at=1,bt=0 --redefineSignalPOIs at,bt -s $RANDOM_SEED --toysFreq
    done
    cd ..
    # combineTool.py -M GoodnessOfFit -d workspace.root -t 500 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
    # --setParameters r_TTW=1,r_TTZ=1,at=1,bt=0 --redefineSignalPOIs at,bt -s 12345 --toysFreq
    cd ..
}

WORKDIR=$PWD

cd /user/nivanden/CMSSW_14_1_0_pre4/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
ulimit -s unlimited

cd $WORKDIR



# Check if it's the first time running by checking if the directory exists or check if a the first argument is "init":
if [ "$1" == "init" ]; then
    combineCards.py y16=DC_SM_2016.txt y17=DC_SM_2017.txt y18=DC_SM_2018.txt > DC_Combined.txt

    t2w DC_Combined.txt workspace.root
    mkdir AsiScan2D AsiScan1D ObsScan2D ObsScan1D
    fitting_asimov workspace.root 40000 500
    fitting_observed workspace.root 40000 500
    impacts_initial
fi

# check if the first argument is fits and that initial fit file exists:
if [ "$1" == "fits" ]; then
    impacts_fits
fi

if [ "$1" == "collect" ]; then
    impacts_collect
fi

if [ "$1" == "scans" ]; then
    postprocess_scans
fi

if [ "$1" == "redo_scans" ]; then
    # clean the directories:
    rm AsiScan2D/*
    rm AsiScan1D/*
    rm ObsScan2D/*
    rm ObsScan1D/*

    fitting_asimov workspace.root $2 $3
    fitting_observed workspace.root $2 $3
fi
