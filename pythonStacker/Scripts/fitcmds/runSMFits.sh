#!/bin/bash

ulimit -s unlimited

# function for text2workspace:
t2w () {
    combineTool.py -M T2W -i */*txt -o workspace.root --parallel 4\
    "-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel  --PO 'map=.*/tttt:r_TTTT[1,-5,5]' --PO 'map=.*/ttW:r_TTW[1,-5,5]' --PO 'map=.*/ttZ:r_TTZ[1,-5,5]' --channel-masks"
}

# Significances:
significance () {
    combineTool.py -M Significance -d */workspace.root --there --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 -m 125 --redefineSignalPOIs r_TTTT\
    -n .ObsSign --job-mode condor --task-name significance_obs --parallel 4
    combineTool.py -M Significance -d */workspace.root --there --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 -m 125 -t -1 --redefineSignalPOIs r_TTTT\
    -n .AsiSign --job-mode condor --task-name significance_asi --parallel 4
}

# Impacts:
impact_initial () {
    # Loop subfolders of SM:
    for i in */; do
        cd $i
        mkdir ImpactsAsi
        mkdir ImpactsObs
        echo $(pwd)
        cp workspace.root ImpactsAsi
        cp workspace.root ImpactsObs
        cd ImpactsAsi
        combineTool.py -M Impacts -d workspace.root --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1\
        --redefineSignalPOIs r_TTTT -t -1 --robustFit 1 --doInitialFit -m 125 --job-mode condor --task-name asi_initial_fit
        cd ..
        cd ImpactsObs
        combineTool.py -M Impacts -d workspace.root --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1\
        --redefineSignalPOIs r_TTTT --robustFit 1 --doInitialFit -m 125 --job-mode condor --task-name obs_initial_fit        
        cd ../..
    done
}

impact_fits () {
    for i in */; do
        cd $i/ImpactsAsi
        combineTool.py -M Impacts -d workspace.root --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1\
        --redefineSignalPOIs r_TTTT -t -1 --robustFit 1 --doFits -m 125 --job-mode condor --task-name asi_fits
        cd ..
        cd ImpactsObs
        combineTool.py -M Impacts -d workspace.root --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1\
        --redefineSignalPOIs r_TTTT --robustFit 1 --doFits -m 125 --job-mode condor --task-name obs_fits
        cd ../..
    done
}

impacts_collect() {
    for i in */; do
        cd $i/ImpactsAsi
        combineTool.py -M Impacts -d workspace.root --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1\
        -t -1 --redefineSignalPOIs r_TTTT -m 125  -o impacts.json
        plotImpacts.py -i impacts.json -o a_impacts
        cd ..
        cd ImpactsObs
        combineTool.py -M Impacts -d workspace.root --redefineSignalPOIs r_TTTT -m 125 -o impacts.json     
        plotImpacts.py -i impacts.json -o a_impacts
        cd ../..
    done
}

goodness_of_fit_simple () {
    cd run2
    mkdir gof 
    cd gof
    cp ../workspace.root ./
    combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --redefineSignalPOIs r_TTTT
    combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --setParameters r_TTTT=1,r_TTW=1.358,r_TTZ=1.099 --redefineSignalPOIs r_TTTT -s 12345 --toysFreq
    # --setParametersForEval r_TTTT=1,r_TTW=1.358,r_TTZ=1.099\


    cd ../../
}


# Goodness of fit tests
goodness_of_fit () {
    for i in */; do
        cd $i
        # if i is run2, continue
        if [ $i == "run2/" ]; then
            cd ..
            continue
        fi
        # skip 2016
        if [ $i != "2016/" ]; then
            cd ..
            continue
        fi
        echo $(pwd)
        cd gof

        rm -rfd gof_cr_only gof_sr2lOnly gof_sr3lOnly gof_sr4lOnly sr2l_sig_only sr2l_ttw_only sr_2l_np_only sr3l_sig_only sr3l_ttw_only sr3l_np_only sr4l_sig_only sr4l_ttw_only gof_cr_2l_only gof_crz_only gof_cr3l_only cr_only_true
        mkdir sr2l_sig_only sr2l_ttw_only sr_2l_np_only sr3l_sig_only sr3l_ttw_only sr3l_np_only sr4l_sig_only sr4l_ttw_only
        mkdir gof_cr_2l_only gof_crz_only gof_cr3l_only
        mkdir gof_cr_only gof_sr2lOnly gof_sr3lOnly gof_sr4lOnly
        cp ../workspace.root ./
        cp workspace.root gof_cr_only/
        cp workspace.root gof_sr2lOnly/
        cp workspace.root gof_sr3lOnly/
        cp workspace.root gof_sr4lOnly/
        cp workspace.root gof_cr_2l_only
        cp workspace.root gof_crz_only
        cp workspace.root gof_cr3l_only
        cp workspace.root sr2l_sig_only
        cp workspace.root sr2l_ttw_only
        cp workspace.root sr_2l_np_only
        cp workspace.root sr3l_sig_only
        cp workspace.root sr3l_ttw_only
        cp workspace.root sr3l_np_only
        cp workspace.root sr4l_sig_only
        cp workspace.root sr4l_ttw_only

        cd sr2l_sig_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..

        cd sr2l_ttw_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..

        cd sr_2l_np_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTZ\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTZ -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..

        cd sr3l_ttw_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_Sig=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_Sig=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..

        cd sr3l_sig_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..

        cd sr3l_np_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_ttw=1,mask_SR3L_Sig=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_ttw=1,mask_SR3L_Sig=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..

        cd sr4l_sig_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_ttw=1,mask_SR3L_Sig=1,mask_SR3L_NP=1,mask_SR4L_ttw=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_ttw=1,mask_SR3L_Sig=1,mask_SR3L_NP=1,mask_SR4L_ttw=1
        cd ..

        cd sr4l_ttw_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_ttw=1,mask_SR3L_Sig=1,mask_SR3L_NP=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTW -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR3L_ttw=1,mask_SR3L_Sig=1,mask_SR3L_NP=1,mask_SR4L_Sig=1
        cd ..

        cd gof_cr_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT\
        --setParametersForEval mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT -s 12345 --toysFreq \
        --setParametersForEval mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..
# 
        cd gof_cr_2l_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT\
        --setParametersForEval mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT -s 12345 --toysFreq \
        --setParametersForEval mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..
# 
        cd gof_crz_only


        cd gof_cr3l_only
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        cd ..


        cd gof_sr2lOnly
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT\
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
# 
        cd ..
        cd gof_sr3lOnly
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
# 
        cd ..
        cd gof_sr4lOnly
        combineTool.py -M GoodnessOfFit -d workspace.root --algo saturated -m 125 --job-mode condor --task-name gof_data\
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1
        combineTool.py -M GoodnessOfFit -d workspace.root -t 20 --algo saturated -m 125 --job-mode condor --task-name gof_asi \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT -s 12345 --toysFreq \
        --setParametersForEval mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR2L_NPClass=1,mask_SR2L_ee_Sig=1,mask_SR2L_em_Sig=1,mask_SR2L_mm_Sig=1,mask_SR2L_ee_ttw=1,mask_SR2L_em_ttw=1,mask_SR2L_mm_ttw=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1
# 
        cd ..
        cd ../../
    done
       


mask_CR2LNP=1,mask_CR2LTTW=1,mask_CR3LZ=1,mask_CR4LZ=1,mask_CR3LNP=1,mask_SR3L_Sig=1,mask_SR3L_ttw=1,mask_SR3L_NP=1,mask_SR4L_ttw=1,mask_SR4L_Sig=1
# Create a loop to run 10 times, each time generate a random number as a seed for the GOF test
        # mkdir highstat
        # cd highstat
        # for nb in {1..100}
        # do
        #     RANDOM_SEED=$(od -An -N2 -tu2 /dev/random | tr -d ' ' | awk '{print $1 % 100000}')
        #     combineTool.py -M GoodnessOfFit -d workspace.root -t 10 --algo saturated -m 125 --job-mode condor --task-name "gof_asi_$RANDOM_SEED" \
        #     --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT -s $RANDOM_SEED --toysFreq
        # done
        # cd ..
}

collect_gof () {
    for i in */; do
        cd $i/gof
        # combineTool.py -M CollectGoodnessOfFit --input higgsCombine*.root -m 125 -o gof.json
        # plotGof.py gof.json --statistic saturated --mass 125.0 -o gof_plot 
        #--title-right="my label"

        # cd highstat
        # hadd -f higgsCombine.Test.GoodnessOfFit.mH125.coll.root higgsCombineTest.GoodnessOfFit.mH125.*.root
        # cp ../higgsCombine.Test.GoodnessOfFit.mH125.root ./
        # combineTool.py -M CollectGoodnessOfFit --input higgsCombine.Test.GoodnessOfFit.mH125.root higgsCombine.Test.GoodnessOfFit.mH125.coll.root -m 125 -o gof.json
        # plotGof.py gof.json --statistic saturated --mass 125.0 -o gof_plot 
        # cd ../..
        # Loop all subdirectories:
        for j in */; do
            cd $j
            combineTool.py -M CollectGoodnessOfFit --input higgsCombine*.root -m 125 -o gof.json
            plotGof.py gof.json --statistic saturated --mass 125.0 -o gof_plot 
            cd ..
        done

        cd ../..
    done
}

fit_diagnostics () {
    for i in */; do
        cd $i
        combineTool.py -M FitDiagnostics -d workspace.root -m 125 --job-mode condor --task-name fitDiagnostics \
        --setParameters r_TTTT=1,r_TTW=1,r_TTZ=1 --redefineSignalPOIs r_TTTT
        cd ..
    done
}


# Check if it's the first time running by checking if the directory exists or check if a the first argument is "init":
if [ "$1" == "init" ]; then
    mkdir -p 2016 2017 2018 run2
    cp DC_SM_2016* 2016/
    cp DC_SM_2017* 2017/
    cp DC_SM_2018* 2018/
    cp DC_SM_201*root run2/
    combineCards.py y16=DC_SM_2016.txt y17=DC_SM_2017.txt y18=DC_SM_2018.txt > run2/DC_Combined.txt

    t2w
    significance
    impact_initial
    fit_diagnostics
fi

# check if the first argument is fits and that initial fit file exists:
if [ "$1" == "fits" ]; then
    impact_fits
fi

if [ "$1" == "collect" ]; then
    impacts_collect
fi
