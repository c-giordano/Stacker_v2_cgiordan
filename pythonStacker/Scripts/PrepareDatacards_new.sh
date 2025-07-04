#!/bin/bash
# Script to generate a single datacard based on arguments from Condor

echo "Sourcing CMSSW environment..."
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /user/cgiordan/CMSSW_13_3_3/src
eval `scram runtime -sh`
export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
echo "Environment setup complete!"

cd /user/cgiordan/CMSSW_13_3_3/src/plots/pythonStacker/
echo "Navigated to $(pwd)..."

MODEL=$1   # First argument (e.g., VectorSinglet)
MASS_POINT=$2 # Second argument (e.g., 0p4)
YEAR=$3       # Third argument (e.g., 2018)
MASS_GEV=$4   # Fourth argument (e.g., 400)
WITH_DATA=$5    # Fifth argument (optional, --data or empty)

CreateProcessFile () {
    local model=$1
    local year=$2
    
    echo "Creating process file for model: ${model}, for the year ${year}"
    cp "settingfiles/Process/BSM_template.json" "settingfiles/Process/bsm_var/${model}.json"
    # replace BSMMODEL with the proper bsmmodel name,
    # replace BSMNAME with the pretty name for this model
    sed -i "s/BSMMODEL/${model}/g" "settingfiles/Process/bsm_var/${model}.json"
    sed -i "s/BSMNAME/${year}/g" "settingfiles/Process/bsm_var/${model}.json"
}


RunDatacardCreation() {
    local model=$1
    local mass=$2
    local year=$3
    local mass_gev=$4
    local with_data=$5

    BSMMODEL="TopPhilic${model}_M${mass}"

    echo "Running datacard creation for BSM model: $BSMMODEL, Year: $year, Mass: $mass_gev GeV"
    
    CreateProcessFile $BSMMODEL $year

    if [[ "$with_data" == "--data" ]]; then
        TMPDIR="output/datacards/tmp/${model}/${mass}/data/"
    else
        TMPDIR="output/datacards/tmp/${model}/${mass}/Asimov/"
    fi

    # TMPDIR="output/datacards/tmp/${model}/${mass}"
    mkdir -p "$TMPDIR"
    echo "Tmp directory for this mass point created : $TMPDIR"

    CMD="python3 buildDatacard.py -vf settingfiles/Variables/base.json \
                                 -sf settingfiles/Uncertainties/${year}_ext_BSM.json \
                                 -pf settingfiles/Process/bsm_var/${BSMMODEL}.json \
                                 -cf settingfiles/Channel/all_channels.json \
                                 -y ${year} \
                                 -dcf settingfiles/Datacards/${year}_full_attempt.json \
                                 -op ${TMPDIR} \
                                 --BSM \
                                 --storage Intermediate_pseudoTotal"

    if [[ -n "$with_data" ]]; then
        CMD+=" --data"
    fi

    echo "Running command: $CMD"
    eval $CMD

    CopyDatacards $model $mass_gev $year $with_data
    
    # mkdir -p "output/bsmlimits/TopPhilic${model}_oldBinning_noSplit_smallDiffs_data/${mass_gev}"

    # ls -l "${TMPDIR}/DC_${year}.txt" "${TMPDIR}/DC_${year}.root"

    # echo "Moving datacard to output directory: output/bsmlimits/TopPhilic${model}_oldBinning_noSplit_smallDiffs_data/${mass_gev}/"
    
    # mv "${TMPDIR}/DC_${year}.txt" "output/bsmlimits/TopPhilic${model}_oldBinning_noSplit_smallDiffs_data/${mass_gev}/"
    # mv "${TMPDIR}/DC_${year}.root" "output/bsmlimits/TopPhilic${model}_oldBinning_noSplit_smallDiffs_data/${mass_gev}/"
    # echo "Datacard for $BSMMODEL (Mass: $mass_gev GeV) created and moved successfully!!!"
}

CopyDatacards() {
    local model=$1
    local mass_gev=$2
    local year=$3
    local with_data=$4

    if [[ "$with_data" == "--data" ]]; then
        OUTDIR="output/bsmlimits/TopPhilic${model}_data_newUnc/${mass_gev}/"
    else
        OUTDIR="output/bsmlimits/TopPhilic${model}_asimov_newUnc/${mass_gev}/"
    fi

    mkdir -p "$OUTDIR"
    echo "Moving datacards to output directory: $OUTDIR"


    mv "${TMPDIR}/DC_${year}.txt" "$OUTDIR"
    mv "${TMPDIR}/DC_${year}.root" "$OUTDIR"

    echo "Datacards for $BSMMODEL (Mass: $mass_gev GeV | Year: $year) moved successfully!"
}


RunDatacardCreation $MODEL $MASS_POINT $YEAR $MASS_GEV $WITH_DATA

# RunDatacardCreation PseudoScalarSinglet 1p2 2018 1200 --data