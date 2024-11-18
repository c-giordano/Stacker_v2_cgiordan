#!/bin/bash
# Script to generate datacards for each BSM model

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /user/nivanden/CMSSW_13_3_3/src
eval `scram runtime -sh`
export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)

cd /user/nivanden/plots/pythonStacker/

CreateProcessFile () {
    cp "settingfiles/Process/BSM_template.json" "settingfiles/Process/bsm_var/${1}.json"
    # replace BSMMODEL with the proper bsmmodel name,
    # replace BSMNAME with the pretty name for this model
    sed -i "s/BSMMODEL/${1}/g" "settingfiles/Process/bsm_var/${1}.json"
    sed -i "s/BSMNAME/${2}/g" "settingfiles/Process/bsm_var/${1}.json"
}

RunDatacardCreation () {
    # args: Base model, Mass point, coupling strength, YEAR, BSMNAME
    BSMMODEL="${1}_M${2}_C${3}"
    CreateProcessFile $BSMMODEL $3
    python3 buildDatacard.py -vf settingfiles/Variables/base.json \
                             -sf settingfiles/Uncertainties/2018.json \
                             -pf settingfiles/Process/bsm_var/${BSMMODEL}.json \
                             -cf settingfiles/Channel/all_channels.json \
                             -y ${2} \
                             -dcf settingfiles/Datacards/${2}_BSM.json \
                             -op output/datacards/tmp/
    # move the datacard to the proper directory
    mkdir -p "bsmlimits/${1}/M${2}_C${3}"
    mv "output/datacards/tmp/DC_${2}_BSM.txt" "bsmlimits/${1}/M${2}_C${3}/"
    mv "output/datacards/tmp/DC_${2}_BSM.root" "bsmlimits/${1}/M${2}_C${3}/"
}

# check if the directory bsm_var exists:
if [ ! -d "settingfiles/Process/bsm_var" ]
then
    mkdir "settingfiles/Process/bsm_var"
fi

# Build directory structure
if [ ! -d "bsmlimits" ]
then
    mkdir "bsmlimits"
fi
if [ ! -d "bsmlimits/TopPhilicScalarSinglet" ]
then
    mkdir "bsmlimits/TopPhilicScalarSinglet"
fi
if [ ! -d "bsmlimits/TopPhilicPseudoScalarSinglet" ]
then
    mkdir "bsmlimits/TopPhilicPseudoScalarSinglet"
fi
if [ ! -d "bsmlimits/TopPhilicVectorSinglet" ]
then
    mkdir "bsmlimits/TopPhilicVectorSinglet"
fi
if [ ! -d "bsmlimits/TopPhilicScalarOctet" ]
then
    mkdir "bsmlimits/TopPhilicScalarOctet"
fi
if [ ! -d "bsmlimits/TopPhilicPseudoScalarOctet" ]
then
    mkdir "bsmlimits/TopPhilicPseudoScalarOctet"
fi
if [ ! -d "bsmlimits/TopPhilicVectorOctet" ]
then
    mkdir "bsmlimits/TopPhilicVectorOctet"
fi

# create datacards for each BSM model
# Output should be specified in a dcf template?