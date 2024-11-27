#!/bin/bash
# Script to generate datacards for each BSM model
echo "Sourcing CMSSW environment..."
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /user/cgiordan/CMSSW_13_3_3/src
eval `scram runtime -sh`
export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
echo "Environment setup complete!"

cd /user/cgiordan/CMSSW_13_3_3/src/plots/pythonStacker/
echo "Navigated to $(pwd)..."

CreateProcessFile () {
    echo "Creating process file for model: ${1}, pretty name: ${2}"
    cp "settingfiles/Process/BSM_template.json" "settingfiles/Process/bsm_var/${1}.json"
    # replace BSMMODEL with the proper bsmmodel name,
    # replace BSMNAME with the pretty name for this model
    sed -i "s/BSMMODEL/${1}/g" "settingfiles/Process/bsm_var/${1}.json"
    sed -i "s/BSMNAME/${2}/g" "settingfiles/Process/bsm_var/${1}.json"
}


RunDatacardCreation () {
    # args: Base model, Mass point, Year
    BSMMODEL="TopPhilic${1}_M${2}"
    echo "Running datacard creation for BSM model: $BSMMODEL, Year: $3, Mass: $4 GeV"
    CreateProcessFile $BSMMODEL $3
    python3 buildDatacard.py -vf settingfiles/Variables/base.json \
                             -sf settingfiles/Uncertainties/${3}.json \
                             -pf settingfiles/Process/bsm_var/${BSMMODEL}.json \
                             -cf settingfiles/Channel/all_channels.json \
                             -y ${3} \
                             -dcf settingfiles/Datacards/${3}_full.json \
                             -op output/datacards/tmp/ \
                             --BSM \
                             --storage Intermediate_VS_v2 
    echo "Moving datacard to output directory: output/bsmlimits/TopPhilic${1}_final/${4}/"
    # move the datacard to the proper directory
    mkdir -p "output/bsmlimits/TopPhilic${1}_final/${4}"
    mv "output/datacards/tmp/DC_${3}.txt" "output/bsmlimits/TopPhilic${1}_final/${4}/"
    mv "output/datacards/tmp/DC_${3}.root" "output/bsmlimits/TopPhilic${1}_final/${4}/"
    echo "Datacard for $BSMMODEL (Mass: $4 GeV) created and moved successfully!!!"
}

# check if the directory bsm_var exists:
if [ ! -d "settingfiles/Process/bsm_var" ]
then
    mkdir "settingfiles/Process/bsm_var"
fi

# Build directory structure
if [ ! -d "output/bsmlimits" ]
then
    mkdir "output/bsmlimits"
fi
# if [ ! -d "output/bsmlimits/TopPhilicScalarSinglet_attempt" ]
# then
#     mkdir "output/bsmlimits/TopPhilicScalarSinglet_attempt"
# fi
# if [ ! -d "output/bsmlimits/TopPhilicPseudoScalarSinglet_final" ]
# then
#     mkdir "output/bsmlimits/TopPhilicPseudoScalarSinglet_final"
# fi
if [ ! -d "output/bsmlimits/TopPhilicVectorSinglet_final" ]
then
    mkdir "output/bsmlimits/TopPhilicVectorSinglet_final"
fi
# if [ ! -d "output/bsmlimits/TopPhilicScalarOctet_final" ]
# then
#     mkdir "output/bsmlimits/TopPhilicScalarOctet_final"
# fi
# if [ ! -d "output/bsmlimits/TopPhilicPseudoScalarOctet_final" ]
# then
#     mkdir "output/bsmlimits/TopPhilicPseudoScalarOctet_final"
# fi
# if [ ! -d "output/bsmlimits/TopPhilicVectorOctet_final" ]
# then
#     mkdir "output/bsmlimits/TopPhilicVectorOctet_final"
# fi

# RunDatacardCreation VectorOctet 0p4 2018 400
# RunDatacardCreation VectorOctet 0p6 2018 600
# RunDatacardCreation VectorOctet 0p8 2018 800
# RunDatacardCreation VectorOctet 1p0 2018 1000
# RunDatacardCreation VectorOctet 1p5 2018 1500

# RunDatacardCreation VectorOctet 0p4 2017 400
# RunDatacardCreation VectorOctet 0p6 2017 600
# RunDatacardCreation VectorOctet 0p8 2017 800
# RunDatacardCreation VectorOctet 1p0 2017 1000
# RunDatacardCreation VectorOctet 1p5 2017 1500

# RunDatacardCreation VectorOctet 0p4 2016 400
# RunDatacardCreation VectorOctet 0p6 2016 600
# RunDatacardCreation VectorOctet 0p8 2016 800
# RunDatacardCreation VectorOctet 1p0 2016 1000
# RunDatacardCreation VectorOctet 1p5 2016 1500

RunDatacardCreation VectorSinglet 0p4 2018 400
RunDatacardCreation VectorSinglet 0p6 2018 600
RunDatacardCreation VectorSinglet 0p8 2018 800
RunDatacardCreation VectorSinglet 1p0 2018 1000
RunDatacardCreation VectorSinglet 1p5 2018 1500

RunDatacardCreation VectorSinglet 0p4 2017 400
RunDatacardCreation VectorSinglet 0p6 2017 600
RunDatacardCreation VectorSinglet 0p8 2017 800
RunDatacardCreation VectorSinglet 1p0 2017 1000
RunDatacardCreation VectorSinglet 1p5 2017 1500

RunDatacardCreation VectorSinglet 0p4 2016 400
RunDatacardCreation VectorSinglet 0p6 2016 600
RunDatacardCreation VectorSinglet 0p8 2016 800
RunDatacardCreation VectorSinglet 1p0 2016 1000
RunDatacardCreation VectorSinglet 1p5 2016 1500

#models_complete=("TopPhilicScalarSinglet" "TopPhilicPseudoScalarSinglet_New" "TopPhilicVectorSinglet_New" \
#        "TopPhilicScalarOctet" "TopPhilicPseudoScalarOctet" "TopPhilicVectorOctet")


#if [ ! -d "output/bsmlimits" ]; then
#    mkdir "output/bsmlimits"
#fi

#for model in "${models_complete[@]}"; do
#    dir="output/bsmlimits/$model"
#    if [ ! -d "$dir" ]; then
#        mkdir -p "$dir"
#        echo "'$dir' created"
#    else
#        echo "'$dir' already exists"
#    fi
#done

#models=("ScalarSinglet")
#mass_points=("0p4" "0p6" "0p8" "1p0" "1p5")
#masses=("400" "600" "800" "1000" "1500")
#years=("2018" "2017" "2016") 

#models=("VectorSinglet")
#mass_points=("0p4")
#masses=("400")
#years=("2018")

#for model in "${models[@]}"; do
#    for year in "${years[@]}"; do
#        for i in "${!mass_points[@]}"; do
#            # this is for debug!
#            echo "RunDatacardCreation $model ${mass_points[$i]} $year ${masses[$i]}"

            # model , mass point, year, mass
#            RunDatacardCreation $model ${mass_points[$i]} $year ${masses[$i]}
#        done
#    done
#done
echo "Datacard creation process complete."
