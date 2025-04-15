

version=$1
dir=/user/nivanden/plots/pythonStacker

set -e

# check if directory exists in output:
if [ -d "${dir}/output/${version}" ]; then
    echo "Output directory already exists. Please remove it before running the script."
    # exit 1
fi

init () {
    mkdir -p "${dir}/output/${version}"
    cd "${dir}/output/${version}"

    mkdir -p SM SM_ttt_tttt Yukawa Yukawa_variations scan_ttt_tttt
    cp ../datacards/DC_SM_2018.* ./
    cp ../datacards/DC_SM_2017.* ./
    cp ../datacards/DC_SM_2016.* ./
}

standardmodel_fit () {
    cp "${dir}/Scripts/fitcmds/runSMFits.sh" SM/
    cp DC_SM* SM/
    # Adjust the cards in the SM dir:
    sed -i '/^cross_section_tttt.*/d' SM/DC_SM_2018.txt
    sed -i '/^cross_section_tttt.*/d' SM/DC_SM_2017.txt
    sed -i '/^cross_section_tttt.*/d' SM/DC_SM_2016.txt

    cd SM
    ./runSMFits.sh init > out.txt &
    cd ..
}

yukawa_central () {
    cp "${dir}/Scripts/fitcmds/runYukawaScan.sh" Yukawa/
    cp DC_SM* Yukawa/
    cd Yukawa
    ./runYukawaScan.sh init > out.txt &
    cd ..
}

yukawa_variations () {
    # Copy the datacards to the Yukawa dirs:
    cp DC_SM* Yukawa_variations/

    cd Yukawa_variations/
    mkdir tth_only tttt_only ttt_only tttt_tth tttt_ttt ttt_tth
    # cp the DCs from Yukawa there
    # Loop dirs:
    for subdir in tth_only tttt_only ttt_only tttt_tth tttt_ttt ttt_tth; do
        cp DC_SM* $subdir/
        cp "${dir}/Scripts/fitcmds/runYukawaScan.sh" $subdir/
    done

    sed -i '/ttH:r_TTH=expr::r_TTH/d' tttt_ttt/runYukawaScan.sh
    sed -i '/tttt:r_TTTT=expr::r_TTTT/d' ttt_tth/runYukawaScan.sh
    sed -i '/tttj:r_TTTJ=expr::r_TTTJ/d' tttt_tth/runYukawaScan.sh
    sed -i '/tttW:r_TTTW=expr::r_TTTW/d' tttt_tth/runYukawaScan.sh

    sed -i '/tttt:r_TTTT=expr::r_TTTT/d' tth_only/runYukawaScan.sh
    sed -i '/tttj:r_TTTJ=expr::r_TTTJ/d' tth_only/runYukawaScan.sh
    sed -i '/tttW:r_TTTW=expr::r_TTTW/d' tth_only/runYukawaScan.sh

    sed -i '/ttH:r_TTH=expr::r_TTH/d' tttt_only/runYukawaScan.sh
    sed -i '/tttj:r_TTTJ=expr::r_TTTJ/d' tttt_only/runYukawaScan.sh
    sed -i '/tttW:r_TTTW=expr::r_TTTW/d' tttt_only/runYukawaScan.sh

    sed -i '/tttt:r_TTTT=expr::r_TTTT/d' ttt_only/runYukawaScan.sh
    sed -i '/ttH:r_TTH=expr::r_TTH/d' ttt_only/runYukawaScan.sh

    for subdir in tth_only tttt_only ttt_only tttt_tth tttt_ttt ttt_tth; do
        cd $subdir
        ./runYukawaScan.sh init > out.txt &
        cd ..
    done
    cd ..
}


sm_ttt_tttt () {
    cp "${dir}/Scripts/fitcmds/runTripleTopFourTopScan.sh" scan_ttt_tttt/
    cp DC_SM_* SM_ttt_tttt/
    sed -i '/^cross_section_ttt.*/d' SM_ttt_tttt/DC_SM_2018.txt
    sed -i '/^cross_section_ttt.*/d' SM_ttt_tttt/DC_SM_2017.txt
    sed -i '/^cross_section_ttt.*/d' SM_ttt_tttt/DC_SM_2016.txt

    cp SM_ttt_tttt/DC_SM_* scan_ttt_tttt/

    cp "${dir}/Scripts/fitcmds/runSM_TripleFourtop.sh" SM_ttt_tttt/
    cp "${dir}/Scripts/fitcmds/runTripleTopFourTopScan.sh" scan_ttt_tttt/

    cd scan_ttt_tttt
    ./runTripleTopFourTopScan.sh init > out.txt &
    cd ../SM_ttt_tttt
    ./runSM_TripleFourtop.sh init > out.txt &
    cd ..
}

if [ "$2" == "all" ]; then
    # note: submits 43k jobs!
    init
    standardmodel_fit
    yukawa_central
    yukawa_variations
    sm_ttt_tttt
fi

if [ "$2" == "init" ]; then
    init
fi

if [ "$2" == "sm" ]; then
    init
    standardmodel_fit
fi

if [ "$2" == "yukawa" ]; then
    init
    yukawa_central
fi

if [ "$2" == "yukawa_var" ]; then
    init
    yukawa_variations
fi

if [ "$2" == "sm_alt" ]; then
    init
    sm_ttt_tttt
fi
