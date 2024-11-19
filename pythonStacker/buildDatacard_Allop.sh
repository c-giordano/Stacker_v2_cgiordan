#!/bin/sh
ulimit -s unlimited
set -e
cd /ada_mnt/ada/user/mshoosht/work/CMSSW_13_3_3/src
export SCRAM_ARCH=slc7_amd64_gcc12
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
cd /user/mshoosht/work/CMSSW_13_3_3/src/plots/pythonStacker/

if [ $1 -eq 0 ]; then
python3 buildDatacard.py -vf settingfiles/Variables/base.json -sf settingfiles/Uncertainties/2016_Signal.json -pf settingfiles/Process/SM.json -cf settingfiles/Channel/all_channels.json -y 2016 -dcf settingfiles/Datacards/2016_BDT.json --EFT --EFTop ctt cQQ1 cQQ8 cQt1 cQt8 ctHIm ctHRe -op ./output/datacards/v15_BDT/Allops/ --All_EFT
fi
if [ $1 -eq 1 ]; then
python3 buildDatacard.py -vf settingfiles/Variables/base.json -sf settingfiles/Uncertainties/2017_Signal.json -pf settingfiles/Process/SM.json -cf settingfiles/Channel/all_channels.json -y 2017 -dcf settingfiles/Datacards/2017_BDT.json --EFT --EFTop ctt cQQ1 cQQ8 cQt1 cQt8 ctHIm ctHRe -op ./output/datacards/v15_BDT/Allops/ --All_EFT
fi
if [ $1 -eq 2 ]; then
python3 buildDatacard.py -vf settingfiles/Variables/base.json -sf settingfiles/Uncertainties/2018_Signal.json -pf settingfiles/Process/SM.json -cf settingfiles/Channel/all_channels.json -y 2018 -dcf settingfiles/Datacards/2018_BDT.json --EFT --EFTop ctt cQQ1 cQQ8 cQt1 cQt8 ctHIm ctHRe -op ./output/datacards/v15_BDT/Allops/ --All_EFT
fi

