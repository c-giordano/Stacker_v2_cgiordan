#!/bin/bash

cd /user/nivanden/CMSSW_13_3_3/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`
export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)
cd /ada_mnt/ada/user/nivanden/plots/pythonStacker

# if no arg is given, build all datacards
if [ $# -eq 0 ]; then
  echo "Building all datacards"
  python3 buildDatacard.py -sf settingfiles/Uncertainties/2016_yukawa.json -pf settingfiles/Process/SM_Yukawa.json -cf settingfiles/Channel/all_channels.json -y 2016 -dcf settingfiles/Datacards/2016_SM.json --data
  python3 buildDatacard.py -sf settingfiles/Uncertainties/2017_yukawa.json -pf settingfiles/Process/SM_Yukawa.json -cf settingfiles/Channel/all_channels.json -y 2017 -dcf settingfiles/Datacards/2017_SM.json --data
  python3 buildDatacard.py -sf settingfiles/Uncertainties/2018_yukawa.json -pf settingfiles/Process/SM_Yukawa.json -cf settingfiles/Channel/all_channels.json -y 2018 -dcf settingfiles/Datacards/2018_SM.json --data
  exit 0
fi

if [ $1 -eq 0 ]; then
  # echo "test12"
  python3 buildDatacard.py -sf settingfiles/Uncertainties/2016_yukawa.json -pf settingfiles/Process/SM_Yukawa.json -cf settingfiles/Channel/all_channels.json -y 2016 -dcf settingfiles/Datacards/2016_SM.json --data
fi
if [ $1 -eq 1 ]; then
  python3 buildDatacard.py -sf settingfiles/Uncertainties/2017_yukawa.json -pf settingfiles/Process/SM_Yukawa.json -cf settingfiles/Channel/all_channels.json -y 2017 -dcf settingfiles/Datacards/2017_SM.json --data
fi
if [ $1 -eq 2 ]; then
  python3 buildDatacard.py -sf settingfiles/Uncertainties/2018_yukawa.json -pf settingfiles/Process/SM_Yukawa.json -cf settingfiles/Channel/all_channels.json -y 2018 -dcf settingfiles/Datacards/2018_SM.json --data
fi



