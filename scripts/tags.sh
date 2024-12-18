# Set Pandora package versions
export PANDORA_PFA_VERSION=v04-09-00
export PANDORA_SDK_VERSION=v03-04-01
export PANDORA_MONITORING_VERSION=v03-06-00
export PANDORA_LAR_CONTENT_VERSION=v04_11_00
export PANDORA_LAR_MLDATA_VERSION=v04-09-00
export PANDORA_LAR_RECO_ND_VERSION=v01-01-03
export EIGEN_VERSION=3.4.0

# Set main working directory by optional run argument
dirName=$1
testArea=$PWD

if [[ $dirName ]]; then
    testArea=$dirName
fi
echo "MY_TEST_AREA is $testArea"
export MY_TEST_AREA=${testArea}

# Set FW_SEARCH_PATH for Pandora xml run files & machine learning data etc
export FW_SEARCH_PATH=$MY_TEST_AREA/LArRecoND/settings
export FW_SEARCH_PATH=$MY_TEST_AREA/LArMachineLearningData:$FW_SEARCH_PATH
