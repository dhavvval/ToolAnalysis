#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <RunNumber>"
    exit 1
fi

RUNNUM=$1
OUTPUTPATH=/exp/annie/app/users/dajana/AmBe_testDump/$RUNNUM
RAWFILE=/exp/annie/app/users/dajana/AmBe_testDump/$RUNNUM/RAWDataR${RUNNUM}S0p0
RAWFILECONF=configfiles/EventBuilderV2_AmBe/my_files.txt
PROCFILE=./ProcessedData_PMT_R${RUNNUM}S0p0
PROCFILECONF=configfiles/BeamClusterAnalysis/my_inputs.txt
WAVEFORMCONF=configfiles/AmBeWaveform/my_inputs.txt

# Create output directory
mkdir -p $OUTPUTPATH

# Check raw file
if [ ! -f "$RAWFILE" ]; then
    echo "Error: Raw file $RAWFILE not found!"
    exit 1
fi

# EventBuilder
echo $RAWFILE > $RAWFILECONF
./Analyse configfiles/EventBuilderV2_AmBe/ToolChainConfig || { echo "EventBuilder failed"; exit 1; }

# BeamClusterAnalysis
echo $PROCFILE > $PROCFILECONF
./Analyse configfiles/BeamClusterAnalysis/ToolChainConfig || { echo "BeamClusterAnalysis failed"; exit 1; }
mv ANNIETree.root $OUTPUTPATH/AmBe_${RUNNUM}_v3.ntuple.root

# AmBeWaveform
echo $PROCFILE > $WAVEFORMCONF
./Analyse configfiles/AmBeWaveform/ToolChainConfig || { echo "AmBeWaveform failed"; exit 1; }

# Move outputs
mv $PROCFILE $OUTPUTPATH
mv RWMBRFWaveforms.root $OUTPUTPATH/AmBeWaveforms_${RUNNUM}.root

echo "Processing for run $RUNNUM completed successfully."

