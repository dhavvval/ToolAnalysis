#!/bin/bash

# Directory containing input ROOT files
INPUT_DIR="/pnfs/annie/persistent/users/dajana/WCSim_files/AmBe_BeamCluster"
# Config file paths
LOAD_WCSIM_CONFIG="/exp/annie/app/users/dajana/EB_BC_TA/configfiles/BeamClusterAnalysisMC/LoadWCSimConfig"
TREEMAKER_CONFIG="/exp/annie/app/users/dajana/EB_BC_TA/configfiles/BeamClusterAnalysisMC/PhaseIITreeMakerConfig"
# ToolChain path
TOOLCHAIN="./configfiles/BeamClusterAnalysisMC/ToolChainConfig"

# Loop over all wcsim_*.root files
for infile in "$INPUT_DIR"/wcsim_*.root; do
    echo "Processing file: $infile"

    # Extract just the file name
    filename=$(basename "$infile")

    # Construct corresponding output file name
    # e.g., wcsim_port4_z0.root -> AmBe_Neutron_port4_z0.root
    suffix=${filename#wcsim_}
    outfile="${INPUT_DIR}/AmBe_Neutron_${suffix}"

    echo "Setting InputFile = $infile"
    echo "Setting OutputFile = $outfile"

    # Update LoadWCSimConfig
    sed -i "/^InputFile /c\InputFile $infile" "$LOAD_WCSIM_CONFIG"

    # Update PhaseIITreemakerConfig
    sed -i "/^OutputFile /c\OutputFile $outfile" "$TREEMAKER_CONFIG"

    # Run compilation and analysis
    make -j4
    ./Analyse "$TOOLCHAIN"

    echo "Done with $filename"
    echo "------------------------"
done
