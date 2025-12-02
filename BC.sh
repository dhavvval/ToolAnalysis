#!/bin/bash

#===============================================================
#Script: batch_beamcluster_processing.sh
#Description:
#For each run number, this script:
#1. Generates my_inputs.txt with all file paths for that run
#2. Updates ANNIEEventTreeMakerConfig to use the correct run number
#3. Runs 'make -j4' and then './Analyse' for processing
#Usage:
#./batch_beamcluster_processing.sh runs_list.txt
#(where runs_list.txt contains one run number per line)
#===============================================================

if [[ -z "$1" ]]; then
echo "Usage: $0 <runs_list_file>"
exit 1
fi

runs_list="$1"

#===============================================================
#Set your base working directory (where make and Analyse are run)
#===============================================================

base_dir="/exp/annie/app/users/dajana/EB_BC_TA"
beamcluster_dir="${base_dir}/configfiles/BeamClusterAnalysis"
output_file="${beamcluster_dir}/my_inputs.txt"
config_file="${beamcluster_dir}/ANNIEEventTreeMakerConfig"

#===============================================================
#Safety checks
#===============================================================

if [[ ! -f "$runs_list" ]]; then
echo "Error: File '$runs_list' not found."
exit 1
fi

if [[ ! -f "$config_file" ]]; then
echo "Error: Config file '$config_file' not found."
exit 1
fi

#===============================================================
#Move to base directory so make/Analyse run correctly
#===============================================================

cd "$base_dir" || { echo "Failed to cd into $base_dir"; exit 1; }

#===============================================================
#Main loop through each run number
#===============================================================

while read -r run; do
# Skip empty or commented lines
[[ -z "$run" || "$run" =~ ^# ]] && continue
echo "=============================="
echo " Processing run $run"
echo "=============================="

processed_dir="/pnfs/annie/persistent/processed/processed_EBV2/R${run}/"

if [[ ! -d "$processed_dir" ]]; then
    echo "Warning: Directory $processed_dir not found. Skipping run $run."
    continue
fi

# Step 1: Create my_inputs.txt at the desired location
find "$processed_dir" -type f -name "*ProcessedData_PMT*" | sort -t'/' -k2V > "$output_file"
echo "Created $output_file for run $run"

# Step 2: Update the config file for the correct output filename
if grep -q "OutputFile" "$config_file"; then
    sed -i "s|OutputFile .*|OutputFile AmBe_${run}_v3.ntuple.root|g" "$config_file"
    echo "Updated OutputFile line in $config_file"
else
    echo "Error: Could not find 'OutputFile' line in $config_file"
    exit 1
fi

# Step 3: Build and run analysis
echo "Running make -j4 in $base_dir"
make -j4 || { echo "make failed for run $run"; exit 1; }

echo "Running Analyse ToolChainConfig"
./Analyse ./configfiles/BeamClusterAnalysis/ToolChainConfig || { echo "Analysis failed for run $run"; exit 1; }

echo "Finished processing run $run"
echo

done < "$runs_list"

echo "All runs processed successfully!"
