#!/bin/bash
# Author: Steven Doran (modified by Dhaval Ajana)

echo ""
echo "Please make sure you have done the following:
  - Toolchain: LoadANNIEEvent, LoadGeometry, FitRWMWaveform
  - Enable 'printToRootFile 1' in FitRWMWaveform
  - Set maxPrintNumber according to number of part files
  - Ensure correct configurations in this script
"
sleep 3

# Check if user provided run list
if [[ -z "$1" ]]; then
    echo ""
    echo "##############################"
    echo "Error: No run list provided."
    echo "Usage: $0 <run_list.txt>"
    echo ""
    exit 1
fi

run_list_file=$1

if [[ ! -f "$run_list_file" ]]; then
    echo "Error: File $run_list_file not found!"
    exit 1
fi

# User-specific parameters
step_size=20
user="dajana"
TA_folder="EB_BC_TA"
toolchain="AmBeWaveform"
offload="/pnfs/annie/persistent/users/${user}/AmBe/"

# Process each run number from list
while read -r run; do
    # Skip empty lines or commented lines
    [[ -z "$run" || "$run" =~ ^# ]] && continue

    echo ""
    echo "=============================="
    echo "Processing Run ${run}"
    echo "=============================="
    echo ""

    pro_dir="/pnfs/annie/persistent/processed/processed_EBV2/R${run}/"

    if [[ ! -d "$pro_dir" ]]; then
        echo "Warning: Directory $pro_dir does not exist. Skipping..."
        continue
    fi

    pro_files=($(ls "$pro_dir" | grep "^ProcessedData_PMT_R${run}S0p*"))
    num_pro_files=${#pro_files[@]}

    if (( num_pro_files == 0 )); then
        echo "No processed files found for run ${run}. Skipping..."
        continue
    fi

    if (( num_pro_files <= step_size )); then
        start_indices=("0")
        end_indices=($(($num_pro_files - 1)))
    else
        start_indices=("0")
        end_indices=()
        for (( i=step_size; i<num_pro_files; i+=step_size )); do
            start_indices+=("$i")
            end_indices+=("$((i - 1))")
        done
        end_indices+=("$((num_pro_files - 1))")
    fi

    echo "Start indices: ${start_indices[@]}"
    echo "End indices: ${end_indices[@]}"

    for (( idx=0; idx<${#start_indices[@]}; idx++ )); do
        p_start=${start_indices[$idx]}
        p_end=${end_indices[$idx]}

        echo ""
        echo "Running over [${p_start},${p_end}]"
        echo ""

        my_files="my_inputs.txt"
        rm -f "$my_files"

        for p in $(seq "$p_start" "$p_end"); do
            echo "${pro_dir}/ProcessedData_PMT_R${run}S0p${p}" >> "$my_files"
        done

        cp "$my_files" /exp/annie/app/users/$user/$TA_folder/configfiles/$toolchain/

        echo ""
        cat /exp/annie/app/users/$user/$TA_folder/configfiles/$toolchain/$my_files
        echo ""

        sleep 3

        singularity shell -B/pnfs:/pnfs,/exp/annie/app/users/$user/temp_directory:/tmp,/exp/annie/data:/exp/annie/data,/exp/annie/app:/exp/annie/app /cvmfs/singularity.opensciencegrid.org/anniesoft/toolanalysis:latest << EOF
        cd /exp/annie/app/users/$user/$TA_folder
        source Setup.sh
        ./Analyse ./configfiles/$toolchain/ToolChainConfig
        exit
EOF

        mkdir -p $offload/$run
        source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
        setup ifdhc v2_5_4

        cp /exp/annie/app/users/$user/$TA_folder/RWMBRFWaveforms.root \
           $offload/$run/AmBeWaveforms_${run}_p${p_start}_p${p_end}.root

        rm -rf /exp/annie/app/users/$user/$TA_folder/RWMBRFWaveforms.root

        echo ""
        ls -lrth $offload/$run
        echo ""
    done

    echo "Finished processing Run ${run}"
    echo ""

done < "$run_list_file"

echo "All runs from $run_list_file processed successfully."

