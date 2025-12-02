#!/bin/bash

# Check if script name is provided
if [[ -z "$1" ]]; then
    echo "Usage: $0 <script_name>"
    exit 1
fi

script_name=$1

echo ""

# Ensure the list file exists
if [[ ! -f "runs.list" ]]; then
    echo "Error: runs.list not found!"
    exit 1
fi

# Ensure the specified script exists and is executable
if [[ ! -x "./${script_name}" ]]; then
    echo "Error: ${script_name} not found or not executable!"
    exit 1
fi

# Process each line in runs.list
while IFS= read -r file; do
    # Skip empty lines or lines with only whitespace
    if [[ -n "$file" ]]; then
        echo "Processing $file..."
        ./"${script_name}" "$file"
        
        # Check for errors in the processing script
        if [[ $? -ne 0 ]]; then
            echo "Error: ${script_name} failed for $file"
            exit 2
        fi
    fi
done < runs.list

echo "All files processed successfully."
echo ""

