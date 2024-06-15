#!/bin/bash

# Define the base directory for easier reference
BASE_DIR="./Data/DPO/USPTO_50K"

# Function to run Python script with parameters
run_script() {
    python scripts.py --log_dir $1 --data_dir $2 --rule_file_path $3 --save_dir $4 --radius $5 ${6:-}
    if [ $? -ne 0 ]; then
        echo "Error executing with parameters: $*"
        exit 1
    fi
}

echo "Starting batch processing..."

# Non-hydrogen cases
run_script "$BASE_DIR/Non_hydrogen/Log/rules_r0_log.txt" "$BASE_DIR/test.json.gz" "$BASE_DIR/Non_hydrogen/R0" "$BASE_DIR/Non_hydrogen/Output" 0
run_script "$BASE_DIR/Non_hydrogen/Log/rules_r1_log.txt" "$BASE_DIR/test.json.gz" "$BASE_DIR/Non_hydrogen/R1" "$BASE_DIR/Non_hydrogen/Output" 1

# Hydrogen cases
run_script "$BASE_DIR/Hydrogen/Log/rules_r0_log.txt" "$BASE_DIR/test.json.gz" "$BASE_DIR/Hydrogen/R0" "$BASE_DIR/Hydrogen/Output" 0
run_script "$BASE_DIR/Hydrogen/Log/rules_r1_log.txt" "$BASE_DIR/test.json.gz" "$BASE_DIR/Hydrogen/R1" "$BASE_DIR/Hydrogen/Output" 1

# Good hydrogen cases
run_script "$BASE_DIR/Good_hydrogen/Log/rules_r0_log.txt" "$BASE_DIR/test.json.gz" "$BASE_DIR/Good_hydrogen/R0" "$BASE_DIR/Good_hydrogen/Output" 0
run_script "$BASE_DIR/Good_hydrogen/Log/rules_r1_log.txt" "$BASE_DIR/test.json.gz" "$BASE_DIR/Good_hydrogen/R1" "$BASE_DIR/Good_hydrogen/Output" 1
run_script "$BASE_DIR/Good_hydrogen/Log/hier_rule_log.txt" "$BASE_DIR/test.json.gz" "$BASE_DIR/Good_hydrogen/R0" "$BASE_DIR/Good_hydrogen/Output" 0 "--hierarchical True"

echo "All processes completed successfully."
