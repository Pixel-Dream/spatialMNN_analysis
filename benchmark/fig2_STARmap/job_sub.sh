#!/bin/bash

# Define the R script to be executed
R_SCRIPT="/users/hzhou1/benchmark/fig2_STARmap/benchmark.R"
# Define the base output path
OUTPUT_BASE_PATH="/users/hzhou1/benchmark/fig2_STARmap/result"

# Define the arguments for each job submission
declare -a R43=(
    "1 spatialMNN,BayesSpace,spatialMNN_par,PRECAST"
    "2 Seurat"
    "3 BASS"
)

declare -a R44=(
    "4 BANKSY"
)


# SLAT, each ref need an independent run
NUM_SLAT_RUNS=3

generate_py_array () {
    local OUTPUT_PATH=$1
    declare -a PY=(
        "5 MENDER"
    )

    for (( i=1; i<=NUM_SLAT_RUNS; i++ )); do
        PY=("${PY[@]}" "$((i+5)) SLAT $i")
    done

    echo "${PY[@]}"
}

submit_job () {
    local ARGS=$1
    local MEM=$2
    local MODULE=$3
    local OUTPUT_PATH=$4
    local OPTIONAL=$5

    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=sim_bench
#SBATCH --output=${OUTPUT_PATH}/log/output_%j.txt
#SBATCH --error=${OUTPUT_PATH}/log/error_%j.txt
#SBATCH --time=1-12:00:00  
#SBATCH --mem=${MEM}         
#SBATCH --cpus-per-task=8 

# Load the required modules (if any)
module load ${MODULE}

# Execute the R script with the current set of arguments
Rscript $R_SCRIPT $ARGS $OUTPUT_PATH $OPTIONAL
EOT
}

# List of output paths to use for submissions
declare -a OUTPUT_PATHS=(
    "$OUTPUT_BASE_PATH/run1"
    "$OUTPUT_BASE_PATH/run2"
    "$OUTPUT_BASE_PATH/run3"
    "$OUTPUT_BASE_PATH/run4"
    "$OUTPUT_BASE_PATH/run5"
)

# Loop through each output path and submit jobs
for OUTPUT_PATH in "${OUTPUT_PATHS[@]}"; do

    mkdir -p "$OUTPUT_PATH/log"
    
    for ARGS in "${R43[@]}"; do
        submit_job "$ARGS" "40G" "R/4.3" "$OUTPUT_PATH"
    done
    
    for ARGS in "${R44[@]}"; do
        submit_job "$ARGS" "40G" "conda_R/4.4" "$OUTPUT_PATH"
    done
    
    # MENDER
    submit_job "5 MENDER" "32G" "R/4.3" "$OUTPUT_PATH"
    # SLAT
    for (( i=1; i<=NUM_SLAT_RUNS; i++ )); do
        ARGS=("$((i+5)) SLAT")
        submit_job "$ARGS" "32G" "R/4.3" "$OUTPUT_PATH" "$i"
    done


done