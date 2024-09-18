#!/bin/bash

# Define the R script to be executed
R_SCRIPT="/users/hzhou1/benchmark/fig5_simulation/benchmark.R"
# Define the base output path
OUTPUT_BASE_PATH="/users/hzhou1/benchmark/fig5_simulation/res"

# Define the arguments for each job submission
declare -a R43=(
    "0 spatialMNN,spatialMNN_par,PRECAST"
    "1 BayesSpace"
    "2 Seurat"
    "3 BASS"
)

declare -a R44=(
    "4 BANKSY"
)


# SLAT, each ref need an independent run
generate_py_array () {
    local OUTPUT_PATH=$1
    local NUM_SLAT_RUNS=$2
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
    local ITER=$5
    local SAMPLENUM=$6
    local OPTIONAL=$7

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
Rscript $R_SCRIPT $ARGS $OUTPUT_PATH $ITER $SAMPLENUM $OPTIONAL
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

declare -a SAMPLE_NUM_ARRAY=(
   2 4 8 16 32 48 64
)

# Loop through each output path and submit jobs
for OUTPUT_PATH in "${OUTPUT_PATHS[@]}"; do
    mkdir -p "$OUTPUT_PATH/log"
    for SAMPLE_NUM in "${SAMPLE_NUM_ARRAY[@]}"; do
        if [[ "$SAMPLE_NUM" == "2" ]]; then
            MEM="16G"
        elif [[ "$SAMPLE_NUM" == "4" ]]; then
            MEM="16G"
        elif [[ "$SAMPLE_NUM" == "8" ]]; then
            MEM="24G"
        elif [[ "$SAMPLE_NUM" == "16" ]]; then
            MEM="48G"
        elif [[ "$SAMPLE_NUM" == "32" ]]; then
            MEM="80G"
        else
            MEM="128G"
        fi

        #MEM="128G"

        for ARGS in "${R43[@]}"; do
            submit_job "$ARGS" "$MEM" "R/4.3" "$OUTPUT_PATH" "1" "$SAMPLE_NUM"
        done
        
        for ARGS in "${R44[@]}"; do
            submit_job "$ARGS" "$MEM" "conda_R/4.4" "$OUTPUT_PATH" "1" "$SAMPLE_NUM"
        done
        
        # MENDER
        submit_job "5 MENDER" "$MEM" "R/4.3" "$OUTPUT_PATH" "1" "$SAMPLE_NUM"
        if [ $SAMPLE_NUM -gt 8 ]; then
            REF_SAMPLE_NUM=8
        else
            REF_SAMPLE_NUM=$SAMPLE_NUM
        fi
        
        # SLAT
        for (( i=1; i<=REF_SAMPLE_NUM; i++ )); do
            ARGS=("$((i+5)) SLAT")
            submit_job "$ARGS" "$MEM" "R/4.3" "$OUTPUT_PATH" "1" "$SAMPLE_NUM" "$i"
        done
    done


done
