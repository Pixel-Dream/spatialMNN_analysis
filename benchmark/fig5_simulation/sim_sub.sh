#!/bin/bash

# Define the R script to be executed
R_SCRIPT="/users/hzhou1/benchmark/fig5_simulation/sim_script.R"
OUTPUT_PATH="/users/hzhou1/benchmark/fig5_simulation/datasets"
# Define the arguments for each job submission
declare -a ARGU1=(
    "1 2"
    "1 4"
    "1 8"
    "1 16"
    "1 32"
)

declare -a ARGU2=(
    "1 48"
    "1 64"
)

submit_job () {
    local ARGS=$1
    local MEM=$2
    local MODULE=$3
    local OUTPUT_PATH=$4

    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=sim_bench
#SBATCH --output=${OUTPUT_PATH}/log/output_%j.txt
#SBATCH --error=${OUTPUT_PATH}/log/error_%j.txt
#SBATCH --time=1-12:00:00  
#SBATCH --mem=${MEM}         
#SBATCH --cpus-per-task=4

# Load the required modules (if any)
module load ${MODULE}

# Execute the R script with the current set of arguments
Rscript $R_SCRIPT $ARGS $OUTPUT_PATH
EOT
}

mkdir -p "$OUTPUT_PATH/log"
    
for ARGS in "${ARGU1[@]}"; do
    submit_job "$ARGS" "32G" "R/4.3" "$OUTPUT_PATH"
done

for ARGS in "${ARGU2[@]}"; do
    submit_job "$ARGS" "80G" "R/4.3" "$OUTPUT_PATH"
done