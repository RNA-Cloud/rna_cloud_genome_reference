#!/bin/bash
set -e # exit on first error

echo "----------------------------------------------------------------------"
echo "Starting rnacloud_genome_reference GRC fixes assessment"
echo "----------------------------------------------------------------------"

echo "Current folder: $(pwd)"

GENOME_AND_ANNOTATION_VERSION=${GENOME_AND_ANNOTATION_VERSION:-"0.0.0"}
echo "Genome and annotation version: $GENOME_AND_ANNOTATION_VERSION"

# Check for -resume flag
RESUME_FLAG=""
if [[ "$1" == "-resume" ]]; then
  RESUME_FLAG="-resume"
  echo "Resume mode enabled"
fi

docker run \
  --rm \
  --name rnacloud_runner \
  --entrypoint nextflow \
  -e GENOME_AND_ANNOTATION_VERSION=$GENOME_AND_ANNOTATION_VERSION \
  -v "$(pwd)"/data:/app/data \
  -v "$(pwd)"/temp:/app/temp \
  -v "$(pwd)"/output:/app/output \
  -v "$(pwd)"/work:/app/work \
  -v "$(pwd)"/reference:/app/reference \
  rnacloud_runner main.nf \
    -with-dag /app/output/dag.mmd \
    -with-report /app/output/report.html \
    -with-timeline /app/output/timeline.html \
    -params-file /app/conf/sources.json \
    $RESUME_FLAG