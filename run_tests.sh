#!/bin/bash
set -e # exit on first error

echo "---------------------------"
echo "🏃 Running tests"
echo "---------------------------"

docker run \
  --rm \
  --name rnacloud_runner \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/temp:/app/temp \
  -v $(pwd)/output:/app/output \
  -v $(pwd)/work:/app/work \
  -v $(pwd)/reference:/app/reference \
  rnacloud_runner -m pytest -v

echo "---------------------------"
echo "🏋️ Preparing gnomAD data for batch run test"
echo "---------------------------"

if [ ! -f "data/gnomad/GRCh38/gnomad_r4_freq.tsv.gz" ]; then
  mkdir -p tests/batch_run/gnomad/GRCh38
  gunzip -c data/gnomad/GRCh38/gnomad_r4_freq.tsv.gz | awk -F"\t" 'NR == 1 || $1=="15" || $1=="22"' | gzip -c > tests/batch_run/gnomad/GRCh38/gnomad_r4_freq.tsv.gz
else
  echo "⏭️ gnomAD data already exists, skipping preparation"
fi

echo "---------------------------"
echo "🏃 Running batch run test"
echo "---------------------------"

BATCH_RUN_BASE_DIR="$(pwd)/tests/batch_run"

echo "Batch run base directory: $BATCH_RUN_BASE_DIR"
echo "Current folder: $(pwd)"

echo "🧹 Cleaning up batch run directories"
rm -rf $BATCH_RUN_BASE_DIR/data
rm -rf $BATCH_RUN_BASE_DIR/temp
rm -rf $BATCH_RUN_BASE_DIR/output
rm -rf $BATCH_RUN_BASE_DIR/work

echo "📂 Creating batch run directories"
mkdir -p $BATCH_RUN_BASE_DIR/data
mkdir -p $BATCH_RUN_BASE_DIR/temp
mkdir -p $BATCH_RUN_BASE_DIR/output
mkdir -p $BATCH_RUN_BASE_DIR/work

echo "📂 Copying gnomAD data to batch run data directory"
cp -r tests/batch_run/gnomad $BATCH_RUN_BASE_DIR/data

GENOME_AND_ANNOTATION_VERSION=${GENOME_AND_ANNOTATION_VERSION:-"0.0.0"}
echo "Genome and annotation version: $GENOME_AND_ANNOTATION_VERSION"

docker run \
  --rm \
  --name rnacloud_runner \
  --entrypoint nextflow \
  -e GENOME_AND_ANNOTATION_VERSION=$GENOME_AND_ANNOTATION_VERSION \
  -v $BATCH_RUN_BASE_DIR/conf:/app/conf \
  -v $BATCH_RUN_BASE_DIR/data:/app/data \
  -v $BATCH_RUN_BASE_DIR/temp:/app/temp \
  -v $BATCH_RUN_BASE_DIR/output:/app/output \
  -v $BATCH_RUN_BASE_DIR/work:/app/work \
  -v "$(pwd)/reference":/app/reference \
  rnacloud_runner main.nf \
    -with-dag /app/output/dag.mmd \
    -with-report /app/output/report.html \
    -with-timeline /app/output/timeline.html \
    -params-file /app/conf/sources.json

echo "---------------------------"
echo "✅ All tests completed successfully"
echo "---------------------------"