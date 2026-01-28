#!/bin/bash

## GITHUB CONFIGURATION
GITHUB_OWNER="kidsneuro-lab"
GITHUB_REPO="rna_cloud_genome_reference"
VERSION="0.0.15"

## GADI CONFIGURATION

# BASE DIRECTORY
BASE_DIR=/g/data/qe93/rna_cloud/reference_data

# ASSEMBLY AND ANNOTATION DIRECTORIES
CUSTOM_ASSEMBLY="${BASE_DIR}/assemblies/GRCh38/custom/rna_cloud/$VERSION"
CUSTOM_ANNOTATION="${BASE_DIR}/annotations/GRCh38/custom/rna_cloud/$VERSION"

rm -rf ${CUSTOM_ASSEMBLY} && mkdir -p ${CUSTOM_ASSEMBLY}
rm -rf ${CUSTOM_ANNOTATION} && mkdir -p ${CUSTOM_ANNOTATION}

echo "🔵 Downloading custom assembly to ${CUSTOM_ASSEMBLY}"
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".fasta.gz"      $CUSTOM_ASSEMBLY
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".fasta.gz.fai"  $CUSTOM_ASSEMBLY
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".fasta.gz.gzi"  $CUSTOM_ASSEMBLY

echo "🔵 Downloading custom annotation to ${CUSTOM_ANNOTATION}"
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".gtf.gz"      $CUSTOM_ANNOTATION
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".gtf.gz.tbi"  $CUSTOM_ANNOTATION
./github-release-downloader.sh $GITHUB_OWNER $GITHUB_REPO $VERSION ".bed"         $CUSTOM_ANNOTATION

echo "🔵 Decompressing custom assembly files"
for f in "$CUSTOM_ASSEMBLY"/*.gz; do
  [ -e "$f" ] || continue
  gunzip -c "$f" > "${f%.gz}"
done

echo "🔵 Decompressing custom annotation files"
for f in "$CUSTOM_ANNOTATION"/*.gz; do
  [ -e "$f" ] || continue
  gunzip -c "$f" > "${f%.gz}"
done

