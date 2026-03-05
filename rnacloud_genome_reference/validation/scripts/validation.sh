#!/bin/bash

set -uo pipefail

if [ $# -ne 9 ]; then
  echo "Usage: $0 <FASTA> <FASTA_INDEX> <GTF> <GTF_INDEX> <MASKED_REGIONS_BED> <UNMASKED_REGIONS_BED> <ANNOTATION_BED> <REFSEQ_MANE_ANNOTATION> <REFSEQ_MANE_BED_FILE>"
  exit 1
fi

FASTA="$1"
FASTA_INDEX="$2"
GTF="$3"
GTF_INDEX="$4"
MASKED_REGIONS_BED="$5"
UNMASKED_REGIONS_BED="$6"
ANNOTATION_BED="$7"
REFSEQ_MANE_ANNOTATION="$8"
REFSEQ_MANE_BED_FILE="$9"

echo "FASTA                  : $FASTA"
echo "FASTA Index            : $FASTA_INDEX"
echo "GTF                    : $GTF"
echo "GTF Index              : $GTF_INDEX"
echo "Masked Regions BED     : $MASKED_REGIONS_BED"
echo "Unmasked Regions BED   : $UNMASKED_REGIONS_BED"
echo "Annotation BED         : $ANNOTATION_BED"
echo "RefSeq MANE Annotation : $REFSEQ_MANE_ANNOTATION"
echo "RefSeq MANE BED File   : $REFSEQ_MANE_BED_FILE"

# Initialize validation status
VALIDATION_STATUS=0
echo 🎬 Starting Genome and Annotation build validation process

echo 🕵️‍♀️ Comparing contigs in GTF and FASTA files...

GTF_CONTIGS=$(mktemp --suffix=.txt)
FASTA_CONTIGS=$(mktemp --suffix=.txt)

tabix -l "$GTF" > "$GTF_CONTIGS"
cut -f 1 "$FASTA_INDEX" > "$FASTA_CONTIGS"

diff "$GTF_CONTIGS" "$FASTA_CONTIGS"

if [ $? -eq 0 ]; then
  echo ✅ Contigs in GTF and FASTA match!
else
  echo ⛔️ Contigs in GTF and FASTA do not match!
  VALIDATION_STATUS=1
fi

echo 🕵️‍♀️ Check that contigs in GTF and FASTA are sorted in the correct order...

GTF_CONTIGS_SORTED=$(mktemp --suffix=.txt)
FASTA_CONTIGS_SORTED=$(mktemp --suffix=.txt)

tabix -l "$GTF" | sort -k1,1V > "$GTF_CONTIGS_SORTED"
cut -f 1 "$FASTA_INDEX" | sort -k1,1V > "$FASTA_CONTIGS_SORTED"

diff "$GTF_CONTIGS_SORTED" "$GTF_CONTIGS" > /dev/null
if [ $? -ne 0 ]; then
  echo ⛔️ Contigs in GTF are not sorted in the correct order!
  VALIDATION_STATUS=1
else
  echo ✅ Contigs in GTF are sorted in the correct order!
fi

diff "$FASTA_CONTIGS_SORTED" "$FASTA_CONTIGS" > /dev/null
if [ $? -ne 0 ]; then
  echo ⛔️ Contigs in FASTA are not sorted in the correct order!
  VALIDATION_STATUS=1
else
  echo ✅ Contigs in FASTA are sorted in the correct order!
fi

rm -f "$GTF_CONTIGS" "$FASTA_CONTIGS" "$GTF_CONTIGS_SORTED" "$FASTA_CONTIGS_SORTED"

echo 🕵️‍♀️ Check that EBV contig is present in the FASTA file...

EBV_FASTA_CONTIG_COUNT=$(cut -f 1 "$FASTA_INDEX" | grep 'chrEBV' | wc -l)
EBV_GTF_CONTIG_COUNT=$(tabix -l "$GTF" | grep 'chrEBV' | wc -l)

if [ "$EBV_FASTA_CONTIG_COUNT" -eq 1 ] && [ "$EBV_GTF_CONTIG_COUNT" -eq 1 ]; then
  echo ✅ EBV contig is present in the FASTA and GTF files!
else
  echo ⛔️ EBV contig is not present in the FASTA or GTF file or contains multiple entries \(FASTA: ${EBV_FASTA_CONTIG_COUNT}, GTF: ${EBV_GTF_CONTIG_COUNT}\)!
  VALIDATION_STATUS=1
fi

echo 🕵️‍♀️ Check that masked regions only contain Ns

bedtools getfasta \
  -fi "$FASTA" \
  -bed "$MASKED_REGIONS_BED" \
  -tab \
  -fo masked_regions.fasta 2>&1 | tee masked_regions.fasta.log

ATCG_COUNT=$(awk -F"\t" '$2~/[ATCG]/' masked_regions.fasta | wc -l)

if [ "$ATCG_COUNT" -eq 0 ]; then
  echo ✅ Masked regions only contain Ns!
else
  echo ⛔️ Masked regions contain non-N bases \(Contigs with ATCG bases: ${ATCG_COUNT}\)!
  VALIDATION_STATUS=1
fi

WARNINGS_COUNT=$(grep 'WARNING' masked_regions.fasta.log | wc -l)

if [ "$WARNINGS_COUNT" -eq 0 ]; then
  echo ✅ No warnings found in masked regions FASTA extraction!
else
  echo ⛔️ Warnings found in masked regions FASTA extraction \(No. of warnings: ${WARNINGS_COUNT}\)!
  VALIDATION_STATUS=1
fi

echo 🕵️‍♀️ Check that unmasked regions only do not contain Ns

bedtools getfasta \
  -fi "$FASTA" \
  -bed "$UNMASKED_REGIONS_BED" \
  -tab \
  -fo unmasked_regions.fasta 2>&1 | tee unmasked_regions.fasta.log

N_COUNT=$(awk -F"\t" '$2~/N/' unmasked_regions.fasta | wc -l)

if [ "$N_COUNT" -eq 0 ]; then
  echo ✅ Unmasked regions only contain ATCGs!
else
  echo ⛔️ Unmasked regions contain N bases \(Regions with N bases: ${N_COUNT}\)!
  VALIDATION_STATUS=1
fi

WARNINGS_COUNT=$(grep 'WARNING' unmasked_regions.fasta.log | wc -l)

if [ "$WARNINGS_COUNT" -eq 0 ]; then
  echo ✅ No warnings found in unmasked regions FASTA extraction!
else
  echo ⛔️ Warnings found in unmasked regions FASTA extraction \(No. of warnings: ${WARNINGS_COUNT}\)!
  VALIDATION_STATUS=1
fi

if [ "$VALIDATION_STATUS" -ne 0 ]; then
  echo ❗ Validation failed!
  exit "$VALIDATION_STATUS"
fi

echo 🕵️‍♀️ Check that annotation bed file transcript names do not start with 'unassigned'
UNASSIGNED_COUNT=$(cat "$ANNOTATION_BED" | cut -f 4 | grep ^unassigned | wc -l)

if [ "$UNASSIGNED_COUNT" -eq 0 ]; then
  echo ✅ No transcript names start with 'unassigned' in annotation bed file!
else
  echo ⛔️ Some transcript names start with 'unassigned' in annotation bed file \(Count: ${UNASSIGNED_COUNT}\)!
  VALIDATION_STATUS=1
fi

echo 🕵️‍♀️ Check that number of lines in RefSeq MANE gtf and bed file are equal
# Note: gtfToGenePred outputs ONE LINE PER TRANSCRIPT model, not per gene.
#
# In the GTF:
#   feature == "gene"        -> one row per gene
#   feature == "transcript"  -> one row per transcript
#
# genePred rows correspond to transcripts (grouped exon structures),
# so the number of lines in the genePred output is expected to match
# the number of transcripts, not the number of genes.
#
# Therefore:
#   wc -l test.genepred
# will often be larger than:
#   awk '$3=="gene"' ...
#
# In the MANE RefSeq GTF most genes have a single transcript, but a small
# number of genes include an additional transcript (e.g. MANE Select +
# MANE Plus Clinical). These extra transcripts produce additional rows
# in the genePred output, explaining why the counts differ slightly.
#
# Example observed counts:
#   genePred rows (transcripts): 19437
#   gene feature rows (genes):   19363
#   difference:                   74 genes with an extra transcript

REFSEQ_MANE_GTF_COUNT=$(gunzip -c "$REFSEQ_MANE_ANNOTATION" | grep -v '^#' | awk -F'\t' '$3=="transcript"' | wc -l)
REFSEQ_MANE_BED_COUNT=$(wc -l < "$REFSEQ_MANE_BED_FILE")

if [ "$REFSEQ_MANE_GTF_COUNT" -eq "$REFSEQ_MANE_BED_COUNT" ]; then
  echo ✅ Number of (MANE Select and Clinical) transcripts in RefSeq MANE gtf and bed file are equal!
else
  echo ⛔️ Number of (MANE Select and Clinical) transcripts in RefSeq MANE gtf and bed file are not equal \(GTF: ${REFSEQ_MANE_GTF_COUNT}, BED: ${REFSEQ_MANE_BED_COUNT}\)!
  VALIDATION_STATUS=1
fi

if [ "$VALIDATION_STATUS" -ne 0 ]; then
  echo ❗ Validation failed!
  exit "$VALIDATION_STATUS"
fi

echo 🎉 All validations passed