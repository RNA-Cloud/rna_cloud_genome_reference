nextflow.enable.dsl=2

process VALIDATE_GENOME_ANNOTATION {
    tag "VALIDATE_GENOME_ANNOTATION"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path fasta
    path fasta_index
    path gtf
    path gtf_index
    path masked_regions_bed
    path unmasked_regions_bed
    path annotation_bed
    val ncbi_assembly_masked_regions_paths
    path refseq_mane_annotation
    path refseq_mane_bed_file

    output:
    path "validation_report.txt", emit: report

    script:
    """
    set -euo pipefail

    echo "Append NCBI assembly masked regions to masked_regions_bed"
    for bed in ${ncbi_assembly_masked_regions_paths.join(' ')}; do
      cat \$bed >> combined_ncbi_assembly_masked_regions.bed
    done
    bedtools subtract -a ${unmasked_regions_bed} -b combined_ncbi_assembly_masked_regions.bed > updated_unmasked_regions.bed

    echo "Running validation script"
    /app/rnacloud_genome_reference/validation/scripts/validation.sh \
      ${fasta} \
      ${fasta_index} \
      ${gtf} \
      ${gtf_index} \
      ${masked_regions_bed} \
      updated_unmasked_regions.bed \
      ${annotation_bed} \
      ${refseq_mane_annotation} \
      ${refseq_mane_bed_file} | tee validation_report.txt
    """
}

process GENOME_AND_ANNOTATION_REPORT {
    tag "GENOME_AND_ANNOTATION_REPORT"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    // Required genome + indices
    path fasta                                   // .fasta.gz
    path fasta_gzi                               // .fasta.gz.gzi
    path fasta_fai                               // .fasta.gz.fai

    // Annotation + index
    path gtf                                     // .gtf.gz
    path gtf_index                               // .gtf.gz.tbi

    // References & metadata
    path assembly_report
    path cen_par_regions
    path grc_fixes_summary
    path config_json                             // conf/sources.json

    // Genome and annotation version
    val version

    output:
    path 'genome_and_annotation_report.md', emit: report

    script:
    def genome_and_annotation_version = version ?: "0.0.0"
    """
    set -euo pipefail

    python -m rnacloud_genome_reference.validation.report \
      --fasta ${fasta} \
      --gtf ${gtf} \
      --assembly-report ${assembly_report} \
      --cen-par-regions ${cen_par_regions} \
      --grc-fixes-summary ${grc_fixes_summary} \
      --config ${config_json} \
      --version ${genome_and_annotation_version} \
      --out genome_and_annotation_report.md
    """
}