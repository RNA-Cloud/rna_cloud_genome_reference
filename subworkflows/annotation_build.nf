nextflow.enable.dsl=2

include { DOWNLOAD_EBV_GTF } from '../modules/download.nf'
include { CONVERT_ANNOTATION_REFSEQ_TO_UCSC } from '../modules/annotation_build.nf'
include { GET_TARGET_CONTIGS } from '../modules/genome_build.nf'
include { SUBSET_GTF } from '../modules/annotation_build.nf'
include { APPEND_GTFS as APPEND_RIBOSOMAL_RNA_GTFS } from '../modules/annotation_build.nf'
include { APPEND_GTFS as APPEND_EBV_GTF } from '../modules/annotation_build.nf'
include { DECOMPRESS_GTF } from '../modules/annotation_build.nf'
include { REMOVE_SECTIONS } from '../modules/annotation_build.nf'
include { SORT_GTF } from '../modules/annotation_build.nf'
include { GENERATE_BED_FILE } from '../modules/annotation_build.nf'
include { GENERATE_REFSEQ_MANE_BED } from '../modules/annotation_build.nf'

workflow BUILD_ANNOTATION_REFERENCE {
    take:
    gtf
    assembly_report
    grc_fixes_assessment
    refseq_mane_annotation

    main:
    DOWNLOAD_EBV_GTF(
        params.genome.ebv_annotation_url
    )

    REMOVE_SECTIONS(
        gtf,
        channel.value([["NC_000021.9", "rRNA"],
                       ["NT_187388.1", "rRNA"],
                       ["NT_167214.1", "rRNA"]])
    )

    def rRNA_gtf_paths = [params.rRNA.NC_000021,
                          params.rRNA.NT_187388,
                          params.rRNA.NT_167214].collect { it -> "${projectDir}/${it}" }

    APPEND_RIBOSOMAL_RNA_GTFS(
        "rRNA",
        REMOVE_SECTIONS.out.gtf,
        channel.value(rRNA_gtf_paths)
    )

    CONVERT_ANNOTATION_REFSEQ_TO_UCSC(
        APPEND_RIBOSOMAL_RNA_GTFS.out.gtf,
        assembly_report
    )

    GET_TARGET_CONTIGS(
        assembly_report,
        grc_fixes_assessment
    )

    def target_contigs = GET_TARGET_CONTIGS.out
    def final_output_prefix = "annotation"

    SUBSET_GTF(
        CONVERT_ANNOTATION_REFSEQ_TO_UCSC.out.gtf,
        CONVERT_ANNOTATION_REFSEQ_TO_UCSC.out.gtf_index,
        target_contigs
    )

    APPEND_EBV_GTF(
        "EBV",
        SUBSET_GTF.out.gtf,
        DOWNLOAD_EBV_GTF.out.gtf.toList()
    )

    SORT_GTF(
        final_output_prefix, // Prefix for output files
        APPEND_EBV_GTF.out.gtf
    )

    GENERATE_BED_FILE(
        final_output_prefix,
        SORT_GTF.out.compressed_gtf
    )

    GENERATE_REFSEQ_MANE_BED(
        refseq_mane_annotation
    )

    emit:
    compressed_gtf                 = SORT_GTF.out.compressed_gtf
    compressed_gtf_index           = SORT_GTF.out.compressed_gtf_index
    uncompressed_gtf               = SORT_GTF.out.uncompressed_gtf
    bed_file                       = GENERATE_BED_FILE.out.bed_file
    refseq_mane_bed_file           = GENERATE_REFSEQ_MANE_BED.out.bed_file
}