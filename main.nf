nextflow.enable.dsl=2

include { DOWNLOAD_GENOME_AND_REFERENCES } from './subworkflows/download_genome_and_references.nf'
include { EXTRACT_GENES } from './modules/grc_fixes.nf'
include { GRC_FIXES_ASSESSMENT } from './subworkflows/grc_fixes_assessment.nf'
include { SPLICE_SITE_GNOMAD_FREQ } from './subworkflows/splice_site_gnomad_freq.nf'
include { BUILD_GENOME_REFERENCE } from './subworkflows/genome_build.nf'
include { BUILD_ANNOTATION_REFERENCE } from './subworkflows/annotation_build.nf'
include { CALCULATE_MD5_SUMMARY } from './modules/common.nf'
include { VALIDATE_GENOME_ANNOTATION } from './modules/validate.nf'
include { GENOME_AND_ANNOTATION_REPORT } from './modules/validate.nf'

workflow {
    DOWNLOAD_GENOME_AND_REFERENCES()

    EXTRACT_GENES(
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf
    )

    GRC_FIXES_ASSESSMENT(
        DOWNLOAD_GENOME_AND_REFERENCES.out.grc_fixes,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_fai_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_gzi_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report,
        EXTRACT_GENES.out.genes,
        DOWNLOAD_GENOME_AND_REFERENCES.out.clinically_relevant
    )

    SPLICE_SITE_GNOMAD_FREQ(
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf_index,
        EXTRACT_GENES.out.genes,
        DOWNLOAD_GENOME_AND_REFERENCES.out.clinically_relevant,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report
    )
    
    BUILD_GENOME_REFERENCE(
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_fai_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.fasta_gzi_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report,
        GRC_FIXES_ASSESSMENT.out.grc_fixes_assessment,
        DOWNLOAD_GENOME_AND_REFERENCES.out.cen_par_mask_regions,
        DOWNLOAD_GENOME_AND_REFERENCES.out.ebv_fasta
    )

    BUILD_ANNOTATION_REFERENCE(
        DOWNLOAD_GENOME_AND_REFERENCES.out.gtf,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report,
        GRC_FIXES_ASSESSMENT.out.grc_fixes_assessment,
        DOWNLOAD_GENOME_AND_REFERENCES.out.refseq_mane_annotation,
        DOWNLOAD_GENOME_AND_REFERENCES.out.ebv_gtf
    )

    // These are regions of the genome that are already masked in the NCBI assembly
    // It is necessary to ensure validation scripts are aware of these regions
    def ncbi_assembly_masked_regions_paths = [params.ncbi_assembly_masked_regions.chr15_KN538374v1_fix].collect { it -> "${projectDir}/${it}" }

    VALIDATE_GENOME_ANNOTATION(
        BUILD_GENOME_REFERENCE.out.compressed_fasta,
        BUILD_GENOME_REFERENCE.out.compressed_fasta_fai_index,
        BUILD_ANNOTATION_REFERENCE.out.compressed_gtf,
        BUILD_ANNOTATION_REFERENCE.out.compressed_gtf_index,
        BUILD_GENOME_REFERENCE.out.mask_regions_bed,
        BUILD_GENOME_REFERENCE.out.unmask_regions_bed,
        BUILD_ANNOTATION_REFERENCE.out.bed_file,
        channel.value(ncbi_assembly_masked_regions_paths),
        DOWNLOAD_GENOME_AND_REFERENCES.out.refseq_mane_annotation,
        BUILD_ANNOTATION_REFERENCE.out.refseq_mane_bed_file
    )

    genome_and_annotation_version = System.getenv("GENOME_AND_ANNOTATION_VERSION") ?: "0.0.0"

    GENOME_AND_ANNOTATION_REPORT(
        BUILD_GENOME_REFERENCE.out.compressed_fasta,
        BUILD_GENOME_REFERENCE.out.compressed_fasta_gzi_index,
        BUILD_GENOME_REFERENCE.out.compressed_fasta_fai_index,
        BUILD_ANNOTATION_REFERENCE.out.compressed_gtf,
        BUILD_ANNOTATION_REFERENCE.out.compressed_gtf_index,
        DOWNLOAD_GENOME_AND_REFERENCES.out.assembly_report,
        DOWNLOAD_GENOME_AND_REFERENCES.out.cen_par_mask_regions,
        GRC_FIXES_ASSESSMENT.out.grc_fixes_assessment,
        "${projectDir}/conf/sources.json",
        genome_and_annotation_version
    )

    def final_outputs = GRC_FIXES_ASSESSMENT.out.grc_fixes_assessment.merge(
        SPLICE_SITE_GNOMAD_FREQ.out.splice_site_pop_freq
        , BUILD_GENOME_REFERENCE.out.compressed_fasta
        , BUILD_GENOME_REFERENCE.out.compressed_fasta_fai_index
        , BUILD_GENOME_REFERENCE.out.compressed_fasta_gzi_index
        , BUILD_GENOME_REFERENCE.out.mask_regions_bed
        , BUILD_GENOME_REFERENCE.out.unmask_regions_bed
        , BUILD_ANNOTATION_REFERENCE.out.compressed_gtf
        , BUILD_ANNOTATION_REFERENCE.out.compressed_gtf_index
        , BUILD_ANNOTATION_REFERENCE.out.uncompressed_gtf
        , BUILD_ANNOTATION_REFERENCE.out.bed_file
        , BUILD_ANNOTATION_REFERENCE.out.refseq_mane_bed_file
    )

    CALCULATE_MD5_SUMMARY(final_outputs)
}