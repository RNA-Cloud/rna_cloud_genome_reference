nextflow.enable.dsl=2

//
// Pull in all the DOWNLOAD_* processes (path is relative to this file)
//
include { DOWNLOAD_AND_INDEX_GENOME }   from '../modules/download.nf'
include { DOWNLOAD_AND_INDEX_GTF }      from '../modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_ASSEMBLY_REPORT }            from '../modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_GRC_FIXES }                  from '../modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_CLINICALLY_RELEVANT_GENES }  from '../modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_CEN_PAR_MASK_REGIONS }       from '../modules/download.nf'
include { DOWNLOAD_FILE as DOWNLOAD_REFSEQ_MANE_ANNOTATION }     from '../modules/download.nf'
include { DOWNLOAD_EBV } from '../modules/download.nf'

//
// A DSL2 sub-workflow that reads params.* directly
//
workflow DOWNLOAD_GENOME_AND_REFERENCES {
    main:
    println "Downloading genome and references with the following parameters:"
    println "Genome FASTA URL              : ${params.genome.fasta_url}"
    println "Genome no-alt FASTA URL       : ${params.genome.no_alt_fasta_url}"
    println "Genome annotation URL         : ${params.genome.annotation_url}"
    println "Genome assembly report URL    : ${params.genome.assembly_report}"
    println "CEN-PAR mask regions URL      : ${params.genome.cen_par_mask_regions}"
    println "Reference GRC fixes URL       : ${params.reference.grc_fixes}"
    println "Clinically relevant genes URL : ${params.reference.clinically_relevant_genes}"
    println "Refseq MANE annotation URL    : ${params.genome.refseq_mane_annotation_url}"

    // build channels from params inside the sub-workflow
    DOWNLOAD_AND_INDEX_GENOME(Channel.from(params.genome.fasta_url))
    DOWNLOAD_AND_INDEX_GTF(Channel.from(params.genome.annotation_url))
    DOWNLOAD_ASSEMBLY_REPORT(Channel.from(params.genome.assembly_report))
    DOWNLOAD_GRC_FIXES(Channel.from(params.reference.grc_fixes))
    DOWNLOAD_CLINICALLY_RELEVANT_GENES(Channel.from(params.reference.clinically_relevant_genes))
    DOWNLOAD_CEN_PAR_MASK_REGIONS(Channel.from(params.genome.cen_par_mask_regions))
    DOWNLOAD_REFSEQ_MANE_ANNOTATION(Channel.from(params.genome.refseq_mane_annotation_url))
    DOWNLOAD_EBV(Channel.from(params.genome.no_alt_fasta_url))
    
    emit:
    // expose the exact same seven channels as before
    fasta                 = DOWNLOAD_AND_INDEX_GENOME.out.fasta
    fasta_fai_index       = DOWNLOAD_AND_INDEX_GENOME.out.fasta_fai_index
    fasta_gzi_index       = DOWNLOAD_AND_INDEX_GENOME.out.fasta_gzi_index
    ebv_fasta             = DOWNLOAD_EBV.out.fasta
    gtf                   = DOWNLOAD_AND_INDEX_GTF.out.gtf
    gtf_index             = DOWNLOAD_AND_INDEX_GTF.out.gtf_index
    assembly_report       = DOWNLOAD_ASSEMBLY_REPORT.out.downloaded_file
    grc_fixes             = DOWNLOAD_GRC_FIXES.out.downloaded_file
    clinically_relevant   = DOWNLOAD_CLINICALLY_RELEVANT_GENES.out.downloaded_file
    cen_par_mask_regions  = DOWNLOAD_CEN_PAR_MASK_REGIONS.out.downloaded_file
    refseq_mane_annotation = DOWNLOAD_REFSEQ_MANE_ANNOTATION.out.downloaded_file
}