nextflow.enable.dsl=2

process DECOMPRESS_GTF {
    tag "DECOMPRESS_GTF"

    input:
    path gtf

    output:
    path "${gtf.simpleName}.gtf", emit: gtf

    script:
    """
    set -euo pipefail

    echo "Decompressing GTF file"
    gunzip -c ${gtf} > ${gtf.simpleName}.gtf
    """
}

process REMOVE_SECTIONS {
    tag "REMOVE_SECTIONS"

    input:
    path gtf   // Compressed GTF
    val pairs  // list of [contig, biotype] pairs

    output:
    path "${gtf.simpleName}.filtered.gtf.gz", emit: gtf

    script:
    """
    set -euo pipefail

    # Print the filtering criteria
    echo "Removing multiple contig-biotype pairs from GTF: ${pairs}"

    # Copy the input GTF to a temporary file for iterative filtering
    cp ${gtf} temp.gtf.gz

    gunzip temp.gtf.gz

    # For each [contig, biotype] pair, filter out lines matching both criteria
    ${pairs.collect { pair -> 
        // Construct an awk command that excludes lines where:
        // - the first column matches the contig
        // - and the attribute field contains the specified biotype
        // Logic: keep line if it's NOT the matching biotype OR NOT the matching contig
        "awk '!(/\\_biotype \\\"${pair[1]}\\\"/) || \$1 != \"${pair[0]}\"' temp.gtf > temp2.gtf && mv temp2.gtf temp.gtf"
    }.join('\n')}

    # Rename the filtered GTF file to the final output
    mv temp.gtf ${gtf.simpleName}.filtered.gtf

    bgzip ${gtf.simpleName}.filtered.gtf
    """
}

process APPEND_GTFS {
    tag "APPEND_GTFS"

    input:
    val appended_contigs  // Plain label
    path gtf              // Compressed GTF
    val additional_gtfs   // List of additional GTF paths as strings

    output:
    path "${appended_contigs}.gtf.gz", emit: gtf

    script:
    """
    set -euo pipefail

    echo "Appending ${additional_gtfs.size()} additional GTF files to ${gtf}"

    gunzip -c ${gtf} > temp.gtf

    # Check that all additional GTF files exist
    for f in ${additional_gtfs.join(' ')}; do
        if [ ! -f "\$f" ]; then
            echo "ERROR: Additional GTF file '\$f' not found." >&2
            exit 1
        fi
    done

    # Concatenate main GTF and all additional GTFs
    cat temp.gtf ${additional_gtfs.join(' ')} > "${appended_contigs}.gtf"

    # Compress appended gtf
    bgzip ${appended_contigs}.gtf

    # Remove temp gtf
    rm temp.gtf
    """
}

process CONVERT_ANNOTATION_REFSEQ_TO_UCSC {
    tag "CONVERT_GENOME_ANNOT_REFSEQ_TO_UCSC"

    input:
    path gtf             // Compressed GTF
    path assembly_report

    output:
    path "${gtf.simpleName}_ucsc.gtf.gz", emit: gtf
    path "${gtf.simpleName}_ucsc.gtf.gz.tbi", emit: gtf_index

    script:
    """
    set -euo pipefail

    gunzip -c ${gtf} > temp.gtf

    awk '
    BEGIN { FS=OFS="\t" }
    NR==FNR {
        # in the assembly report, skip headers
        if (\$0 ~ /^#/) next
        refseq = \$7
        ucsc   = \$NF
        sub(/\\r\$/, "", refseq)
        sub(/\\r\$/, "", ucsc)
        map[refseq] = ucsc
        next
    }
    {
        # now in the GTF: if the seqname is in our map, swap it
        if (\$1 in map) \$1 = map[\$1]
        print
    }
    ' ${assembly_report} temp.gtf | bedtools sort > ${gtf.simpleName}_ucsc.gtf

    echo "Removing temp file"
    rm temp.gtf

    echo "Compressing GTF file"
    bgzip ${gtf.simpleName}_ucsc.gtf

    echo "Creating index"
    tabix ${gtf.simpleName}_ucsc.gtf.gz
    """
}

process SUBSET_GTF {
    tag "SUBSET_GTF"

    input:
    path gtf           // Compressed GTF
    path gtf_index
    val target_contigs

    output:
    path "subset.gtf.gz", emit: gtf
    path "subset.gtf.gz.tbi", emit: gtf_index

    script:
    """
    set -euo pipefail

    echo "Filtering GTF for target contigs"
    tabix ${gtf} ${target_contigs} > subset.gtf

    echo "Compressing subset GTF"
    bgzip subset.gtf

    echo "Index GTF"
    tabix subset.gtf.gz
    """
}

process SORT_GTF {
    tag "SORT_GTF"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val final_output_prefix
    path gtf                // Compressed GTF

    output:
    path "${final_output_prefix}.gtf.gz", emit: compressed_gtf
    path "${final_output_prefix}.gtf.gz.tbi", emit: compressed_gtf_index
    path "${final_output_prefix}.gtf", emit: uncompressed_gtf

    script:
    """
    set -euo pipefail

    echo "Sorting GTF file"
    /app/rnacloud_genome_reference/genome_build/scripts/gtf_sort.sh \
      ${gtf} ${final_output_prefix}.gtf

    echo "Compressing GTF file"
    bgzip -c ${final_output_prefix}.gtf > ${final_output_prefix}.gtf.gz

    echo "Indexing GTF file"
    tabix ${final_output_prefix}.gtf.gz
    """
}

process GENERATE_REFSEQ_MANE_BED {
    tag "GENERATE_REFSEQ_MANE_BED"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path refseq_mane_annotation

    output:
    path "annotation_mane.bed", emit: bed_file

    script:
    """
    set -euo pipefail

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
    #   awk '\$3=="gene"' ...
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

    echo "Generating BED file from RefSeq MANE annotation GTF"
    gtfToGenePred -ignoreGroupsWithoutExons ${refseq_mane_annotation} refseq_mane_annotation.genepred

    echo "Converting genePred to BED file"
    genePredToBed refseq_mane_annotation.genepred annotation_mane.bed
    """
}

process GENERATE_BED_FILE {
    tag "GENERATE_BED_FILE"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val final_output_prefix
    path compressed_gtf   // Compressed GTF

    output:
    path "${final_output_prefix}.bed", emit: bed_file

    script:
    """
    set -euo pipefail

    # A. BEGIN sets the field separators
    # B. The first block finds lines with string 'transcript_id "unassigned_transcript'.
    # C. Saves the gene_id and original transcript_id values from the matches in B.
    # D. Builds a new, unique ID.
    # E. We then use sub() to replace the old value with the new one.
    # F. We print the modified line, then skip to the next.
    # G. The final {print} deals with the lines that did not match B.

    awk 'BEGIN {OFS=FS="\t"}
    /transcript_id "unassigned_transcript/ {
        gene_id = ""
        trans_id = ""

        # extract gene_id
        if (match(\$9, /gene_id "[^"]+"/)) {
            tmp = substr(\$9, RSTART, RLENGTH)
            sub(/^gene_id "/, "", tmp)
            sub(/"\$/, "", tmp)
            gene_id = tmp
        }

        # extract transcript_id
        if (match(\$9, /transcript_id "unassigned_transcript_[^"]+"/)) {
            tmp = substr(\$9, RSTART, RLENGTH)
            sub(/^transcript_id "/, "", tmp)
            sub(/"\$/, "", tmp)
            trans_id = tmp
        }

        if (gene_id != "" && trans_id != "") {
            new_tid = \$1 "_" gene_id "_" trans_id
            sub(trans_id, new_tid, \$9)
        }
        print
        next
    }
    { print }' <(gunzip -c ${compressed_gtf}) > temp.gtf

    echo "Compressing GTF file"
    bgzip temp.gtf
    mv temp.gtf.gz ${final_output_prefix}.renamed_unassigned_tx.gtf.gz

    echo "Generating BED file from GTF"
    gtfToGenePred -ignoreGroupsWithoutExons ${final_output_prefix}.renamed_unassigned_tx.gtf.gz ${final_output_prefix}.renamed_unassigned_tx.genepred

    echo "Converting genePred to BED file"
    genePredToBed ${final_output_prefix}.renamed_unassigned_tx.genepred ${final_output_prefix}.bed
    """
}
