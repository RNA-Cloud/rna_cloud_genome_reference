import argparse
import json
import re
import logging
import os
from contextlib import redirect_stdout

import numpy as np
import pandas as pd
import pysam

from rnacloud_genome_reference.common.gtf import GTFHandler
from rnacloud_genome_reference.common.utils import AssemblyReportParser
from rnacloud_genome_reference.genome_build.common import GRC_FIXES_QUERY

logger = logging.getLogger(__name__)

class AssemblyStatsProvider():
    def __init__(self, fasta: str):
        self._fasta = fasta
        self._fasta_loaded =  pysam.FastaFile(fasta)

    def get_contigs_count(self) -> int:
        return len(self._fasta_loaded.references)

    def get_file_size(self) -> int:
        return os.path.getsize(self._fasta)

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return json.load(f)

class ContigReportBuilder:
    def __init__(self):
        self._fasta: pysam.FastaFile
        self._assembly_report: AssemblyReportParser
        self._grc: pd.DataFrame
        self._cen_par_regions: pd.DataFrame
        self._contigs_report: pd.DataFrame
        self._rRNA_contigs: list[str]

    def with_fasta(self, fasta_path: str):
        self._fasta = pysam.FastaFile(fasta_path)
        return self

    def with_assembly_report(self, assembly_report_path: str):
        self._assembly_report = AssemblyReportParser(assembly_report_path)
        return self

    def with_grc_fixes_summary(self, grc_fixes_summary_path: str):
        self._grc = pd.read_csv(grc_fixes_summary_path, sep='\t', low_memory=False)
        return self

    def with_cen_par_regions(self, cen_par_regions_path: str):
        self._cen_par_regions = pd.read_csv(cen_par_regions_path, sep='\t', low_memory=False)
        return self
    
    def with_rRNA_contigs(self, rRNA_contigs: list[str]):
        self._rRNA_contigs = rRNA_contigs
        return self

    def build_base_report(self):
        self._contigs_report = (
            self._assembly_report.regions
            .query('`UCSC-style-name`.isin(@self._fasta.references)')
            [['UCSC-style-name', 'RefSeq-Accn', 'Sequence-Role', 'Sequence-Length']]
            .rename(columns={
                'UCSC-style-name': 'chr_ucsc',
                'RefSeq-Accn': 'chr_refseq',
                'Sequence-Role': 'role',
                'Sequence-Length': 'length'
            })
        )
        return self

    def merge_grc_fixes(self):
        grc_filtered = self._grc.query(GRC_FIXES_QUERY)

        primary_contigs_genes = (
            grc_filtered.groupby('chr_ucsc')['gene_name']
            .agg(lambda x: ', '.join(x.drop_duplicates()))
            .reset_index()
            .rename(columns={'gene_name': 'masked_genes'})
        )

        fix_contigs_genes = (
            grc_filtered.groupby('alt_chr_ucsc')['gene_name']
            .agg(lambda x: ', '.join(x.drop_duplicates()))
            .reset_index()
            .rename(columns={'alt_chr_ucsc': 'chr_ucsc', 'gene_name': 'genes_with_grc_fixes'})
        )

        self._contigs_report = self._contigs_report.merge(primary_contigs_genes, on='chr_ucsc', how='left')
        self._contigs_report = self._contigs_report.merge(fix_contigs_genes, on='chr_ucsc', how='left')
        return self

    def add_ebv_contig(self):
        ebv_contig = pd.DataFrame({
            'chr_ucsc': ['chrEBV'],
            'chr_refseq': [None],
            'role': ['EBV contig'],
            'length': [self._fasta.get_reference_length('chrEBV')]
        })
        self._contigs_report = pd.concat([self._contigs_report, ebv_contig], ignore_index=True)
        return self

    def annotate_cen_par(self):
        self._contigs_report['cen_par'] = np.where(
            self._contigs_report['chr_ucsc'].isin(self._cen_par_regions['masked_copy_chr_name']),
            'Y', None
        )
        return self
    
    def annotate_rRNA_regions(self):
        self._contigs_report['rRNA'] = np.where(
            self._contigs_report['chr_refseq'].replace(r'\.\d+', '', regex=True).isin(self._rRNA_contigs),
            'Y', None
        )
        return self

    def annotation_redundant_5S_regions(self):
        # TODO: Remove hard coded redundant 5S regions
        self._contigs_report['masked_redundant_rRNA_regions'] = np.where(
            self._contigs_report['chr_refseq'].isin(['NC_000001.11']),
            'Y', None
        )
        return self

    def build(self) -> pd.DataFrame:
        return self._contigs_report

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build genome reference report and write to Markdown.")
    # Primary inputs (now CLI params)
    parser.add_argument("--fasta", required=False, default="output/GCF_000001405.40_GRCh38.p14_rna_cloud.fasta.gz")
    parser.add_argument("--gtf", required=False, default="output/GCF_000001405.40_GRCh38.p14_rna_cloud.gtf.gz")
    parser.add_argument("--assembly-report", required=False, default="data/GCF_000001405.40_GRCh38.p14_assembly_report.txt")
    parser.add_argument("--cen-par-regions", required=False, default="data/unmasked_cognates_of_masked_CEN_PAR.txt")
    parser.add_argument("--grc-fixes-summary", required=False, default="output/grc_fixes_assessment.tsv")
    parser.add_argument("--config", required=False, default="conf/sources.json", help="Path to sources.json")
    parser.add_argument("--version", required=False, default="0.0.0", help="Version of the genome and annotation build")
    # Output
    parser.add_argument("--out", required=True, help="Output Markdown filename (e.g., report.md)")
    return parser.parse_args(argv)

def main(argv: list[str] | None = None):
    args = parse_args(argv)

    fasta = args.fasta
    gtf = args.gtf
    assembly_report = args.assembly_report
    cen_par_regions = args.cen_par_regions
    grc_fixes_summary = args.grc_fixes_summary
    config_path = args.config

    # Load config from CLI-provided path
    config = load_config(config_path)

    with open(args.out, "w", encoding="utf-8") as f, redirect_stdout(f):
        assembly_report_parser = AssemblyReportParser(assembly_report)

        assembly_info = assembly_report_parser.assembly_info
        assembly = assembly_info.get('Assembly name', 'Unknown Assembly')
        submitter = assembly_info.get('Submitter', 'Unknown Submitter')
        genbank_accession = assembly_info.get('GenBank assembly accession', 'Unknown Accession')
        refseq_accession = assembly_info.get('RefSeq assembly accession', 'Unknown Accession')

        print("# Reference Genome and Annotation")
        print(f"- Assembly: {assembly}")
        print(f"- Submitter: {submitter}")
        print(f"- GenBank Accession: {genbank_accession}")
        print(f"- RefSeq Accession: {refseq_accession}")
        print()

        print("# File manifest")
        print(f"- Version: {args.version}")
        print(f"- FASTA: {os.path.basename(fasta)}")
        print(f"- GTF: {os.path.basename(gtf)}")
        print()

        print("# Sources")
        print("## Human genome")
        print(f"- FASTA: {config['genome']['fasta_url']}")
        print(f"- GTF: {config['genome']['annotation_url']}")
        print(f"- CEN PAR Mask regions: {config['genome']['cen_par_mask_regions']}")
        print(f"- MANE Refseq annotation: {config['genome']['refseq_mane_annotation_url']}")
        print("## EBV genome")
        print(f"- FASTA: {config['genome']['ebv_fasta_url']}")
        print(f"- GTF: {config['genome']['ebv_annotation_url']}")
        print("## Additional references")
        print(f"- GRC fixes: {config['reference']['grc_fixes']}")
        print(f"- Clinically relevant genes: {config['reference']['clinically_relevant_genes']}")
        print()

        print("# Genome assembly statistics")
        assembly_stats = AssemblyStatsProvider(fasta)
        file_size = assembly_stats.get_file_size()
        file_size_mb = file_size / (1024 * 1024)

        print(f"- Contigs: {assembly_stats.get_contigs_count()}")
        print(f"- File Size: {file_size_mb:.2f} MB")
        print()

        print("# Genome contigs report")
        contig_report_builder = ContigReportBuilder()
        contigs_report = (contig_report_builder
                          .with_fasta(fasta)
                          .with_assembly_report(assembly_report)
                          .with_grc_fixes_summary(grc_fixes_summary)
                          .with_cen_par_regions(cen_par_regions)
                          .with_rRNA_contigs(list(config['rRNA'].keys()))
                          .build_base_report()
                          .merge_grc_fixes()
                          .add_ebv_contig()
                          .annotate_cen_par()
                          .annotate_rRNA_regions()
                          .annotation_redundant_5S_regions()
                          .build())
        
        print(contigs_report.fillna('').to_markdown(index=False, tablefmt="github"))
        print()

        print("# Feature count")
        gtf_handler = GTFHandler(gtf)
        feature_counts = gtf_handler.get_feature_counts()
        print(pd.DataFrame(feature_counts.items(), columns=['feature', 'count'])
              .sort_values(by='count', ascending=False)
              .to_markdown(index=False, tablefmt="github"))
        print()

        print("# Gene biotypes count")
        gene_biotype_counts = gtf_handler.get_gene_biotype_counts()
        print(pd.DataFrame(gene_biotype_counts.items(), columns=['gene_biotype', 'count'])
              .sort_values(by='count', ascending=False)
              .to_markdown(index=False, tablefmt="github"))

if __name__ == "__main__":
    main()