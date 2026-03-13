import logging

import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def simplify_and_annotate_grc_fixes(grc_fixes_path: str, assembly_report_path: str, output_path: str) -> None:
    grc_fixes_temp = pd.read_csv(grc_fixes_path, sep='\t',
                                 low_memory=False,
                                 dtype={'parent_name': str})
    logger.info(f"Loaded GRC fixes from {grc_fixes_path} with {len(grc_fixes_temp)} entries.")

    regions = pd.read_csv(assembly_report_path, 
                      sep='\t',
                      comment='#',
                      low_memory=False,
                      header=None,
                      names=['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name'])
    logger.info("Loaded assembly report with {} entries.".format(len(regions)))

    grc_fixes = grc_fixes_temp.groupby(['parent_name', 'parent_start', 'parent_stop', 'ori', 'alt_scaf_acc', 'alt_scaf_start', 'alt_scaf_stop']).agg({
                    'issue_id': lambda x: ';'.join(x.astype(str)),
                    'type': lambda x: ';'.join(x.astype(str)),
                    'summary': lambda x: ';'.join(x.astype(str)),
                    'description': lambda x: ';'.join(x.astype(str))
                }).reset_index()
    logger.info(f"Grouped GRC fixes to {len(grc_fixes)} entries.")
    
    # Merge with assembly report to get UCSC-style names
    regions['parent_name'] = regions['Sequence-Name'].where(regions['Sequence-Name'] != regions['Assigned-Molecule'], regions['Assigned-Molecule'])
    
    # Merge GRC fixes with regions to get UCSC-style names
    grc_fixes = grc_fixes.merge(regions[['parent_name','UCSC-style-name','RefSeq-Accn']], on='parent_name', how='inner')

    # Rename columns for clarity
    grc_fixes.rename(columns={'UCSC-style-name': 'chr_ucsc',
                              'RefSeq-Accn': 'chr_refseq'}, inplace=True)
    
    # Rename alt_scaf_acc to GenBank-Accn for consistency
    grc_fixes.rename(columns={'alt_scaf_acc': 'GenBank-Accn'}, inplace=True)

    # Merge with regions again to get alternative chromosome names
    grc_fixes = grc_fixes.merge(regions[['GenBank-Accn','UCSC-style-name','RefSeq-Accn']], on='GenBank-Accn', how='inner')

    # Rename columns for alternative chromosome names
    grc_fixes.rename(columns={'UCSC-style-name': 'alt_chr_ucsc',
                              'RefSeq-Accn': 'alt_chr_refseq'}, inplace=True)
    
    # Check if any columns are NA
    if grc_fixes.isna().any().any():
        logger.error("NA values found in the grc_fixes DataFrame.")
        raise ValueError("There are NA values in the grc_fixes DataFrame. Please check the data.")

    grc_fixes.to_csv(output_path, sep='\t', index=False, header=True)
    logger.info(f"Simplified GRC fixes saved to {output_path}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Simplify and annotate GRC fixes.")
    parser.add_argument("grc_fixes_path", help="Path to the GRC fixes TSV file.")
    parser.add_argument("assembly_report_path", help="Path to the assembly report TSV file.")
    parser.add_argument("output_path", help="Path to save the simplified GRC fixes TSV file.")

    args = parser.parse_args()

    simplify_and_annotate_grc_fixes(args.grc_fixes_path, args.assembly_report_path, args.output_path)