import logging
import pandas as pd

logger = logging.getLogger(__name__)

class AssemblyReportParser:
    def __init__(self, assembly_report: str):
        self.assembly_report = assembly_report

        self.regions = self._load_assembly_report()
        self.contig_lengths = self.regions.set_index('UCSC-style-name')['Sequence-Length'].to_dict()
        self.refseq_to_ucsc_map = self._load_refseq_to_ucsc_map()
        self.ucsc_to_refseq_map = self._load_ucsc_to_refseq_map()
        self.assembly_info = self._load_assembly_info()

    def _load_assembly_info(self) -> dict[str, str]:
        logger.debug(f"Loading assembly info from {self.assembly_report}")
        assembly_info = {}

        with open(self.assembly_report, "r", encoding="utf-8") as f:
            for line in f:
                if not line.startswith("#"):
                    break

                if line.startswith("##"):
                    continue

                line = line[1:].strip()
                if not line or ":" not in line:
                    continue

                key, value = line.split(":", 1)
                assembly_info[key.strip()] = value.strip()

        logger.debug(f"Loaded assembly info: {assembly_info}")
        return assembly_info        

    def _load_assembly_report(self) -> pd.DataFrame:
        # Load the assembly report into a DataFrame
        logger.debug(f"Loading assembly report from {self.assembly_report}")
        regions = pd.read_csv(self.assembly_report, 
                    sep='\t',
                    comment='#',
                    low_memory=False,
                    header=None,
                    names=['Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name'])
        return regions

    def get_contig_range(self, ucsc_contig_name: str) -> tuple[int, int]:
        # Get the start and end positions of a contig based on its UCSC style name
        logger.debug(f"Getting contig ranges for UCSC contig name {ucsc_contig_name}")

        if ucsc_contig_name not in self.contig_lengths:
            logger.error(f"UCSC contig name {ucsc_contig_name} not found in assembly report")
            raise ValueError(f"UCSC contig name {ucsc_contig_name} not found in assembly report")

        start = 1 # 1-based
        end = self.contig_lengths[ucsc_contig_name]
        return start, end

    def _load_refseq_to_ucsc_map(self) -> dict:
        # Load the RefSeq to UCSC map from the assembly report
        logger.debug(f"Loading RefSeq to UCSC map")
        chromosome_map = self.regions.set_index('RefSeq-Accn')['UCSC-style-name'].to_dict()

        logger.debug(f"Loaded chromosome map: {chromosome_map}")
        return chromosome_map
    
    def _load_ucsc_to_refseq_map(self) -> dict:
        # Create a reverse map from UCSC to RefSeq
        logger.debug("Creating UCSC to RefSeq map")
        ucsc_to_refseq_map = {v: k for k, v in self.refseq_to_ucsc_map.items()}
        logger.debug(f"Created UCSC to RefSeq map: {ucsc_to_refseq_map}")
        return ucsc_to_refseq_map
    
    def refseq_to_ucsc(self, refseq_id: str) -> str:
        # Convert RefSeq ID to UCSC style name
        logger.debug(f"Converting RefSeq ID {refseq_id} to UCSC style name")
        ucsc_name = self.refseq_to_ucsc_map.get(refseq_id, None)
        if ucsc_name is None:
            logger.error(f"RefSeq ID {refseq_id} not found in RefSeq to UCSC map")
            raise ValueError(f"RefSeq ID {refseq_id} not found in RefSeq to UCSC map")
        return ucsc_name
    
    def ucsc_to_refseq(self, ucsc_name: str) -> str:        # Convert UCSC style name to RefSeq ID
        logger.debug(f"Converting UCSC style name {ucsc_name} to RefSeq ID")
        refseq_id = self.ucsc_to_refseq_map.get(ucsc_name, None)
        if refseq_id is None:
            logger.error(f"UCSC style name {ucsc_name} not found in UCSC to RefSeq map")
            raise ValueError(f"UCSC style name {ucsc_name} not found in UCSC to RefSeq map")
        return refseq_id
