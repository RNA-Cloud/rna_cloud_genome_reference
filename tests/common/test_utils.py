import pytest

from rnacloud_genome_reference.common.utils import AssemblyReportParser

@pytest.fixture
def chromosome_converter():
    # Create a ChromosomeConverter instance with a mock assembly report path
    return AssemblyReportParser(assembly_report='tests/fixtures/GCF_000001405.40_GRCh38.p14_assembly_report.txt')

def test_load_refseq_to_ucsc_map(chromosome_converter):
    # Test if the RefSeq to UCSC map is loaded correctly
    refseq_to_ucsc_map = chromosome_converter.refseq_to_ucsc_map
    assert isinstance(refseq_to_ucsc_map, dict)
    assert len(refseq_to_ucsc_map) > 0  # Ensure the map is not empty

def test_load_ucsc_to_refseq_map(chromosome_converter):
    # Test if the UCSC to RefSeq map is loaded correctly
    ucsc_to_refseq_map = chromosome_converter.ucsc_to_refseq_map
    assert isinstance(ucsc_to_refseq_map, dict)
    assert len(ucsc_to_refseq_map) > 0  # Ensure the map is not empty

def test_refseq_to_ucsc_conversion(chromosome_converter):
    # Test conversion from RefSeq ID to UCSC style name
    refseq_id = 'NC_000001.11'
    ucsc_name = chromosome_converter.refseq_to_ucsc(refseq_id)
    assert isinstance(ucsc_name, str)
    assert ucsc_name == 'chr1'

def test_ucsc_to_refseq_conversion(chromosome_converter):
    # Test conversion from UCSC style name to RefSeq ID
    ucsc_name = 'chr1'
    refseq_id = chromosome_converter.ucsc_to_refseq(ucsc_name)
    assert isinstance(refseq_id, str)
    assert refseq_id == 'NC_000001.11'

def test_invalid_refseq_to_ucsc(chromosome_converter):
    # Test conversion with an invalid RefSeq ID
    with pytest.raises(ValueError, match="RefSeq ID invalid_id not found in RefSeq to UCSC map"):
        chromosome_converter.refseq_to_ucsc('invalid_id')

def test_get_contig_range(chromosome_converter):
    # Test getting contig range for a valid UCSC contig name
    ucsc_name = 'chr1'
    start, end = chromosome_converter.get_contig_range(ucsc_name)
    assert start == 1
    assert end == 248956422

def test_get_contig_range_invalid(chromosome_converter):
    # Test getting contig range for an invalid UCSC contig name
    with pytest.raises(ValueError, match="UCSC contig name invalid_chr not found in assembly report"):
        chromosome_converter.get_contig_range('invalid_chr')  # Replace with an actual invalid UCSC contig name

def test_assembly_info(chromosome_converter):
    # Test if the assembly info is loaded correctly
    assembly_info = chromosome_converter.assembly_info
    assert isinstance(assembly_info, dict)
    assert 'Assembly name' in assembly_info
    assert assembly_info['Assembly name'] == 'GRCh38.p14'