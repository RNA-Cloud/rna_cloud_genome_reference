import pandas as pd
import pytest
from unittest.mock import patch, Mock
from rnacloud_genome_reference.common.gnomad import GnomadProvider  # replace 'your_module' with actual filename

class TestGnomadProvider:
    def test_query_gnomad_real_call(self):
        """
        This test will hit the real gnomAD GraphQL API.
        """
        provider = GnomadProvider()
        # Use a region with known data (e.g. chromosome 1, a short region)
        chrom = "X"
        start = 154380901
        stop = 154380901
        
        results = provider.query_gnomad(chrom, start, stop)
        
        assert len(results) == 3, f"Expected 3 results, got {len(results)}"

        # Expected result
        # {
        #   "variant_id": "X-154380901-CG-C",
        #   "chrom": "X",
        #   "pos": 154380901,
        #   "ref": "CG",
        #   "alt": "C",
        #   "lof_filter": "END_TRUNC",
        #   "joint": {
        #     "ac": 1,
        #     "an": 1211812,
        #     "hemizygote_count": 0,
        #     "homozygote_count": 0,
        #     "filters": []
        #   }
        # }

        assert results[2].chrom == "X", f"Expected chromosome X, got {results[2].chrom}"
        assert results[2].pos == 154380901, f"Expected position 154380901, got {results[2].pos}"
        assert results[2].ref == "CG", f"Expected reference CG, got {results[2].ref}"
        assert results[2].alt == "C", f"Expected alternate C, got {results[2].alt}"
        assert results[2].lof_filter == "END_TRUNC", f"Expected lof_filter 'END_TRUNC', got {results[2].lof_filter}"
        assert results[2].ac == 1, f"Expected allele count 1, got {results[2].ac}"
        assert results[2].an == 1211812, f"Expected allele number 1211812, got {results[2].an}"
        assert results[2].hemizygote_count == 0, f"Expected hemizygote count 0, got {results[2].hemizygote_count}"
        assert results[2].homozygote_count == 0, f"Expected homozygote count 0, got {results[2].homozygote_count}"
        assert results[2].filters == [], f"Expected no filters, got {results[2].filters}"
        assert results[2].clinvar_variation_id is None, f"Expected no ClinVar variation ID, got {results[2].clinvar_variation_id}"
        assert results[2].clinical_significance is None, f"Expected no clinical significance, got {results[2].clinical_significance}"
        assert results[2].review_status is None, f"Expected no review status, got {results[2].review_status}"
        assert results[2].filters_count == 0, f"Expected filters count 0, got {results[2].filters_count}"

        # {
        #   "variant_id": "X-154380901-C-T",
        #   "chrom": "X",
        #   "pos": 154380901,
        #   "ref": "C",
        #   "alt": "T",
        #   "lof_filter": null,
        #   "joint": {
        #     "ac": 10,
        #     "an": 1210315,
        #     "hemizygote_count": 5,
        #     "homozygote_count": 0,
        #     "filters": []
        #   }
        # },
    
        # "clinvar_variants": [
        #     {
        #     "variant_id": "X-154380901-C-T",
        #     "clinvar_variation_id": "1597495",
        #     "clinical_significance": "Likely benign",
        #     "review_status": "criteria provided, single submitter"
        #     }
        # ],

        assert results[0].chrom == "X", f"Expected chromosome X, got {results[0].chrom}"
        assert results[0].pos == 154380901, f"Expected position 154380901, got {results[0].pos}"
        assert results[0].ref == "C", f"Expected reference C, got {results[0].ref}"
        assert results[0].alt == "T", f"Expected alternate T, got {results[0].alt}"
        assert results[0].lof_filter is None, f"Expected lof_filter None, got {results[0].lof_filter}"
        assert results[0].ac == 10, f"Expected allele count 10, got {results[0].ac}"
        assert results[0].an == 1210315, f"Expected allele number 1210315, got {results[0].an}"
        assert results[0].hemizygote_count == 5, f"Expected hemizygote count 5, got {results[0].hemizygote_count}"
        assert results[0].homozygote_count == 0, f"Expected homozygote count 0, got {results[0].homozygote_count}"
        assert results[0].filters == [], f"Expected no filters, got {results[0].filters}"
        assert results[0].clinvar_variation_id == "1597495", f"Expected ClinVar variation ID 1597495, got {results[0].clinvar_variation_id}"
        assert results[0].clinical_significance == "Likely benign", f"Expected clinical significance 'Likely benign', got {results[0].clinical_significance}"
        assert results[0].review_status == "criteria provided, single submitter", f"Expected review status 'criteria provided, single submitter', got {results[0].review_status}"
        assert results[0].filters_count == 0, f"Expected filters count 0, got {results[0].filters_count}"

    def test_transform_clinvar_variants(self):
        input_data = [
                {
                    "variant_id": "5-13867994-T-TA",
                    "clinvar_variation_id": "454771",
                    "clinical_significance": "Likely benign",
                    "review_status": "criteria provided, single submitter"
                },
                {
                    "variant_id": "5-13867994-TA-T",
                    "clinvar_variation_id": "215494",
                    "clinical_significance": "Benign/Likely benign",
                    "review_status": "criteria provided, multiple submitters, no conflicts"
                }
            ]

        expected_output = {
                "5-13867994-T-TA": {
                    "clinvar_variation_id": "454771",
                    "clinical_significance": "Likely benign",
                    "review_status": "criteria provided, single submitter"
                },
                "5-13867994-TA-T": {
                    "clinvar_variation_id": "215494",
                    "clinical_significance": "Benign/Likely benign",
                    "review_status": "criteria provided, multiple submitters, no conflicts"
                }
        }

        response = GnomadProvider._transform_clinvar_variants(input_data)

        assert len(response) == len(expected_output), f"Expected {len(expected_output)} items, got {len(response)}"

        for a, b in zip(response, expected_output):
            assert a == b, f"Expected {b}, got {a}"

    def test_query_gnomad_mocked_single_variant(self):
        """
        Unit test: mock requests.post to test logic in isolation.
        """
        provider = GnomadProvider()

        mock_response = Mock()
        mock_response.json.return_value = {
            "data": {
                "region": {
                "clinvar_variants": [
                    {
                    "variant_id": "5-13867994-T-TA",
                    "clinvar_variation_id": "454771",
                    "clinical_significance": "Likely benign",
                    "review_status": "criteria provided, single submitter"
                    },
                    {
                    "variant_id": "5-13867994-TA-T",
                    "clinvar_variation_id": "215494",
                    "clinical_significance": "Benign/Likely benign",
                    "review_status": "criteria provided, multiple submitters, no conflicts"
                    },
                    {
                    "variant_id": "5-13867994-TAA-T",
                    "clinvar_variation_id": "1554235",
                    "clinical_significance": "Benign",
                    "review_status": "criteria provided, single submitter"
                    },
                    {
                    "variant_id": "5-13867994-T-TAA",
                    "clinvar_variation_id": "1533409",
                    "clinical_significance": "Benign",
                    "review_status": "criteria provided, single submitter"
                    }
                ],
                "variants": [
                    {
                    "variant_id": "5-13867994-T-TAA",
                    "chrom": "5",
                    "pos": 13867994,
                    "ref": "T",
                    "alt": "TAA",
                    "lof_filter": None,
                    "joint": {
                        "ac": 29,
                        "an": 1601142,
                        "hemizygote_count": 0,
                        "homozygote_count": 0,
                        "filters": []
                    }
                    },
                    {
                    "variant_id": "5-13867994-T-TAAAA",
                    "chrom": "5",
                    "pos": 13867994,
                    "ref": "T",
                    "alt": "TAAAA",
                    "lof_filter": None,
                    "joint": {
                        "ac": 0,
                        "an": 1597786,
                        "hemizygote_count": 0,
                        "homozygote_count": 0,
                        "filters": []
                    }
                    },
                    {
                    "variant_id": "5-13867994-TAA-T",
                    "chrom": "5",
                    "pos": 13867994,
                    "ref": "TAA",
                    "alt": "T",
                    "lof_filter": None,
                    "joint": {
                        "ac": 1,
                        "an": 1601066,
                        "hemizygote_count": 0,
                        "homozygote_count": 0,
                        "filters": []
                    }
                    },
                    {
                    "variant_id": "5-13867994-TA-T",
                    "chrom": "5",
                    "pos": 13867994,
                    "ref": "TA",
                    "alt": "T",
                    "lof_filter": None,
                    "joint": {
                        "ac": 688509,
                        "an": 1599786,
                        "hemizygote_count": 0,
                        "homozygote_count": 149799,
                        "filters": []
                    }
                    },
                    {
                    "variant_id": "5-13867994-T-C",
                    "chrom": "5",
                    "pos": 13867994,
                    "ref": "T",
                    "alt": "C",
                    "lof_filter": None,
                    "joint": {
                        "ac": 2,
                        "an": 1601180,
                        "hemizygote_count": 0,
                        "homozygote_count": 0,
                        "filters": ['test']
                    }
                    },
                    {
                    "variant_id": "5-13867994-T-TA",
                    "chrom": "5",
                    "pos": 13867994,
                    "ref": "T",
                    "alt": "TA",
                    "lof_filter": 'test',
                    "joint": {
                        "ac": 390,
                        "an": 1599388,
                        "hemizygote_count": 0,
                        "homozygote_count": 0,
                        "filters": []
                    }
                    }
                ]
                }
            }
        }
        mock_response.raise_for_status = Mock()
        
        def mock_post(*args, **kwargs):
            return mock_response
        
        with patch("requests.post", mock_post):
            results = provider.query_gnomad("5", 13867994, 13867994)

            assert len(results) == 6, f"Expected 6 result, got {len(results)}"
            assert results[0].chrom == "5", f"Expected chromosome 5, got {results[0].chrom}"
            assert results[0].pos == 13867994, f"Expected position 13867994, got {results[0].pos}"
            assert results[0].ref == "T", f"Expected reference T, got {results[0].ref}"
            assert results[0].alt == "TAA", f"Expected alternate TAA, got {results[0].alt}"
            assert results[0].ac == 29, f"Expected allele count 29, got {results[0].ac}"
            assert results[0].an == 1601142, f"Expected allele number 1601142, got {results[0].an}"
            assert results[0].hemizygote_count == 0, f"Expected hemizygote count 0, got {results[0].hemizygote_count}"
            assert results[0].homozygote_count == 0, f"Expected homozygote count 0, got {results[0].homozygote_count}"
            assert results[0].filters == [], f"Expected no filters, got {results[0].filters}"
            assert results[0].clinvar_variation_id == "1533409", f"Expected ClinVar variation ID 1533409, got {results[0].clinvar_variation_id}"
            assert results[0].clinical_significance == "Benign", f"Expected clinical significance 'Benign', got {results[0].clinical_significance}"
            assert results[0].review_status == "criteria provided, single submitter", f"Expected review status 'criteria provided, single submitter', got {results[0].review_status}"
            assert results[0].filters_count == 0, f"Expected filters count 0, got {results[0].filters_count}"

            assert results[1].chrom == "5", f"Expected chromosome 5, got {results[1].chrom}"
            assert results[1].pos == 13867994, f"Expected position 13867994, got {results[1].pos}"
            assert results[1].ref == "T", f"Expected reference T, got {results[1].ref}"
            assert results[1].alt == "TAAAA", f"Expected alternate TAAAA, got {results[1].alt}"
            assert results[1].ac == 0, f"Expected allele count 0, got {results[1].ac}"
            assert results[1].an == 1597786, f"Expected allele number 1597786, got {results[1].an}"
            assert results[1].hemizygote_count == 0, f"Expected hemizygote count 0, got {results[1].hemizygote_count}"
            assert results[1].homozygote_count == 0, f"Expected homozygote count 0, got {results[1].homozygote_count}"
            assert results[1].filters == [], f"Expected no filters, got {results[1].filters}"
            assert results[1].clinvar_variation_id is None, f"Expected no ClinVar variation ID, got {results[1].clinvar_variation_id}"
            assert results[1].clinical_significance is None, f"Expected no clinical significance, got {results[1].clinical_significance}"
            assert results[1].review_status is None, f"Expected no review status, got {results[1].review_status}"
            assert results[1].filters_count == 0, f"Expected filters count 0, got {results[1].filters_count}"

            assert results[4].filters == ['test'], f"Expected filters ['test'], got {results[4].filters}"
            assert results[4].filters_count == 1, f"Expected filters count 1, got {results[4].filters_count}"

            assert results[5].lof_filter == 'test', f"Expected lof_filter 'test', got {results[5].lof_filter}"
    
    @pytest.mark.parametrize(
        "start, stop, max_range, expected",
        [
            # Single base range
            (1, 1, 50000, [(1, 1)]),
            # Exactly max_range
            (1, 50000, 50000, [(1, 50000)]),
            # One over max_range
            (1, 50001, 50000, [(1, 50000), (50001, 50001)]),
            # Multiple full chunks plus remainder
            (100, 160000, 50000, [(100, 50099), (50100, 100099), (100100, 150099), (150100, 160000)]),
            # start == stop edge case at high value
            (123456, 123456, 50000, [(123456, 123456)]),
            # start greater than stop should produce empty list
            (10, 5, 50000, []),
        ]
    )
    def test_split_ranges(self, start: int, stop: int, max_range: int, expected: list[tuple[int, int]]) -> None:
        result = GnomadProvider._split_ranges(start, stop, max_range)
        assert result == expected, f"Expected {expected} but got {result} for range {start}-{stop} with max_range {max_range}"

    def test_fetch_gnomad_stats_for_region(self):
        provider = GnomadProvider()

        results = provider.fetch_gnomad_stats_for_region("5", 13867994, 13867994)

        if results is not None:
            assert results.shape[0] == 6, f"Expected 6 rows, got {results.shape[0]}"

            assert 'chrom' in results.columns, "Expected 'chrom' column in results"
            assert 'pos' in results.columns, "Expected 'pos' column in results"
            assert 'ref' in results.columns, "Expected 'ref' column in results"
            assert 'alt' in results.columns, "Expected 'alt' column in results"
            assert 'lof_filter' in results.columns, "Expected 'lof_filter' column in results"
            assert 'ac' in results.columns, "Expected 'ac' column in results"
            assert 'an' in results.columns, "Expected 'an' column in results"
            assert 'hemizygote_count' in results.columns, "Expected 'hemizygote_count' column in results"
            assert 'homozygote_count' in results.columns, "Expected 'homozygote_count' column in results"
            assert 'filters' in results.columns, "Expected 'filters' column in results"
            assert 'clinvar_variation_id' in results.columns, "Expected 'clinvar_variation_id' column in results"
            assert 'clinical_significance' in results.columns, "Expected 'clinical_significance' column in results"
            assert 'review_status' in results.columns, "Expected 'review_status' column in results"
            assert 'filters_count' in results.columns, "Expected 'filters_count' column in results"

    def test_fetch_gnomad_stats_for_region_no_results(self):
        provider = GnomadProvider()

        results = provider.fetch_gnomad_stats_for_region("5", 13867997, 13867997)

        assert results is None, "Expected no results for non-existent region"

        