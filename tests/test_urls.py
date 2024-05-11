import pytest
import requests
import responses
from requests import HTTPError

from gene_validator import GeneValidator


@responses.activate
@pytest.mark.parametrize("url, expected_status_code", [
    ("All_Mammalia.gene_info.gz", 200),
    ("Bos_taurus.gene_info.gz", 200),
    ("Canis_familiaris.gene_info.gz", 200),
    ("Homo_sapiens.gene_info.gz", 200),
    ("Mus_musculus.gene_info.gz", 200),
    ("Pan_troglodytes.gene_info.gz", 200),
    ("Rattus_norvegicus.gene_info.gz", 200),
    ("Sus_scrofa.gene_info.gz", 200),
])
def test_good_url(url, expected_status_code):
    base_url = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/'
    responses.add(
        method=responses.GET,
        url=base_url + url,
        status=requests.codes.ok,
        match_querystring=True
    )
    gene_validator = GeneValidator(url)

    assert gene_validator._response.status_code == expected_status_code


@pytest.mark.parametrize("url, expected_exception", [
    ('Bos_taurus.gene_info', HTTPError),
    ('wrong_url', HTTPError),
    (None, TypeError),
    (123, TypeError)

])
def test_bad_url(url, expected_exception):
    with pytest.raises(expected_exception):
        GeneValidator(url)
