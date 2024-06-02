import tempfile
from pathlib import Path

import pytest
import requests_mock
from requests.compat import urljoin

from gene_validator import GeneValidator


BASE_URL = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/'
RAW_FILE_PATH = 'data/raw.txt'
HOMO_SAPIENS_VALIDATOR = GeneValidator("Homo sapiens")


@pytest.mark.slow
@pytest.mark.parametrize("name, expected_url", [
    ("All Mammalia", urljoin(BASE_URL, "All_Mammalia.gene_info.gz")),
    ("Bos taurus", urljoin(BASE_URL, "Bos_taurus.gene_info.gz")),
    ("Canis familiaris", urljoin(BASE_URL, "Canis_familiaris.gene_info.gz")),
    ("Homo sapiens", urljoin(BASE_URL, "Homo_sapiens.gene_info.gz")),
    ("Mus musculus", urljoin(BASE_URL, "Mus_musculus.gene_info.gz")),
    ("Pan troglodytes", urljoin(BASE_URL, "Pan_troglodytes.gene_info.gz")),
    ("Rattus norvegicus", urljoin(BASE_URL, "Rattus_norvegicus.gene_info.gz")),
    ("Sus scrofa", urljoin(BASE_URL, "Sus_scrofa.gene_info.gz")),
])
def test_good_init_name(name, expected_url):
    with requests_mock.Mocker() as m:
        m.get(expected_url, text='data')
        gene_validator = GeneValidator(name)

    assert gene_validator.url == expected_url


@pytest.mark.parametrize("name, expected_exception", [
    ('Bos_taurus.gene_info', ValueError),
    ('wrong_url', ValueError),
    (None, ValueError),
])
def test_bad_init_name(name, expected_exception):
    with pytest.raises(expected_exception):
        GeneValidator(name)


@pytest.mark.parametrize("name, expected_status_code", [
    ("Homo sapiens", 200),
    ("Mus musculus", 200),
    ("Pan troglodytes", 200),
])
def test_good_get_response(name, expected_status_code):
    gene_validator = GeneValidator(name)
    assert gene_validator.get_response().status_code == expected_status_code


def test_get_text_data():
    file = Path(RAW_FILE_PATH)
    assert file.is_file()


def test_get_text_data_is_not_empty():
    with open(RAW_FILE_PATH, 'r') as file:
        content = file.readlines()
    assert content


@pytest.mark.parametrize("name", ["Homo sapiens", "Mus musculus"])
def test_get_ref_data(name):
    gene_validator = GeneValidator(name)
    content = gene_validator.get_ref_data()
    assert content


REF_DATA = {
    'COL3A1': 'COL3A1',
    'COL5A2': 'COL5A2',
    'GCP-2': 'GCP-2'
}


def test_positive_check_values():
    valid_genes = ['COL3A1', 'COL5A2']
    empty_genes = []
    assert HOMO_SAPIENS_VALIDATOR.check_values(valid_genes, REF_DATA)
    assert HOMO_SAPIENS_VALIDATOR.check_values(empty_genes, REF_DATA)


def test_negative_check_values():
    invalid_genes = ['COL3A1', 'INV']
    assert not HOMO_SAPIENS_VALIDATOR.check_values(invalid_genes, REF_DATA)


@pytest.mark.parametrize("clean, bad, expected", [
    (["ACHE AMOT CDK5R1 CDK6 CELSR1 CNTFR CRMP1 DPYSL2 ETS2"],
     ["10090 Irebp Chrna7 Agrp"],
     "CLEAN\nACHE AMOT CDK5R1 CDK6 CELSR1 CNTFR CRMP1 DPYSL2 ETS2\nBAD\n10090 Irebp Chrna7 Agrp"),
    (["APOH APP CCND2 COL3A1 COL5A2 TUBGCP2"],
     [],
     "CLEAN\nAPOH APP CCND2 COL3A1 COL5A2 TUBGCP2\nBAD\n")
])
def test_write_to_file(clean, bad, expected):
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file_path = temp_file.name

    HOMO_SAPIENS_VALIDATOR.write_to_file(clean, bad, temp_file_path)

    with open(temp_file_path, 'r') as file:
        content = file.read()
        assert content == expected
