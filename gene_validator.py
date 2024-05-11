import csv
import gzip
import re
import sys
import tempfile
from typing import Optional

import requests
from requests import HTTPError


class GeneValidator:
    raw_data_path: str = './data/raw.txt'
    chunk_size: int = 512

    def __init__(self, path: str, base_url: str = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/',
                 result_path: str = './data/result.txt'):
        self.base_url = base_url
        self.result_path = result_path
        self.url = self.base_url + path
        self._response = self._get_response(self.url)

    def _get_response(self, url: str) -> Optional[requests.Response]:
        """
        Receives a response from the server at the given URL.
        :param url: url to get data from
        :return: bytes
        """
        try:
            response: requests.Response = requests.get(self.url)
            response.raise_for_status()
        except HTTPError:
            print(f'{url} is not a valid')
            raise

        return response

    @staticmethod
    def get_text_data(path: str) -> list[str]:
        """
        Reads raw data from a text file
        :param path: path to file
        :return: A list of strings containing the raw data from the file
        """
        with open(path, 'r') as file:
            return file.readlines()

    def get_ref_data(self) -> dict[str, str]:
        """
        Extracts reference data from responce to temp_file from a gzip-compressed file.
        Loads data to a dictionary

        :return: dict synonim gene: main gene
        """
        ref: dict[str, str] = {}
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            with self._response as fd:
                # for chunk in fd.iter_content(chunk_size=self.chunk_size):
                #     temp_file.write(chunk)
                temp_file.write(fd.content)

            with gzip.open(temp_file.name, 'rt', encoding='utf-8') as file:
                reader: csv.DictReader = csv.DictReader(file, delimiter='\t')

                for gene_info in reader:
                    main_symbol: str = gene_info['Symbol']
                    ref.update({main_symbol: main_symbol})
                    synonyms: str = gene_info['Synonyms']
                    if synonyms != '-':
                        ref.update({synonym: main_symbol for synonym in synonyms.split("|")})

            return ref

    @staticmethod
    def check_values(line: list[str], reference: dict[str, str]) -> bool:
        """
        Checks if all elements in a line are present in the reference dictionary
        :param line: a string of genes
        :param reference: dict with synonyms of genes as key
        :return: bool
        """
        return all(gene in reference.keys() for gene in line)

    @staticmethod
    def write_to_file(clean: list[str], bad: list[str], path: str) -> None:
        """
        Writes the contents of the clean and bad lists to a file at the specified path.
        At the beginning of the file, the label "CLEAN" is recorded,
        followed by the contents of the clean list,
        then the label "BAD" and the contents of the bad list are recorded.

        Example:
        CLEAN
        APOH APP CCND2 COL3A1 COL5A2 TUBGCP2
        BAD
        2 10090 Irebp Chrna7 Agrp
        3 9606 A1BG AACT ADAR3 ADDB fffe

        :param clean: list of strings with clean data
        :param bad: list of strings with bad data
        :param path: path to save the result file
        :return: None
        """
        with open(path, 'w') as file:
            file.write("CLEAN\n")
            file.write("\n".join(clean))
            file.write("\nBAD\n")
            file.write("\n".join(bad))

    def validate(self):
        """
        Reads the raw file at the specified path
        and compares it with the reference specified during class initialization
        checks the data in the raw file for correctness
        saves the data to a new file at the specified path
        :return: None
        """
        clear_data: list[str] = []
        bad_data: list[str] = []

        ref_data: dict[str, str] = self.get_ref_data()
        raw_data: list[str] = self.get_text_data(self.raw_data_path)

        for idx, row in enumerate(raw_data):
            clear_line: list[str] = re.sub(r'[^\w-]', ' ', row).split()

            if not self.check_values(clear_line, ref_data):
                bad_data.append(str(idx) + ' ' + ' '.join(clear_line))
            else:
                clear_data.append(' '.join([ref_data[gene] for gene in clear_line]))

        self.write_to_file(clear_data, bad_data, self.result_path)


if __name__ == '__main__':
    filename: str = sys.stdin.readline().strip()
    gene_validator: GeneValidator = GeneValidator(filename)
    gene_validator.validate()
