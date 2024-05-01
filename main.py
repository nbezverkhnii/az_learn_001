import csv
import gzip
import re
from typing import List, Dict


def get_ref_data(path: str) -> Dict[str, str]:
    """
    Extracts reference data from a gzip-compressed file.

    :param path: path to file
    :return: dict synonim gene: main gene
    """
    ref: Dict[str, str] = {}

    with gzip.open(path, 'rt', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')

        for gene_info in reader:
            main_symbol: str = gene_info['Symbol']
            ref.update({main_symbol: main_symbol})
            synonyms: str = gene_info['Synonyms']
            if synonyms != '-':
                ref.update({synonym: main_symbol for synonym in synonyms.split("|")})

        return ref


def get_raw_data(path: str) -> List[str]:
    """
    Reads raw data from a file
    :param path: path to file
    :return: A list of strings containing the raw data from the file
    """
    with open(path, 'r') as file:
        return file.readlines()


def check_values(line: List[str], reference: Dict[str, str]) -> bool:
    """
    Checks if all elements in a line are present in the reference dictionary
    :param line: a string of genes
    :param reference: dict with synonyms of genes as key
    :return: bool
    """
    return all(gene in reference.keys() for gene in line)


def write_to_file(clean: List[str], bad: List[str], path: str) -> None:
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


def main() -> None:
    """
    Main function
    Reads raw and reference files along the specified paths
    checks the data in the raw file for correctness
    saves the data to a new file at the specified path
    :return: None
    """
    RAW_DATA_PATH = './data/raw.txt'
    REF_DATA_PATH = './data/Homo_sapiens.gene_info.gz'
    RESULT_PATH = './result.txt'
    clear_data: List[str] = []
    bad_data: List[str] = []

    ref_data: Dict[str, str] = get_ref_data(REF_DATA_PATH)
    raw_data: List[str] = get_raw_data(RAW_DATA_PATH)

    for idx, row in enumerate(raw_data):
        clear_line: List[str] = re.sub(r'[^\w-]', ' ', row).split()

        if not check_values(clear_line, ref_data):
            bad_data.append(str(idx) + ' ' + ' '.join(clear_line))
        else:
            clear_data.append(' '.join([ref_data[gene] for gene in clear_line]))

    write_to_file(clear_data, bad_data, RESULT_PATH)


if __name__ == '__main__':
    main()
