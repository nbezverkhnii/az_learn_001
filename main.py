from typing import List, Dict
import gzip
import re
import csv


def get_ref_data(path: str) -> Dict[str, str]:
    ref: Dict[str, str] = {}

    with gzip.open(path, 'rt', encoding='utf-8') as f:  # Открываем файл в текстовом режиме
        reader = csv.DictReader(f, delimiter='\t')

        for gene_info in reader:
            main_symbol: str = gene_info['Symbol']
            ref.update({main_symbol: main_symbol})
            synonyms: str = gene_info['Synonyms']
            if synonyms != '-':
                ref.update({
                    synonym: main_symbol
                    for synonym in synonyms.split("|")
                })

        return ref


def get_raw_data(path: str) -> List[str]:
    with open(path, 'r') as f:
        return f.readlines()


def check_values(line: List[str], ref: Dict[str, str]) -> bool:
    return all(el in ref.keys() for el in line)


def write_to_file(clean, bad) -> None:
    with open(RESULT_PATH, 'w') as f:
        f.write("CLEAN\n")
        f.write("\n".join(clean))
        f.write("\nBAD\n")
        f.write("\n".join(bad))


if __name__ == '__main__':
    RAW_DATA_PATH = './data/raw.txt'
    REF_DATA_PATH = './data/Homo_sapiens.gene_info.gz'
    RESULT_PATH = './result.txt'
    clear_data = []
    bad_data = []

    ref_data: Dict[str, str] = get_ref_data(REF_DATA_PATH)
    raw_data: List[str] = get_raw_data(RAW_DATA_PATH)

    for idx, row in enumerate(raw_data):
        clear_line: List[str] = re.sub(r'[^\w-]', ' ', row).split()

        if not check_values(clear_line, ref_data):
            bad_data.append(str(idx) + ' ' + ' '.join(clear_line))
        else:
            clear_data.append(' '.join([ref_data[gene] for gene in clear_line]))

    write_to_file(clear_data, bad_data)
