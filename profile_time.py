from utils import time_of_function
from gene_validator import GeneValidator


@time_of_function
def profile_time(url: str) -> None:
    gene_validator = GeneValidator(url)
    gene_validator.validate()


profile_time('All_Mammalia.gene_info.gz')
profile_time('Homo_sapiens.gene_info.gz')

# Время выполнения: 304.701 сек.
# Время выполнения: 8.075 сек.

# coverage 49%
