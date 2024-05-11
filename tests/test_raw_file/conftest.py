import pytest


@pytest.fixture(autouse=True)
def clean_raw_text_file():
    with open("tests/data/raw_text.txt", "w"):
        pass
