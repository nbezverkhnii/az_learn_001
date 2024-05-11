from pathlib import Path

import pytest


def test_raw_file():
    raw_file_path = 'data/raw.txt'
    file = Path(raw_file_path)
    assert file.is_file()


def test_raw_file_is_not_empty():
    raw_file_path = 'data/raw.txt'
    with open(raw_file_path, 'r') as file:
        content = file.readlines()
    assert content