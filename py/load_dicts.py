"""
Tong Shu Li
Last updated: 2016-02-25

Dictionary loaders for gene-gene relation extractor.
"""
import os

from file_util import read_file

def get_gold_std_docids(BASE_FOLDER):
    """Load gold standard document identifiers."""

    LOC_GS_SKIP = os.path.join(BASE_FOLDER, "data", "plos_journals_dip_mint_pmids.txt")
    LOC_10K_SKIP = os.path.join(BASE_FOLDER, "data", "plos_docids_sample_10000.txt")

    id_set_A = [line.split("\t")[0] for line in read_file(LOC_GS_SKIP)]
    id_set_B = [line.rstrip(".pdf") for line in read_file(LOC_10K_SKIP)]

    return set(id_set_A + id_set_B)
