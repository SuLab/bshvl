"""
Tong Shu Li
Last updated: 2016-02-25

Simple file utilities.
"""

def read_file(fname):
    with open(fname, "r") as fin:
        for line in fin:
            yield line.rstrip('\n')
