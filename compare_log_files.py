#!/usr/bin/env python3

# quick script to compare DeepDive logging output for gene_relations.py refactor
# run from $DEEPDIVE_HOME/app/deepdive_genegene/ with `python3 test_gene_relations_log_output.py original_log_file.txt`

# sandip chatterjee, 2/23/16

import os
import re
import sys
import tempfile
import subprocess

def main():

    original_log_file = sys.argv[1]
    deepdive_home_dir = os.getenv('DEEPDIVE_HOME')
    new_log_file = '/'.join((deepdive_home_dir,
                             'app',
                             'bshvl',
                             'run',
                             'LATEST',
                             'log.txt'
                             ))

    diff_filename = 'results.txt'
    with tempfile.NamedTemporaryFile() as t1, tempfile.NamedTemporaryFile() as t2, open(diff_filename, 'w') as f:

        orig_lines = get_relevant_lines(original_log_file)
        new_lines = get_relevant_lines(new_log_file)
        t1.write(''.join(orig_lines).encode('utf-8'))
        t1.flush()
        t2.write(''.join(new_lines).encode('utf-8'))
        t2.flush()

        try:
            diff_output = subprocess.check_output(['diff', t1.name, t2.name])
            print('Files are the same')
        except subprocess.CalledProcessError as e:
            diff_output = e.output

            print('Files are different, diff saved in file '+diff_filename)
            f.write(diff_output.decode('utf-8'))

def get_relevant_lines(logfile_path, regex_prefix='&&&&&'):

    ''' return lines that contain regex_prefix at the beginning of
        the line, after all of the timestamp/process executor info
    '''

    useful_file_lines = []

    with open(logfile_path) as f:
        for line in f:
            matched_index = line.find(regex_prefix)
            if matched_index != -1:
                useful_file_lines.append(line[matched_index-1:].replace(regex_prefix, ''))

    return useful_file_lines

if __name__ == '__main__':
    main()