#!/usr/bin/env python
"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Loads docids into database.

Refactored by: Tong Shu Li
Last updated: 2016-02-25
"""
import json
import os

from helper.easierlife import log
from helper.easierlife import BASE_FOLDER

# location of the text corpus.
INPUT_FOLDER = BASE_FOLDER + "/data/plos_full_conll"
INPUT_NEG_FOLDER = BASE_FOLDER + "/data/NEGARTICLE"

log("~~~~~~~~~~~~~~~~~" + INPUT_FOLDER)
for docid in os.listdir(INPUT_FOLDER):
    if (docid.startswith("journal.pbio.")
        or docid.startswith("journal.pgen.")
        or docid.startswith("journal.pone.")):

        file_name = os.path.join(INPUT_FOLDER, docid)

        print json.dumps({"docid": docid.rstrip(".nlp"), "folder": file_name})

log("~~~~~~~~~~~~~~~~~" + INPUT_NEG_FOLDER)
for docid in os.listdir(INPUT_NEG_FOLDER):
    if not docid.startswith("."):
        file_name = os.path.join(INPUT_NEG_FOLDER, docid)
        print json.dumps({"docid": docid, "folder": file_name})
