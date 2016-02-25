#!/usr/bin/env python
"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Loads sentences into database.
"""
import json

#from extractor.EntityExtractor_Drug import *

from helper.easierlife import log
from helper.easierlife import serialize
from helper.easierlife import deserialize
from helper.easierlife import get_inputs

for row in get_inputs():
    doc = deserialize(row["document"])
    log(doc.docid)

    for sent in doc.sents:
        sentence_text = " ".join(word.word for word in sent.words)

        print json.dumps({
            "docid": doc.docid,
            "sentid": sent.sentid,
            "sentence": serialize(sent),
            "text": sentence_text
        })
