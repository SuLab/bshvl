#!/usr/bin/env python
"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Loads documents into database.
"""
import fileinput

from helper.easierlife import log
from helper.easierlife import serialize

from dstruct.Document import Document
from dstruct.Word import Word

for row in fileinput.input():
    docid, folder = row.rstrip("\n").split("\t")

    log(docid)

    doc = Document(docid)

    for line in open(folder):
        ss = line.rstrip().split("\t")

        if len(ss) < 3:
            continue

        insent_id, word, pos, ner, lemma, deppath, deppar, sentid, box = ss

        doc.push_word(Word(insent_id, word, pos, ner, lemma, deppath, deppar, sentid, box))

    print "\t".join(["\\N", docid, serialize(doc)])
