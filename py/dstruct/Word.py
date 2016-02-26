#!/usr/bin/env python
"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Word class.
"""

class Word(object):
    def __init__(self, insent_id, text, pos, ner, lemma, dep_label, dep_par,
        sentid, box):

        self.insent_id = int(insent_id) - 1

        self.word = text
        self.pos = pos
        self.ner = ner
        self.lemma = lemma
        self.dep_label = dep_label

        self.dep_par = int(dep_par) - 1
        self.sentid = int(sentid.split("_")[-1]) - 1

        # If do not do this, outputing an Array in the language will crash
        self.lemma = self.lemma.replace('"', "''")
        self.lemma = self.lemma.replace('\\', "_")

    def __repr__(self):
        return self.word

    def get_feature(self):
        return self.lemma if self.ner == 'O' else self.ner
