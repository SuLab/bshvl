#!/usr/bin/env python
"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Document class.

Refactored by: Tong Shu Li
Last updated: 2016-02-25
"""
from dstruct.Sentence import Sentence

class Document(object):
    def __init__(self, doc_id):
        self.docid = doc_id
        self.sents = [Sentence()]

    def push_word(self, word):
        if not self.sents[-1].push_word(word):
            self.sents.append(Sentence())
            self.sents[-1].push_word(word)
