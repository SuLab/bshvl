#!/usr/bin/env python
"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Document class.

Refactored by: Tong Shu Li
Last updated: 2016-02-25
"""
import math
import copy
import re

from helper.easierlife import *
from dstruct.Sentence import *
from dstruct.Word import *
from dstruct.Box import *

class Document(object):
    def __init__(self, doc_id):
        self.docid = doc_id
        self.sents = []
        self.sents.append(Sentence())

    def push_word(self, word):
        if not self.sents[-1].push_word(word):
            self.sents.append(Sentence())
            self.sents[-1].push_word(word)
