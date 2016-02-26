#!/usr/bin/env python
"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Sentence class.

Refactored by: Tong Shu Li
Last updated: 2016-02-25
"""
from dstruct.Word import Word

class Sentence(object):

    def __init__(self):
        self.sentid = None
        self.words = []

        # to avoid bad parse tree that have self-recursion
        self.MAXLEN = 1000

    def __repr__(self):
        return " ".join(word.word for word in self.words)

    def push_word(self, word):
        assert isinstance(word, Word)

        # time to create a new Sentence object
        if self.sentid is not None and self.sentid != word.sentid:
            return False

        if self.sentid is None:
            self.sentid = word.sentid

        self.words.append(word)
        return True

#-------------------------------------------------------------------------------

    def get_path_till_root(self, word_index):
        """Get all words in the dependency tree from word_index to the root."""
        path = []
        c = word_index
        MAXLEN = self.MAXLEN
        while MAXLEN > 0:
            MAXLEN -= 1
            try:
                if c == -1:
                    break

                path.append(c)
                c = self.words[c].dep_par
            except:
                break

        return path


    def get_common_ancestor(self, pathA, pathB):
        """This function finds the least common ancestor of two paths returned by
        get_path_till_root.

        This function by Mallory does __not__ have the off-by-one error that the
        get_common_root() function in dependency_path.py has.

        Therefore if possible use this version, and not that version.

        Tong Shu Li
        2016-02-25
        """
        pos = 1
        while (pos <= len(pathA) and pos <= len(pathB)
            and pathA[-pos] == pathB[-pos]):
            pos += 1

        return pathA[-pos+1] if pos > 1 else None


    def get_direct_dependency_path_between_words(self, idx1, idx2):
        """Given two word idx1 and idx2, where idx2 is the parent of idx1, return the
        words on the dependency path."""

        words_on_path = []
        c = idx1
        MAXLEN = self.MAXLEN

        while MAXLEN > 0:
            MAXLEN -= 1
            try:
                if c == -1:
                    break

                if c == idx2:
                    break

                if c == idx1:
                    # we do not include the word of idx1
                    words_on_path.append(str(self.words[c].dep_label))
                else:
                    words_on_path.append(str(self.words[c].dep_label) + "|" + self.words[c].get_feature())

                c = self.words[c].dep_par
            except:
                break

        return words_on_path


    def get_word_dep_path(self, idx1, idx2):
        """Given two word idx1 and idx2, return the dependency path feature between them."""
        path1 = self.get_path_till_root(idx1)
        path2 = self.get_path_till_root(idx2)

        parent = self.get_common_ancestor(path1, path2)

        words_from_idx1_to_parents = self.get_direct_dependency_path_between_words(idx1, parent)
        words_from_idx2_to_parents = self.get_direct_dependency_path_between_words(idx2, parent)

        return "-".join(words_from_idx1_to_parents) + "@" + "-".join(words_from_idx2_to_parents)


    def get_prev_wordobject(self, mention):
        pos = mention.prov_words[0].insent_id
        return self.words[pos - 1] if pos > 0 else None


    def dep_parent(self, mention):
        begin = mention.prov_words[0].insent_id
        end = mention.prov_words[-1].insent_id

        left_path = [
            self.get_word_dep_path(i, j)
            for j in range(begin)
                for i in range(begin, end+1)
        ]

        right_path = [
            self.get_word_dep_path(i, j)
            for j in range(end+1, len(self.words))
                for i in range(begin, end+1)
        ]

        return sorted(left_path + right_path, key = len)[:5]


    def dep_path(self, entity1, entity2):
        begin1 = entity1.prov_words[0].insent_id
        end1 = entity1.prov_words[-1].insent_id
        begin2 = entity2.prov_words[0].insent_id
        end2 = entity2.prov_words[-1].insent_id

        paths = []
        for idx1 in range(begin1, end1+1):
            for idx2 in range(begin2, end2+1):
                paths.append(self.get_word_dep_path(idx1, idx2))

        # we pick the one that is shortest
        path = ""
        ll = 100000000
        for p in paths:
            if len(p) < ll:
                path = p
                ll = len(p)

        return path
