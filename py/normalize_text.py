"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Refactored by: Tong Shu Li (tongli@scripps.edu)
Last updated: 2016-02-25
"""
import re

def normalize(word):
    """
    Name: normalize
    Input: word
    Return: normalized word

    Normalization for word characters
    """
    return word.encode("ascii", "ignore").replace("'", '_').replace('{', '-_-').replace('}','-__-').replace('"', '-___-').replace(', ,', ',__')

def normalize_utf(word):
    """
    Name: normalize_utf
    Input: text word
    Return: normalized word

    Replaces common UTF codes with appropriate characters.
    """
    word = re.sub('\\xe2\\x80\\x94', '-', word)
    word = re.sub('\\xef\\xac\\x81', 'fi', word)
    word = re.sub('\\xc2\\xb0', "DEGREE", word)
    word = re.sub('\\xe2\\x80\\x99', "\'", word)
    word = re.sub('\\xef\\xac\\x82', "fl", word)
    word = re.sub('\\xc2\\xa3', 'POUND', word)
    word = re.sub('\\xe2\\x80\\x98', "\'", word)
    word = re.sub('\\xe2\\x80\\x9c', "\"", word)
    word = re.sub('\\xe2\\x80\\x9d', "\"", word)
    word = re.sub('\\xe2\\x80\\x93', "-", word)
    word = re.sub('\\xe2\\x80\\x94', "--", word)
    word = re.sub('\\xe2\\x80\\xa6', "...", word)
    word = re.sub('\\xc2\\x82', ',', word)       # High code comma
    word = re.sub('\\xc2\\x84', ',,', word)       # High code double comma
    word = re.sub('\\xc2\\x85', '...', word)     # Tripple dot
    word = re.sub('\\xc2\\x88', '^', word)       # High carat
    word = re.sub('\\xc2\\x91', '\'', word)    # Forward single quote
    word = re.sub('\\xc2\\x92', '\'', word)    # Reverse single quote
    word = re.sub('\\xc2\\x93', '\"', word)    # Forward double quote
    word = re.sub('\\xc2\\x94', '\"', word)    # Reverse double quote
    word = re.sub('\\xc2\\x95', '_', word)
    word = re.sub('\\xc2\\x96', '-', word)       # High hyphen
    word = re.sub('\\xc2\\x97', '--', word)      # Double hyphen
    word = re.sub('\\xc2\\x99', '_', word)
    word = re.sub('\\xc2\\xa0', '_', word)
    word = re.sub('\\xc2\\xa6', '|', word)       # Split vertical bar
    word = re.sub('\\xc2\\xab', '<<', word)      # Double less than
    word = re.sub('\\xc2\\xbb', '>>', word)      # Double greater than
    word = re.sub('\\xc2\\xbc', '1/4', word)     # one quarter
    word = re.sub('\\xc2\\xbd', '1/2', word)     # one half
    word = re.sub('\\xc2\\xbe', '3/4', word)     # three quarters
    word = re.sub('\\xca\\xbf', '\'', word)    # c-single quote
    word = re.sub('\\xcc\\xa8', '', word)        # modifier - under curve
    word = re.sub('\\xcc\\xb1', '', word)         # modifier - under line
    word = re.sub('\\xc2\\xa7', 'CODE', word)

    return word
