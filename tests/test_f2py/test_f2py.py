# -*- coding: utf-8 -*-
"""
Created on Sun Dec  7 18:32:46 2014

@author: prowe
"""


# Note that the following import works even if "times3.so" does not 
# exist. It imports anything named times3.*.so
from tests.test_f2py import times3


def test():
    m = times3.times3(5)

    assert m == 15