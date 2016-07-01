# -*- coding: utf-8 -*-
# Written by Guillaume Schworer, 2015

class FakeSampler(object):
    def __init__(self, flatchain, flatlnprobability):
        self.flatchain = flatchain
        self.flatlnprobability = flatlnprobability
