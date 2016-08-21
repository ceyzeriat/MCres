#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#  
#  MCres - Results post-processor for emcee MCMC simulations
#  Copyright (C) 2016  Guillaume Schworer
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  
#  For any information, bug report, idea, donation, hug, beer, please contact
#    guillaume.schworer@obspm.fr
#
###############################################################################
"""
Easy stuff to process MCMC emcee results
"""

#__all__ = ['MCres', 'version', 'FakeSampler','load']

from .MCres import MCres
from .FakeSampler import FakeSampler
from .load import load
from ._version import __version__, __major__, __minor__, __micro__

