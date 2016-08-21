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

from .FakeSampler import FakeSampler
from .MCres import MCres
try:
    import astropy.io.fits as pf
except:
    import pyfits as pf


def load(filename):
    """
    Loads a MCres result file saved with :func:`MCres.save`.

    Args:
      * filename (str): the file to open, which must have fits format. The loader will check for `EXT` header cards with `MCRRES` or `MCRESFLT` values. and return a :class:`MCres` object accordingly.

    Returns:
      * A :class:`MCres` object

    Raises:
      IOError: if both `filename` and `filename` with the `.mcres.fits` extension fail to be found.
    """
    ext = '.mcres.fits'
    card = 'EXT'
    cardres = 'MCRES'
    cardresfilter = 'MCRESFLT'
    cardnfilter = 'NFILTERS'
    nwalker = 'NWALKER'
    niter = 'NITER'
    nburn = 'NBURN'
    nparam = 'NPARAM'
    try:
        hdulist = pf.open(filename)
    except IOError:
        if filename.find(ext)==-1: filename += ext
        hdulist = pf.open(filename)
    for item in hdulist:
        if item.header.get(card)==cardres:
            paramstr = [str(item.header.get('PARAM'+str(i))) for i in range(item.header.get(nparam, 0))]
            sampler = FakeSampler(item.data[:,:-1], item.data[:,-1])
            res = MCres(sampler=sampler, paramstr=paramstr, nwalker=int(item.header.get(nwalker, 0)), burnInIts=int(item.header.get(nburn, 0)), niters=int(item.header.get(niter, 0)))
            break
    for item in hdulist:
        if item.header.get(card)==cardresfilter:
            for i in range(item.header.get(cardnfilter, 0)):
                res._addfilter(param=str(item.header.get('PARAM'+str(i))), v_min=float(item.header.get('MIN'+str(i))), v_max=float(item.header.get('MAX'+str(i))))
    return res
