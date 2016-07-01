# -*- coding: utf-8 -*-
# Written by Guillaume Schworer, 2015

from FakeSampler import FakeSampler
from MCres import MCres
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
