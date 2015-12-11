# -*- coding: utf-8 -*-
# Written by Guillaume Schworer, 2015

from FakeSampler import FakeSampler
from MCres import MCres
try:
    import astropy.io.fits as pf
except:
    import pyfits as pf


def load(filename):
        ext = '.mcres.fits'
        if filename.find(ext)==-1: filename += ext
        hdulist = pf.open(filename)
        for item in hdulist:
            if item.header.get('EXT')=='MCRES':
                paramstr = [str(item.header.get('PARAM'+str(i))) for i in range(item.header.get('NPARAM', 0))]
                sampler = FakeSampler(item.data[:,:-1], item.data[:,-1])
                res = MCres(sampler=sampler, paramstr=paramstr, nwalker=int(item.header.get('NWALKER', 0)), burnInIts=int(item.header.get('NBURN', 0)), niters=int(item.header.get('NITER', 0)))
        for item in hdulist:
            if item.header.get('EXT')=='MCRESFLT':
                for i in range(item.header.get('NFILTERS', 0)):
                    res._addfilter(param=str(item.header.get('PARAM'+str(i))), v_min=float(item.header.get('MIN'+str(i))), v_max=float(item.header.get('MAX'+str(i))))
        return res

