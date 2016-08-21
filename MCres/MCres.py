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

import numpy as np
import matplotlib.pyplot as plt
from corner import corner as cornercorner
from matplotlib.image import NonUniformImage as mplimageNonUniformImage
from scipy.stats import binned_statistic_2d
try:
    import astropy.io.fits as pf
except:
    import pyfits as pf
from time import strftime

from . import _core


class MCres(object):
    """
    Sampler can be the emcee sampler or a path to a '*.res.fits' file
    paramstr is a list of (unique) string keywords that relate to each free parameter fitted in the emcee simulation. If left None, it will be set automatically to p1, p2, etc.
    nwalkers, nwalkers, burnInIts are the emcee MCMC simulation parameters
    """
    def __init__(self, sampler, paramstr=None, nwalkers=None, niters=None, burnInIts=None):
        self._init(sampler=sampler, paramstr=paramstr, nwalkers=nwalkers, niters=niters, burnInIts=burnInIts)


    def _init(self, sampler=None, paramstr=None, nwalkers=None, niters=None, burnInIts=None):
        if paramstr is None:
            paramstr = [p+str(i) for i in range(sampler.flatchain.shape[1])]
        elif sampler.flatchain.shape[1]!=len(paramstr) or len(set(paramstr))!=len(paramstr):
            raise Exception("paramstr must be same size as parameter count and each element must be unique")
        self.nwalkers = int(0 if nwalkers is None else nwalkers)
        self.niters = int(0 if niters is None else niters)
        self.burnInIts = int(0 if burnInIts is None else burnInIts)
        self._chaindtype = [(str(item), np.float64) for item in paramstr]
        dum = np.zeros(sampler.flatchain.shape[0], dtype=self._chaindtype)
        for i, item in enumerate(paramstr):
            dum[item] = sampler.flatchain[:,i].copy()
        if hasattr(self, 'chainraw'): # if already has some iterations
            self.chainraw = np.hstack((self.chainraw, dum))
        else:
            self.chainraw = dum
        if hasattr(self, 'lnprob'): # if already has some iterations
            self.lnprob = np.ma.array(np.hstack((self.lnprob.data, sampler.flatlnprobability)))
        else: 
            self.lnprob = np.ma.array(sampler.flatlnprobability)
        self.chain = np.ma.array(self.chainraw)
        self._infmask = np.logical_not(self.lnprob.data > -np.inf)
        self.paramstr = paramstr
        self.sampler = sampler
        self.delfilters()


    def _info(self):
        return _core.font.blue+"<MCMC Fit>%s\n %s steps\n params: %s"%(_core.font.normal, self.chainraw.shape[0], self.paramstr)
    def __repr__(self):
        return self._info()
    def __str__(self):
        return self._info()


    def addfilter(self, param, v_min=None, v_max=None):
        """
        Apply a filter on a parameter
        """
        if param not in self.paramstr:
            print(_core.font.red+"ERROR: Parameter '%s' not found." % (param)+_core.font.normal)
            return
        vals = self.chain[param]
        if v_max is None: v_max = vals.max()
        if v_min is None: v_min = vals.min()
        self._addfilter(param=param, v_min=v_min, v_max=v_max)        
        print(_core.font.blue+"Applied filter on '%s'.\nMin: %f, Max: %f. %i data points left." % (param, v_min, v_max, (self.chain.mask[self.paramstr[0]]==False).sum())+_core.font.normal)
        

    def _addfilter(self, param, v_min, v_max):
        self._filtermask = ((self._filtermask) | (self.chain[param]<v_min) | (self.chain[param]>v_max)) # True === hide
        newmask = (self._filtermask | self._infmask)
        self.chain.mask = newmask
        self.lnprob.mask = newmask
        self._filters.append((param, v_min, v_max))
        if (self.chain.mask[param]==False).sum()==0:
            print(_core.font.red+"ERROR: Your filter removed all data."+_core.font.normal)


    def remfilter(self, ind):
        """
        Remove one visualization filter, of index ind in the filters list
        """
        self._filters.pop(ind)
        savefilter = self._filters
        self.delfilters()
        for param, v_min, v_max in savefilter:
            self._addfilter(param=param, v_min=v_min, v_max=v_max)


    def delfilters(self):
        """
        Remove all visualization filters applied to the data
        """
        self.chain.mask = self._infmask
        self.lnprob.mask = self._infmask
        self._filtermask = np.zeros(self.chainraw.shape, dtype=bool)
        self._filters = []


    def wrap(self, param, center=np.pi, cycle=2*np.pi):
        if param not in self.paramstr:
            print(_core.font.red+"ERROR: Parameter '%s' not found." % (param)+_core.font.normal)
            return
        cycle = float(cycle)
        at = float(center)%cycle-cycle/2.
        self.chain[param] = (self.chain[param]-at)%(cycle)+at
        print(_core.font.blue+"Wrapping applied. Window for '%s' is between %f and %f." % (param, at, at+cycle)+_core.font.normal)


    @property
    def filters(self):
        return self._filters
    @filters.setter
    def filters(self, value):
        raise Exception("Read-only. Use addfilter method")


    @property
    def chain2D(self):
        return self.chain.view(np.float64).reshape(self.chain.shape + (-1,))
    @chain2D.setter
    def chain2D(self, value):
        raise Exception("Read-only.")


    @property
    def chainraw2D(self):
        return self.chainraw.view(np.float64).reshape(self.chainraw.shape + (-1,)).data
    @chain2D.setter
    def chainraw2D(self, value):
        raise Exception("Read-only.")
    
    
    def corner(self, raw=False, bins=50, quantiles=[0.16, 0.5, 0.84], **kwargs):
        """
        Plot pretty corners for the whole simulation
        """
        if raw:
            cornercorner(self.chainraw2D, labels=self.paramstr, quantiles=quantiles, bins=bins, **kwargs)
        else:
            cornercorner(self.chain2D, labels=self.paramstr, quantiles=quantiles, bins=bins, **kwargs)


    def fitparam(self, plot=False, params=[], perbin=131, q=[0.16, 0.84], best=0.5):
        """
        Fit all or specific parameter with a gaussian, because everything is gaussian
        """
        if params==[]:
            paramiter = enumerate(self.paramstr)
            params = self.paramstr
        else:
            paramiter = []
            if not hasattr(params, '__iter__'): params = [params]
            for item in params:
                if item not in self.paramstr:
                    print("Error with %s, you must choose each element of params among %s" % (item, self.paramstr))
                else:
                    paramiter.append((self.paramstr.index(item), item))
        res = []
        q = [float(best)] + list(map(float, q))
        if plot:
            plt.figure()
            plt_x = int(np.sqrt(len(params)))
            plt_y = int(np.ceil(len(params)*1./plt_x))
        for i, (cnt, item) in enumerate(paramiter):
            data = self.chain[item].compressed()
            res.append(_core.quantile(x=data, q=q))
            if plot:
                y, x = np.histogram(data, bins=int((self.chain.mask[item]==False).sum()/float(perbin)))
                xarrayed = _core.bins_to_array(x)
                plt.subplot(plt_y, plt_x, i+1)
                plt.plot(xarrayed, y, 'b-')
                plt.xlabel(self.paramstr[cnt])
                plt.plot([res[-1][0],res[-1][0]], [0,y.max()], 'r-')
                for ind in range(len(q[1:])):
                    plt.plot([res[-1][ind+1], res[-1][ind+1]], [0, y.max()], 'b--')
        res = _core.np.asarray(res).T
        res[1:] -= res[0]
        return res[0].tolist(), [item.tolist() for item in res[1:]]


    @property
    def best(self):
        return np.asarray(self.chain[self.lnprob.argmax()]).tolist()


    @property
    def bounds(self):
        return [(self.chain[par].min(), self.chain[par].max()) for par in self.paramstr]


    def _map(self, param_x, param_y, data=None, method=np.mean, noplot=False, fig=None, ax=None, cax=None, bin_x=50, bin_y=50, cmap="jet", cm_min=None, cm_max=None, axescolor='w', polar=False, showmax=True, **kwargs):
        """
        Return a 2D histogram of the MC chain, showing the walker density per bin
        """
        if param_x not in self.paramstr or param_y not in self.paramstr:
            print("You must choose param_x and param_y among %s" % self.paramstr)
            return
        x = self.chain[param_x]
        if hasattr(x, 'compressed'): x = x.compressed()
        y = self.chain[param_y]
        if hasattr(y, 'compressed'): y = y.compressed()
        if data is not None:
            H, bin_y, bin_x = binned_statistic_2d(y, x, data, method, bins=(bin_y, bin_x))[:3]
        else:
            H, bin_y, bin_x = np.histogram2d(y, x, bins=(bin_y, bin_x))
        H[np.isnan(H)] = np.nanmin(H)
        maxd = np.unravel_index(np.argmax(H), H.shape)
        if cm_min is None: cm_min = np.min(H)
        if cm_max is None: cm_max = np.max(H)
        cmap, norm, mappable = _core.colorbar(cmap=cmap, cm_min=cm_min, cm_max=cm_max)
        if noplot: return H, bin_y, bin_x, cmap, norm, mappable
        if not polar:
            if fig is None: fig = plt.figure()
            if ax is None: ax = fig.add_subplot(111)
            im = mplimageNonUniformImage(ax, cmap=cmap, norm=norm, interpolation='bilinear')
            arrbin_x = _core.bins_to_array(bin_x)
            arrbin_y = _core.bins_to_array(bin_y)
            im.set_data(arrbin_x, arrbin_y, H)
            ax.images.append(im)
            if showmax is True: ax.plot(arrbin_x[maxd[1]], arrbin_y[maxd[0]], '^w', ms=7)
            ax.set_xlim(bin_x[0], bin_x[-1])
            ax.set_xlabel(param_x)
            ax.set_ylim(bin_y[0], bin_y[-1])
            ax.set_ylabel(param_y)
        else:
            if ax is None: fig, ax = plt.subplots(subplot_kw={'projection':'polar'})
            T,R = np.meshgrid(bin_x,bin_y)
            pax = ax.pcolormesh(T, R, H, cmap=cmap, norm=norm)
            ax.grid(True)
            if showmax is True: plt.polar(T[maxd], R[maxd], '^w', ms=7)
            ax.set_ylim(0, bin_y[-1])
            ax.set_title(param_x+' vs '+param_y)
        ax.grid(True, color=axescolor)
        ax.tick_params(axis='both', colors=axescolor)
        fig.colorbar(mappable, cax=cax)


    def MCmap(self, param_x, param_y, bin_x=50, bin_y=50, cmap="jet", cm_min=None, cm_max=None, axescolor='w', polar=False, showmax=True, **kwargs):
        """
        Return a 2D histogram of the MC chain, showing the walker density per bin
        """
        self._map(param_x=param_x, param_y=param_y, bin_x=bin_x, bin_y=bin_y, cmap=cmap, cm_min=cm_min, cm_max=cm_max, axescolor=axescolor, polar=polar, showmax=showmax, **kwargs)


    def Pbmap(self, param_x, param_y, bin_x=50, bin_y=50, cmap="jet", cm_min=None, cm_max=None, axescolor='w', polar=False, showmax=True, **kwargs):
        """
        Return a 2D histogram of the MC chain, showing the best loglikelihood per bin
        """
        self._map(param_x=param_x, param_y=param_y, data=self.lnprob.compressed(), method=np.max, bin_x=bin_x, bin_y=bin_y, cmap=cmap, cm_min=cm_min, cm_max=cm_max, axescolor=axescolor, polar=polar, showmax=showmax, **kwargs)


    def save(self, name, clobber=False, append=False):
        ext = '.mcres.fits'
        if not append:
            hdulist = pf.HDUList()
            if name.find(ext)==-1: name += ext
        else:
            hdulist = pf.open(name, mode='append')

        hdu = pf.PrimaryHDU()
        hdu.header.set('EXT', 'MCRES', comment='MCres object')
        hdu.header.set('DATE', strftime('%Y%m%dT%H%M%S'), comment='Creation Date')
        hdu.header.set('NWALKER', int(self.nwalkers), comment='Number of walkers used')
        hdu.header.set('NBURN', int(self.burnInIts), comment='Number of burn-in interation used')
        hdu.header.set('NITER', int(self.niters), comment='Number of interations used')
        hdu.header.set('NPARAM', len(self.paramstr), comment='Number of parameters')
        for i, par in enumerate(self.paramstr):
            hdu.header.set('PARAM'+str(i), par, comment='Parameter '+str(i))

        hdu.data = np.concatenate((self.chainraw2D.data, self.lnprob.data.reshape((-1, 1))), axis=1)

        hdu.header.add_comment('Data shape is (n, p+1). n is the number of steps (nwalkers*niter), p is the number of free parameters. data[0:-1] contains the values of all steps. data[-1] is the calculated log-likelihood.')
        hdu.header.add_comment('Written by Guillaume SCHWORER')
        hdulist.append(hdu)

        if len(self._filters)>0:
            hdu = pf.PrimaryHDU()
            hdu.header.set('EXT', 'MCRESFLT', comment='MCres filters')
            hdu.header.set('DATE', strftime('%Y%m%dT%H%M%S'), comment='Creation Date')
            hdu.header.set('NFILTERS', len(self._filters), comment='Creation Date')
            for i, (param, v_min, v_max) in enumerate(self._filters):
                hdu.header.set('PARAM'+str(i), param, comment='Parameter to filter')
                hdu.header.set('MIN'+str(i), v_min, comment='Lower bound of filter')
                hdu.header.set('MAX'+str(i), v_max, comment='Lower bound of filter')
            hdu.header.add_comment('Written by Guillaume SCHWORER')
            hdulist.append(hdu)
        if append:
            hdulist.flush()
            hdulist.close()
        else:
            hdulist.writeto(name, clobber=clobber)
