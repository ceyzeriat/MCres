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
from matplotlib.pyplot import Normalize as matplotlibpyplotNormalize
import matplotlib.cm as cm

class font:
    white = '\033[97m'
    black = '\033[38;5;16m'
    gray = '\033[90m'
    red = '\033[31m'
    green = '\033[32m'
    yellow = '\033[33m'
    orange = '\033[38;5;166m'
    blue = '\033[34m'
    magenta = '\033[35m'

    nocolor = '\033[39m'

    bold = '\033[1m'
    nobold = '\033[21m'

    underlined = '\033[4m'
    nounderlined = '\033[24m'

    dim = '\033[2m'
    nodim = '\033[22m'

    normal = nodim + nobold + nobold + nocolor

    clear = chr(27)+"[2J"

    @staticmethod
    def pos(line, col):
        return "\033["+str(int(line))+";"+str(int(col))+"H"



def quantile(x, q, weights=None):
    if weights is None:
        return np.percentile(x, [100. * qi for qi in q])
    else:
        idx = np.argsort(x)
        xsorted = x[idx]
        cdf = np.add.accumulate(weights[idx])
        cdf /= cdf[-1]
        return np.interp(q, cdf, xsorted).tolist()


def bins_to_array(arr):
    arr = np.asarray(arr)
    if arr.size==1: raise Exception("Size must be greater than one")
    if arr.ndim!=1: raise Exception("Dimension must be one")
    arr = np.r_[(arr[1]+arr[0])/2. ,arr[1:-1].copy()*2]
    arr[0]
    for i in np.arange(arr[1:].size)+1:
        arr[i] -= arr[i-1]
    return arr


def colorbar(cmap="jet", cm_min=0, cm_max=1):
    cmap = cm.get_cmap(cmap)
    norm = matplotlibpyplotNormalize(cm_min, cm_max)
    mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
    mappable._A = []
    return cmap, norm, mappable
