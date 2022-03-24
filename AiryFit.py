#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 15:11:30 2021

@author: Eliana

Airy fit for light sheet data, normalised and centred
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from scipy import optimize
#%%
# define Airy function according to Siviloglou's paper OSA 2007
def airy_beam(x, x0, x_shift, a, scale):
        
    s = (-x + x_shift) / x0 #x0 arbitrary transverse scale
    
    z = 488*1e-3
    w_0 = 488 * 1e-3 # wavelength
    n = 1.51
    k = 2 * np.pi * n / w_0
    ksi = z/(k * x0**2) #normalised propagation distance
    
    t = s - (ksi/2)**2 + 1J*a*ksi
    
    ai, aip, bi, bip = special.airy(t)

    phi =  (ai * np.exp(a*s - (a/2*ksi**2) - 1J * (ksi**3/12) + 1J * (a**2 * ksi/2) + 1J * (s*ksi/2)))

    return scale * ((phi.apply(lambda r : r.real))**2 + (phi.apply(lambda r : r.real))**2)

#%%
file = 'folder/lineprofile.csv'

df = pd.read_csv(file)
#df.head()

# save the first and second columns as x and y
x = df['X']
y = df['Y']

# subtract the background
y_norm = (y - y[:50].mean())/(y.max()- y[:50].mean())

# centre the peak on zero
x_norm = x - x[np.argmax(y)]

# parameters for Airy function optimisation
x0 = 6
x_shift = -5
a = 0.7 #0.07
scala = 2
# uncomment this to trial the parameters 
# y_airy = airy_beam_normalised(x_norm, x0, x_shift, a, scala)

popt, _ = optimize.curve_fit(airy_beam, x_norm, y_norm, p0 = (x0, x_shift, a, scala))

#%%
width = [a for a in range(len(y_norm)) if y_norm[a] > np.max(y_norm)/2.]

w = []
for i in range(1,len(width)):
    if width[i-1]== width[i]-1:
        w.append(width[i-1])
        i += 1
    else:
        break
   
w.append(width[i-1])
del i
fwhm = np.absolute(x_norm[w[0]]-x_norm[w[-1]])

#%%
fig, ax = plt.subplots(figsize=(8,6))

ax.plot(x_norm, y_norm , 'k', label = 'Experimental data')
ax.plot(x_norm, airy_beam(x_norm, *popt), 'darkorange', label = 'Airy function, fit')
ax.hlines(np.max(y_norm)/2., x_norm[w[0]],x_norm[w[-1]], color = 'orangered', label = 'FWHM$_{exp}$ = %.1f $\mu m$' %fwhm)

xlabel = 'x ($\mu m$)'
ylabel = 'Normalised intensity ($a.u.$)'
title = 'Light sheet profile at y = 3.0 $mm$'

ax.set_title(title)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.legend()

ax.set_xlim((-40,100))
fig.tight_layout()


#%%
from matplotlib.font_manager import FontProperties

font = FontProperties()
font.set_family('sans-serif')

font.set_size(16)

fig, ax = plt.subplots(figsize=(4,4*1.1))

ax.plot(x_norm, y_norm , 'k', label = 'exp')
ax.plot(x_norm, airy_beam(x_norm, *popt), 'darkorange', label = 'fit')
#ax.hlines(np.max(y_norm)/2., x_norm[w[0]],x_norm[w[-1]], color = 'orangered', label = 'FWHM$_{exp}$ = %.1f $\mu m$' %fwhm)

xlabel = 'x (μm)'
ylabel = 'Normalised intensity (a.u.)'
title = 'Profile at y = 3.0 mm'

ax.set_title(title, font = font)
ax.set_xlabel(xlabel, font = font)
ax.set_ylabel(ylabel, font = font)
ax.legend(fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=16)


ax.set_xlim((-40,100))
fig.tight_layout()