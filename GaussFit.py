#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 15:11:30 2021

@author: Eliana

Airy fit for light sheet data, normalised
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from scipy import optimize


#%%
file = 'folder/lineprofile.csv'

df = pd.read_csv(file)#, names=column_names, skiprows=[0])
df.head()
# save the first and second columns as x and y
x = df['X']
y = df['Y']

#convert px to um
#px_size = 4.456 #um/px
#x = x * 4.456 
#convert to mm
#x = x* 1e3

# subtract the background
y_norm = (y - y[:50].mean())/(y.max()- y[:50].mean())

x_norm = x - x[np.argmax(y)]

# define gaussian function for optimisation fit
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(x - mean)**2 / 2 / stddev**2)

#popt, _ = optimize.curve_fit(gaussian, x_norm, y_norm, p0 = (100, 50*1e3, 150*1e3))
popt, _ = optimize.curve_fit(gaussian, x_norm, y_norm, p0 = (100, 50, 150))


plt.figure(1)
plt.plot(x_norm, y_norm)
plt.plot(x_norm, gaussian(x_norm, *popt))

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
from matplotlib.font_manager import FontProperties

font = FontProperties()
font.set_family('sans-serif')
font.set_size(16)


fig, ax = plt.subplots(figsize=(4,4*1.1))

ax.plot(x_norm, y_norm , 'k', label = 'exp')
ax.plot(x_norm, gaussian(x_norm, *popt), 'darkorange', label = 'fit')
#ax.hlines(np.max(y_norm)/2., x_norm[w[0]],x_norm[w[-1]], color = 'orangered', label = 'FWHM$_{exp}$ = %.1f $\mu m$' %fwhm)

xlabel = 'x (μm)'
ylabel = 'Normalised intensity (a.u.)'
title = 'Profile at y = 0 mm'

ax.set_title(title, font = font)
ax.set_xlabel(xlabel, font = font)
ax.set_ylabel(ylabel, font = font)
ax.legend(fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=16)


ax.set_xlim((-70,70))
fig.tight_layout()