#!/usr/bin/env python
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
import os
import csv
import sys
import scipy
import numpy as numpy

# -------------------------------------------------------------------
import matplotlib as MPL
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as PyPlot
from matplotlib.font_manager import FontProperties
import scipy.interpolate as interp

# ------------------------------------------------------------------
# Read csv data
# ------------------------------------------------------------------
resultsOP = []              
#resultsOP = numpy.empty(101)
with open('LineStresses.dat', 'r') as csvfile:
    csvreader = csv.reader(csvfile , delimiter = ',')
    next(csvreader)

    for ii in csvreader:
    	resultsOP.append(ii)
#----------------------------------------------------------------------
#coordinatesOPnp = numpy.arange(cubesizeOP)*dx
#print resultsOP
resultsOPnp = numpy.array(resultsOP)
resultsOPnp = resultsOPnp.astype(float)
coordinatesOPnp = resultsOPnp[:,0]
resultsOPnp = resultsOPnp[:,1:5]
cubesizeOP = numpy.shape(resultsOPnp)[1]

#print(coordinatesOPnp)
fig1 = PyPlot.figure()

# Settings
PyPlot.xticks(fontsize=14)
PyPlot.yticks(fontsize=14)
#PyPlot.xlim(-1.5,1.5)
#PyPlot.ylim(-500,500)
# Labels Definiton

# Set Grid
#PyPlot.grid(True, color='0.75')
  
PyPlot.xlabel(r'$x/L_x$', fontsize=18)
PyPlot.ylabel(r'Stress $\sigma$ (Pa)', fontsize=18)
colors = ['r','g','b','y']
markers = ['o','x','o','x']
lstyle =['solid', 'solid', 'dashed', 'dashed']

for index in [0, 1, 2, 3]:
    PyPlot.plot(coordinatesOPnp, resultsOPnp[:,index], '--', linewidth=2, color = colors[index],
     linestyle = lstyle[index], marker = markers[index], markevery=5)

#xformatter = FuncFormatter(my_formatter)
#xformatter = PyPlot.FormatStrFormatter('%.1e')
#yformatter = PyPlot.FormatStrFormatter('%1d')

ax = PyPlot.gca()

# Hide these grid behind plot objects
ax.set_axisbelow(True)

# Legend Definition
PyPlot.legend((r'$\sigma_{r,OpenPhase}$', r'$\sigma_{t,OpenPhase}$',r'$\sigma_{r,analytic}$', r'$\sigma_{t,analytic}$'),
			prop = FontProperties(size=16), loc=0)

# Save Figure
PyPlot.savefig('result.pdf', format='pdf',dpi=300, bbox_inches='tight')
PyPlot.show()

