# -------------------------------------------------------------
# -------------------------------------------------------------
#        	  OpenPhaseSnapshot
#					 Version: 17.11.2014
#		
# 	2D plotting tool for OpenPhase
# Call function with the following parameters:
# - Name of vtk file created by OpenPhase (including .vtk-ending)
# - (Optional) [x,y,z]  Cutting plane normal (x,y,z)
# - (Optional) -s[how] 	   Direct display of figure
# - (Optional) -o[verview] 	   Contour plots (overview for tensors)
# - (Optional) -singlecontour Contour plot
# - (Optional) 	--sc#n 	Contour plot for tensorial value n (1-6)
# - (Optional) -dx#n 	   Give distance between data points
#
# 	OpenPhase, ICAMS, Ruhr-Universitaet Bochum
#		philipp.s.engels@rub.de
# -------------------------------------------------------------
# -------------------------------------------------------------

#import matplotlib
import sys, os
import csv # csv module  
import re # String manipulation
import vtk # VTK support
import numpy as np # Numpy
import matplotlib.pyplot as mpl # Matplotlib main 
import matplotlib.cm as cm # Matplot contour map
import matplotlib.ticker as ticker  # Module for tick positioning
import matplotlib.colors as mplcolors # Module to convert numbers to RGB colors

from matplotlib import rc # LaTex Interpreter
from vtk.util import numpy_support as npvtk

# -------------------------------------------------------------
# -------------------------------------------------------------

class Snapshot(object):

	def __init__(self):	
		self.timestep = 0
		self.cpnormal = ''	
		self.filename = ''
		self.plottype = 'singlecontour'
		self.varname = ''
		self.dx = 1.0e-3
		self.stressunitcorrection = 1.e-6  #m to mm (MPa)
		self.showonscreen = False
		self.istensorvariable = 999

	def getCall(self, call):
		calllength = np.size(call)
		#if calllength < 2: # 1 + 2 necessary + 2 optional
		#	sys.exit("Error: Parameter missing. Quit!")
		for index, item in enumerate(call):
			#if item in Varnamearray:
			#	varname = str(call[index])
			if item in ['-overview', '-o']:
				self.plottype = 'overview'
			if item in ['-l', '-lines']:
				self.plottype = 'lines'
			if (item in ['-show', '-s']):
				self.showonscreen = True#str(call[index])
			if item in ['-sc#1','-sc#2','-sc#3','-sc#4','-sc#5','-sc#6']:
				self.istensorvariable =  int(re.split('[#]',str(call[index]))[1])-1
				self.plottype = 'singlecontour'

			#if item in ['-dx#']:
			#	print (item)
			#	dx =  int(re.split('[#]',str(call[index]))[1])-1
			#	print(dx)
			if item in ['x','y','z']:
				self.cpnormal = str(call[index])
			if '.vts' in item:
				self.filename = str(call[index])
				self.varname =  re.split('[_.]',self.filename)[0]
				self.timestep = str(call[index])
				self.timestep = int(re.split('[_.]',self.timestep)[1])				
	
		if os.path.exists(self.filename) == False:
			sys.exit("Error: File not found. Exit!")
		if self.varname == '':
			sys.exit('Error: No output variable given. Exit!')
		if self.cpnormal == '':
			self.cpnormal = 'z'
			print('Warning: No cutting plane given. Set normal to <0,0,1>')

#	return [filename, varname, options.cpnormal, dx,
#		showonscreen, stressunitcorrection, timestep, plottype, istensorvariable]


def setprintoptions(): # Define Plot Options
	mpl.rc('text', usetex=True)
	mpl.rc('font', family='serif')
	mpl.rc('font', size=16)
	numcontourline = 800
	figdpi = 150
	figformat = 'png' # eps, pdf -> bad results for contour plots
	return [figdpi, figformat, numcontourline]


def getcuttingplanenormal(): # return xi,yi,zi where zi is set to options.cpnormal
	if options.cpnormal == 'x':
		temp = [2,0,1]
	if options.cpnormal == 'y':
		temp = [0,2,1]
	if options.cpnormal == 'z':
		temp = [0,1,2]

	for i in temp:
		if temp[i] == 0:
			xi = i
		if temp[i] == 1:
			yi = i 
		if temp[i] == 2:
			zi = i
	
	return [xi, yi, zi]

def standardplotoptions(istensorvariable, subplotnames, plottype):
	subfig.set_xlabel(axesnames[xi])
	
	subfig.grid(True)

	if plottype == 'singlecontour':
			subfig.set_ylabel(unitname)
			subfig.set_ylabel(axesnames[yi])

	if plottype == 'lines':
			subfig.legend()
			subfig.set_ylabel(unitname)

def setsubplotoptions(istensorvariable, subplotnames):
	
	axesnames = ['x','y','z']
	
	subfig.set_xlabel(axesnames[xi])
	subfig.set_ylabel(axesnames[yi])

	if multiplesubplotsatonce == True:
		subfig.set_title(subplotnames[i] + unitname)
	elif istensorvariable != 999:
		subfig.set_title(subplotnames[istensorvariable] + unitname)
	else:
		subfig.set_title(subplotnames + unitname)

# -------------------------------------------------------------
#			Start
# -------------------------------------------------------------
# Read call parameters

print('                     ')
print('                     ')

options = Snapshot()
multiplesubplotsatonce = False
options.getCall(sys.argv)

#[filename, varname, options.cpnormal, dx, showonscreen, stressunitcorrection,
#	 timestep, plottype,istensorvariable] = initiate(sys.argv)

# -------------------------------------------------------------
# Open VTK-File

reader = vtk.vtkXMLStructuredGridReader()
reader.SetFileName(options.filename)
reader.UpdateInformation()
#reader.ReadAllScalarsOn()
reader.Update()
dim3 = reader.GetOutput().GetDimensions()
ngridpoints = dim3[0]*dim3[1]*dim3[2]

# -------------------------------------------------------------
# Transfer VTK-Data to Array 

if (options.varname == 'AccumulatedStresses'):
	StressNames = ['SigmaAcc_xx','SigmaAcc_yy','SigmaAcc_zz','SigmaAcc_yz','SigmaAcc_xz','SigmaAcc_xy']
	array = np.arange(6*ngridpoints, dtype='f').reshape(6, dim3[0],dim3[1],dim3[2])
	for i in range(6):
		if not (reader.GetOutput().GetPointData().GetArray(StressNames[i]) == None):
			arrayVTK = reader.GetOutput().GetPointData().GetArray(StressNames[i])
			temparray = npvtk.vtk_to_numpy(arrayVTK)
			array[i,:,:,:] = temparray.reshape(dim3[0],dim3[1],dim3[2])*options.stressunitcorrection
		else:
			array[i,:,:,:] =  np.zeros(ngridpoints, dtype='f').reshape(dim3[0],dim3[1],dim3[2])

if (options.varname == 'Stresses'):
	StressNames = ['Sigma_xx','Sigma_yy','Sigma_zz','Sigma_yz','Sigma_xz','Sigma_xy']
	array = np.arange(6*ngridpoints, dtype='f').reshape(6, dim3[0],dim3[1],dim3[2])
	for i in range(6):
		if not (reader.GetOutput().GetPointData().GetArray(StressNames[i]) == None):
			print(reader.GetOutput().GetPointData().GetArray(StressNames[i]))
			arrayVTK = reader.GetOutput().GetPointData().GetArray(StressNames[i])
			temparray = npvtk.vtk_to_numpy(arrayVTK)
			array[i,:,:,:] = temparray.reshape(dim3[0],dim3[1],dim3[2])*options.stressunitcorrection
		else:
			array[i,:,:,:] =  np.zeros(ngridpoints, dtype='f').reshape(dim3[0],dim3[1],dim3[2])

if options.varname == 'Strains':
	StrainNames = ['epsilon_xx','epsilon_yy','epsilon_zz','gamma_yz','gamma_xz','gamma_xy']
	array = np.arange(6*ngridpoints, dtype='f').reshape(6, dim3[0],dim3[1],dim3[2])
	for i in range(6):
		if not (reader.GetOutput().GetPointData().GetArray(StrainNames[i]) == None):
			arrayVTK = reader.GetOutput().GetPointData().GetArray(StrainNames[i])
			temparray = npvtk.vtk_to_numpy(arrayVTK)
			array[i,:,:,:] = temparray.reshape(dim3[0],dim3[1],dim3[2])
		else:
			array[i,:,:,:] =  np.zeros(ngridpoints, dtype='f').reshape(dim3[0],dim3[1],dim3[2])

if options.varname == 'PlasticStrains':
	StrainNames = ['PEpsilon_xx','PEpsilon_yy','PEpsilon_zz','PEpsilon_yz','PEpsilon_xz','PEpsilon_xy']
	array = np.arange(6*ngridpoints, dtype='f').reshape(6, dim3[0],dim3[1],dim3[2])
	for i in range(6):
		if not (reader.GetOutput().GetPointData().GetArray(StrainNames[i]) == None):
			arrayVTK = reader.GetOutput().GetPointData().GetArray(StrainNames[i])
			temparray = npvtk.vtk_to_numpy(arrayVTK)
			array[i,:,:,:] = temparray.reshape(dim3[0],dim3[1],dim3[2])
		else:
			array[i,:,:,:] =  np.zeros(ngridpoints, dtype='f').reshape(dim3[0],dim3[1],dim3[2])

if (options.varname == 'EffectivePlasticStrains'):
	StrainNames = ['PEpsilon_xx','PEpsilon_yy','PEpsilon_zz','PEpsilon_yz','PEpsilon_xz','PEpsilon_xy']
	array = np.arange(6*ngridpoints, dtype='f').reshape(6, dim3[0],dim3[1],dim3[2])
	for i in range(6):
		if not (reader.GetOutput().GetPointData().GetArray(StrainNames[i]) == None):
			arrayVTK = reader.GetOutput().GetPointData().GetArray(StrainNames[i])
			temparray = npvtk.vtk_to_numpy(arrayVTK)
			array[i,:,:,:] = temparray.reshape(dim3[0],dim3[1],dim3[2])
		else:
			array[i,:,:,:] =  np.zeros(ngridpoints, dtype='f').reshape(dim3[0],dim3[1],dim3[2])

if options.varname == 'PhaseField':
	PhaseFieldNames = ['PhaseFields','Interfaces']
	array = np.arange(2*ngridpoints).reshape(2, dim3[0],dim3[1],dim3[2]) #, dtype='f'
	
	for i in range(2):
		arrayVTK = reader.GetOutput().GetPointData().GetArray(PhaseFieldNames[i])
		temparray = npvtk.vtk_to_numpy(arrayVTK)
		array[i,:,:,:] = temparray.reshape(dim3[0],dim3[1],dim3[2])

'''
	shufflearray = np.random.randint(np.max(array), size=dim3[2]*dim3[1]*dim3[0])
	
	# <Experimental> Shuffle routine for phase field indices.   
	for i in range(dim3[2]):
		for j in range(dim3[1]): 	
			for k in range(dim3[0]): 
				array[0,i,j,k] = shufflearray[array[0,i,j,k]]	
'''

VarNames = ['Temperature', 'Concentration', 'DislocationDensity']
arraynames = ['T', 'c(tot)', 'DD_ave']	
for i in range(3):
	if options.varname == VarNames[i]:
		arrayVTK = reader.GetOutput().GetPointData().GetArray(arraynames[i])
		temparray = npvtk.vtk_to_numpy(arrayVTK)
		array = np.arange(ngridpoints, dtype='f').reshape(dim3[0],dim3[1],dim3[2])	
		#array = temparray.reshape(dim3[2],dim3[1],dim3[0])
		#array = array.transpose(2,1,0)
VarNames = 'MechResiduum'
arraynames = 'residuum'	
if options.varname == VarNames:
	arrayVTK = reader.GetOutput().GetPointData().GetArray(arraynames)
	temparray = npvtk.vtk_to_numpy(arrayVTK)
	array = np.arange(ngridpoints, dtype='f').reshape(dim3[0],dim3[1],dim3[2])	
	#array = temparray.reshape(dim3[2],dim3[1],dim3[0])
# -------------------------------------------------------------
# Get Cutting Plane Normal

[xi, yi, zi] = plotcoordinates = getcuttingplanenormal()

if dim3[zi] == 1: 
	cpindex = 0 # pure 2D
else:
	cpindex = dim3[zi]/2
# -------------------------------------------------------------
# Plot

axesnames = ['x','y','z']
subplotindex = [231,232,233,234,235,236]

x = np.arange(0,dim3[xi]) * options.dx
y = np.arange(0,dim3[yi]) * options.dx
X,Y = np.meshgrid(x,y)

if options.varname in ['AccumulatedStresses','Stresses', 'Strains', 'PlasticStrains', 'EffectivePlasticStrains', 'PhaseField']: # Multiple variables in VTK
	if options.cpnormal == 'x':
		z = array[:,cpindex,:,:] 
	if options.cpnormal == 'y':
		z = array[:,:,cpindex,:]
	if options.cpnormal == 'z':
		z = array[:,:,:,cpindex]
if options.varname in ['Temperature', 'Concentration', 'DislocationDensity', 'MechResiduum']:  # Single variable in VTK
	
	if options.cpnormal == 'x':
		z = array[cpindex,:,:] 
	if options.cpnormal == 'y':
		z = array[:,cpindex,:]
	if options.cpnormal == 'z':
		z = array[:,:,cpindex]

[figdpi, figformat, numcontourline] = setprintoptions()

myfigure = mpl.figure(figsize=(10.0, 10.0))
#myfigure.suptitle(varname + ' - Time Step: ' + str(timestep), fontsize=24, fontweight='bold')

mainaxes = mpl.axes()
mainaxes.axes.set_visible(False)

if (options.varname in ['Strains', 'PlasticStrains', 'EffectivePlasticStrains', 'Stresses', 'AccumulatedStresses']):
	getlength = np.shape(z)[2]
	#print(np.shape(z))
	if (getlength % 2 == 0):
	    print ('Warning: Non-odd cube length. Average coordinates: ', getlength/2, 'and ', getlength/2-1)
	    resreturn = (z[:,:,getlength/2]+z[:,:,getlength/2-1])/2.0
	else:
	    resreturn = z[:,:,getlength/2]
		
	if (options.varname == 'Stresses' or options.varname == 'AccumulatedStresses'):
		unitname = ' (MPa)'
		formatstring = '%0.2e'
		subplotnames = ['$\sigma_{xx}$','$\sigma_{yy}$','$\sigma_{zz}$',
				'$\sigma_{yz}$','$\sigma_{xz}$','$\sigma_{xy}$']
	if options.varname == 'Strains':
		unitname = ' (-)'
		formatstring = '%0.2e'
		subplotnames = ['$\epsilon_{xx}$','$\epsilon_{yy}$','$\epsilon_{zz}$',
				'$\epsilon_{yz}$','$\epsilon_{xz}$','$\epsilon_{xy}$']
	if options.varname == 'PlasticStrains' or options.varname == 'EffectivePlasticStrains':
		unitname = ' (-)'
		formatstring = '%0.2e'
		subplotnames = ['$\epsilon_{xx}^p$','$\epsilon_{yy}^p$','$\epsilon_{zz}^p$',
				'$\epsilon_{yz}^p$','$\epsilon_{xz}^p$','$\epsilon_{xy}^p$']

	if (options.plottype == 'lines'): # Lineplot
		subfig = myfigure.add_subplot(111,aspect='auto')
		
		for i in np.arange(6):
			subfig.plot(x,resreturn[i,:],'-o', label=subplotnames[i])
		
			#mpl.ylim(np.min(z[i,:,getlength/2]),np.max(z[i,:,getlength/2]))
		standardplotoptions(options.istensorvariable, subplotnames, options.plottype)
		
		# -------------------------------------------------------------
		# CSV Output of line data
		with open(options.varname + '.csv', 'wb') as csvfile:
			mycsv = csv.writer(csvfile, delimiter=';')
			#for i in range(6):
			#	mycsv.writerow(z[i,:,getlength/2]) #[subplotnames[i]] + 
			mycsv.writerow(x) #[subplotnames[i]] + 
			mycsv.writerows(resreturn[:,:]) #[subplotnames[i]] + 

		# -------------------------------------------------------------
	if (options.plottype == 'singlecontour'): #Contourplot of single tensor variable
		if (options.istensorvariable == 999):
			print('No tensor variable given! Enter: ')
			options.istensorvariable = input() - 1
		
		subfig = myfigure.add_subplot(111,aspect='equal')
		#getlength = np.shape(z)[2]
		subplot = subfig.contourf(X,Y,np.transpose(z[options.istensorvariable,:,:]),
				numcontourline,interpolation='nearest', cmap=cm.coolwarm)
		cbar = myfigure.colorbar(subplot, orientation='vertical',shrink=0.8, format=formatstring)
		setsubplotoptions(options.istensorvariable, subplotnames)
		#subfig.set_aspect('equal')
		#subfig.set_title(subplotnames[int(istensorvariable)] + unitname)

	if (options.plottype == 'overview'): # Generate overview contour plots
		multiplesubplotsatonce = True
		myfigure.set_size_inches(20,15)
		for i in range(6):
			subfig = myfigure.add_subplot(subplotindex[i],aspect='equal')
			subplot = subfig.contourf(X,Y,np.transpose(z[i,:,:]),
					numcontourline,interpolation='nearest', cmap=cm.coolwarm)
			cbar = myfigure.colorbar(subplot, orientation='vertical',shrink=0.8, format=formatstring)
		setsubplotoptions(options.istensorvariable, subplotnames)

if options.varname == 'PhaseField':
	
	#multiplesubplotsatonce = False
	unitname = ''
	subfig = myfigure.add_subplot(111, aspect='auto')
	subplotnames = r'$\phi$ '
	getlength = np.shape(z)[2]

	formatstring = '%0.2e'

	if (options.plottype == 'lines'):
		subfig.plot(x,np.transpose(z[0,:,getlength/2]),'-o', label=subplotnames) #, label=subplotnames[i] , format=formatstring
		standardplotoptions(options.istensorvariable, subplotnames, options.plottype)
		# -------------------------------------------------------------
		# CSV Output of line data
		with open(options.varname + '.csv', 'wb') as csvfile:
			mycsv = csv.writer(csvfile, delimiter=';')
			#for i in range(6):
			#	mycsv.writerow(z[i,:,getlength/2]) #[subplotnames[i]] + 
			mycsv.writerow(x) #[subplotnames[i]] + 
			mycsv.writerows(z[:,:,getlength/2]) #[subplotnames[i]] + 

		# -------------------------------------------------------------

	if (options.plottype == 'singlecontour'): # Generate overview contour plots
		subplot = subfig.contourf(X,Y,np.transpose(z[0,:,:]),numcontourline,
				interpolation='nearest', cmap=cm.coolwarm, format=formatstring)
		#subfig.contour(X,Y,np.transpose(z[0,:,:]), levels=np.arange(99), #,levels=np.remainder(npz[0,:,:]
		#	linewidths=1,colors='blue') #  !!noch keine Idee!!
		cbar = mpl.colorbar(subplot, orientation='vertical',shrink=0.8) #,	
		#ticks=np.arange(np.min(z), np.max(z),4)
		cbar.ax.set_ylabel(unitname)
		#subplotnames = '-'
		standardplotoptions(options.istensorvariable, subplotnames, options.plottype)

if options.varname in ['Concentration', 'Temperature', 'DislocationDensity', 'MechResiduum']:
	
	#multiplesubplotsatonce = False
	
	if options.varname == 'Concentration':
		subplotnames = 'c'
		unitname = r' $(\text{mol}/m^3)$'
		formatstring = '%3.0f'
	if options.varname == 'Temperature':
		subplotnames = 'T'
		unitname = ' (K)'
		formatstring = '%3.0e'
	if options.varname == 'DislocationDensity':
		subplotnames = r'$\rho$'
		unitname = r' $(1/m^2)$'
		formatstring = '%3.0e'
	if options.varname == 'MechResiduum':
		subplotnames = r''
		unitname = r''
		formatstring = '%3.0e'

	subfig = myfigure.add_subplot(111,aspect='auto')
	getlength = np.shape(z)[1]
	
	if (options.plottype == 'lines'):
		subfig.plot(x,z[:,getlength/2],'-o',format=formatstring) #, label=subplotnames[i]
		standardplotoptions(options.istensorvariable, subplotnames)
		# -------------------------------------------------------------
		# CSV Output of line data
		with open(options.varname + '.csv', 'wb') as csvfile:
			mycsv = csv.writer(csvfile, delimiter=';')
			#for i in range(6):
			#	mycsv.writerow(z[i,:,getlength/2]) #[subplotnames[i]] + 
			mycsv.writerow(x) #[subplotnames[i]] + 
			mycsv.writerows(z[:,:,getlength/2]) #[subplotnames[i]] + 

		# -------------------------------------------------------------

	if (options.plottype == 'singlecontour' and options.varname not in ['MechResiduum']): # Generate contour plots

		if options.varname in ['DislocationDensity']: # LogPlot
				#z2 = np.ma.masked_where(z<= 0, z+1)
				#lev_exp = np.arange(np.floor(np.log10(z2.min())+2),
						#       np.ceil(np.log10(z2.max())+1))			
				#levs = np.power(10, lev_exp)
				levels = np.linspace(1.0e10, z.max(), numcontourline)
				subplot = subfig.contourf(X,Y,np.transpose(z), levels=levels, interpolation='nearest', cmap=cm.coolwarm) 
		else: 
				subplot = subfig.contourf(X,Y,np.transpose(z),numcontourline,
						ticks=np.arange(np.min(z), np.max(z),4),interpolation='nearest', cmap=cm.coolwarm)
		#bar = mpl.colorbar(subplot, orientation='vertical',shrink=0.8, format=formatstring) 
		#cbar.ax.set_ylabel(unitname)
		setsubplotoptions(options.istensorvariable, subplotnames)

	if options.varname in ['MechResiduum']:
		subplot = subfig.hist(np.log10(z.flatten()), bins= 200, color='b', cumulative = False, log=False, histtype='stepfilled')
		#setsubplotoptions(options.istensorvariable, subplotnames)
		
# Print Figure
print(options.filename[:-4])
myfigure.savefig(options.filename[:-4] + '.' + figformat,dpi=figdpi, bbox_inches='tight')

# Show Figure
if options.showonscreen == True:
	mpl.show()

# -------------------------------------------------------------
#			 End
# -------------------------------------------------------------
