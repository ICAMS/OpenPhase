#!/usr/bin/env pvbatch
# -----------------------------------------------------------------------------
# Script: post-render.py
# Purpose: Render VTK simulation data to PNG images
# Creator: Raphael Schiedung
# -----------------------------------------------------------------------------
# -----------------------------
# User-configurable parameters
# -----------------------------
# Define image resolution here (width x height)
IMAGE_WIDTH  = 1920
IMAGE_HEIGHT = 1080
INPUT_DIR='VTK'
OUTPUT_DIR='Images'
OUTPUT_NAME ='SolidPhaseFraction.png'
MAX_VTK=1000000000
# -----------------------------

#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
import glob
import os
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Define directories relative to the script location
script_dir = os.path.dirname(os.path.abspath(__file__))
vtk_dir    = os.path.join(script_dir, INPUT_DIR)
image_dir  = os.path.join(script_dir, OUTPUT_DIR)

# Create the folder if it doesn't exist
os.makedirs(image_dir, exist_ok=True)

# File patterns and output file
pattern = 'PhaseField_*.vts'
png_dir = os.path.join(image_dir, OUTPUT_NAME) 

# Get the sorted list of files matching the pattern
full_pattern = os.path.join(vtk_dir, pattern)
vtk_files = sorted(glob.glob(full_pattern))
num_vtk_files = len(vtk_files)
if num_vtk_files == 0:
    raise RuntimeError("No VTK files found — please check your path or pattern.")
else:
    print(f"Found {num_vtk_files} VTK files.")

# Create the reader with the list
phaseField_reader = XMLStructuredGridReader(
    registrationName='PhaseField',
    FileName=vtk_files
)
phaseField_reader.PointArrayStatus = ['Interfaces', 'Flags', 'PhaseFields', 'PhaseFraction_0', 'PhaseFraction_1', 'Junctions', 'Variants', 'ParentGrain']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on phaseField_reader
phaseField_reader.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
phaseField_Display = Show(phaseField_reader, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
phaseField_Display.Representation = 'Surface'
phaseField_Display.ColorArrayName = [None, '']
phaseField_Display.SelectTCoordArray = 'None'
phaseField_Display.SelectNormalArray = 'None'
phaseField_Display.SelectTangentArray = 'None'
phaseField_Display.OSPRayScaleArray = 'Flags'
phaseField_Display.OSPRayScaleFunction = 'PiecewiseFunction'
phaseField_Display.SelectOrientationVectors = 'None'
phaseField_Display.ScaleFactor = 44.900000000000006
phaseField_Display.SelectScaleArray = 'Flags'
phaseField_Display.GlyphType = 'Arrow'
phaseField_Display.GlyphTableIndexArray = 'Flags'
phaseField_Display.GaussianRadius = 2.245
phaseField_Display.SetScaleArray = ['POINTS', 'Flags']
phaseField_Display.ScaleTransferFunction = 'PiecewiseFunction'
phaseField_Display.OpacityArray = ['POINTS', 'Flags']
phaseField_Display.OpacityTransferFunction = 'PiecewiseFunction'
phaseField_Display.DataAxesGrid = 'GridAxesRepresentation'
phaseField_Display.PolarAxes = 'PolarAxesRepresentation'
phaseField_Display.ScalarOpacityUnitDistance = 14.463821502707397
phaseField_Display.SelectInputVectors = [None, '']
phaseField_Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
phaseField_Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
phaseField_Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [10000.0, 34.5, 224.5]
renderView1.CameraFocalPoint = [0.0, 34.5, 224.5]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(phaseField_Display, ('POINTS', 'PhaseFraction_1'))

# rescale color and/or opacity maps used to include current data range
phaseField_Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
phaseField_Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'PhaseFraction_1'
phaseFraction_1LUT = GetColorTransferFunction('PhaseFraction_1')

# get color legend/bar for phaseFraction_1LUT in view renderView1
phaseFraction_1LUTColorBar = GetScalarBar(phaseFraction_1LUT, renderView1)

# Properties modified on phaseFraction_1LUTColorBar
phaseFraction_1LUTColorBar.Title = 'Solid Phase Fracton'
phaseFraction_1LUTColorBar.RangeLabelFormat = '%-#6.1f'

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

renderView1.ApplyIsometricView()

# reset view to fit data
renderView1.ResetCamera(False)
renderView1.ResetActiveCameraToPositiveZ()

# reset view to fit data
renderView1.ResetCamera(False)
renderView1.ResetActiveCameraToPositiveY()


# reset view to fit data
renderView1.ResetCamera(False)
renderView1.ResetActiveCameraToPositiveX()

# rescale color and/or opacity maps used to exactly fit the current data range
animationScene1.GoToLast()
#phaseField_Display.RescaleTransferFunctionToDataRange(False, True)

# change scalar bar placement
phaseFraction_1LUTColorBar.Orientation = 'Horizontal'
phaseFraction_1LUTColorBar.Position = [0.43931050314457115, 0.8980762564991334]
phaseFraction_1LUTColorBar.ScalarBarLength = 0.1026655326945149

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(IMAGE_WIDTH, IMAGE_HEIGHT)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1127.2113276303223, 34.5, 225.01272104325156]
renderView1.CameraFocalPoint = [0.0, 34.5, 225.01272104325156]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 295.840041956145

# save animation
SaveAnimation(png_dir, renderView1, ImageResolution=[IMAGE_WIDTH, IMAGE_HEIGHT],
    TransparentBackground=1,
    FrameWindow=[0, min(num_vtk_files-1,MAX_VTK)])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(IMAGE_WIDTH, IMAGE_HEIGHT)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1127.2113276303223, 34.5, 225.01272104325156]
renderView1.CameraFocalPoint = [0.0, 34.5, 225.01272104325156]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 295.840041956145

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
