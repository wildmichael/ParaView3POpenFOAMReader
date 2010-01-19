#!/usr/bin/env python
#/*=========================================================================
#
#    Copyright (c) 2010 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
#    All rights reserved.
#    See License_v1.2.txt for details.
#
#    This software is distributed WITHOUT ANY WARRANTY; without even
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#    PURPOSE.  See the above copyright notice for more information.
#
#=========================================================================*/

# A simple VTK-Python scripting example for use with vtkPOpenFOAMReader.
#
# Set current directory to the case directory and animate through all
# available time steps by coloring the mesh with the specified field.
#
# Usage: vtkPythonExample1.py fieldName caseType
#   where
#     fieldName: name of the field used for color mapping
#     caseType: 0 = decomposed case, 1 = reconstructed case
# Ex.:
#   cd run/tutorials/multiphase/interFoam/laminar/damBreakFine
#   vtkPythonExample1.py alpha1 0

import sys, vtk

if (len(sys.argv) < 3):
  print 'Usage: %s fieldName caseType' % sys.argv[0]
  print '  fieldName: name of the field used for color mapping'
  print '  caseType: 0 = decomposed case, 1 = reconstructed case'
  quit()

fieldName = sys.argv[1] # Field name used for color mapping
caseType = int(sys.argv[2]) # Case type

# Create and setup the reader
try:
  reader = vtk.vtkNewPOpenFOAMReader() # This package
except AttributeError:
  reader = vtk.vtkPOpenFOAMReader() # VTK 5.5-cvs now has the reader
# SetFileName(): Basically, only the path part of the given file name
# string matters to the reader, i.e.:
#  * The specified file name does _not_ have to really exist. Thus the
#    file name can even be a whitespace (as shown in this example).
#  * If no prepending path is present in the given string (as shown in
#    this example), the current directory is regarded as the case
#    directory.
# However there is an exception: when the file is named
# 'controlDict.foam', the reader assumes the file is in the system/
# subdirectory of the intended case directory.
reader.SetFileName(' ')
reader.SetCaseType(caseType) # 0/1 = Decomposed/Reconstructed case
reader.SetCreateCellToPoint(0) # 0/1 = Turn off/on point filtered arrays

# Create the geometry filter that converts the multiblock dataset to a
# polydata that can be handled by vtkPolyDataMapper
gf = vtk.vtkCompositeDataGeometryFilter()
gf.SetInputConnection(0, reader.GetOutputPort(0))

# Create and setup the mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(gf.GetOutputPort())
mapper.CreateDefaultLookupTable() # Use default color lookup table
mapper.SetScalarModeToUseCellFieldData() # Use cell data for color mapping
mapper.SelectColorArray(fieldName) # Specify array name used for color mapping

# Create the actor, renderer and render window
actor = vtk.vtkActor()
actor.SetMapper(mapper)
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
ren.AddActor(actor)
ren.SetBackground(0.1, 0.2, 0.4) # Background color
renWin.SetSize(1024, 768) # Window size

# Scan time steps and get the number of time steps
reader.UpdateInformation() # Scan time steps and create metadata
exe = reader.GetExecutive()
outInfo = exe.GetOutputInformation(0)
timeStepsKey = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()
nTimeSteps = outInfo.Length(timeStepsKey) # Get the number of time steps
print 'Number of time steps = ', nTimeSteps

# Step through the time steps
for stepI in range(nTimeSteps):
  timeValue = outInfo.Get(timeStepsKey, stepI) # Get the stepI-th time value
  print 'Time step = ', stepI, ', Time value = ', timeValue
  exe.SetUpdateTimeStep(0, timeValue) # Set the reader time to timeValue
  reader.Modified() # Reading is triggered upon mapper.Update()
  # Force reading now. If automatic ranging is not necessary, this can
  # be omitted as well.
  mapper.Update()
  # Automatically set the color range. Note that the array is not
  # available until mapper.Update() is executed.
  mapper.SetScalarRange(
      gf.GetOutput().GetCellData().GetArray(fieldName).GetRange())
  renWin.Render() # Render
