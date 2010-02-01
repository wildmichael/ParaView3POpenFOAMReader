#!/usr/bin/env pvpython
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

# A simple ParaView-Python scripting example for use with vtkPOpenFOAMReader.
#
# Synopsis: Set current directory to the case directory and animate
#   through all available time steps by coloring the mesh with the
#   specified field.
#
# Usage: pvPythonExample1.py fieldName caseType
#   where
#     fieldName: name of the field used for color mapping
#     caseType: 0 = decomposed case, 1 = reconstructed case
# Ex.:
#   cd run/tutorials/multiphase/interFoam/laminar/damBreakFine
#   pvPythonExample1.py alpha1 0

from paraview.simple import *
import sys, time

if (len(sys.argv) < 3):
  print 'Usage: %s fieldName caseType' % sys.argv[0]
  print '  fieldName: name of the field used for color mapping'
  print '  caseType: 0 = decomposed case, 1 = reconstructed case'
  quit()

fieldName = sys.argv[1] # Field name used for color mapping
caseType = int(sys.argv[2]) # Case type

# Create and setup the reader

# Unlike vtkPythonExample1.py, a try-except clause is not necessary to
# access the fallback builtin reader that comes with ParaView 3.7-cvs,
# since vtkNewOpenFOAMReader (a plugin version) and the builtin
# version shares a common proxy name of OpenFOAMReader. If both the
# plugin and the builtin versions are present in a single ParaView
# installation, the plugin version is used by ParaView's design.

# FileName: Note that the file name should be specified at the reader
# instantiation. As To the file name, basically only the path part of
# the given file name string matters to the reader, i.e.:
#  * The specified file does _not_ have to really exist, nor does its
#    name have to end with the '.foam' extension. Thus the file name
#    can even be a whitespace (as shown in this example).
#  * If no prepending path is present in the given string (as shown in
#    this example), the current directory is regarded as the case
#    directory.
# However there is an exception: when the file is named
# 'controlDict.foam', the reader assumes the file is in the system/
# subdirectory of the intended case directory.

reader = OpenFOAMReader(FileName = ' ')
reader.CaseType = caseType # 0/1 = Decomposed/Reconstructed case
# 0/1 = Turn off/on point filtered arrays
reader.Createcelltopointfiltereddata = 0

# Mark the reader as visible
Show()

# Set color mapping
dp = GetDisplayProperties()
dp.LookupTable = MakeBlueToRedLT(0.0, 1.0) # Use blue-to-red color lookup table
dp.ColorAttributeType = 'CELL_DATA' # Use cell data for color mapping
dp.ColorArrayName = fieldName # Specify array name used for color mapping

# Get the active view and set the view size
view = GetActiveView()
view.Background = [0.1, 0.2, 0.4] # Background color
view.ViewSize = [1024, 768] # View size

# We are going a hard way hereafter by not using the following shortcut
# Render()
# AnimateReader(reader)

# Scan time steps and get the number of time steps
reader.UpdatePipelineInformation() # Scan time steps and create metadata
nTimeSteps = len(reader.TimestepValues) # Get the number of time steps
print 'Number of time steps = ', nTimeSteps

# Step through the time steps
for stepI in range(nTimeSteps):
  timeValue = reader.TimestepValues[stepI] # Get the stepI-th time value
  print 'Time step = ', stepI, ', Time value = ', timeValue
  view.ViewTime = timeValue # Set the view time to timeValue
  # Automatically set the color range
  ra = reader.CellData[fieldName].GetRange()
  dp.LookupTable = MakeBlueToRedLT(ra[0], ra[1])
  Render() # Render
