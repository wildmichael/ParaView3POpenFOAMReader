## Copyright (c) 2008 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>
##
##     This software is distributed WITHOUT ANY WARRANTY; without even
##     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##     PURPOSE.

## A simple example of Python scription with the vtkOpenFOAMReader.
## For general information on Python scripting, see
## http://paraview.org/Wiki/images/f/f9/Servermanager2.pdf
##
## This script is a part of the vtkOpenFOAMReader distribution.

## The path to the icoFoam/cavity case. Modify as appropriate.
cavityPath = 'D:/ohshima/OpenFOAM/ohshima-1.4/run/tutorials/icoFoam/cavity/system/controlDict'

## Importing servermanager is not necessary if running from Tools->Python Shell 
#from paraview import servermanager

## Connect to a new built-in server
#connection = servermanager.Connect()

## Set the file name to be read and create the reader object. The file
## name doesn't have to end with the .foam suffix since we can
## explicitly specify the OpenFOAM reader.
reader = servermanager.sources.OpenFOAMReader(FileName = cavityPath)

## Force reader to get array, time and mesh information from the case directory
reader.UpdatePipelineInformation()

## We can get some information at this stage
#reader.TimestepValues
#reader.CellArrayInfo
#reader.PointArrayInfo
#reader.PatchArrayInfo

## Force reading data
#reader.UpdatePipeline()

## Get the first render view
view = servermanager.GetRenderView()
## To create a new view, use CreateRenderView()
#view = servermanager.CreateRenderView()
## Create a new representation object
rep = servermanager.CreateRepresentation(reader, view)

## Get timestep values in the case
tsteps = reader.TimestepValues
## Create time annotation filter using the reader as input
annTime = servermanager.filters.TimeToTextConvertor(Input=reader)
## Create a text representation from the time annotation
timeRep = servermanager.rendering.TextSourceRepresentation(Input=annTime)
## Append the text representation to the view
view.Representations.append(timeRep)
## Set the timestep value to the third timestep value
view.ViewTime = tsteps[2]

## Set the representation to Surface
## 0 = Points, 1 = Wireframe, 2 = Surface, 3 = Outline
rep.Representation = 2

## To color the geometry, we first get the lookup table object 
lt = servermanager.rendering.PVLookupTable()
## Set the object to the representation
rep.LookupTable = lt
## Set the coloring attribute type
## 0 = Point data, 1 = Cell data
rep.ColorAttributeType = 0
## Color by CellToPoint[p]
rep.ColorArrayName = 'CellToPoint[p]'
## The lookup table is used to map scalar values of the array to
## colors. The values in the lookup table have 4 numbers: 1 scalar
## value and 3 color values (R, G, B). This particular table has 2
## values, -5 = (0, 0, 1) (blue) and 5 = (1, 0, 0) (red). The colors in
## between are interpolated.
lt.RGBPoints = [-5, 0, 0, 1, 5, 1, 0, 0]
# Use the HSV color space to interpolate
lt.ColorSpace = 1 # HSV

## Reset the viewing camera
view.ResetCamera()
## Get the active camera
#cam = view.GetActiveCamera()
## Set the elevation
#cam.Elevation(45)

## Render the view
view.StillRender()

## Set the timestep value to the last timestep
view.ViewTime = tsteps[5]
## Render again
view.StillRender()
