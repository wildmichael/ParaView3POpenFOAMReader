#/*=========================================================================
#
#    Copyright (c) 2008-2009 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
#    All rights reserved.
#    See License_v1.2.txt for details.
#
#    This software is distributed WITHOUT ANY WARRANTY; without even
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#    PURPOSE.  See the above copyright notice for more information.
#
#=========================================================================*/

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

INCLUDE(${QT_USE_FILE})

#SET(Project_CMAKE_DIR "${Project_SOURCE_DIR}/CMake")
#INCLUDE(${Project_CMAKE_DIR}/FindRegex.cmake)
INCLUDE("CMake/FindRegex.cmake")

# Server side plugin

ADD_PARAVIEW_PLUGIN(
  POpenFOAMReaderPlugin "1.0" 
  SERVER_MANAGER_XML POpenFOAMReaderSM.xml
  SERVER_MANAGER_SOURCES vtkNewOpenFOAMReader.cxx vtkPOpenFOAMReader.cxx
  )

