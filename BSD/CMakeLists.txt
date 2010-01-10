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

# Determine if in-source-tree plugin build

IF(DEFINED INSOURCE_PLUGIN)
  IF(NOT INSOURCE_PLUGIN)
    INCLUDE(${PARAVIEW_USE_FILE})
  ENDIF()
ELSE()
  IF(ParaView_SOURCE_DIR)
    SET(INSOURCE_PLUGIN TRUE)
  ELSE()
    SET(INSOURCE_PLUGIN FALSE)

    FIND_PACKAGE(ParaView REQUIRED)
    INCLUDE(${PARAVIEW_USE_FILE})
  ENDIF()
ENDIF()

# Construct a version number string

IF(NOT DEFINED PARAVIEW_VERSION)
  SET(PARAVIEW_VERSION "${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}")
ENDIF()

# Add include directories

IF(INSOURCE_PLUGIN)
  INCLUDE_DIRECTORIES(${PARAVIEW_INCLUDE_DIRS} ${VTK_INCLUDE_DIR})
ENDIF()

# Find a regular expression library

#SET(Project_CMAKE_DIR "${Project_SOURCE_DIR}/CMake")
#INCLUDE(${Project_CMAKE_DIR}/FindRegex.cmake)
INCLUDE("CMake/FindRegex.cmake")

# Server side plugin

IF(PARAVIEW_VERSION VERSION_GREATER "3.6")
  ADD_PARAVIEW_PLUGIN(
    POpenFOAMReaderPlugin "1.0"
    REQUIRED_ON_SERVER
    SERVER_MANAGER_XML POpenFOAMReaderSM.xml
    SERVER_MANAGER_SOURCES vtkNewOpenFOAMReader.cxx vtkPOpenFOAMReader.cxx
    )
ELSE()
  ADD_PARAVIEW_PLUGIN(
    POpenFOAMReaderPlugin "1.0"
    SERVER_MANAGER_XML POpenFOAMReaderSM.xml
    SERVER_MANAGER_SOURCES vtkNewOpenFOAMReader.cxx vtkPOpenFOAMReader.cxx
    )
ENDIF()

# Set installation destination

IF(NOT INSOURCE_PLUGIN)
  # ParaView 3.6.2 has been modified after release so that plugins are
  # built without manifests if MSVC is used. With the modification
  # there seems to be no way to tell exactly whether plugins are built
  # as SHARED or MODULE so there seems to be no way to tell whether
  # RUNTIME or LIBRARY target should be specified either.
  INSTALL(TARGETS POpenFOAMReaderPlugin DESTINATION "${CMAKE_INSTALL_PREFIX}")
ENDIF()
