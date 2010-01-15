#/*=========================================================================
#
#    This file is part of pqPOpenFOAMPanel.
#
#    vtkPOpenFOAMReader is free software: you can redistribute it
#    and/or modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    vtkPOpenFOAMReader is distributed in the hope that it will be
#    useful, but WITHOUT ANY WARRANTY; without even the implied
#    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#    See the GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with vtkPOpenFOAMReader.  If not, see
#    <http://www.gnu.org/licenses/>.
#
#=========================================================================*/
# Copyright (c) 2008-2010 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
# All rights reserved.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Determine if in-source-tree plugin build

IF(DEFINED INSOURCE_PLUGIN)
  # VTK_DIR is defined in ${PARAVIEW_USE_FILE}
  IF(NOT INSOURCE_PLUGIN AND NOT VTK_DIR)
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

IF(PARAVIEW_BUILD_QT_GUI)

  # Construct a version number string

  IF(NOT DEFINED PARAVIEW_VERSION)
    SET(PARAVIEW_VERSION "${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}")
  ENDIF()

  # Add include directories

  IF(INSOURCE_PLUGIN)
    INCLUDE_DIRECTORIES(
      ${PARAVIEW_INCLUDE_DIRS}
      ${PARAVIEW_GUI_INCLUDE_DIRS}
      ${VTK_INCLUDE_DIR}
    )
    INCLUDE(${QT_USE_FILE})
  ENDIF()

  # Additional UI elements for the POpenFOAMReader panel

  QT4_WRAP_CPP(
    MOC_SRCS pqPOpenFOAMPanel.h
    )

  ADD_PARAVIEW_OBJECT_PANEL(
    IFACES IFACE_SRCS
    CLASS_NAME pqPOpenFOAMPanel
    XML_NAME OpenFOAMReader
    XML_GROUP sources
    )

  IF(PARAVIEW_VERSION VERSION_GREATER "3.6")
    ADD_PARAVIEW_PLUGIN(
      POpenFOAMPanelPlugin "1.0"
      REQUIRED_ON_CLIENT
      GUI_INTERFACES ${IFACES}
      GUI_SOURCES pqPOpenFOAMPanel.cxx
        ${MOC_SRCS} ${RCS_SRCS} ${IFACE_SRCS}
      )
  ELSE()
    ADD_PARAVIEW_PLUGIN(
      POpenFOAMPanelPlugin "1.0"
      GUI_INTERFACES ${IFACES}
      SOURCES pqPOpenFOAMPanel.cxx
        ${MOC_SRCS} ${RCS_SRCS} ${IFACE_SRCS}
    )
  ENDIF()

  # Set installation destination

  IF(NOT INSOURCE_PLUGIN)
    # ParaView 3.6.2 has been modified after release so that plugins
    # are built without manifests if MSVC is used. With the
    # modification there seems to be no way to tell exactly whether
    # plugins are built as SHARED or MODULE so there seems to be no
    # way to tell whether RUNTIME or LIBRARY target should be
    # specified either.
    INSTALL(TARGETS POpenFOAMPanelPlugin DESTINATION "${CMAKE_INSTALL_PREFIX}")
  ENDIF()
ENDIF()
