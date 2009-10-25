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
# Copyright (c) 2008-2009 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
# All rights reserved.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

INCLUDE(${QT_USE_FILE})

IF(PARAVIEW_BUILD_QT_GUI)

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

  ADD_PARAVIEW_PLUGIN(
    POpenFOAMPanelPlugin "1.0" 
    GUI_INTERFACES ${IFACES}
    SOURCES pqPOpenFOAMPanel.cxx
      ${MOC_SRCS} ${RCS_SRCS} ${IFACE_SRCS}
    )

ENDIF()
