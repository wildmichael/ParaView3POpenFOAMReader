#/*=========================================================================
#
#    Copyright (c) 2009 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
#    All rights reserved.
#
#    This software is distributed WITHOUT ANY WARRANTY; without even
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#    PURPOSE.
#
#=========================================================================*/

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# Find a ParaView build first in order to get its configuration

FIND_PACKAGE(ParaView REQUIRED)

# Set default installation directory

IF(PARAVIEW_LIBRARY_DIRS)
  IF(CMAKE_INSTALL_PREFIX STREQUAL "/usr/local"
      OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files/Project")
    IF(APPLE)
      SET(CMAKE_INSTALL_PREFIX "${PARAVIEW_LIBRARY_DIRS}/Plugins")
    ELSEIF(DEFINED CMAKE_BUILD_TYPE)
      SET(CMAKE_INSTALL_PREFIX "${PARAVIEW_LIBRARY_DIRS}/plugins")
    ELSE()
      SET(CMAKE_INSTALL_PREFIX "${PARAVIEW_LIBRARY_DIRS}/Release/plugins")
    ENDIF()
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE PATH
      "Install path prefix, prepended onto install directories." FORCE)
  ENDIF()
ENDIF()

# Set default build type to Release if Makefile-based build system is used

IF(DEFINED CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF()

ADD_SUBDIRECTORY(BSD)
ADD_SUBDIRECTORY(GPL)
