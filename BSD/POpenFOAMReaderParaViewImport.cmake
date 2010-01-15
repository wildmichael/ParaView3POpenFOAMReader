#/*=========================================================================
#
#    Copyright (c) 2009-2010 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
#    All rights reserved.
#
#    This software is distributed WITHOUT ANY WARRANTY; without even
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#    PURPOSE.  See the above copyright notice for more information.
#
#=========================================================================*/

INCLUDE(${POpenFOAMReader_SOURCE_DIR}/CMake/FindRegex.cmake)

SET(POpenFOAMReader_SRCS 
  ${POpenFOAMReader_SOURCE_DIR}/vtkNewOpenFOAMReader.cxx
  ${POpenFOAMReader_SOURCE_DIR}/vtkNewPOpenFOAMReader.cxx
  )

INCLUDE_DIRECTORIES(${POpenFOAMReader_SOURCE_DIR})
#INCLUDE_DIRECTORIES(${POpenFOAMReader_SOURCE_DIR}/..)

PARAVIEW_INCLUDE_WRAPPED_SOURCES("${POpenFOAMReader_SRCS}")

PARAVIEW_INCLUDE_SERVERMANAGER_RESOURCES(
  "${POpenFOAMReader_SOURCE_DIR}/POpenFOAMReaderSM.xml"
  )
