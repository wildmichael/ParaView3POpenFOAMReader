#/*=========================================================================
#
#    Copyright (c) 2009 Takuya OSHIMA <oshima@eng.niigata-u.ac.jp>.
#    All rights reserved.
#
#    This software is distributed WITHOUT ANY WARRANTY; without even
#    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#    PURPOSE.  See the above copyright notice for more information.
#
#=========================================================================*/

# Find Regex

MARK_AS_ADVANCED(REGEX_INCLUDE_DIR)

# prevent FIND_PATH from finding regex.h under
# /System/Library/Frameworks/Ruby.framework/Headers on OS X
SET(CMAKE_FIND_FRAMEWORK_OLD "${CMAKE_FIND_FRAMEWORK}")
SET(CMAKE_FIND_FRAMEWORK "NEVER")

FIND_PATH(REGEX_INCLUDE_DIR regex.h)

SET(CMAKE_FIND_FRAMEWORK "${CMAKE_FIND_FRAMEWORK_OLD}")

IF(REGEX_INCLUDE_DIR)
  INCLUDE_DIRECTORIES(${REGEX_INCLUDE_DIR})
  ADD_DEFINITIONS(-DVTK_FOAMFILE_HAVE_REGEX)
ENDIF(REGEX_INCLUDE_DIR)
