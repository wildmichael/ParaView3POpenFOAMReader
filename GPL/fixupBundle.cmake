# run as follows under the plugin build directory:
# cmake -DParaView_DIR=<full path to ParaView build directory> \
#   -DParaView_App=<full path to ParaView application bundle (ParaView.app)> \
#   -P fixupBundle.cmake

# Make sure this CMake has BundleUtilities.cmake:
#
if(NOT EXISTS "${CMAKE_ROOT}/Modules/BundleUtilities.cmake")
  message(FATAL_ERROR "error: BundleUtilities.cmake not found. Use CMake 2.6.2 or later.")
endif(NOT EXISTS "${CMAKE_ROOT}/Modules/BundleUtilities.cmake")


# Avoid following symlinks encountered during FILE GLOB_RECURSE calls:
#
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0009 NEW)
endif(COMMAND CMAKE_POLICY)


# gp_item_default_embedded_path_override item default_embedded_path_var
#
# Return the path that others should refer to the item by when the item
# is embedded inside a bundle.
#
# This is a project-specific override of BundleUtilities.cmake's
# gp_item_default_embedded_path
#
# Added *Panel.dylib pattern from the original one that comes with ParaView.
#
function(gp_item_default_embedded_path_override item default_embedded_path_var)

  # By default, embed items as set by gp_item_default_embedded_path:
  #
  set(path "${${default_embedded_path_var}}")

  # But for ParaView...
  #
  # ...embed *.dylib in the Libraries folder:
  #
  if(item MATCHES "\\.dylib$")
    set(path "@executable_path/../Libraries")
  endif(item MATCHES "\\.dylib$")

  # ...embed *Plugin.dylib in the Plugins folder:
  #
  if(item MATCHES "Plugin\\.dylib$")
    set(path "@executable_path/../Plugins")
  endif(item MATCHES "Plugin\\.dylib$")

  # ...embed *Panel.dylib in the Plugins folder:
  #
  if(item MATCHES "Panel\\.dylib$")
    set(path "@executable_path/../Plugins")
  endif(item MATCHES "Panel\\.dylib$")

  # ...embed *.so in the Plugins folder:
  #
  if(item MATCHES "\\.so$")
    set(path "@executable_path/../Plugins")
  endif(item MATCHES "\\.so$")

  # ...embed *Python.so in the Libraries folder:
  #
  if(item MATCHES "Python\\.so$")
    set(path "@executable_path/../Libraries")
  endif(item MATCHES "Python\\.so$")

  set(${default_embedded_path_var} "${path}" PARENT_SCOPE)
endfunction(gp_item_default_embedded_path_override)

include(BundleUtilities)

# script mode does not accept load_cache()
if(ParaView_App)
  set(app ${ParaView_App})
else(ParaView_App)
  message(FATAL_ERROR "ParaView_App is not set.")
endif(ParaView_App)

if(ParaView_DIR)
  set(dirs "${ParaView_DIR}/bin")
else(ParaView_DIR)
  message(FATAL_ERROR "ParaView_DIR is not set.")
endif(ParaView_DIR)

# plugins that have to be fixed up
set(plugins "libPOpenFOAMReaderSMPlugin.dylib;libPOpenFOAMPanel.dylib")

foreach(plugin ${plugins})
  set(libs "${libs}$ENV{PWD}/${plugin};")
endforeach(plugin)

get_bundle_and_executable("${app}" bundle executable valid)
if(valid)
  # resolve dependencies
  get_bundle_keys("${app}" "${libs}" "${dirs}" keys)

  get_filename_component(exepath "${executable}" PATH)

  foreach(plugin ${plugins})
    message(STATUS "fixing ${plugin}")

    # get variable prefix for the corresponding plugin
    get_item_key("${plugin}" key)

    # copy plugin into the bundle Plugin directory
    copy_resolved_item_into_bundle(
      "${${key}_RESOLVED_ITEM}"
      "${${key}_RESOLVED_EMBEDDED_ITEM}"
      )

    # fix the dependencies of the copied plugin inside the bundle
    fixup_bundle_item("${${key}_RESOLVED_EMBEDDED_ITEM}" "${exepath}" "${dirs}")
  endforeach(plugin)
else(valid)
  message(STATUS "error: ${app} is not a valid application bundle.")
endif(valid)
