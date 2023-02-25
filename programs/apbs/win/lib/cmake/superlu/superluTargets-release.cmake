#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "superlu::superlu" for configuration "Release"
set_property(TARGET superlu::superlu APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(superlu::superlu PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/superlu.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS superlu::superlu )
list(APPEND _IMPORT_CHECK_FILES_FOR_superlu::superlu "${_IMPORT_PREFIX}/lib/superlu.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
