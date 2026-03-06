#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "cp2k::cp2k" for configuration "Release"
set_property(TARGET cp2k::cp2k APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cp2k::cp2k PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libcp2k.so.2026.1"
  IMPORTED_SONAME_RELEASE "libcp2k.so.2026.1"
  )

list(APPEND _cmake_import_check_targets cp2k::cp2k )
list(APPEND _cmake_import_check_files_for_cp2k::cp2k "${_IMPORT_PREFIX}/lib64/libcp2k.so.2026.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
