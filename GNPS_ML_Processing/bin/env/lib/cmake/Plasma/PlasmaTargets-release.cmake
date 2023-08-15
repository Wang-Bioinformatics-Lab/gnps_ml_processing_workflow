#----------------------------------------------------------------
# Generated CMake target import file for configuration "RELEASE".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Plasma::plasma_shared" for configuration "RELEASE"
set_property(TARGET Plasma::plasma_shared APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Plasma::plasma_shared PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libplasma.so.1000.1.0"
  IMPORTED_SONAME_RELEASE "libplasma.so.1000"
  )

list(APPEND _cmake_import_check_targets Plasma::plasma_shared )
list(APPEND _cmake_import_check_files_for_Plasma::plasma_shared "${_IMPORT_PREFIX}/lib/libplasma.so.1000.1.0" )

# Import target "Plasma::plasma-store-server" for configuration "RELEASE"
set_property(TARGET Plasma::plasma-store-server APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Plasma::plasma-store-server PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/plasma-store-server"
  )

list(APPEND _cmake_import_check_targets Plasma::plasma-store-server )
list(APPEND _cmake_import_check_files_for_Plasma::plasma-store-server "${_IMPORT_PREFIX}/bin/plasma-store-server" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
