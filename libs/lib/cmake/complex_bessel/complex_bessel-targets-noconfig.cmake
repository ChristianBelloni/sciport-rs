#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "complex_bessel::complex_bessel" for configuration ""
set_property(TARGET complex_bessel::complex_bessel APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(complex_bessel::complex_bessel PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libcomplex_bessel.dylib"
  IMPORTED_SONAME_NOCONFIG "@rpath/libcomplex_bessel.dylib"
  )

list(APPEND _cmake_import_check_targets complex_bessel::complex_bessel )
list(APPEND _cmake_import_check_files_for_complex_bessel::complex_bessel "${_IMPORT_PREFIX}/lib/libcomplex_bessel.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
