# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-src"
  "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-build"
  "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-subbuild/googletest-populate-prefix"
  "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-subbuild/googletest-populate-prefix/tmp"
  "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
  "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-subbuild/googletest-populate-prefix/src"
  "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/jamessquires/Music Dev/test_sc_plugin_cookie/hardclipadaa/_deps/googletest-subbuild/googletest-populate-prefix/src/googletest-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
