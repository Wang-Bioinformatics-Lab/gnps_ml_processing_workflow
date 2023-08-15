# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
#
# This config sets the following variables in your project::
#
#   Gandiva_FOUND - true if Gandiva found on the system
#
# This config sets the following targets in your project::
#
#   Gandiva::gandiva_shared - for linked as shared library if shared library is built
#   Gandiva::gandiva_static - for linked as static library if static library is built


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was GandivaConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set(ARROW_LLVM_VERSIONS "15.0;14.0;13.0;12.0;11.1;11.0;10;9;8;7")

include(CMakeFindDependencyMacro)
find_dependency(Arrow)
if(DEFINED CMAKE_MODULE_PATH)
  set(GANDIVA_CMAKE_MODULE_PATH_OLD ${CMAKE_MODULE_PATH})
else()
  unset(GANDIVA_CMAKE_MODULE_PATH_OLD)
endif()
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
find_dependency(LLVMAlt)
if(DEFINED GANDIVA_CMAKE_MODULE_PATH_OLD)
  set(CMAKE_MODULE_PATH ${GANDIVA_CMAKE_MODULE_PATH_OLD})
  unset(GANDIVA_CMAKE_MODULE_PATH_OLD)
else()
  unset(CMAKE_MODULE_PATH)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/GandivaTargets.cmake")

arrow_keep_backward_compatibility(Gandiva gandiva)

check_required_components(Gandiva)

arrow_show_details(Gandiva GANDIVA)
