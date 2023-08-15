// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#define ARROW_VERSION_MAJOR 10
#define ARROW_VERSION_MINOR 0
#define ARROW_VERSION_PATCH 1
#define ARROW_VERSION ((ARROW_VERSION_MAJOR * 1000) + ARROW_VERSION_MINOR) * 1000 + ARROW_VERSION_PATCH

#define ARROW_VERSION_STRING "10.0.1"

#define ARROW_SO_VERSION "1000"
#define ARROW_FULL_SO_VERSION "1000.1.0"

#define ARROW_CXX_COMPILER_ID "GNU"
#define ARROW_CXX_COMPILER_VERSION "11.3.0"
#define ARROW_CXX_COMPILER_FLAGS "-fvisibility-inlines-hidden -fmessage-length=0 -march=nocona -mtune=haswell -ftree-vectorize -fPIC -fstack-protector-strong -fno-plt -O2 -ffunction-sections -pipe -isystem /home/user/SourceCode/GNPS_ML_Processing_Workflow/bin/env/include -fdebug-prefix-map=/home/conda/feedstock_root/build_artifacts/apache-arrow_1673819175290/work=/usr/local/src/conda/libarrow-10.0.1 -fdebug-prefix-map=/home/user/SourceCode/GNPS_ML_Processing_Workflow/bin/env=/usr/local/src/conda-prefix -fdiagnostics-color=always -fuse-ld=gold -O2 -DNDEBUG -ftree-vectorize"

#define ARROW_BUILD_TYPE "RELEASE"

#define ARROW_GIT_ID "86b1fc2b89e0c8d3b3f02f67b56337be743dad8f"
#define ARROW_GIT_DESCRIPTION ""

#define ARROW_PACKAGE_KIND ""

#define ARROW_COMPUTE
#define ARROW_CSV
/* #undef ARROW_CUDA */
#define ARROW_DATASET
#define ARROW_FILESYSTEM
#define ARROW_FLIGHT
#define ARROW_FLIGHT_SQL
#define ARROW_IPC
#define ARROW_JEMALLOC
#define ARROW_JEMALLOC_VENDORED
#define ARROW_JSON
#define ARROW_ORC
#define ARROW_PARQUET
#define ARROW_SUBSTRAIT

#define ARROW_GCS
#define ARROW_S3
#define ARROW_USE_NATIVE_INT128
/* #undef ARROW_WITH_MUSL */
/* #undef ARROW_WITH_OPENTELEMETRY */
/* #undef ARROW_WITH_UCX */

#define GRPCPP_PP_INCLUDE
