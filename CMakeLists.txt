cmake_minimum_required(VERSION 3.16)
project(moresampler C)

# Use the compiler settings from ciglet's Makefile.
set(MORESAMPLER_ARCH "64" CACHE STRING "Target architecture (32 or 64)")
set_property(CACHE MORESAMPLER_ARCH PROPERTY STRINGS 32 64)

if(NOT MORESAMPLER_ARCH STREQUAL "32" AND
   NOT MORESAMPLER_ARCH STREQUAL "64")
    message(FATAL_ERROR "MORESAMPLER_ARCH must be 32 or 64")
endif()

set(FP_TYPE "float" CACHE STRING "Floating-point type used by the C sources")
set_property(CACHE FP_TYPE PROPERTY STRINGS float double)
if(NOT FP_TYPE)
    message(FATAL_ERROR "FP_TYPE must not be empty")
endif()

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Match original default build type when user did not pass one.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Detect whether the selected compiler already targets 32-bit.
set(MORESAMPLER_IS_I686 FALSE)

execute_process(
    COMMAND "${CMAKE_C_COMPILER}" -dumpmachine
    OUTPUT_VARIABLE MORESAMPLER_DUMPMACHINE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(MORESAMPLER_DUMPMACHINE MATCHES "^i[3-6]86" OR
   CMAKE_SIZEOF_VOID_P EQUAL 4 OR
   CMAKE_C_COMPILER_ARCHITECTURE_ID MATCHES "X86")
    set(MORESAMPLER_IS_I686 TRUE)
endif()

include(CheckCCompilerFlag)

if(MORESAMPLER_ARCH STREQUAL "32")
    if(NOT MORESAMPLER_IS_I686)
        if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            message(FATAL_ERROR "A 32-bit GNU or Clang compiler is required")
        endif()
        check_c_compiler_flag("-m32" MORESAMPLER_FLAG_M32)
        if(NOT MORESAMPLER_FLAG_M32)
            message(FATAL_ERROR
                "The selected compiler is not 32-bit and does not support -m32")
        endif()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m32")
    endif()
else()
    if(MORESAMPLER_IS_I686)
        if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            message(FATAL_ERROR "A 64-bit GNU or Clang compiler is required")
        endif()
        check_c_compiler_flag("-m64" MORESAMPLER_FLAG_M64)
        if(NOT MORESAMPLER_FLAG_M64)
            message(FATAL_ERROR
                "The selected compiler is 32-bit and does not support -m64")
        endif()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m64")
    endif()
endif()

# Common flags matching ciglet/Makefile.
set(MORESAMPLER_C_FLAGS
    -I${CMAKE_CURRENT_SOURCE_DIR}/external
    -DFP_TYPE=${FP_TYPE}
    -std=c99
    -Wall
    -fPIC
    -fopenmp
)
# ciglet Makefile uses _POSIX_C_SOURCE=2. Keep it on 32-bit; GCC 16
# treats the resulting undeclared pipe/fork in ciglet.h as errors.
if(MORESAMPLER_ARCH STREQUAL "32")
    list(APPEND MORESAMPLER_C_FLAGS -D_POSIX_C_SOURCE=2)
endif()

# ciglet Makefile: Release=-Ofast, Debug=-Og -g.
set(CMAKE_C_FLAGS_RELEASE "-Ofast"
    CACHE STRING "Release C flags" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-Ofast -g"
    CACHE STRING "RelWithDebInfo C flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-Og -g"
    CACHE STRING "Debug C flags" FORCE)
set(CMAKE_C_FLAGS_MINSIZEREL "-Os"
    CACHE STRING "MinSizeRel C flags" FORCE)

# OpenMP is required for every build.
find_package(OpenMP REQUIRED)

# ---- source files ----------------------------------------------------------
set(MORESAMPLER_SOURCES
    main.c
    logging.c
    utils.c
    config.c
    flag.c
    # lmix module
    lmix/lmix.c
    lmix/lmix_xchg.c
    lmix/lmix_analysis.c
    lmix/lmix_synth.c
    lmix/lmix_utils.c
    # xchg module
    xchg.c
    # resampler
    resampler.c
    resampler_args.c
    # wavtool
    wavtool.c
    # external: ciglet
    external/ciglet/ciglet.c
    external/ciglet/external/fftsg_h.c
    external/ciglet/external/fast_median.c
    external/ciglet/external/wavfile.c
    # external: libaiff 5.0 (AIFF/AIFC reader for audioread)
    external/libaiff/libaiff.c
    external/libaiff/aifx.c
    external/libaiff/extended.c
    external/libaiff/float32.c
    external/libaiff/g711.c
    external/libaiff/iff.c
    external/libaiff/lpcm.c
    external/libaiff/pascal.c
    # external: libgvps (Viterbi for nebula_est / phase unwrap)
    external/libgvps/gvps_obsrv.c
    external/libgvps/gvps_full.c
    external/libgvps/gvps_sampled.c
    external/libgvps/gvps_variable.c
    # external: libllsm
    external/libllsm/llsm-layer0.c
    external/libllsm/llsm-layer1.c
    external/libllsm/envelope.c
    external/libllsm/math-funcs.c
    # external: mrq
    external/mrq/mrq.c
    # nebula (model, IF/SNR features, F0 estimation/refinement)
    nebula/nebula_model.c
    nebula/nebula_features.c
    nebula/nebula_f0.c
    nebula/nebula_wave.c
    # llsmsynth
    llsmsynth.c
    # otogen (Lua 5.3 host + HSMM alignment)
    otogen/otogen.c
    otogen/otogen_hash.c
    otogen/otogen_dict.c
    otogen/otogen_analyze.c
    # external: FastLZ (otogen G2P blobs)
    external/FastLZ/fastlz.c
    # external: cJSON
    external/cJSON/cJSON.c
    # external: liblrhsmm
    external/liblrhsmm/common.c
    external/liblrhsmm/math-funcs.c
    external/liblrhsmm/mempool.c
    external/liblrhsmm/model.c
    external/liblrhsmm/data.c
    external/liblrhsmm/generate.c
    external/liblrhsmm/inference.c
    external/liblrhsmm/estimate.c
)

# Lua 5.3.3 runtime core (no lua.c or luac.c).
set(LUA_CORE_SOURCES
    external/lua/lapi.c
    external/lua/lauxlib.c
    external/lua/lbaselib.c
    external/lua/lbitlib.c
    external/lua/lcode.c
    external/lua/lcorolib.c
    external/lua/lctype.c
    external/lua/ldblib.c
    external/lua/ldebug.c
    external/lua/ldo.c
    external/lua/ldump.c
    external/lua/lfunc.c
    external/lua/lgc.c
    external/lua/linit.c
    external/lua/liolib.c
    external/lua/llex.c
    external/lua/lmathlib.c
    external/lua/lmem.c
    external/lua/loadlib.c
    external/lua/lobject.c
    external/lua/lopcodes.c
    external/lua/loslib.c
    external/lua/lparser.c
    external/lua/lstate.c
    external/lua/lstring.c
    external/lua/lstrlib.c
    external/lua/ltable.c
    external/lua/ltablib.c
    external/lua/ltm.c
    external/lua/lundump.c
    external/lua/lutf8lib.c
    external/lua/lvm.c
    external/lua/lzio.c
)
list(APPEND MORESAMPLER_SOURCES ${LUA_CORE_SOURCES})

add_executable(moresampler ${MORESAMPLER_SOURCES})

target_include_directories(moresampler PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/lmix
    ${CMAKE_CURRENT_SOURCE_DIR}/nebula
    ${CMAKE_CURRENT_SOURCE_DIR}/external/ciglet
    ${CMAKE_CURRENT_SOURCE_DIR}/external/ciglet/external
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libgvps
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libllsm
    ${CMAKE_CURRENT_SOURCE_DIR}/external/mrq
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libaiff
    ${CMAKE_CURRENT_SOURCE_DIR}/external/cJSON
    ${CMAKE_CURRENT_SOURCE_DIR}/external/liblrhsmm
    ${CMAKE_CURRENT_SOURCE_DIR}/external/lua
    ${CMAKE_CURRENT_SOURCE_DIR}/external/FastLZ
    ${CMAKE_CURRENT_SOURCE_DIR}/otogen
)

# Keep cJSON symbols private to this executable.
target_compile_definitions(moresampler PRIVATE
    CJSON_HIDE_SYMBOLS
    MORESAMPLER_ARCH=${MORESAMPLER_ARCH}
)

target_compile_options(moresampler PRIVATE ${MORESAMPLER_C_FLAGS})

target_link_libraries(moresampler PRIVATE OpenMP::OpenMP_C)

# GCC's OpenMP imported target may expose the compile flag but omit it from
# the executable link command when using the MinGW toolchain bundled with
# CLion.  Pass the detected OpenMP flags explicitly so libgomp is linked.
if(OpenMP_C_FLAGS)
    target_link_options(moresampler PRIVATE ${OpenMP_C_FLAGS})
endif()

if(WIN32)
    target_link_libraries(moresampler PRIVATE shlwapi)
    # Keep the standard Windows security flags for the executable.
    if(MINGW)
        target_link_options(moresampler PRIVATE
            "LINKER:--dynamicbase"
            "LINKER:--nxcompat"
        )
        # Strip release binaries.
        if(CMAKE_BUILD_TYPE MATCHES "^(Release|MinSizeRel)$")
            target_link_options(moresampler PRIVATE "LINKER:-s")
        endif()
    endif()
endif()

# Status summary
message(STATUS "==== moresampler build config ====")
message(STATUS "  C compiler : ${CMAKE_C_COMPILER}")
message(STATUS "  target     : ${CMAKE_C_COMPILER_TARGET}")
message(STATUS "  build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  arch       : ${MORESAMPLER_ARCH}-bit")
message(STATUS "  OpenMP     : ${OpenMP_C_FOUND}")
message(STATUS "  C flags    : ${MORESAMPLER_C_FLAGS}")
message(STATUS "==================================")
