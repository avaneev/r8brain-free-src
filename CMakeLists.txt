cmake_minimum_required(VERSION 3.15)
project(r8brain-free-src C CXX)

option(R8B_DLL "Build r8brain-free-src as shared library" ON)
option(R8B_IPP "Use IPP FFT instead of Ooura FFT" OFF)
option(R8B_PFFFT "Use PFFFT instead of Ooura FFT" OFF)
option(R8B_PFFFT_DOUBLE "Use PFFFT instead of Ooura FFT" OFF)
option(R8B_STATIC_CRT "Link MSVC CRT statically" OFF)
option(R8B_FASTTIMING "Enable fast timing in r8brain-free-src" OFF)

if (R8B_PFFFT_DOUBLE AND R8B_PFFFT)
	message(FATAL_ERROR "If you want to use PFFFT with double, you should first of all disable R8B_PFFFT.")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND WIN32)
	message(FATAL_ERROR "ClangCL is not supported yet.")
endif()

file(GLOB R8B_SRC
	"CDSPBlockConvolver.h"
	"CDSPFIRFilter.h"
	"CDSPFracInterpolator.h"
	"CDSPHBDownsampler.h"
	"CDSPHBUpsampler.h"
	"CDSPHBUpsampler.inc"
	"CDSPProcessor.h"
	"CDSPRealFFT.h"
	"CDSPResampler.h"
	"CDSPSincFilterGen.h"

	"r8bbase.cpp"
	"r8bbase.h"
	
	"r8bconf.h"
	"r8butil.h"
)

file(GLOB R8B_DLL_SRC
	"DLL/r8bsrc.cpp"
	"DLL/r8bsrc.h"
)

if (R8B_PFFFT)
	file(GLOB FFT_SRC
		"pffft.cpp"
		"pffft.h"
	)
else()
	if (R8B_PFFFT_DOUBLE)
		file(GLOB FFT_SRC
			"pffft_double/pffft_double.c"
			"pffft_double/pffft_double.h"
		)
	else()
		file(GLOB FFT_SRC "fft4g.h")
	endif()
endif()

if (R8B_DLL)
	add_library(r8brain SHARED ${R8B_SRC} ${FFT_SRC} ${R8B_DLL_SRC})
else()
	add_library(r8brain STATIC ${R8B_SRC} ${FFT_SRC})
endif()
	
target_include_directories(r8brain PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) # public include for usage
target_include_directories(r8brain PRIVATE "DLL")
target_include_directories(r8brain PRIVATE "pffft_double")

if (R8B_FASTTIMING)
	target_compile_definitions(r8brain PUBLIC "R8B_FASTTIMING=1")
endif()

if (R8B_IPP)
	target_compile_definitions(r8brain PUBLIC "R8B_IPP=1")
endif()

if (R8B_PFFFT)
	target_compile_definitions(r8brain PUBLIC "R8B_PFFFT=1")
endif()

if (R8B_PFFFT_DOUBLE)
	target_compile_definitions(r8brain PUBLIC "R8B_PFFFT_DOUBLE=1")
endif()

if (R8B_DLL)
	target_compile_definitions(r8brain PUBLIC R8B_DLL)
	target_compile_definitions(r8brain PRIVATE R8B_DLL_EXPORT)
endif()
	
set_target_properties(r8brain
	PROPERTIES
		LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/"
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/"
)

if (R8B_STATIC_CRT)
	set_property(TARGET r8brain PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()
