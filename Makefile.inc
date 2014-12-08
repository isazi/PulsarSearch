
PROJ_BASE	:= $(HOME)/src/PulsarSearch

# https://github.com/isazi/utils
ISA_UTILS_INCLUDE := $(HOME)/src/utils/include
# https://github.com/isazi/OpenCL
ISA_OPENCL_INCLUDE := $(HOME)/src/OpenCL/include
# https://github.com/isazi/AstroData
ISA_ASTRO_INCLUDE := $(HOME)/src/AstroData/include
# https://github.com/isazi/Dedispersion
DEDISPERSION_INCLUDE := $(HOME)/src/Dedispersion/include
# https://github.com/isazi/Transpose
TRANSPOSE_INCLUDE := $(HOME)/src/Transpose/include
# https://github.com/isazi/Folding
FOLDING_INCLUDE := $(HOME)/src/Folding/include
# https://github.com/isazi/SNR
SNR_INCLUDE := $(HOME)/src/SNR/include
# Boost
BOOST_INCLUDE := $(HOME)/src/boost/include
BOOST_LIB := $(HOME)/src/boost/lib
# CImg
CIMG_INCLUDE := $(HOME)/src/CImg

INCLUDES := -I"$(PROJ_BASE)/include" -I"$(OPENCL_INCLUDE)" -I"$(ISA_UTILS_INCLUDE)" -I"$(ISA_OPENCL_INCLUDE)" -I"$(ISA_ASTRO_INCLUDE)" -I"$(DEDISPERSION_INCLUDE)" -I"$(TRANSPOSE_INCLUDE)" -I"$(FOLDING_INCLUDE)" -I"$(SNR_INCLUDE)" -I"$(BOOST_INCLUDE)"
CPU_INCLUDES := -I"$(PROJ_BASE)/include" -I"$(ISA_UTILS_INCLUDE)" -I"$(ISA_ASTRO_INCLUDE)" -I"$(DEDISPERSION_INCLUDE)" -I"$(TRANSPOSE_INCLUDE)" -I"$(FOLDING_INCLUDE)" -I"$(SNR_INCLUDE)" -I"$(BOOST_INCLUDE)"
IMG_INCLUDES := -I"$(PROJ_BASE)/include" -I"$(ISA_UTILS_INCLUDE)" -I"$(ISA_ASTRO_INCLUDE)" -I"$(CIMG_INCLUDE)"
LIBS := -L"$(OPENCL_LIB)" -L"$(BOOST_LIB)"

ifneq ($(debug), 1)
	CFLAGS := -std=c++11 -O3 -g0 -Wall -fopenmp
else
	CFLAGS := -std=c++11 -O0 -g3 -Wall -fopenmp
endif

LDFLAGS := -lm -lOpenCL -lhdf5 -lhdf5_cpp -lboost_mpi -lboost_serialization
CPU_LDFLAGS := -lm
IMG_LDFLAGS := -lm -lX11

CC := g++
MPICXX := mpicxx

