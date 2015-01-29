
# https://github.com/isazi/utils
UTILS := $(HOME)/src/utils
# https://github.com/isazi/OpenCL
OPENCL := $(HOME)/src/OpenCL
# https://github.com/isazi/AstroData
ASTRODATA := $(HOME)/src/AstroData
# https://github.com/isazi/Dedispersion
DEDISPERSION := $(HOME)/src/Dedispersion
# https://github.com/isazi/Transpose
TRANSPOSE := $(HOME)/src/Transpose
# https://github.com/isazi/Folding
FOLDING := $(HOME)/src/Folding
# https://github.com/isazi/SNR
SNR := $(HOME)/src/SNR
# Boost
BOOST := $(HOME)/src/boost
# CImg
CIMG := $(HOME)/src/CImg

INCLUDES := -I"include" -I"$(ASTRODATA)/include" -I"$(UTILS)/include" -I"$(DEDISPERSION)/include" -I"$(TRANSPOSE)/include" -I"$(FOLDING)/include" -I"$(SNR)/include"
CL_INCLUDES := $(INCLUDES) -I"$(OPENCL)/include"
CL_LIBS := -L"$(OPENCL_LIB)"
BOOST_LIBS := -L"$(BOOST)/lib" 

CFLAGS := -std=c++11 -Wall
ifneq ($(debug), 1)
	CFLAGS += -O3 -g0
else
	CFLAGS += -O0 -g3
endif

LDFLAGS := -lm
CL_LDFLAGS := $(LDFLAGS) -lOpenCL
BOOST_LDFLAGS := -lboost_mpi -lboost_serialization 
HDF5_LDFLAGS := -lhdf5_cpp
CIMG_LDFLAGS := $(LDFLAGS) -lX11 -fopenmp

CC := g++
MPI := mpicxx

# Dependencies
KERNELS := $(DEDISPERSION)/bin/Dedispersion.o $(TRANSPOSE)/bin/Transpose.o $(FOLDING)/bin/Folding.o $(SNR)/bin/SNR.o
DEPS := $(ASTRODATA)/bin/Observation.o $(ASTRODATA)/bin/ColorMap.o $(UTILS)/bin/ArgumentList.o $(UTILS)/bin/Timer.o $(UTILS)/bin/utils.o
CL_DEPS := $(DEPS) $(OPENCL)/bin/Exceptions.o $(OPENCL)/bin/InitializeOpenCL.o $(OPENCL)/bin/Kernel.o 


all: bin/readConfiguration.o bin/PulsarSearch bin/searchImage bin/searchPercentile

bin/readConfiguration.o: $(DEPS) $(KERNELS) include/readConfiguration.hpp src/readConfiguration.cpp
	$(CC) -o bin/readConfiguration.o -c src/readConfiguration.cpp $(INCLUDES) $(CFLAGS)

bin/PulsarSearch: $(DEPS) $(KERNELS) $(ASTRODATA)/include/ReadData.hpp $(ASTRODATA)/include/Generator.hpp bin/readConfiguration.o
	$(MPI) -o bin/PulsarSearch src/PulsarSearch.cpp bin/readConfiguration.o $(KERNELS) $(CL_DEPS) $(CL_INCLUDES) $(CL_LIBS) $(BOOST_LIBS) $(HDF5_LDFLAGS) $(CL_LDFLAGS) $(CFLAGS)

bin/searchImage: $(DEPS) src/searchImage.cpp
	$(CC) -o bin/searchImage src/searchImage.cpp $(DEPS) $(INCLUDES) -I"$(CIMG)/include" $(CIMG_LDFLAGS) $(CFLAGS)

bin/searchPercentile: $(DEPS) src/searchPercentile.cpp
	$(CC) -o bin/searchPercentile src/searchPercentile.cpp $(DEPS) $(INCLUDES) $(LDFLAGS) $(CFLAGS)

clean:
	-@rm bin/*

