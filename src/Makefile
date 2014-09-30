
debug := 0
include ../Makefile.inc

all: clean PulsarSearch SearchImage

PulsarSearch: PulsarSearch.cpp
	$(MPICXX) -o $(PROJ_BASE)/bin/PulsarSearch PulsarSearch.cpp readConfiguration.cpp $(INCLUDES) $(LIBS) $(CFLAGS) $(LDFLAGS)

PulsarSearchSequential: PulsarSearchSequential.cpp
	$(MPICXX) -o $(PROJ_BASE)/bin/PulsarSearchSequential PulsarSearchSequential.cpp $(CPU_INCLUDES) $(LIBS) $(CFLAGS) $(LDFLAGS)

SearchImage: searchImage.cpp
	$(INTELCXX) -o $(PROJ_BASE)/bin/SearchImage searchImage.cpp $(IMG_INCLUDES) $(CFLAGS) $(IMG_LDFLAGS)

clean:
	rm -f $(PROJ_BASE)/bin/PulsarSearch
	rm -f $(PROJ_BASE)/bin/SearchImage
