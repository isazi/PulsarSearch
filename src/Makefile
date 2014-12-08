
debug := 0
include ../Makefile.inc

all: clean PulsarSearch SearchImage SearchPercentile

PulsarSearch: PulsarSearch.cpp
	$(MPICXX) -o $(PROJ_BASE)/bin/PulsarSearch PulsarSearch.cpp readConfiguration.cpp $(INCLUDES) $(LIBS) $(CFLAGS) $(LDFLAGS)

SearchImage: searchImage.cpp
	$(CC) -o $(PROJ_BASE)/bin/SearchImage searchImage.cpp $(IMG_INCLUDES) $(CFLAGS) $(IMG_LDFLAGS)

SearchPercentile: searchPercentile.cpp
	$(CC) -o $(PROJ_BASE)/bin/SearchPercentile searchPercentile.cpp $(CPU_INCLUDES) $(CFLAGS) $(CPU_LDFLAGS)

clean:
	rm -f $(PROJ_BASE)/bin/PulsarSearch
	rm -f $(PROJ_BASE)/bin/SearchImage
	rm -f $(PROJ_BASE)/bin/SearchPercentile
