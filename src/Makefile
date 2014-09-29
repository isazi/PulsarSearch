
debug		:= 0
include		../Makefile.inc

all: clean PulsarSearch SearchImage

PulsarSearch: main.cpp
	$(MPICXX) -o $(PROJ_BASE)/bin/PulsarSearch main.cpp readConfiguration.cpp $(INCLUDES) $(LIBS) $(CFLAGS) $(LDFLAGS)

SearchImage: searchImage.cpp
	$(INTELCXX) -o $(PROJ_BASE)/bin/SearchImage searchImage.cpp $(IMG_INCLUDES) $(CFLAGS) $(IMG_LDFLAGS)

clean:
	rm -f $(PROJ_BASE)/bin/PulsarSearch
	rm -f $(PROJ_BASE)/bin/SearchImage
