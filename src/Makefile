
debug		:= 0
include		../Makefile.inc


all: clean PulsarSearch SearchImage
 

PulsarSearch: main.cpp
	$(MPICXX) -o $(PROJ_BASE)/bin/PulsarSearch main.cpp readConfiguration.cpp $(UTILS_SRC) $(INCLUDES) $(LIBS) $(CFLAGS) $(LDFLAGS)

SearchImage: searchImage.cpp
	$(CC) -o $(PROJ_BASE)/bin/SearchImage searchImage.cpp $(UTILS_SRC) $(IMG_INCLUDES) $(CFLAGS) $(IMG_LDFLAGS)

clean:
	rm -f $(PROJ_BASE)/bin/PulsarSearch
	rm -f $(PROJ_BASE)/bin/SearchImage
