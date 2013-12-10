
debug		:= 1
include		../Makefile.inc


all: clean PulsarSearchTest
 

PulsarSearchTest: main.cpp
	$(MPICXX) -o $(PROJ_BASE)/bin/PulsarSearchTest main.cpp readConfiguration.cpp $(UTILS_SRC) $(INCLUDES) $(LIBS) $(CFLAGS) $(LDFLAGS) 

clean:
	rm -f $(PROJ_BASE)/bin/PulsarSearchTest
