
debug		:= 0
include		../Makefile.inc


all: clean PulsarSearch
 

PulsarSearch: main.cpp
	$(LINKER) -o $(PROJ_BASE)/bin/PulsarSearch main.cpp $(UTILS_SRC) $(INCLUDES) $(LIBS) $(CFLAGS) $(LDFLAGS) 

clean:
	rm -f $(PROJ_BASE)/bin/PulsarSearch
