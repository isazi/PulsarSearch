
debug		:= 1
include		../Makefile.inc


all: clean PulsarSearchTest


PulsarSearchTest: main.cpp
	$(MPICXX) -o $(PROJ_BASE)/bin/PulsarSearchTest main.cpp ../src/readConfiguration.cpp $(UTILS_SRC) $(INCLUDES) $(LIBS) $(CFLAGS) $(LDFLAGS)

readConfigurationTest: readConfiguration.cpp
	$(CC) -o $(PROJ_BASE)/bin/readConfigurationTest readConfiguration.cpp ../src/readConfiguration.cpp $(UTILS_SRC) $(CPU_INCLUDES) $(CPU_CFLAGS) $(CPU_LDFLAGS)

clean:
	rm -f $(PROJ_BASE)/bin/PulsarSearchTest
	rm -f $(PROJ_BASE)/bin/readConfigurationTest
