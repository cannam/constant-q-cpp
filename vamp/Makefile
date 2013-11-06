
DEFINES := -DUSE_PTHREADS

CFLAGS := -I../.. $(CFLAGS) $(DEFINES)

#CXXFLAGS := -I../.. -Wall -fPIC -g $(CXXFLAGS) $(DEFINES)
CXXFLAGS := -I../.. -Wall -fPIC -O3 -ffast-math -ftree-vectorize $(CXXFLAGS) $(DEFINES)

LDFLAGS := $(LDFLAGS)

#VG	:= valgrind

LIBS 		:= ../../qm-dsp/libqm-dsp.a ../../vamp-plugin-sdk/libvamp-sdk.a -lpthread 

SOURCES	:= CQVamp.cpp libmain.cpp

OBJECTS := $(SOURCES:.cpp=.o)
OBJECTS := $(OBJECTS:.c=.o)

EXTRA_OBJECTS	:= ../cpp-qm-dsp/ConstantQ.o ../cpp-qm-dsp/CQKernel.o 

PLUGIN_EXT := so
PLUGIN = cqvamp.$(PLUGIN_EXT)

all: $(PLUGIN)

$(PLUGIN):	$(OBJECTS)
	$(CXX) -shared -o $@ $^ $(EXTRA_OBJECTS) $(LDFLAGS) $(LIBS)

clean: 
	rm -f *.o 

depend:
	makedepend -Y $(SOURCES)

# DO NOT DELETE

CQVamp.o: CQVamp.h ../cpp-qm-dsp/ConstantQ.h ../cpp-qm-dsp/CQKernel.h
libmain.o: CQVamp.h
