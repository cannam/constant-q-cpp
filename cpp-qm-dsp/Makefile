
DEFINES := -DUSE_PTHREADS

CFLAGS := -I../.. $(CFLAGS) $(DEFINES)

CXXFLAGS := -I../.. -Wall -g $(CXXFLAGS) $(DEFINES)
#CXXFLAGS := -I../.. -Wall -O3 -ffast-math -ftree-vectorize $(CXXFLAGS) $(DEFINES)

LDFLAGS := $(LDFLAGS)

#VG	:= valgrind

LIBS 		:= ../../qm-dsp/libqm-dsp.a -lpthread
PROGRAM_LIBS	:= -lsndfile

SOURCES	:= CQKernel.cpp test.cpp

OBJECTS := $(SOURCES:.cpp=.o)
OBJECTS := $(OBJECTS:.c=.o)

PROGRAM	:= test

all: $(PROGRAM)

test:	$(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS) $(PROGRAM_LIBS)

clean: 
	rm -f *.o 

