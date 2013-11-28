
DEFINES := -DUSE_PTHREADS

CFLAGS := -I../.. $(CFLAGS) $(DEFINES)

#CXXFLAGS := -I../.. -fPIC -Wall -g $(CXXFLAGS) $(DEFINES)
CXXFLAGS := -I../.. -fPIC -Wall -O3 -ffast-math -ftree-vectorize $(CXXFLAGS) $(DEFINES)

LDFLAGS := $(LDFLAGS)

#VG	:= valgrind

LIBS 		:= ../../qm-dsp/libqm-dsp.a -lpthread
PROGRAM_LIBS	:= -lsndfile

SOURCES	:= CQKernel.cpp ConstantQ.cpp test.cpp

OBJECTS := $(SOURCES:.cpp=.o)
OBJECTS := $(OBJECTS:.c=.o)

PROGRAM	:= test

all: $(PROGRAM)

test:	$(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS) $(PROGRAM_LIBS)

clean: 
	rm -f *.o 

depend:
	makedepend -Y $(SOURCES)

# DO NOT DELETE

CQKernel.o: CQKernel.h
ConstantQ.o: ConstantQ.h CQKernel.h
test.o: ConstantQ.h CQKernel.h
