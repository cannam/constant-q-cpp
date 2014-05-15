
LIB_DIR	  := src
INC_DIR	  := cq
TEST_DIR  := test
VAMP_DIR  := vamp

VAMPSDK_DIR  ?= ../vamp-plugin-sdk

PLUGIN_EXT	?= .so

CXX	?= g++
CC	?= gcc

CFLAGS := $(CFLAGS) 
CXXFLAGS := -I. -I$(VAMPSDK_DIR) -I$(INC_DIR) -I$(LIB_DIR) $(CXXFLAGS)

LDFLAGS := $(LDFLAGS) 
PLUGIN_LDFLAGS := $(LDFLAGS) $(PLUGIN_LDFLAGS)
TEST_LDFLAGS := $(LDFLAGS) -lboost_unit_test_framework
PF_LDFLAGS := $(LDFLAGS) -lsndfile

PLUGIN	:= cqvamp$(PLUGIN_EXT)
TEST	:= $(TEST_DIR)/test
PF	:= $(TEST_DIR)/processfile

LIB_HEADERS := $(INC_DIR)/CQBase.h $(LIB_DIR)/CQKernel.h $(INC_DIR)/ConstantQ.h $(INC_DIR)/CQSpectrogram.h $(INC_DIR)/CQInverse.h
LIB_SOURCES := $(LIB_DIR)/CQKernel.cpp $(LIB_DIR)/ConstantQ.cpp $(LIB_DIR)/CQSpectrogram.cpp $(LIB_DIR)/CQInverse.cpp

VAMP_HEADERS := $(VAMP_DIR)/CQVamp.h $(VAMP_DIR)/CQChromaVamp.h $(VAMP_DIR)/Pitch.h
VAMP_SOURCES := $(VAMP_DIR)/CQVamp.cpp $(VAMP_DIR)/CQChromaVamp.cpp $(VAMP_DIR)/libmain.cpp $(VAMP_DIR)/Pitch.cpp

HEADERS	     := $(LIB_HEADERS) $(VAMP_HEADERS)
SOURCES	     := $(LIB_SOURCES) $(VAMP_SOURCES)
OBJECTS	     := $(SOURCES:.cpp=.o)

TEST_SOURCES := $(TEST_DIR)/test.cpp 
TEST_OBJECTS := $(TEST_SOURCES:.cpp=.o) $(OBJECTS)

PF_SOURCES := $(TEST_DIR)/processfile.cpp
PF_OBJECTS := $(PF_SOURCES:.cpp=.o) $(OBJECTS)

LIBS	:= $(VAMPSDK_DIR)/libvamp-sdk.a -lpthread

all: $(PLUGIN) $(TEST) $(PF)

$(PLUGIN):	$(OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(PLUGIN_LDFLAGS)

$(TEST):	$(TEST_OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(TEST_LDFLAGS)

$(PF):	$(PF_OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(PF_LDFLAGS)

clean:		
	rm -f $(OBJECTS) $(TEST_OBJECTS) $(PF_OBJECTS)

distclean:	clean
	rm -f $(PLUGIN)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(TEST_SOURCES) $(PF_SOURCES) $(HEADERS)

# DO NOT DELETE

src/CQKernel.o: src/CQKernel.h src/dsp/MathUtilities.h src/dsp/nan-inf.h
src/CQKernel.o: src/dsp/FFT.h src/dsp/Window.h
src/ConstantQ.o: src/CQKernel.h src/dsp/Resampler.h src/dsp/MathUtilities.h
src/ConstantQ.o: src/dsp/nan-inf.h src/dsp/FFT.h
src/CQInverse.o: src/dsp/Resampler.h src/dsp/MathUtilities.h
src/CQInverse.o: src/dsp/nan-inf.h src/dsp/FFT.h
vamp/CQVamp.o: vamp/CQVamp.h
vamp/CQChromaVamp.o: vamp/CQChromaVamp.h
vamp/libmain.o: vamp/CQVamp.h vamp/CQChromaVamp.h
cq/ConstantQ.o: cq/CQBase.h src/CQKernel.h
cq/CQSpectrogram.o: cq/ConstantQ.h cq/CQBase.h src/CQKernel.h
cq/CQInverse.o: cq/CQBase.h src/CQKernel.h
