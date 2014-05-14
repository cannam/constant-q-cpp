
LIB_DIR	  := cpp-qm-dsp
VAMP_DIR  := vamp

QMDSP_DIR ?= ../qm-dsp
VAMPSDK_DIR  ?= ../vamp-plugin-sdk

PLUGIN_EXT	?= .so

CXX	?= g++
CC	?= gcc

CFLAGS := $(CFLAGS) 
CXXFLAGS := -I. -I$(VAMPSDK_DIR) -I$(QMDSP_DIR) $(CXXFLAGS)

LDFLAGS := $(LDFLAGS) 
PLUGIN_LDFLAGS := $(LDFLAGS) $(PLUGIN_LDFLAGS)
TEST_LDFLAGS := $(LDFLAGS) -lboost_unit_test_framework
PF_LDFLAGS := $(LDFLAGS) -lsndfile

PLUGIN	:= cqvamp$(PLUGIN_EXT)
TEST	:= $(LIB_DIR)/test
PF	:= $(LIB_DIR)/processfile

LIB_HEADERS := $(LIB_DIR)/CQBase.h $(LIB_DIR)/CQKernel.h $(LIB_DIR)/ConstantQ.h $(LIB_DIR)/CQSpectrogram.h $(LIB_DIR)/CQInverse.h
LIB_SOURCES := $(LIB_DIR)/CQKernel.cpp $(LIB_DIR)/ConstantQ.cpp $(LIB_DIR)/CQSpectrogram.cpp $(LIB_DIR)/CQInverse.cpp

VAMP_HEADERS := $(VAMP_DIR)/CQVamp.h $(VAMP_DIR)/CQChromaVamp.h
VAMP_SOURCES := $(VAMP_DIR)/CQVamp.cpp $(VAMP_DIR)/CQChromaVamp.cpp $(VAMP_DIR)/libmain.cpp

HEADERS	     := $(LIB_HEADERS) $(VAMP_HEADERS)
SOURCES	     := $(LIB_SOURCES) $(VAMP_SOURCES)
OBJECTS	     := $(SOURCES:.cpp=.o)

TEST_SOURCES := $(LIB_DIR)/test.cpp 
TEST_OBJECTS := $(TEST_SOURCES:.cpp=.o) $(OBJECTS)

PF_SOURCES := $(LIB_DIR)/processfile.cpp
PF_OBJECTS := $(PF_SOURCES:.cpp=.o) $(OBJECTS)

LIBS	:= $(QMDSP_DIR)/libqm-dsp.a $(VAMPSDK_DIR)/libvamp-sdk.a -lpthread

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

cpp-qm-dsp/CQKernel.o: cpp-qm-dsp/CQKernel.h
cpp-qm-dsp/ConstantQ.o: cpp-qm-dsp/ConstantQ.h cpp-qm-dsp/CQBase.h
cpp-qm-dsp/ConstantQ.o: cpp-qm-dsp/CQKernel.h
cpp-qm-dsp/CQSpectrogram.o: cpp-qm-dsp/CQSpectrogram.h cpp-qm-dsp/ConstantQ.h
cpp-qm-dsp/CQSpectrogram.o: cpp-qm-dsp/CQBase.h cpp-qm-dsp/CQKernel.h
cpp-qm-dsp/CQInverse.o: cpp-qm-dsp/CQInverse.h cpp-qm-dsp/CQBase.h
cpp-qm-dsp/CQInverse.o: cpp-qm-dsp/CQKernel.h
vamp/CQVamp.o: vamp/CQVamp.h cpp-qm-dsp/CQSpectrogram.h
vamp/CQVamp.o: cpp-qm-dsp/ConstantQ.h cpp-qm-dsp/CQBase.h
vamp/CQVamp.o: cpp-qm-dsp/CQKernel.h
vamp/CQChromaVamp.o: vamp/CQChromaVamp.h cpp-qm-dsp/CQSpectrogram.h
vamp/CQChromaVamp.o: cpp-qm-dsp/ConstantQ.h cpp-qm-dsp/CQBase.h
vamp/CQChromaVamp.o: cpp-qm-dsp/CQKernel.h
vamp/libmain.o: vamp/CQVamp.h cpp-qm-dsp/CQSpectrogram.h
vamp/libmain.o: cpp-qm-dsp/ConstantQ.h cpp-qm-dsp/CQBase.h
vamp/libmain.o: cpp-qm-dsp/CQKernel.h
cpp-qm-dsp/test.o: cpp-qm-dsp/CQSpectrogram.h cpp-qm-dsp/ConstantQ.h
cpp-qm-dsp/test.o: cpp-qm-dsp/CQBase.h cpp-qm-dsp/CQKernel.h
cpp-qm-dsp/processfile.o: cpp-qm-dsp/ConstantQ.h cpp-qm-dsp/CQBase.h
cpp-qm-dsp/processfile.o: cpp-qm-dsp/CQKernel.h cpp-qm-dsp/CQInverse.h
cpp-qm-dsp/ConstantQ.o: cpp-qm-dsp/CQBase.h cpp-qm-dsp/CQKernel.h
cpp-qm-dsp/CQSpectrogram.o: cpp-qm-dsp/ConstantQ.h cpp-qm-dsp/CQBase.h
cpp-qm-dsp/CQSpectrogram.o: cpp-qm-dsp/CQKernel.h
cpp-qm-dsp/CQInverse.o: cpp-qm-dsp/CQBase.h cpp-qm-dsp/CQKernel.h
vamp/CQVamp.o: cpp-qm-dsp/CQSpectrogram.h cpp-qm-dsp/ConstantQ.h
vamp/CQVamp.o: cpp-qm-dsp/CQBase.h cpp-qm-dsp/CQKernel.h
