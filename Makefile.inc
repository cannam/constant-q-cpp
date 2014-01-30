
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

PLUGIN	:= cqvamp$(PLUGIN_EXT)
TEST	:= $(LIB_DIR)/test

LIB_HEADERS := $(LIB_DIR)/CQKernel.h $(LIB_DIR)/ConstantQ.h
LIB_SOURCES := $(LIB_DIR)/CQKernel.cpp $(LIB_DIR)/ConstantQ.cpp

VAMP_HEADERS := $(VAMP_DIR)/CQVamp.h
VAMP_SOURCES := $(VAMP_DIR)/CQVamp.cpp $(VAMP_DIR)/libmain.cpp

HEADERS	     := $(LIB_HEADERS) $(VAMP_HEADERS)
SOURCES	     := $(LIB_SOURCES) $(VAMP_SOURCES)
OBJECTS	     := $(SOURCES:.cpp=.o)

TEST_SOURCES := $(LIB_DIR)/test.cpp $(SOURCES)
TEST_OBJECTS := $(TEST_SOURCES:.cpp=.o)

LIBS	:= $(QMDSP_DIR)/libqm-dsp.a $(VAMPSDK_DIR)/libvamp-sdk.a -lpthread

all: $(PLUGIN) $(TEST)

$(PLUGIN):	$(OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(PLUGIN_LDFLAGS)

$(TEST):	$(TEST_OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(TEST_LDFLAGS)

clean:		
	rm -f $(OBJECTS)

distclean:	clean
	rm -f $(PLUGIN)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(HEADERS)

# DO NOT DELETE

cpp-qm-dsp/CQKernel.o: cpp-qm-dsp/CQKernel.h
cpp-qm-dsp/ConstantQ.o: cpp-qm-dsp/ConstantQ.h cpp-qm-dsp/CQKernel.h
vamp/CQVamp.o: vamp/CQVamp.h cpp-qm-dsp/ConstantQ.h cpp-qm-dsp/CQKernel.h
vamp/libmain.o: vamp/CQVamp.h
cpp-qm-dsp/ConstantQ.o: cpp-qm-dsp/CQKernel.h
