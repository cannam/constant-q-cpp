
LIB_DIR	  := src
INC_DIR	  := cq
TEST_DIR  := test
VAMP_DIR  := vamp
KFFT_DIR  := src/ext/kissfft

VAMPSDK_DIR  ?= ../vamp-plugin-sdk

PLUGIN_EXT	?= .so

CXX	?= g++
CC	?= gcc
#VALGRIND	?= valgrind -q

GENERAL_FLAGS	:= -I. -I$(VAMPSDK_DIR) -I$(INC_DIR) -I$(LIB_DIR) -I$(KFFT_DIR) -I$(KFFT_DIR)/tools -Dkiss_fft_scalar=double

CFLAGS := $(GENERAL_FLAGS) $(CFLAGS) 
CXXFLAGS := $(GENERAL_FLAGS) $(CXXFLAGS)

LDFLAGS := $(LDFLAGS) 
PLUGIN_LDFLAGS := $(LDFLAGS) $(PLUGIN_LDFLAGS)
TEST_LDFLAGS := $(LDFLAGS) -lboost_unit_test_framework
PF_LDFLAGS := $(LDFLAGS) -lsndfile

LIB	:= libcq.a
PLUGIN	:= cqvamp$(PLUGIN_EXT)
PF	:= $(TEST_DIR)/processfile

LIB_HEADERS	:= \
	$(INC_DIR)/CQBase.h \
	$(INC_DIR)/CQKernel.h \
	$(INC_DIR)/ConstantQ.h \
	$(INC_DIR)/CQSpectrogram.h \
	$(INC_DIR)/CQInverse.h \
	$(LIB_DIR)/dsp/FFT.h \
	$(LIB_DIR)/dsp/KaiserWindow.h \
	$(LIB_DIR)/dsp/MathUtilities.h \
	$(LIB_DIR)/dsp/nan-inf.h \
	$(LIB_DIR)/dsp/Resampler.h \
	$(LIB_DIR)/dsp/SincWindow.h \
	$(LIB_DIR)/dsp/Window.h \
	$(KFFT_DIR)/kiss_fft.h \
	$(KFFT_DIR)/tools/kiss_fftr.h

LIB_SOURCES	:= \
	$(LIB_DIR)/CQKernel.cpp \
	$(LIB_DIR)/ConstantQ.cpp \
	$(LIB_DIR)/CQSpectrogram.cpp \
	$(LIB_DIR)/CQInverse.cpp \
	$(LIB_DIR)/dsp/FFT.cpp \
	$(LIB_DIR)/dsp/KaiserWindow.cpp \
	$(LIB_DIR)/dsp/MathUtilities.cpp \
	$(LIB_DIR)/dsp/Resampler.cpp \
	$(LIB_DIR)/dsp/SincWindow.cpp \
	$(KFFT_DIR)/kiss_fft.c \
	$(KFFT_DIR)/tools/kiss_fftr.c

VAMP_HEADERS	:= \
	$(VAMP_DIR)/CQVamp.h \
	$(VAMP_DIR)/CQChromaVamp.h \
	$(VAMP_DIR)/Pitch.h

VAMP_SOURCES	:= \
	$(VAMP_DIR)/CQVamp.cpp \
	$(VAMP_DIR)/CQChromaVamp.cpp \
	$(VAMP_DIR)/libmain.cpp \
	$(VAMP_DIR)/Pitch.cpp

TEST_SOURCES	:= \
	$(TEST_DIR)/TestFFT.cpp \
	$(TEST_DIR)/TestMathUtilities.cpp \
	$(TEST_DIR)/TestResampler.cpp \
	$(TEST_DIR)/TestWindow.cpp \
	$(TEST_DIR)/TestCQKernel.cpp \
	$(TEST_DIR)/TestCQFrequency.cpp \
	$(TEST_DIR)/TestCQTime.cpp

HEADERS	     := $(LIB_HEADERS) $(VAMP_HEADERS)
SOURCES	     := $(LIB_SOURCES) $(VAMP_SOURCES)

LIB_OBJECTS  := $(LIB_SOURCES:.cpp=.o)
LIB_OBJECTS  := $(LIB_OBJECTS:.c=.o)

OBJECTS	     := $(SOURCES:.cpp=.o)
OBJECTS	     := $(OBJECTS:.c=.o)

TEST_OBJECTS	:= $(TEST_SOURCES:.cpp=.o)
TEST_TARGETS	:= $(TEST_SOURCES:.cpp=)

PF_SOURCES := $(TEST_DIR)/processfile.cpp
PF_OBJECTS := $(PF_SOURCES:.cpp=.o) $(OBJECTS)

LIBS	:= $(VAMPSDK_DIR)/libvamp-sdk.a -lpthread

all: $(LIB) $(PLUGIN) $(TEST_TARGETS) $(PF)
	for t in $(TEST_TARGETS); do echo; echo "Running $$t"; $(VALGRIND) ./"$$t" || exit 1; done && echo && $(VALGRIND) "./test/test-inverse.sh" && echo 'Tests complete'

$(PLUGIN):	$(OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(PLUGIN_LDFLAGS)

$(PF):	$(PF_OBJECTS)
	$(CXX) -o $@ $^ $(LIBS) $(PF_LDFLAGS)

$(LIB):	$(LIB_OBJECTS)
	ar cr $@ $^
	ranlib $@

Test%:	Test%.o
	$(CXX) -o $@ $^ $(LIB) $(LIBS) $(TEST_LDFLAGS)

clean:		
	rm -f $(OBJECTS) $(TEST_OBJECTS) $(PF_OBJECTS)

distclean:	clean
	rm -f $(PLUGIN) $(TEST_TARGETS)

depend:
	makedepend -Y -fMakefile.inc $(SOURCES) $(TEST_SOURCES) $(PF_SOURCES) $(HEADERS)

# DO NOT DELETE

src/CQKernel.o: src/dsp/MathUtilities.h src/dsp/nan-inf.h src/dsp/FFT.h
src/CQKernel.o: src/dsp/Window.h
src/ConstantQ.o: src/dsp/Resampler.h src/dsp/MathUtilities.h
src/ConstantQ.o: src/dsp/nan-inf.h src/dsp/FFT.h
src/CQInverse.o: src/dsp/Resampler.h src/dsp/MathUtilities.h
src/CQInverse.o: src/dsp/nan-inf.h src/dsp/FFT.h
src/dsp/FFT.o: src/dsp/FFT.h src/dsp/MathUtilities.h src/dsp/nan-inf.h
src/dsp/KaiserWindow.o: src/dsp/KaiserWindow.h src/dsp/MathUtilities.h
src/dsp/KaiserWindow.o: src/dsp/nan-inf.h
src/dsp/MathUtilities.o: src/dsp/MathUtilities.h src/dsp/nan-inf.h
src/dsp/Resampler.o: src/dsp/Resampler.h src/dsp/MathUtilities.h
src/dsp/Resampler.o: src/dsp/nan-inf.h src/dsp/KaiserWindow.h
src/dsp/Resampler.o: src/dsp/SincWindow.h
src/dsp/SincWindow.o: src/dsp/SincWindow.h
src/ext/kissfft/kiss_fft.o: src/ext/kissfft/_kiss_fft_guts.h
src/ext/kissfft/kiss_fft.o: src/ext/kissfft/kiss_fft.h
src/ext/kissfft/tools/kiss_fftr.o: src/ext/kissfft/tools/kiss_fftr.h
src/ext/kissfft/tools/kiss_fftr.o: src/ext/kissfft/kiss_fft.h
src/ext/kissfft/tools/kiss_fftr.o: src/ext/kissfft/_kiss_fft_guts.h
vamp/CQVamp.o: vamp/CQVamp.h cq/CQSpectrogram.h cq/ConstantQ.h cq/CQBase.h
vamp/CQVamp.o: cq/CQParameters.h cq/CQKernel.h vamp/Pitch.h
vamp/CQChromaVamp.o: vamp/CQChromaVamp.h cq/CQSpectrogram.h cq/ConstantQ.h
vamp/CQChromaVamp.o: cq/CQBase.h cq/CQParameters.h cq/CQKernel.h vamp/Pitch.h
vamp/libmain.o: vamp/CQVamp.h cq/CQSpectrogram.h cq/ConstantQ.h cq/CQBase.h
vamp/libmain.o: cq/CQParameters.h cq/CQKernel.h vamp/CQChromaVamp.h
vamp/Pitch.o: vamp/Pitch.h
test/TestCQFrequency.o: cq/CQSpectrogram.h cq/ConstantQ.h cq/CQBase.h
test/TestCQFrequency.o: cq/CQParameters.h cq/CQKernel.h src/dsp/Window.h
test/TestCQTime.o: cq/CQSpectrogram.h cq/ConstantQ.h cq/CQBase.h
test/TestCQTime.o: cq/CQParameters.h cq/CQKernel.h src/dsp/Window.h
test/TestFFT.o: src/dsp/FFT.h
test/TestMathUtilities.o: src/dsp/MathUtilities.h src/dsp/nan-inf.h
test/TestResampler.o: src/dsp/Resampler.h src/dsp/Window.h src/dsp/FFT.h
test/TestWindow.o: src/dsp/Window.h
test/processfile.o: cq/ConstantQ.h cq/CQBase.h cq/CQParameters.h
test/processfile.o: cq/CQKernel.h
cq/CQKernel.o: cq/CQParameters.h
cq/ConstantQ.o: cq/CQBase.h cq/CQParameters.h cq/CQKernel.h
cq/CQSpectrogram.o: cq/ConstantQ.h cq/CQBase.h cq/CQParameters.h
cq/CQSpectrogram.o: cq/CQKernel.h
cq/CQInverse.o: cq/CQBase.h cq/CQKernel.h cq/CQParameters.h
src/dsp/MathUtilities.o: src/dsp/nan-inf.h
src/ext/kissfft/tools/kiss_fftr.o: src/ext/kissfft/kiss_fft.h
vamp/CQVamp.o: cq/CQSpectrogram.h cq/ConstantQ.h cq/CQBase.h
vamp/CQVamp.o: cq/CQParameters.h cq/CQKernel.h
