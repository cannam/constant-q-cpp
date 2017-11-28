
Constant-Q Library
==================

A C++ library and Vamp plugin implementing the Constant-Q transform
of a time-domain signal.

https://code.soundsoftware.ac.uk/projects/constant-q-cpp

The Constant-Q transform is a time-to-frequency-domain transform
related to the short-time Fourier transform, but with output bins
spaced logarithmically in frequency, rather than linearly. The output
bins are therefore linearly spaced in terms of musical pitch.

This implementation is reasonably fast and is causal, operating
block-by-block on the input, though with quite high latency depending
on the frequency range specified. By default it produces output at a
higher time resolution than some other implementations, using multiple
kernel atoms per time block. The inverse transform is approximate
rather than exact (see the paper cited below for details).


About this library
------------------

This library provides:

 * Forward transform: time-domain to complex Constant-Q bins
 * Forward spectrogram: time-domain to interpolated Constant-Q
   magnitude spectrogram
 * Inverse transform: complex Constant-Q bins to time domain

The Vamp plugin provides:

 * Constant-Q magnitude spectrogram with high and low frequency
   extents defined in Hz
 * Constant-Q magnitude spectrogram with high and low frequency
   extents defined as MIDI pitch values
 * Pitch chromagram obtained by folding a Constant-Q spectrogram
   around into a single-octave range


Building the library and plugin
-------------------------------

To compile this code, use "make -f <file>" where <file> is one of
Makefile.linux, Makefile.osx, Makefile.mingw32 depending on
platform. These files set up some flags and include Makefile.inc,
which defines the input files and so forth. You could equally write
your own Makefile which does the same.

The Vamp plugin part of the build expects to find a compiled version
of the Vamp Plugin SDK in a neighbouring folder (../vamp-plugin-sdk).
The unit tests also require Boost headers to be available.

Note that this code uses the KissFFT library (compiled from source
bundled here) in its non-default double-precision mode via the use of
the flag -Dkiss_fft_scalar=double in Makefile.inc. If you want to
build this code using an external KissFFT library or custom build
scripts, you may need to do some work to ensure you have the right
configuration.


Credits
-------

The method is drawn from Christian Schörkhuber and Anssi Klapuri,
"Constant-Q transform toolbox for music processing", SMC 2010. See the
file CITATION for details. If you use this code in research work,
please cite this paper.

The C++ implementation is by Chris Cannam, Copyright 2014-2015 Queen
Mary, University of London.

The library is provided under a liberal BSD/MIT-style open source
licence. See the file COPYING for more information.

