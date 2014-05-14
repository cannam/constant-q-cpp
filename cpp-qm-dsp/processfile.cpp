
#include "ConstantQ.h"
#include "CQInverse.h"

#include <sndfile.h>

#include <iostream>

using std::vector;
using std::cerr;
using std::endl;

#include <cstring>

#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <cstdlib>

int main(int argc, char **argv)
{
    double maxFreq = 0;
    double minFreq = 0;
    int bpo = 0;
    bool help = false;
    
    int c;

    while (1) {
	int optionIndex = 0;
	
	static struct option longOpts[] = {
	    { "help", 0, 0, 'h', },
	    { "maxfreq", 1, 0, 'x', },
	    { "minfreq", 1, 0, 'n', },
	    { "bpo", 1, 0, 'b' },
	    { 0, 0, 0, 0 },
	};

	c = getopt_long(argc, argv,
			"hx:n:b:",
			longOpts, &optionIndex);
	if (c == -1) break;

	switch (c) {
	case 'h': help = true; break;
	case 'x': maxFreq = atof(optarg); break;
	case 'n': minFreq = atof(optarg); break;
	case 'b': bpo = atoi(optarg); break;
	default: help = true; break;
	}
    }

    if (help || (optind + 2 != argc && optind + 3 != argc)) {
	cerr << endl;
	cerr << "Usage: " << argv[0] << " [options] infile.wav outfile.wav [differencefile.wav]" << endl;
	cerr << endl;
	cerr << "Options:" << endl;
	cerr << "  -x<X>, --maxfreq <X>  Maximum frequency (default = sample rate / 3)" << endl;
	cerr << "  -n<X>, --minfreq <X>  Minimum frequency (default = 100, actual min may vary)" << endl;
	cerr << "  -b<X>, --bpo <X>      Bins per octave   (default = 60)" << endl;
	cerr << "  -h, --help            Print this help" << endl;
	cerr << endl;
	cerr << "This rather useless program simply performs a forward Constant-Q transform with" << endl;
	cerr << "the requested parameters, followed by its inverse, and writes the result to the" << endl;
	cerr << "output file. If a diff file name is provided, the difference between input and" << endl;
	cerr << "output signals is also written to that. All this accomplishes is to produce a" << endl;
	cerr << "signal that approximates the input: it's intended for test purposes only." << endl;
	cerr << endl;
	cerr << "(Want to calculate and obtain a Constant-Q spectrogram? Use the CQVamp plugin" << endl;
	cerr << "in a Vamp plugin host.)" << endl;
	cerr << endl;
	return 2;
    }

    char *fileName = strdup(argv[optind++]);
    char *fileNameOut = strdup(argv[optind++]);
    char *diffFileName = (optind < argc ? strdup(argv[optind++]) : 0);
    bool doDiff = (diffFileName != 0);

    SNDFILE *sndfile;
    SNDFILE *sndfileOut;
    SNDFILE *sndDiffFile = 0;
    SF_INFO sfinfo;
    SF_INFO sfinfoOut;
    SF_INFO sfinfoDiff;
    memset(&sfinfo, 0, sizeof(SF_INFO));

    sndfile = sf_open(fileName, SFM_READ, &sfinfo);
    if (!sndfile) {
	cerr << "ERROR: Failed to open input file \"" << fileName << "\": "
	     << sf_strerror(sndfile) << endl;
	return 1;
    }

    sfinfoOut.channels = 1;
    sfinfoOut.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    sfinfoOut.frames = sfinfo.frames;
    sfinfoOut.samplerate = sfinfo.samplerate;
    sfinfoOut.sections = sfinfo.sections;
    sfinfoOut.seekable = sfinfo.seekable;

    sndfileOut = sf_open(fileNameOut, SFM_WRITE, &sfinfoOut) ;
    if (!sndfileOut) {
	cerr << "ERROR: Failed to open output file \"" << fileNameOut << "\" for writing: "
	     << sf_strerror(sndfileOut) << endl;
	return 1;
    }

    if (doDiff) {
	sfinfoDiff = sfinfoOut;
	sndDiffFile = sf_open(diffFileName, SFM_WRITE, &sfinfoDiff);
	if (!sndDiffFile) {
	    cerr << "ERROR: Failed to open diff output file \"" << diffFileName << "\" for writing: "
		 << sf_strerror(sndDiffFile) << endl;
	    return 1;
	}
    }

    int ibs = 1024;
    int channels = sfinfo.channels;
    float *fbuf = new float[channels * ibs];

    if (maxFreq == 0.0) maxFreq = sfinfo.samplerate / 3;
    if (minFreq == 0.0) minFreq = 100;
    if (bpo == 0) bpo = 60;

    ConstantQ cq(sfinfo.samplerate, minFreq, maxFreq, bpo);
    CQInverse cqi(sfinfo.samplerate, minFreq, maxFreq, bpo);

    cerr << "max freq = " << cq.getMaxFrequency() << ", min freq = "
	 << cq.getMinFrequency() << ", octaves = " << cq.getOctaves() << endl;

    cerr << "octave boundaries: ";
    for (int i = 0; i < cq.getOctaves(); ++i) {
	cerr << cq.getMaxFrequency() / pow(2, i) << " ";
    }
    cerr << endl;

    int inframe = 0;
    int outframe = 0;
    int latency = cq.getLatency() + cqi.getLatency();

    vector<double> buffer;

    double maxdiff = 0.0;
    int maxdiffidx = 0;

    cerr << "forward latency = " << cq.getLatency() << ", inverse latency = " 
	 << cqi.getLatency() << ", total = " << latency << endl;

    timeval tv;
    (void)gettimeofday(&tv, 0);

    while (inframe < sfinfo.frames) {

        int count = -1;
	
	if ((count = sf_readf_float(sndfile, fbuf, ibs)) < 0) {
	    break;
	}

	vector<double> cqin;
	for (int i = 0; i < count; ++i) {
	    double v = fbuf[i * channels];
	    if (channels > 1) {
		for (int c = 1; c < channels; ++c) {
		    v += fbuf[i * channels + c];
		}
		v /= channels;
	    }
	    cqin.push_back(v);
	}

	if (doDiff) {
	    buffer.insert(buffer.end(), cqin.begin(), cqin.end());
	}
	
	vector<double> cqout = cqi.process(cq.process(cqin));

	for (int i = 0; i < int(cqout.size()); ++i) {
	    if (cqout[i] > 1.0) cqout[i] = 1.0;
	    if (cqout[i] < -1.0) cqout[i] = -1.0;
	}

	if (outframe >= latency) {

	    sf_writef_double(sndfileOut, 
			     cqout.data(), 
			     cqout.size());

	} else if (outframe + (int)cqout.size() >= latency) {

	    int offset = latency - outframe;
	    sf_writef_double(sndfileOut, 
			     cqout.data() + offset,
			     cqout.size() - offset);
	}

	if (doDiff) {
	    for (int i = 0; i < (int)cqout.size(); ++i) {
		if (outframe + i >= latency) {
		    int dframe = outframe + i - latency;
		    if (dframe >= (int)buffer.size()) cqout[i] = 0;
		    else cqout[i] -= buffer[dframe];
		    if (fabs(cqout[i]) > maxdiff &&
			dframe > sfinfo.samplerate && // ignore first/last sec
			dframe + sfinfo.samplerate < sfinfo.frames) {
			maxdiff = fabs(cqout[i]);
			maxdiffidx = dframe;
		    }
		}
	    }
	    
	    if (outframe >= latency) {

		sf_writef_double(sndDiffFile, 
				 cqout.data(), 
				 cqout.size());

	    } else if (outframe + (int)cqout.size() >= latency) {

		int offset = latency - outframe;
		sf_writef_double(sndDiffFile, 
				 cqout.data() + offset,
				 cqout.size() - offset);
	    }
	}

	inframe += count;
	outframe += cqout.size();
    }

    vector<double> r = cqi.process(cq.getRemainingOutput());
    vector<double> r2 = cqi.getRemainingOutput();

    r.insert(r.end(), r2.begin(), r2.end());

    for (int i = 0; i < int(r.size()); ++i) {
	if (r[i] > 1.0) r[i] = 1.0;
	if (r[i] < -1.0) r[i] = -1.0;
    }

    sf_writef_double(sndfileOut, r.data(), r.size());
    if (doDiff) {
	for (int i = 0; i < (int)r.size(); ++i) {
	    if (outframe + i >= latency) {
		int dframe = outframe + i - latency;
		if (dframe >= (int)buffer.size()) r[i] = 0;
		else r[i] -= buffer[dframe];
		if (fabs(r[i]) > maxdiff &&
		    dframe > sfinfo.samplerate && // ignore first/last sec
		    dframe + sfinfo.samplerate < sfinfo.frames) {
		    maxdiff = fabs(r[i]);
		    maxdiffidx = dframe;
		}
	    }
	}
	sf_writef_double(sndDiffFile, r.data(), r.size());
    }
    outframe += r.size();

    sf_close(sndfile);
    sf_close(sndfileOut);

    if (doDiff) {
	sf_close(sndDiffFile);
    }

    cerr << "in: " << inframe << ", out: " << outframe - latency << endl;

    if (doDiff) {
	double db = 10 * log10(maxdiff);
	cerr << "max diff [excluding first and last second of audio] is "
	     << maxdiff << " (" << db << " dBFS)"
	     << " at sample index " << maxdiffidx << endl;
    }
    
    timeval etv;
    (void)gettimeofday(&etv, 0);
        
    etv.tv_sec -= tv.tv_sec;
    if (etv.tv_usec < tv.tv_usec) {
	etv.tv_usec += 1000000;
	etv.tv_sec -= 1;
    }
    etv.tv_usec -= tv.tv_usec;
        
    double sec = double(etv.tv_sec) + (double(etv.tv_usec) / 1000000.0);
    cerr << "elapsed time (not counting init): " << sec << " sec, frames/sec at input: " << inframe/sec << endl;

    return 0;
}

