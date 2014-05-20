#!/bin/bash

# Test that the forward-inverse CQ transform produces the same output
# (to a given noise level, and within the appropriate frequency range)
# as its input.
#
# This requires the program "processfile" to be compiled and in the
# same directory, and an audio file filtered-whitenoise-480-14600.wav
# to be in the subdir "data". The audio file contains white noise
# band-limited to the 480-14600Hz range. This is fed as input to a
# forward-inverse CQ chain restricted to the range 465-14700 Hz (5
# octaves); the processfile program calculates the output and performs
# a sample-by-sample diff against the input. We then check that the
# diff is below a suitable noise floor.

mydir=`dirname "$0"`

process="$mydir/processfile"
if [ ! -x "$process" ]; then
    echo "ERROR: $mydir/processfile not found or not executable"
    exit 1
fi
infile="$mydir/data/filtered-whitenoise-480-14600.wav"
if [ ! -f "$infile" ]; then
    echo "ERROR: Test file $infile not found"
    exit 1
fi

outfile="/tmp/$$.out.wav"
difffile="/tmp/$$.diff.wav"
logfile="/tmp/$$.log.txt"
trap "rm -f ${outfile} ${difffile} ${logfile}" 0

"$process" -x 14700 -n 465 -b 36 "$infile" "$outfile" "$difffile" 2>&1 | tee "$logfile" || exit 1

int_db=`grep 'max diff' "$logfile" | sed 's/^[^(]*(//' | sed 's/[^0-9-].*//'`
good=`expr "$int_db" "<" "-20"`
if [ "$good" == "1" ]; then
    echo "Forward-inverse process is satisfactory"
    exit 0
else
    echo "Forward-inverse not OK: Rounded dB value $int_db is too high -- should be < -20"
    exit 1
fi

