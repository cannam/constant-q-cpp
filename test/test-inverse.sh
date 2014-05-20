#!/bin/bash
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

