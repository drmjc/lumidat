#!/bin/bash
# lumidat: Extract data from Illumina Gene Expression iDAT files
# how:
# - uses lumidat-1.2.jar
# You can pass in the idat files 3 ways:
# 1) as a zip file, using the "-inputfile mydata.zip" parameter
# 2) as trailing command line arguments
# 3) via stdin, by specifying '-' as the only commandline parameter (in addition to the options).
# 
# In stdin mode, make sure you've already escaped spaces and other odd chars
# (hint: see escapepath within idat.finder.sh)
#
# usage:
# lumidat.sh -manifestfile /path/to/manifest.txt
# cat batch1.ids | lumidat.sh -manifestfile /misc/ICGCPancreas/icgc_data/icgc_gex/2010-11-11/HumanHT-12_V4_0_R2_15002873_B.txt -bg false -collapse none -probeID ProbeID -prefix batch1 -
#
# CHANGELOG
# 2012-05-30: doesn't work if your prefix has spaces/hyphens in it...
# 2012-10-02: dropped readarray, which is non-portable; increased default memory
# 2013-01-08: updated to use lumidat-1.2, which natively handles the trailing '-' parameter.
# 
# Mark Cowley
################################################################################

# cause execution to terminate, with an optional message,
# and optional error code.
# $1: optional error message
# $2: optional error code. default = 1
function die {
	msg="######### ABORTING: ${1-Unexpected error}"
	exitcode=${2-1}
	echo >&2 "$msg"
	exit $exitcode
}

################################################################################

MEMORY=-Xmx4000M

#
# Where's the JAR?
# 
if [[ `uname` == "Darwin" ]]; then
    JAR=~/bin/lumidat-1.2.jar
elif [[ `uname` == "Linux" ]]; then 
    JAR=/misc/ICGCPancreas/bin/lumidat-1.2.jar
else
    JAR=./lumidat-1.2.jar
fi
[[ -e $JAR ]] || die "I can't find the JAR file.  Aborting."

#
# Run lumidat, accepting all command line arguments, and, if specified,
# the iDAT array names from stdin
#
# cat - | java $MEMORY -jar $JAR "$@"
cmd="cat - | java $MEMORY -jar $JAR "$@""
# echo $cmd
eval $cmd || die "This command failed: $cmd"
