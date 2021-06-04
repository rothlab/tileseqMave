#!/bin/bash

usage () {
  cat << EOF

checkFASTQ.sh v0.0.1

Compares the R1 and R2 FASTQ files for each sample in a given directory.
Usage: checkFASTQ <DIR> [-v|--verbose] [-l|--log <LOGFILE>]
  -v | --verbose : Print more details about progress
  -l | --log : The log file to which file differences will be written.
               Defaults to checkFASTQ.log
  -h | --help : Show this help text
EOF
 exit $1
}

#Parse Arguments
PARAMS=""
VERBOSE=1 #OFF by default
LOGFILE="checkFASTQ.log"
while (( "$#" )); do
  case "$1" in
    -v|--verbose)
      VERBOSE=0
      shift
      ;;
    -h|--help)
      usage 0
      shift
      ;;
    -l|--log)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        LOGFILE=$2
        echo "Setting log file $LOGFILE"
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      usage 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
# set positional arguments in their proper place
eval set -- "$PARAMS"

#The first positional argument is the work dir (defaults to "./")
WORKDIR=${1:-"."}

#create temporary files
T1=$(mktemp)
T2=$(mktemp)
TDIFF=$(mktemp)

#overwrite existing logfile
rm -f $LOGFILE

#Find all samples in the work directory
SAMPLES=$(basename -a `ls ${WORKDIR}/*R1*fastq.gz`|cut -f 1 -d_) 

if [ -z "$SAMPLES" ]; then
  echo "ERROR: No matching files found!"
  exit 1
fi

#final status of all checks
STATUS=0

for SAMPLE in $SAMPLES; do

  if [ $VERBOSE -eq 0 ]; then
  	echo "Checking $SAMPLE"
  fi

  #find the R1 and R2 files
  R1=$(ls ${WORKDIR}/${SAMPLE}_*R1*.fastq.gz)
  R2=$(ls ${WORKDIR}/${SAMPLE}_*R2*.fastq.gz)

  # Check that all files are present and unique
  N1=$(echo R1|wc -w)
  N2=$(echo R2|wc -w)
  if [ -z $R1 ]; then
    echo "ERROR: No R1 file found for sample $SAMPLE"|tee -a $LOGFILE
    exit 1
  elif [ -z $R2 ]; then
    echo "ERROR: No R2 file found for sample $SAMPLE"|tee -a $LOGFILE
    exit 1
  elif (( ($N1 > 1) || ($N2 > 1) )); then
    echo "ERROR: No exact matches for sample $SAMPLE"|tee -a $LOGFILE
    exit 1
  fi

  #copy the FASTQ headers of each file into temp files
  # zgrep -P "^@" $R1|cut -f1 -d" ">$T1
  # zgrep -P "^@" $R2|cut -f1 -d" ">$T2
  zcat $R1|sed -n '1~4p'|cut -f1 -d" ">$T1
  zcat $R2|sed -n '1~4p'|cut -f1 -d" ">$T2

  #check the header files for differences
  diff $T1 $T2>$TDIFF
  NUMDIFF=$(cat $TDIFF|wc -l)
  if [ "$NUMDIFF" -ne 0 ]; then
  	echo "FAILED: $NUMDIFF reads in sample $SAMPLE differ between R1 and R2"|tee -a $LOGFILE
  	cat $TDIFF>>$LOGFILE
  	printf "\n\n">>$LOGFILE
  elif [ $VERBOSE -eq 0 ]; then
  	echo "OK"
  fi

  #if the check failed, update status 
  #unintuitively we have to use OR instead of AND because booleans
  #work backwards in BASH
  (( STATUS=$STATUS||$NUMDIFF ))

done

#delete temp files
rm $T1 $T2 

if [ "$STATUS" -eq 0 ]; then
  	echo "All files match!"|tee -a $LOGFILE
fi

exit $STATUS
