#!/usr/bin/bash

#fail on error, even within pipes; require variable declarations, disable history tracking
set -euo pipefail +H

#helper function to print usage information
usage () {
  cat << EOF

mutcount.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2021

A job manager for tileseq_mutcount

Usage: mutcount.sh [-b|--blacklist <BLACKLIST>] [-c|--conda <ENV>] <INDIR> <OUTDIR> <PARAMS>

<INDIR>        : The input directory containing the fastq.gz files
<OUTDIR>       : The output directory
<PARAMS>       : A parameter sheet JSON file
-b|--blacklist : An optional comma-separated blacklist of nodes to avoid
-c|--conda     : Conda environment to activate for jobs
-q|--queue     : Submit jobs to given HPC queue
--usePhiX      : calibrate Q-scores based on phiX reads instead of WT controls

EOF
 exit $1
}

#Parse Arguments
PARAMS=""
BLACKLIST=""
CONDAARG=""
QUEUEARG=""
USEPHIX=0
while (( "$#" )); do
  case "$1" in
    -h|--help)
      usage 0
      shift
      ;;
    -b|--blacklist)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        BLACKLIST=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -c|--conda)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        CONDAARG="--conda $2"
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -q|--queue)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        QUEUEARG="--queue $2"
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    --usePhiX)
      USEPHIX=1
      shift
      ;;
    --) # end of options indicates that the main command follows
      shift
      PARAMS="$PARAMS $@"
      eval set -- ""
      ;;
    -*|--*=) # unsupported flags
      echo "ERROR: Unsupported flag $1" >&2
      usage 1
      ;;
    *) # positional parameter
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
#reset command arguments as only positional parameters
eval set -- "$PARAMS"

INDIR="$1"
if [[ -z "$INDIR" ]]; then
  echo "No input directory provided!"
elif ! [[ -d "$INDIR" ]]; then
  echo "$INDIR is not a valid directory!">&2
  exit 1
fi

OUTDIR="$2"
if [[ -z "$OUTDIR" ]]; then
  echo "No output directory provided!"
elif ! [[ -d "$OUTDIR" ]]; then
  echo "$OUTDIR is not a valid directory!">&2
  exit 1
fi

#parameter file
PARAMETERS=${3:-parameters.json}
if ! [[ -r "$PARAMETERS" ]]; then
  if [[ "$PARAMETERS" != *.json ]]; then
    echo "Parameter file $PARAMETERS must be a .json file!">&2
    exit 1
  fi
  echo "$PARAMETERS not found or unreadable!">&2
  exit 1
fi

if [[ -z $BLACKLIST ]]; then
  BLARG=""
else
  BLARG="--blacklist $BLACKLIST"
fi

TMPDIR=$(mktemp -d -p ./)


#extract an element from a json file
# first argument: json file name
# second argument /-separated path to element
# for example: if json content is {'foo': 1, 'bar': {'baz': 'hello', 'buz': 4}}
# then the path "bar/baz" would extract element "hello"
extractJSON() {
  python3 -c '
import sys
import json
args = sys.argv[1:]
infile=args[0]
path=args[1].split("/")
with open(infile,"r") as stream:
  data = json.load(stream)
elem=data
for p in path:
  elem=elem[p]
if (type(elem) is dict):
  firstkey=next(iter(elem))
  if type(elem[firstkey]) is list:
    for i in range(len(elem[firstkey])):
      print(" ".join([str(elem[col][i]) for col in elem]))
  else:
    for k in elem:
      print("%s %s"%(k,str(elem[k])))
elif (type(elem) is list):
  if type(elem[0]) is dict:
    for i in range(len(elem)):
      print(" ".join([str(elem[i][key]) for key in elem[0]]))
  else:
    print(" ".join(elem))
else:
  print(elem)
' "$1" "$2"
}

#helper function to extract a table column
extractCol() {
  echo "$1"|awk "{print \$$2}"
}

#helper function to concatenate strings with a separator
#parameters: item, list, separator
cons() {
  SEP="${3:-;}"
  if [[ -z $2 ]]; then
    echo "$1"
  else
    echo "$2$SEP$1"
  fi
}

#helper function to find matching array entries
#parameters: arrray, string to match
matches() {
  #load the array with the given name via eval
  eval ARR=\( \${${1}[@]} \)
  for ((i=0;i<${#ARR[@]};i++)); do
    [[ "${ARR[$i]}" == $2 ]] && echo "$i"
  done
}

#get the intersection of lists of indices 
intersect() {
  LISTA=$1
  shift
  while (( "$#" )); do
    LISTB=$1
    shift
    LISTA="$(comm -12 <(echo "$LISTA"|sort) <(echo "$LISTB"|sort) |sort -n)"
  done
  echo "$LISTA"
}

#PHASE 0: Validate parameters and files, setup template libraries

#Read sample table
SAMPLETABLE="$(extractJSON "$PARAMETERS" samples)"
mapfile -t SIDS < <(extractCol "$SAMPLETABLE" 1)
mapfile -t STILES < <(extractCol "$SAMPLETABLE" 2)
mapfile -t SCONDS < <(extractCol "$SAMPLETABLE" 3)
mapfile -t STPS < <(extractCol "$SAMPLETABLE" 4)
mapfile -t SREPS < <(extractCol "$SAMPLETABLE" 5)

# mapfile SAMPLELINES < <(extractJSON "$PARAMETERS" samples)

#Check if R1 and R2 FASTQ files can be found and read for each sample
MISSING=""; declare -a R1S; declare -a R2S
for ((i=0; i<${#SIDS[@]}; i++)); do
  SID="${SIDS[$i]}"
  R1S[$i]=$(ls "${INDIR}/${SID}_"*R1*fastq.gz)
  if [[ ! -r  "${R1S[$i]}" ]]; then
    MISSING=$(cons "${SID}_R1" $MISSING)
  fi
  R2S[$i]=$(ls "${INDIR}/${SID}_"*R2*fastq.gz)
  if [[ ! -r "${R2S[$i]}" ]]; then
    MISSING=$(cons "${SID}_R2" $MISSING)
  fi
done

if [[ -n "$MISSING" ]]; then
  echo "ERROR: No FASTQ files found for the following samples: $MISSING">&2
  exit 1
fi

#Find the appropriate WT sample for each non-WT sample
CONDEFS="$(extractJSON $PARAMETERS conditions/definitions)"
declare -a SWTCOND
for ((i=0; i<${#SIDS[@]}; i++)); do
  WTCONDNAME=$(grep "is_wt_control_for ${SCONDS[$i]}" < <(echo "$CONDEFS")|awk '{print $1}')
  if [[ -n $WTCONDNAME ]]; then
    j=$(intersect "$(matches SCONDS $WTCONDNAME)" "$(matches STILES ${STILES[$i]})" "$(matches STPS "${STPS[$i]}")" "$(matches SREPS ${SREPS[$i]})")
    SWTCOND[$i]=${SIDS[$j]}
  fi
done

# Create reference libraries
REFDIR="${OUTDIR}/ref"
mkdir -p "$REFDIR"

PROJECT=$(extractJSON "$PARAMETERS" project)
REFSEQ=$(extractJSON "$PARAMETERS" template/seq)
REFFASTA="${REFDIR}/${PROJECT}.fasta"
echo ">$PROJECT">"$REFFASTA"
echo "$REFSEQ">>"$REFFASTA"

echo "Building reference library"
bowtie2-build -f "$REFFASTA" "${REFFASTA%.fasta}"

if [[ "$USEPHIX" == "1" ]]; then
  #download phix genome reference
  echo "Downloading PhiX reference genome..."
  PHIX_REFSEQID="NC_001422.1"
  PHIX_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${PHIX_REFSEQ}&rettype=fasta&retmode=text"
  PHIXFASTA="${REFDIR}/phix.fasta"
  curl "$PHIX_URL">"$PHIXFASTA"

  echo "Building PhiX reference library"
  bowtie2-build -f "$PHIXFASTA" "${PHIXFASTA%.fasta}"
fi

#PHASE 1: Run alignment jobs and monitor

extractJobID() {
  RETVAL=${1//$'\n'/ }
  echo "${RETVAL##* }"
}

SAMDIR="${OUTDIR}/sam_files"
mkdir -p ${SAMDIR}
JOBS=""; JOBTAGS=""
for ((i=0; i<${#SIDS[@]}; i++)); do
  LOGFILE="${SAMDIR}/bowtie_${SIDS[$i]}_R1.log"
  SAMFILE="${SAMDIR}/${SIDS[$i]}_R1.sam"
  RETVAL=$(submitjob.sh -n "bowtie${SIDS[$i]}" \
    -c ${CPUS} -m 1G -l $LOGFILE -e $LOGFILE $CONDAARG $QUEUEARG -- \
    bowtie2 --no-head --norc --no-sq --rdg 12,1 --rfg 12,1 --local --threads ${CPUS}\
    -x "${REFFASTA%.fasta}" -U "${R1S[$i]}" -S "$SAMFILE" )
  cons $(extractJobID "$RETVAL") "$JOBS" ","
  cons "${SIDS[$i]}_R1" "$JOBTAGS" " "

  LOGFILE="${SAMDIR}/bowtie_${SIDS[$i]}_R2.log"
  SAMFILE="${SAMDIR}/${SIDS[$i]}_R2.sam"
  RETVAL=$(submitjob.sh -n "bowtie${SIDS[$i]}" \
    -c ${CPUS} -m 1G -l $LOGFILE -e $LOGFILE $CONDAARG $QUEUEARG -- \
    bowtie2 --no-head --nofw --no-sq --rdg 12,1 --rfg 12,1 --local --threads ${CPUS}\
    -x "${REFFASTA%.fasta}" -U "${R2S[$i]}" -S "$SAMFILE" )
  cons $(extractJobID "$RETVAL") "$JOBS" ","
  cons "${SIDS[$i]}_R2" "$JOBTAGS" " "
done
if [[ "$USEPHIX" == "1" ]]; then

fi


# bowtie2 --no-head --norc --no-sq --rdg 12,1 --rfg 12,1 --local -x {ref} -U {os.path.abspath(r1)} -S {r1_sam_file}


#PHASE 2: Run calibrateQC jobs
echo "Running QC calibrations."

CALIBDIR="${OUTDIR}/calibrations/"
mkdir -p "$CALIBDIR"

JOBS=""

if [[ "$USEPHIX" == "0" ]]; then
  WTCONDNAMES=$(grep "is_wt_control_for" < <(echo "$CONDEFS")|awk '{print $1}'|sort|uniq)
  for WTCONDNAME in $WTCONDNAMES; do
    for i in $(matches SCONDS $WTCONDNAME); do
      # echo "${R1S[$i]} ${R2S[$i]}"
      #tsm calibratePhred {wt_r1_sam} -p {self._param} -o {phred_output_r1} -l {log_f} --silent --cores {self._cores}
      JOBOUT="${CALIBDIR}/${SIDS[$i]}_R1.csv"
      JOBLOG="${CALIBDIR}/${SIDS[$i]}_R1.log"
      JOBLOGINT="${CALIBDIR}/${SIDS[$i]}_R1_internal.log"
      RETVAL=$(submitjob.sh -n "calibrate${SIDS[$i]}R1" -c 8 -m 1G \
        -l "$JOBLOG" -e "$JOBLOG" $CONDAARG $QUEUEARG -- \
        tsm calibratePhred "${R1S[$i]}" -p "$PARAMETERS" \
        -o "$JOBOUT" -l "$JOBLOGINT" --cores 8)
      #extract job ID from return value and add to jobs list
      cons $(extractJobID "$RETVAL") $JOBS ","

      JOBOUT="${CALIBDIR}/${SIDS[$i]}_R2.csv"
      JOBLOG="${CALIBDIR}/${SIDS[$i]}_R2.log"
      JOBLOGINT="${CALIBDIR}/${SIDS[$i]}_R2_internal.log"
      RETVAL=$(submitjob.sh -n "calibrate${SIDS[$i]}R2" -c 8 -m 1G \
        -l "$JOBLOG" -e "$JOBLOG" $CONDAARG $QUEUEARG -- \
        tsm calibratePhred "${R2S[$i]}" -p "$PARAMETERS" \
        -o "$JOBOUT" -l "$JOBLOGINT" --cores 8)
      #extract job ID from return value and add to jobs list
      cons $(extractJobID "$RETVAL") $JOBS ","
    done
  done
else # if USEPHIX==1
  ###########################
  # TODO: Implement PhiX here
  ###########################
fi
waitForJobs.sh -v "$JOBS"


#PHASE 3: Run variant calling jobs and monitor

#PHASE 4: Validate results