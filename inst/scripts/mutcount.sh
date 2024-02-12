#!/usr/bin/bash

#fail on error, even within pipes; require variable declarations, disable history tracking
set -eEuo pipefail +H
# set -eEo pipefail +H

log() {
  printf '%s : %s\n' "$(date +%Y/%m/%d-%H:%M:%S)" "$*"
}
logWarn() {
  printf '%s \033[33mWARNING: %s\033[m\n' "$(date +%Y/%m/%d-%H:%M:%S)" "$*"
}
logErr() {
  printf '%s \033[31mERROR: %s\033[m\n' "$(date +%Y/%m/%d-%H:%M:%S)" "$*">&2
}

#if an error occurred, print where it happened
function handle_error {
  local retval=$?
  local line=$1
  logErr "Mutcount failed in line $line: $BASH_COMMAND"
  exit $retval
}
trap 'handle_error $LINENO' ERR

BLARG=""
CONDAARG=""
QUEUEARG=""
USEPHIX=0
CPUS=6
VARCALLMEM=15G
INDIR=""
WORKDIR=$(pwd)
PARAMETERS=""
VALIDATEFASTQ=1
SAMDIR=""
SCRATCH="/scratch"

#Testing parameters
# INDIR=CBS_FASTQ_fixed/
# PARAMETERS=parameters.json
# CONDAARG='--conda pacybara'

#Display Logo and version number
Rscript -e 'library(tileseqMave)'

#helper function to print usage information
usage () {
  cat << EOF

mutcount.sh v0.0.1 

by Jochen Weile <jochenweile@gmail.com> 2023

A new job manager for tileseq_mutcount

Usage: mutcount.sh [-f|--fastq <DIRECTORY>] [-p|--parameters <JSONFILE>]
       [-o|--output <DIRECTORY>]  [--usePhiX] [--skipValidation]
       [-e|--existingAlignments <DIRECTORY> ]
       [-t|--threads <INTEGER>] [-m|--memory <INTEGER>G]  
       [-b|--blacklist <BLACKLIST>] [-c|--conda <ENV>] [-q|--queue <QUEUENAME>]

-f|--fastq              : The input directory containing the fastq.gz files. Required.
-p|--parameters         : A parameter sheet JSON file. Required.
-o|--output             : The working directory to which results will be written. Defaults to current directory.
--usePhiX               : calibrate Q-scores based on phiX reads instead of WT controls
--skipValidation        : skip validation of FASTQ files
-e|--existingAlignments : Use existing alignments in the given folder
-t|--threads            : Number of threads to use per job for alignment and variant calling. Default: $CPUS
-m|--memory             : Amount of memory to request for each variant calling job. Default: $VARCALLMEM
-c|--conda              : Conda environment to activate for jobs. (Optional)
-b|--blacklist          : A comma-separated blacklist of nodes to avoid. (Optional)
-q|--queue              : Submit jobs to given HPC queue. (Optional)
--scratch               : Path to cluster's local scratch folder. Default: $SCRATCH

EOF
 exit $1
}

#Parse Arguments
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      usage 0
      shift
      ;;
    -f|--fastq)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        INDIR=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -o|--output)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        WORKDIR=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -p|--parameters)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        PARAMETERS=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -e|--existingAlignments)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        SAMDIR=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -b|--blacklist)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        BLARG="--blacklist $2"
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
    -t|--threads)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        CPUS=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    -m|--memory)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        VARCALLMEM=$2
        shift 2
      else
        echo "ERROR: Argument for $1 is missing" >&2
        usage 1
      fi
      ;;
    --scratch)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        SCRATCH=$2
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
    --skipValidation)
      VALIDATEFASTQ=0
      shift
      ;;
    --) # end of options indicates that the main command follows
      shift
      PARAMS="$PARAMS $*"
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

# INDIR="$1"
if [[ -z "$INDIR" ]]; then
  echo "No fastq directory provided! Use -f or --fastq to specify.">&2
  usage 1
elif ! [[ -d "$INDIR" ]]; then
  echo "$INDIR is not a valid directory!">&2
  usage 1
fi

# WORKDIR="$2"
if [[ -z "$WORKDIR" ]]; then
  logErr "No output directory provided!"
  usage 1
elif ! [[ -d "$WORKDIR" ]]; then
  logErr "$WORKDIR is not a valid directory!"
  usage 1
fi

#TODO: Make a result folder in the Output directory named:
# <tag>_YYYY-MM-DD-hh-mm-ss/

#parameter file
# PARAMETERS=${3:-parameters.json}
if [[ -z "$PARAMETERS" ]]; then
  logErr "No parameter sheet provided! Use -p or --parameters to specify."
  usage 1
elif ! [[ -r "$PARAMETERS" ]]; then
  #TODO: Auto-convert csv to json if provided?
  if [[ "$PARAMETERS" != *.json ]]; then
    logErr "Parameter file $PARAMETERS must be a .json file!"
    usage 1
  fi
  logErr "Parameter sheet $PARAMETERS doesn't exist or is unreadable!"
  exit 1
fi


#Helper function to extract an element from a json file
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

#helper function to find matching entries in an array.
#Note that the array is passed *by reference*, not by value!
#Arguments: (1) name of array, (2) string to match
matches() {
  #load the array with the given name via eval
  eval ARR=\( \${${1}[@]} \)
  for ((i=0;i<${#ARR[@]};i++)); do
    [[ "${ARR[$i]}" == $2 ]] && echo "$i"
  done
}

# #get the intersection of lists of indices 
# intersect() {
#   LISTA=$1
#   shift
#   while (( "$#" )); do
#     LISTB=$1
#     shift
#     LISTA="$(comm -12 <(echo "$LISTA"|sort) <(echo "$LISTB"|sort) |sort -n)"
#   done
#   echo "$LISTA"
# }

#PHASE 0: Validate parameters and files, setup template libraries

#Extract project name
PROJECT=$(extractJSON "$PARAMETERS" project)
#Replace space with dashes in project name
PROJECT="${PROJECT// /-}"
PROJECT="${PROJECT//_/-}"

#Create or find output directory
#if SAMDIR is defined, then were using existing alignments,
#so our ourput folder is the parent folder of SAMDIR
if [[ -n "$SAMDIR" ]]; then
  OUTDIR=$(dirname "$SAMDIR")
else
  #otherwise we make a new outdir from scratch
  TIMESTAMP=$(date +%Y-%m-%d-%H-%M-%S)
  OUTDIR="${WORKDIR}/${PROJECT}_${TIMESTAMP}"
  mkdir -p "$OUTDIR"
fi


#Read sample table
SAMPLETABLE="$(extractJSON "$PARAMETERS" samples)"
#Extract columns into arrays: SIDS=SampleIDs, STILES=SampleTiles,
# SCONDS=SampleConditions, STPS, SampleTimePoints, SREPS=SampleReplicates
mapfile -t SIDS < <(extractCol "$SAMPLETABLE" 1)
mapfile -t STILES < <(extractCol "$SAMPLETABLE" 2)
mapfile -t SCONDS < <(extractCol "$SAMPLETABLE" 3)
# mapfile -t STPS < <(extractCol "$SAMPLETABLE" 4)
# mapfile -t SREPS < <(extractCol "$SAMPLETABLE" 5)

# mapfile SAMPLELINES < <(extractJSON "$PARAMETERS" samples)

#Check if R1 and R2 FASTQ files can be found and read for each sample
MISSING=""; declare -a R1S; declare -a R2S
for ((i=0; i<${#SIDS[@]}; i++)); do
  SID="${SIDS[$i]}"
  R1S[$i]=$(echo "${INDIR}/${SID}_"*R1*fastq.gz)
  if [[ ! -r  "${R1S[$i]}" ]]; then
    MISSING=$(cons "${SID}_R1" "$MISSING" "; ")
  fi
  R2S[$i]=$(echo "${INDIR}/${SID}_"*R2*fastq.gz)
  if [[ ! -r "${R2S[$i]}" ]]; then
    MISSING=$(cons "${SID}_R2" "$MISSING" "; ")
  fi
done
if [[ "$USEPHIX" == "1" ]]; then
  PHIXR1=$(echo "${INDIR}/Undetermined_"*R1*fastq.gz)
  if [[ ! -r  "$PHIXR1" ]]; then
    MISSING=$(cons "Undetermined_R1" "$MISSING" "; ")
  fi
  PHIXR2=$(echo "${INDIR}/Undetermined_"*R2*fastq.gz)
  if [[ ! -r  "$PHIXR2" ]]; then
    MISSING=$(cons "Undetermined_R2" "$MISSING" "; ")
  fi
fi

if [[ -n "$MISSING" ]]; then
  logErr "No FASTQ files found for the following samples: $MISSING"
  exit 1
fi

#If we're using existing alignments then the whole first section can be skipped!
if [[ -z "$SAMDIR" ]]; then

  #Check that R1 and R2 files are not corrupted
  if [[ "$VALIDATEFASTQ" == 1 ]]; then
    printReadIDS() {
      FASTQ=$1
      zcat "$FASTQ"|awk '{if (NR % 4 == 1) {print $1}}'
    }

    BROKEN=""
    for ((i=0; i<${#SIDS[@]}; i++)); do
      printf "Checking FASTQ file integrity for sample %s ... " "${SIDS[$i]}"
      ISDIFF=$(diff -q <(printReadIDS "${R1S[$i]}") <(printReadIDS "${R2S[$i]}")||true)
      if [[ "$ISDIFF" == Files*differ ]]; then
        echo "❌"
        BROKEN=$(cons "${SID}" "$BROKEN" "; ")
      else
        echo "✅"
      fi
    done

    if [[ -n "$BROKEN" ]]; then
      logErr "R1/R2 FASTQ files for the following samples have mismatching cluster IDs: $BROKEN"
      exit 1
    fi
  fi

  # Create reference libraries
  REFDIR="${OUTDIR}/ref"
  mkdir -p "$REFDIR"

  GENE=$(extractJSON "$PARAMETERS" template/geneName)
  REFSEQ=$(extractJSON "$PARAMETERS" template/seq)
  REFFASTA="${REFDIR}/${GENE}.fasta"
  echo ">$GENE">"$REFFASTA"
  echo "$REFSEQ">>"$REFFASTA"

  log "Building reference library"
  bowtie2-build -f --quiet "$REFFASTA" "${REFFASTA%.fasta}"

  if [[ "$USEPHIX" == "1" ]]; then
    #download phix genome reference
    log "Downloading PhiX reference genome..."
    PHIX_REFSEQID="NC_001422.1"
    EFETCH_BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    PHIX_URL="${EFETCH_BASE}?db=nuccore&id=${PHIX_REFSEQID}&rettype=fasta&retmode=text"
    PHIXFASTA="${REFDIR}/phix.fasta"
    curl "$PHIX_URL">"$PHIXFASTA"

    log "Building PhiX reference library"
    bowtie2-build -f --quiet "$PHIXFASTA" "${PHIXFASTA%.fasta}"
  fi

  #PHASE 1: Run alignment jobs and monitor
  echo "
  ###########################################
  # PHASE 1: Performing sequence alignments #
  ###########################################
  "

  extractJobID() {
    RETVAL=${1//$'\n'/ }
    echo "${RETVAL##* }"
  }


  #helper function to submit alignment jobs
  submitAlignments() {
    JOBTAGS="$1"
    SAMDIR="$2"
    JOBS="" 
    # for ((i=0; i<${#SIDS[@]}; i++)); do
    for JOBTAG in $JOBTAGS; do
      #generate a random alphanumeric ID for the scratch file to avoid name collisions
      ALPHATAG=$(tr -dc A-Za-z0-9 </dev/urandom | head -c 16||true)
      SCRATCHFILE="${SCRATCH}/${JOBTAG}_${ALPHATAG}.sam"
      if [[ "$JOBTAG" == *phix* ]]; then
        if [[ "$JOBTAG" == *R1 ]]; then
          LOGFILE="${SAMDIR}/bowtie_phix_R1.log"
          SAMFILE="${SAMDIR}/phix_R1.sam"
          FQ="$PHIXR1"
        elif [[ "$JOBTAG" == *R2 ]]; then
          LOGFILE="${SAMDIR}/bowtie_phix_R2.log"
          SAMFILE="${SAMDIR}/phix_R2.sam"
          FQ="$PHIXR2"
        else 
          logErr "Invalid job tag. Report this as a bug!"
          exit 100
        fi
        RETVAL=$(submitjob.sh -n "bowtiePhix" -c "${CPUS}" -m 1G \
          -l "$LOGFILE" -e "$LOGFILE" $CONDAARG $QUEUEARG $BLARG \
          --report --skipValidation -- \
          bowtie2 --no-unal --no-head --no-sq --rdg 12,1 --rfg 12,1 \
          --local --threads "${CPUS}" --reorder -x "${PHIXFASTA%.fasta}" \
          -U "$FQ" -S "$SCRATCHFILE" '&&' mv "$SCRATCHFILE" "$SAMFILE")
        JOBS=$(cons "$(extractJobID "$RETVAL")" "$JOBS" ",")
      else
        i="${JOBTAG%%_*}"
        if [[ "$JOBTAG" == *R1 ]]; then
          LOGFILE="${SAMDIR}/bowtie_${SIDS[$i]}_R1.log"
          SAMFILE="${SAMDIR}/${SIDS[$i]}_R1.sam"
          FQ="${R1S[$i]}"
          FWREV="--norc"
        elif [[ "$JOBTAG" == *R2 ]]; then
          LOGFILE="${SAMDIR}/bowtie_${SIDS[$i]}_R2.log"
          SAMFILE="${SAMDIR}/${SIDS[$i]}_R2.sam"
          FQ="${R2S[$i]}"
          FWREV="--nofw"
        else 
          logErr "Invalid job tag. Report this as a bug!"
          exit 100
        fi
        #profiling showed that bowtie uses < 1GB RAM even with 8 threads
        # RETVAL=$(submitjob.sh -n "bowtie${SIDS[$i]}" -c "${CPUS}" -m 2G \
        RETVAL=$(submitjob.sh -n "bowtie${SIDS[$i]}" -c "${CPUS}" -m 2G \
          -l "$LOGFILE" -e "$LOGFILE" $CONDAARG $QUEUEARG $BLARG \
          --report --skipValidation -- \
          bowtie2 --no-head ${FWREV} --no-sq --rdg 12,1 --rfg 12,1 --local \
          --threads "${CPUS}" --reorder \
          -x "${REFFASTA%.fasta}" -U "$FQ" -S "$SCRATCHFILE" \
          '&&' mv "$SCRATCHFILE" "$SAMFILE" )
          # bowtie2 ${FWREV} --rdg 12,1 --rfg 12,1 --local \
          # --threads "${CPUS}" -x "${REFFASTA%.fasta}" -U "$FQ" \
          # '|samtools sort' -@ "$CPUS" -n -u - \
          # '|samtools view' --no-header -o "$SAMFILE" -O SAM -)
        JOBS=$(cons "$(extractJobID "$RETVAL")" "$JOBS" ",")
      fi
    done
    echo "$JOBS"
  }

  #helper function to check if any jobs have failed
  findFailedAlignments() {
    TAGS="$1"
    FAILEDJOBTAGS=""
    for JOBTAG in $TAGS; do
      LOGFILE="${SAMDIR}/bowtie_${JOBTAG#*_}.log"
      if [[ -e "$LOGFILE" ]]; then
        #FIXME: This won't work on PBS systems, as it appends extra lines with job stats to the end of logs
        # LASTLOGLINE=$(tail -1 "$LOGFILE")
        # if [[ "$LASTLOGLINE" != "Job completed successfully." ]]; then
        #   FAILEDJOBTAGS="$(cons "$JOBTAG" "$FAILEDJOBTAGS" ' ')"
        # fi
        if ! tail -4 "$LOGFILE"|grep -q "Job completed successfully"; then
          FAILEDJOBTAGS="$(cons "$JOBTAG" "$FAILEDJOBTAGS" ' ')"
        fi
      else
        logErr "Log file for $JOBTAG is missing!"
        FAILEDJOBTAGS="$(cons "$JOBTAG" "$FAILEDJOBTAGS" ' ')"
      fi
    done
    echo "$FAILEDJOBTAGS"
  }

  #Assemble list of alignment jobs to submit
  JOBTAGS=""
  for ((i=0; i<${#SIDS[@]}; i++)); do
    JOBTAGS="$(cons "${i}_${SIDS[$i]}_R1 ${i}_${SIDS[$i]}_R2" "$JOBTAGS" ' ')"
  done
  if [[ "$USEPHIX" == "1" ]]; then
    JOBTAGS="$(cons "NA_phix_R1 NA_phix_R2" "$JOBTAGS" ' ')"
  fi

  #submit alignment jobs
  SAMDIR="${OUTDIR}/sam_files"
  mkdir -p "${SAMDIR}"
  JOBS="$(submitAlignments "$JOBTAGS" "$SAMDIR")"
  waitForJobs.sh -v "$JOBS"

  #find any failed jobs and re-run them up to 3 times
  FAILEDJOBTAGS="$(findFailedAlignments "$JOBTAGS")"
  RETRIES=0
  while [[ -n "$FAILEDJOBTAGS" && "$RETRIES" -lt 3 ]]; do
    ((RETRIES++))||true
    logWarn "$(echo "$FAILEDJOBTAGS"|wc -w) alignment jobs failed and will be re-submitted."
    JOBS="$(submitAlignments "$FAILEDJOBTAGS" "$SAMDIR")"
    waitForJobs.sh -v "$JOBS"
    FAILEDJOBTAGS="$(findFailedAlignments "$FAILEDJOBTAGS")"
  done
  #if there are still failed jobs left, throw an error
  if [[ -n "$FAILEDJOBTAGS" ]]; then
    logErr "The following alignment jobs failed and have exhausted all re-tries: $FAILEDJOBTAGS"
    exit 1
  fi

  log "All $(echo "$JOBTAGS"|wc -w) alignment jobs reported success!"

fi #this is the end of "if [[ -z $SAMDIR]] "

log "Validating alignments..."

# VALIDATE SAM FILES
ISBROKEN=0
for ((i=0; i<${#SIDS[@]}; i++)); do
  SAMR1="${SAMDIR}/${SIDS[$i]}_R1.sam"
  SAMR2="${SAMDIR}/${SIDS[$i]}_R2.sam"
  if [[ -e "$SAMR1" && -e "$SAMR2" ]]; then
    ISDIFF=$(diff -q <(awk '{print $1}' "$SAMR1") <(awk '{print $1}' "$SAMR2")||true)
    if [[ "$ISDIFF" == Files*differ ]]; then
      log "Alignments for sample ${SIDS[$i]} ❌"
      ISBROKEN=1
    else
      log "Alignments for sample ${SIDS[$i]} ✅"
    fi
  else
    log "Alignments for sample ${SIDS[$i]} missing!"
    ISBROKEN=1
  fi
done
if [[ "$ISBROKEN" == 1 ]]; then
  logErr "Alignment (SAM) files failed validation."
  exit 1
fi

log "All alignments appear healthy."

#PHASE 2: Run calibrateQC jobs
echo "
#################################
# PHASE 2: Calibrating Q-scores #
#################################
"

#output file naming scheme required by tilseqMut: 
# f"{wt_id}_T{self._sample_tile}_R1_calibrate_phred.csv"

submitCalibrations() {
  TAGS="$1"
  OUTDIR="$2"
  JOBS=""
  for JOBTAG in $TAGS; do
    PREFIX="${JOBTAG#*_}"
    i="${JOBTAG%%_*}"
    R12="${PREFIX#*_}"
    # OUT="${OUTDIR}/${PREFIX}.csv"
    OUT="${OUTDIR}/${SIDS[$i]}_T${STILES[$i]}_${R12}_calibrate_phred.csv"
    LOG="${OUTDIR}/${PREFIX}.log"
    LOG2="${OUTDIR}/${PREFIX}_internal.log"
    SAMFILE="${SAMDIR}/${PREFIX}.sam"
    RETVAL=$(submitjob.sh -n "calibrate${SIDS[$i]}R1" -c "$CPUS" -m 1G \
      -l "$LOG" -e "$LOG" $CONDAARG $QUEUEARG $BLARG \
      --report --skipValidation -- \
      tsm calibratePhred "$SAMFILE" -p "$PARAMETERS" \
      -o "$OUT" -l "$LOG2" --cores "$CPUS")
    # echo "$RETVAL">&2
    #extract job ID from return value and add to jobs list
    JOBS=$(cons "$(extractJobID "$RETVAL")" "$JOBS" ',')
  done
  echo "$JOBS"
}

#helper function to check if any jobs have failed
findFailedCalibrations() {
  TAGS="$1"
  LOGDIR="$2"
  FAILEDJOBTAGS=""
  for JOBTAG in $TAGS; do
    LOGFILE="${LOGDIR}/${JOBTAG#*_}.log"
    # LASTLOGLINE=$(tail -1 "$LOGFILE")
    # if [[ "$LASTLOGLINE" != "Job completed successfully." ]]; then
    #   FAILEDJOBTAGS="$(cons "$JOBTAG" "$FAILEDJOBTAGS" ' ')"
    # fi
    if ! tail -4 "$LOGFILE"|grep -q "Job completed successfully"; then
      FAILEDJOBTAGS="$(cons "$JOBTAG" "$FAILEDJOBTAGS" ' ')"
    fi
  done
  echo "$FAILEDJOBTAGS"
}

#Find the appropriate WT sample for each non-WT sample
CONDEFS="$(extractJSON "$PARAMETERS" "conditions/definitions")"
# declare -a SWTCOND
# for ((i=0; i<${#SIDS[@]}; i++)); do
#   WTCONDNAME=$(echo "$CONDEFS" | grep "is_wt_control_for ${SCONDS[$i]}" | awk '{print $1}')
#   if [[ -n $WTCONDNAME ]]; then
#     #find the row indices that match all the following fields:
#     j=$(intersect "$(matches "SCONDS" $WTCONDNAME)" "$(matches "STILES" ${STILES[$i]})" "$(matches "STPS" "${STPS[$i]}")" "$(matches "SREPS" ${SREPS[$i]})")
#     SWTCOND[$i]=${SIDS[$j]}
#   fi
# done

#Assemble list of calibration jobs to be run later
JOBTAGS=""
if [[ "$USEPHIX" == 1 ]]; then
  JOBTAGS="NA_phix_R1 NA_phix_R2"
else # i.e. if USEPHIX=0
  WTCONDNAMES=$(grep "is_wt_control_for" < <(echo "$CONDEFS")|awk '{print $1}'|sort|uniq)
  for WTCONDNAME in $WTCONDNAMES; do
    for i in $(matches SCONDS "$WTCONDNAME"); do
      SID=${SIDS[$i]}
      JOBTAGS="$(cons "${i}_${SID}_R1 ${i}_${SID}_R2" "$JOBTAGS" ' ')"
    done
  done
fi

#make a directory for the calibration output
CALIBDIR="${OUTDIR}/calibrations/"
mkdir -p "$CALIBDIR"
JOBS="$(submitCalibrations "$JOBTAGS" "$CALIBDIR")"
waitForJobs.sh -v "$JOBS"

#find any failed jobs and re-run them up to 3 times
FAILEDJOBTAGS="$(findFailedCalibrations "$JOBTAGS" "$CALIBDIR")"
RETRIES=0
while [[ -n "$FAILEDJOBTAGS" && "$RETRIES" -lt 3 ]]; do
  ((RETRIES++))||true
  logWarn "$(echo "$FAILEDJOBTAGS"|wc -w) calibration jobs failed and will be re-submitted."
  JOBS="$(submitCalibrations "$FAILEDJOBTAGS" "$CALIBDIR")"
  waitForJobs.sh -v "$JOBS"
  FAILEDJOBTAGS="$(findFailedAlignments "$FAILEDJOBTAGS" "$CALIBDIR")"
done
#if there are still failed jobs left, throw an error
if [[ -n "$FAILEDJOBTAGS" ]]; then
  logErr "The following calibration jobs failed and have exhausted all re-tries: $FAILEDJOBTAGS"
  exit 1
fi

log "All $(echo "$JOBTAGS"|wc -w) calibration jobs reported success!"

#PHASE 3: Run variant calling jobs and monitor
echo "
#############################
# PHASE 3: Calling variants #
#############################
"

submitVarcalls() {
  TAGS="$1"
  VCDIR="$2"
  if [[ "$USEPHIX" == 1 ]]; then
    CALIBARG="--calibratePhredPhix"
  else
    CALIBARG="--calibratePhredWT"
  fi
  JOBS=""
  for JOBTAG in $TAGS; do
    SID="${JOBTAG#*_}"
    VCLOG="${VCDIR}/${SID}_varcall.log"
    R1SAM="${SAMDIR}/${SID}_R1.sam"
    R2SAM="${SAMDIR}/${SID}_R2.sam"
    RETVAL=$(submitjob.sh -n "PTM-${SID}" -c "$CPUS" -m "$VARCALLMEM" \
      -t "36:00:00" -l "$VCLOG" -e "$VCLOG" $CONDAARG $QUEUEARG $BLARG \
      --report --skipValidation -- \
      tileseq_mut -n "$PROJECT" -r1 "$R1SAM" -r2 "$R2SAM" -o "$VCDIR" \
      -p "$PARAMETERS" --skip_alignment -log info -c "$CPUS" \
      $CALIBARG \
    )
    # echo "$RETVAL">&2
    #extract job ID from return value and add to jobs list
    JOBS="$(cons "$(extractJobID "$RETVAL")" "$JOBS" ',')"
  done
  echo "$JOBS"
}

findFailedVarcalls() {
  TAGS="$1"
  VCDIR="$2"
  FAILEDJOBTAGS=""
  for JOBTAG in $TAGS; do
    SID="${JOBTAG#*_}"
    VCLOG="${VCDIR}/${SID}_varcall.log"
    # LASTLOGLINE=$(tail -1 "$VCLOG")
    # if [[ "$LASTLOGLINE" != "Job completed successfully." ]]; then
    #   FAILEDJOBTAGS="$(cons "$JOBTAG" "$FAILEDJOBTAGS" ' ')"
    # fi
    if ! tail -4 "$VCLOG"|grep -q "Job completed successfully"; then
      FAILEDJOBTAGS="$(cons "$JOBTAG" "$FAILEDJOBTAGS" ' ')"
    fi
  done
  echo "$FAILEDJOBTAGS"
}


#Assemble list of varcall jobs to submit
JOBTAGS=""
for ((i=0; i<${#SIDS[@]}; i++)); do
  JOBTAGS="$(cons "${i}_${SIDS[$i]}" "$JOBTAGS" ' ')"
done

TIMESTAMP=$(date +%Y-%m-%d-%H-%M-%S)
VARCALLDIR="${OUTDIR}/${PROJECT}_${TIMESTAMP}_mut_count"
mkdir -p "$VARCALLDIR"
#Move calibration results to varcall output folder (so that 
# tileseqMut knows where to find them)
mv "${CALIBDIR}/"*_calibrate_phred.csv "${VARCALLDIR}/"
JOBS=$(submitVarcalls "$JOBTAGS" "$VARCALLDIR")
waitForJobs.sh -v "$JOBS"


#find any failed jobs and re-run them up to 3 times
FAILEDJOBTAGS="$(findFailedVarcalls "$JOBTAGS" "$VARCALLDIR")"
RETRIES=0
while [[ -n "$FAILEDJOBTAGS" && "$RETRIES" -lt 3 ]]; do
  ((RETRIES++))||true
  logWarn "$(echo "$FAILEDJOBTAGS"|wc -w) variant calling jobs failed and will be re-submitted."
  JOBS="$(submitVarcalls "$FAILEDJOBTAGS" "$VARCALLDIR")"
  waitForJobs.sh -v "$JOBS"
  FAILEDJOBTAGS="$(findFailedVarcalls "$FAILEDJOBTAGS" "$VARCALLDIR")"
done
#if there are still failed jobs left, throw an error
if [[ -n "$FAILEDJOBTAGS" ]]; then
  logErr "The following variant calling jobs failed and have exhausted all re-tries: $FAILEDJOBTAGS"
  exit 1
fi

#PHASE 4: Validate results

FAILEDJOBTAGS=""
for JOBTAG in $JOBTAGS; do
  SID="${JOBTAG#*_}"
  # VCLOG="${VARCALLDIR}/${SID}_varcall.log"
  COUNTFILE="${VARCALLDIR}/counts_sample_${SID}.csv"
  NVARS=$(grep -cP '^c\.' "$COUNTFILE")
  if (( NVARS > 0 )); then
    log "Sample $SID has $NVARS unique variant calls."
  else
    lowWarn "No variants found in sample ${SID}!"
  fi
done

#CLEANUP
log "Cleaning up and compressing files..."

for SAM in "${SAMDIR}"/*.sam; do
  submitjob.sh -n cleanup -c 1 -m 1G $CONDAARG $BLARG -- \
  gzip "$SAM"
done
tar czvf "${SAMDIR}/logs.tgz" "${SAMDIR}/"*.log && rm "${SAMDIR}/"*.log
tar czvf "${VARCALLDIR}/calibrationLogs.tgz" "${CALIBDIR}/"* && rm -r "${CALIBDIR}"

log "Mutcount is done!"


# submitjob.sh -n mutcount -c 2 -m 2G -l mutcount.log -e mutcount.log \
#   --conda pacybara -- \
#   ~/projects/tileseqMave/inst/scripts/mutcount.sh \
#   -f ~/tileseq/SUMO1/SUMO1_FASTQ_joint \
#   -p ~/tileseq/SUMO1/parameters.json --conda pacybara
