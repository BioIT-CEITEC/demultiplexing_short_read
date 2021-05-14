#!/bin/bash

source /opt/anaconda3/etc/profile.d/conda.sh

STAGED_DIR=/mnt/ssd/ssd_1/snakemake

ORIGINAL_BCL_DIR=${1}; shift
TEMP_RAW_BCL_DIR=${1}; shift
BCL2FASTQ_LOGS_OUTPUT=${1}; shift
LOG=${1}; shift
RUNDIRNAME=${1}; shift

ALL_ARGS=""
for var in "$@"
do
    ALL_ARGS="${ALL_ARGS} ${var}"
done

echo "## COMMAND: rsync -rt ${ORIGINAL_BCL_DIR} ${TEMP_RAW_BCL_DIR}"
rsync -rt ${ORIGINAL_BCL_DIR} ${TEMP_RAW_BCL_DIR}

echo "## COMMAND: mkdir ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}"
mkdir ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}

echo "## LOG FILE: $LOG"

STLOG="${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/${RUNDIRNAME}.analysis.log"
echo "## LOG FILE (TEMP): $STLOG"


echo "## COMMAND: cat $LOG > $STLOG"
cat $LOG > $STLOG

#bcl2fastq command (does not need redirection of stdout and stderr into $LOG because this entire script is redirected into $LOG)
conda activate bcl2fastq
echo "## COMMAND: conda activate bcl2fastq"
conda env list
echo "## VERSION: $(bcl2fastq --version 2>&1 | grep 'bcl2fastq')"
echo "## COMMAND: bcl2fastq $ALL_ARGS 2>&1"
bcl2fastq $ALL_ARGS 2>&1


if [ $? -eq 0 ]
then
  cp -rt ${BCL2FASTQ_LOGS_OUTPUT}/${RUNDIRNAME} ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/*/
  LOG_DIR=$(dirname "${LOG}")
  echo "## COMMAND: cp -r ${LOG_DIR} ${BCL2FASTQ_LOGS_OUTPUT}"
  echo "## COMMAND: cp -rt ${BCL2FASTQ_LOGS_OUTPUT}/${RUNDIRNAME} ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/*/"
  cp -r ${LOG_DIR} ${BCL2FASTQ_LOGS_OUTPUT}
  cp -rt ${BCL2FASTQ_LOGS_OUTPUT}/${RUNDIRNAME} ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/*/
  ORIG_RUN_NAME=$(basename "${ORIGINAL_BCL_DIR}")
  rm -rf ${TEMP_RAW_BCL_DIR}${ORIG_RUN_NAME}
  rm -rf ${LOG_DIR}
  touch ${STAGED_DIR}/raw_fastq_temp/${RUNDIRNAME}/task_finished_OK
else
  rm -f $STLOG
  cp $LOG $STLOG
fi

echo "## COMMAND: conda deactivate"
conda deactivate


# conda deactivate
