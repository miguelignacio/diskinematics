#!/bin/bash
#
# ------------------------------------------------ #
#
# Submit the entire analysis to the bird HTC cluster
# 
# 1) generate a new directory to 
#    + preserve the code/steerings
#    + hold the output files
# 2) run all jobs, by executing at a bird node: run_HTC.sh
#
# ------------------------------------------------ #


if [[ -z $2 ]]; then
    echo "Usage: "
    echo "$0 <jobname> <QUEUE-file>"
    echo ""
    echo "Example:"
    echo "submit_HTC.sh jobname QUEUE_all.txt"
    exit 0
fi

QUEUEFILE=$2
if [[ -f $QUEUEFILE ]] ; then
    echo "OK! QUEUE-file found. Staring with:"
    head -2 $QUEUEFILE
    echo "..."
    tail -2 $QUEUEFILE
else
    echo "ERROR! QUEUE-file NOT found!"
    exit 1
fi

# making job directory
JOBNAME=`date +%y-%m-%d-%H-%M`-$1
echo $JOBNAME
if [[ -d $JOBNAME ]]; then
    echo "Output directory for the given job name already exists. Exiting"
    exit 0
fi


# compile executable and library
make all

# copy all files to job directory
echo "Now copying files to work directory $JOBNAME"
ln -s ../LumiFiles .
mkdir -p $JOBNAME/lib
cp *.C *.h *.steer $JOBNAME/.
cp -r ../lib/$PLATFORM/ $JOBNAME/lib/.
cp ../bin/$PLATFORM/EventShapes $JOBNAME/.
cp -r Steering $JOBNAME/.
cp QUEUE* $JOBNAME/.
cp run_HTC.sh $JOBNAME/.

mkdir -p $JOBNAME/output
mkdir -p $JOBNAME/log

# go to job directory and make CONDOR submit file:
cd $JOBNAME
cat <<EOF > HTC.condor
executable               = run_HTC.sh
# # switch to either transfer the executable per job to the node
# # or from the shared storage for all jobs (don't touch during your job upstart!)
transfer_executable      = True
universe                 = vanilla
arguments                = \$(STEERING) \$(CHAIN) \$(OUT) ${PWD}
output                   = log/\$(CHAIN).out
error                    = log/\$(CHAIN).error
log                      = log/\$(CHAIN).log
RequestMemory            = 2048
should_transfer_files    = Yes
getenv                   = False
when_to_transfer_output  = ON_EXIT
requirements             = OpSysAndVer == "CentOS7"
+RequestRuntime          = 3600
#+RequestRuntime         = 86399
#+RequestRuntime         = 604800
#+RequestRuntime         = 259200
+MyProject               = "h1"
#+MyProject              = "atlas"
queue STEERING  CHAIN OUT from ${QUEUEFILE}
EOF

echo "Submitting jobs for runperiod: all"
condor_submit HTC.condor
