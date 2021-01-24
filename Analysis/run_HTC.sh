#!/bin/zsh

echo "Running EventShape analysis"
echo "Repeat execution, by running:"
echo " ------------------------------------------------------------------------ "
echo "run_HTC.sh  $@"
echo " ------------------------------------------------------------------------ "
echo "  STEERING:   $1"
echo "  CHAIN:      $2"
echo "  OUTPUTDIR:  $3  (optional - default: 'output')"
echo "  PWD:        $4  (optional - default '.')"
echo " ------------------------------------------------------------------------ "
echo " "
echo "Setting up environment"
if [[ -n $4 ]]; then
    cd $4
fi


echo " ------------------------------------------------------------------------ "
echo "This node:"
uname -a
echo "HOST=$HOST"
echo "HOSTNAME=$HOSTNAME"
echo " ------------------------------------------------------------------------ "
env
echo " ------------------------------------------------------------------------ "

THISDIR=$PWD
export H1ANALYSISIDR=$PWD
cd /afs/desy.de/group/h1/root/checkout2/oo-releases/relvol10/releases/4.1.1
source env_2020.sh
source thish1.sh
# go back
cd $THISDIR

OUTDIR=output
if [[ -n $3 ]]; then
    OUTDIR=$3
    mkdir -p $OUTDIR
fi

echo ""
#echo "Now running analysis: $PWD/EventShapes"
#./EventShapes -t -f $1 -c $2 -o $OUTDIR/$2.nominal.root
#echo "Now running analysis: $PWD/EventShapes. SYSTEMATIC VARIATION 0 (Jet energy scale)"
#./EventShapes -t -f $1 -c $2 -o $OUTDIR/$2.sys_0.root -s 0
#echo "Now running analysis: $PWD/EventShapes. SYSTEMATIC VARIATION 1 (cluster energy scale)"
#./EventShapes -t -f $1 -c $2 -o $OUTDIR/$2.sys_1.root -s 1
#echo "Now running analysis: $PWD/EventShapes. SYSTEMATIC VARIATION 6 (electron energy)"
#./EventShapes -t -f $1 -c $2 -o $OUTDIR/$2.sys_6.root -s 6
echo "Now running analysis: $PWD/EventShapes. SYSTEMATIC VARIATION 9 (electron angle)"
./EventShapes -t -f $1 -c $2 -o $OUTDIR/$2.sys_9.root -s 9

