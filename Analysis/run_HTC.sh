#!/bin/zsh

echo "Running EventShape analysis"
echo "Repeat execution, by running:
echo " ------------------------------------------------------------------------ "
echo "$0 $@"
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
echo "Now running analysis: $PWD/EventShapes"
./EventShapes -f $1 -c $2 -o $OUTDIR/$2.root


