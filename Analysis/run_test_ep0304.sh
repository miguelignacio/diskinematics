#!/bin/bash

OUTDIR=Test0304Out
mkdir -p $OUTDIR

for chain in DataEplus0304  DataEplus0304_1 DjangoEplus0304_test  RapgapEplus0304_test 
do
    echo Running: $chain
    nohup nice ../bin/x86_64-centos7-gcc9-opt/EventShapes -f Steering/es.test.ep0304dis.steer -c $chain -o $OUTDIR/$chain.root &> $OUTDIR/log.$chain.txt &
done

echo ""
echo "--------------------------------------------------"
echo "  All jobs running in background."
echo "  This takes about ~2min."
echo "  Please check! (use $> ps l;  or $> htop)"
echo ""
echo "  Once done. merge with:"
echo "  $> htop -f merged.root $OUTDIR/*root"
echo ""
echo "  Then plot. using e.g. $> root -b -l plot_all.cxx"
echo "--------------------------------------------------"
echo ""
