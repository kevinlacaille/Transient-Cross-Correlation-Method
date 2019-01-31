#!/bin/bash

export STARLINK_DIR=/Users/kevinlacaille/star-2015B
#source $STARLINK_DIR/etc/cshrc
source $STARLINK_DIR/etc/login
source $STARLINK_DIR/etc/profile

export ADAM_NOPROMPT=1
export ADAM_EXIT=1
#export SMURF_THREADS=1

source $SMURF_DIR/smurf.sh > /dev/null
source $KAPPA_DIR/kappa.sh > /dev/null

set -e

smurf

fcf = 2294.2398997

set dirlist = "20160112_00056  20160112_00068  20160124_00056  20160125_00064  20160521_00012  20160112_00058  20160112_00069  20160125_00043  20160427_00026  20160521_00013  20160112_00059  20160118_00057  20160125_00050  20160428_00017  20160604_00037 20160112_00065  20160118_00064  20160125_00054  20160507_00024  20160604_00038"

for dir in $dirlist
do
   echo "doing $dir"

   cd $dir

   ndf2fits in=qmap.sdf out=qmap.fits
   ndf2fits in=umap.sdf out=umap.fits

   cd ..
done



maths "x*ia-ib" ia="20160115_shifted.sdf" ib="../../850_R1_crop/IC348_20151222_850_R1_f200_mjya2_crop.sdf"
