#!/bin/bash

export STARLINK_DIR=/stardev
source $STARLINK_DIR/etc/cshrc
source $STARLINK_DIR/etc/login
source $STARLINK_DIR/etc/profile

export ADAM_NOPROMPT=1
export ADAM_EXIT=1
#export SMURF_THREADS=1

source $SMURF_DIR/smurf.sh > /dev/null
source $KAPPA_DIR/kappa.sh > /dev/null

set -e

smurf

makemap in="../../raw_data/DATE/*.sdf" out="DATE_shifted.sdf" config=^../dimmconfig_R1_f200.lis pixsize=3 pointing="DATE_offset"

