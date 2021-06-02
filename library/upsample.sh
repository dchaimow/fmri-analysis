#!/bin/bash
#
# upsample.sh <filename> <factor> <method>
#
# -upsamples in all 3 dimensions using AFNI 3dresample

mapfile=$1
factor=$2
method=$3

scaled_mapfile=$(dirname ${mapfile})/scaled_$(basename ${mapfile})

delta_x=$(3dinfo -di ${mapfile})
delta_y=$(3dinfo -dj ${mapfile})
delta_z=$(3dinfo -dk ${mapfile})
sdelta_x=$(echo "(($delta_x / ${factor}))"|bc -l)
sdelta_y=$(echo "(($delta_x / ${factor}))"|bc -l)
sdelta_z=$(echo "(($delta_z / ${factor}))"|bc -l)
3dresample -dxyz $sdelta_x $sdelta_y $sdelta_z \
           -rmode ${method} \
           -overwrite \
           -prefix ${scaled_mapfile} \
           -input ${mapfile}
