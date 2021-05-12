#!/bin/bash
#
# boldcorrect.sh <basename> <TR> <shiftFraction>
#
# - does bold correction of VASO by dividing nulled volumes by time shifted notnulled volumes
# - shiftFraction is shift relative to TR, positive value shifts back in time

fBaseName=$1
shiftFraction=$2

# shift notnulled (BOLD)
slicetimer -i ${fBaseName}_notnulled -o ${fBaseName}_notnulled_tshift.nii \
           --tglobal=-${shiftFraction}
# BOLD correction
fslmaths  ${fBaseName}_nulled -div ${fBaseName}_notnulled_tshift.nii \
          -max 0 -min 5 ${fBaseName}_vaso.nii
