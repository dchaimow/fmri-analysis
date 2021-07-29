#!/bin/bash
#
# upsample_for-layer-sampling.sh <statsdir> <rim> [<additional> ...]
#
for filename in ${fBase}_${fBaseExt}_*_cope* ${fBase}_${fBaseExt}_*_tstat*
 do
     f_upsample $filename NN
 done    


