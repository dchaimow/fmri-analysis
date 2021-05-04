#!/bin/bash
#
# avgruns.sh <avgname> <run1_file> <run2_file> ...
#

avgName=$1
fileNames=${@:2}

3dMean -prefix ${avgName} ${fileNames}
