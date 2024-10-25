#!/usr/bin/env python3

"""
command line tool to convert voxel coordinates to world coordinates

usage: voxel-to-world.py <nifti_file> <x> <y> <z>
"""

import sys
import nibabel as nib

if len(sys.argv) != 5:
    print("Usage: voxel-to-world.py <nifti_file> <x> <y> <z>")
    sys.exit(1)

nifti_file = sys.argv[1]
x = int(sys.argv[2])
y = int(sys.argv[3])
z = int(sys.argv[4])

affine = nib.load(nifti_file).affine
world = affine.dot([x, y, z, 1])

print(world[0], world[1], world[2])