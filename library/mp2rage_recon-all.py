#! /usr/bin/env python3
import anatomy
import sys

if __name__ == "__main__":
    if len(sys.argv)==3:
        anatomy.mp2rage_recon_all(sys.argv[1],sys.argv[2])
    elif len(sys.argv)==4:
        anatomy.set_spm_path(sys.argv[3])
        anatomy.mp2rage_recon_all(sys.argv[1],sys.argv[2])
    else:
        print('Usage: mp2rage_recon-all.py inv2_file uni_file [spm_path]\n' +
              '(requires freesurfer, intended for version 7.1)')
        
