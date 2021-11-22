#! /usr/bin/env python3
import anatomy
import sys

if __name__ == "__main__":
    if len(sys.argv)==3:
        anatomy.mp2rage_recon_all(sys.argv[1],sys.argv[2])
    elif len(sys.argv)==4:
        anatomy.mp2rage_recon_all(sys.argv[1],sys.argv[2],output_fs_dir=sys.argv[3])
    elif len(sys.argv)==5:
        anatomy.set_spm_path(sys.argv[4])
        anatomy.mp2rage_recon_all(sys.argv[1],sys.argv[2],output_fs_dir=sys.argv[3])
    else:
        print('Usage: mp2rage_recon-all.py inv2_file uni_file [output_fs_dir] [spm_path]\n' +
              '(requires freesurfer, intended for version 7.1)')
        
