#! /usr/bin/env python3
import sys
import os
from nipype.interfaces import spm,matlab

spm_path = '/data/pt_02389/Software/spm12'

def spm_bias_correction():
    # bias correct INV2
    seg = spm.NewSegment()
    seg.inputs.channel_info = (0.001, 20, (False, True))
    seg.inputs.affine_regularization = 'mni'
    seg.inputs.sampling_distance = 3
    seg.inputs.warping_regularization = [0, 0.001, 0.5, 0.05, 0.2]
    seg.inputs.write_deformation_fields = [False, False]
    return seg
                          
if __name__ == "__main__":
    if len(sys.argv) == 3:
        matlab.MatlabCommand.set_default_paths(sys.argv[2])
    else:
        matlab.MatlabCommand.set_default_paths(spm_path)
        
    if len(sys.argv) in [2,3]:
        infile = sys.argv[1]
        bias_correct = spm_bias_correction()        
        bias_correct.inputs.channel_files = infile
        bias_correct.run(cwd=os.path.dirname(os.path.abspath(infile)))
    else:
        print('Usage: spm_bias-correct.py infile [spm_path]')
    
    
