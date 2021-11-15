#from nipype.interfaces import spm
from nipype.interfaces import fsl
from nipype.interfaces import matlab
from nipype.interfaces import cat12
import nipype.pipeline.engine as pe
import nibabel as nib
import numpy as np
import os

spm_path = '/Users/denis/Software/spm12'
matlab.MatlabCommand.set_default_matlab_cmd(
    "/Applications/MATLAB_R2020b.app/bin/matlab -nodesktop -nosplash")

def load_niimg(niimg):
    if type(niimg) is str:
        return nib.load(niimg)
    else:
        return niimg

def normalize(niimg_in,out_file=None):
    niimg_in = load_niimg(niimg_in)
    data = niimg_in.get_fdata()
    data_norm = (data-np.min(data))/(np.max(data)-np.min(data))
    niimg_out = nib.Nifti1Image(data_norm,niimg_in.affine,niimg_in.header)
    if out_file:
        nib.save(niimg_out,out_file)
    return niimg_out

def multiply(niimg_in1, niimg_in2, out_file=None):
    niimg_in1 = load_niimg(niimg_in1)
    niimg_in2 = load_niimg(niimg_in2)
    data1 = niimg_in1.get_fdata()
    data2 = niimg_in2.get_fdata()
    data_mult = data1 * data2
    niimg_out = nib.Nifti1Image(data_mult,niimg_in1.affine,niimg_in1.header)
    if out_file:
        nib.save(niimg_out,out_file)
    return niimg_out

def mprageize(inv2_file, uni_file):
    """ 
    Based on Sri Kashyap (https://github.com/srikash/presurfer/blob/main/func/presurf_MPRAGEise.m)
    """
    #mprageize_wf = pe.Workflow(name='mprageize')
    
    # bias correct INV2
    seg = spm.NewSegment()
    seg.inputs.channel_files = inv2_file
    seg.inputs.channel_info = (0.001, 30, (False, True))
    tissue1 = ((os.path.join(spm_path,'tpm','TPM.nii'), 1), 2, (False,False), (False, False))
    tissue2 = ((os.path.join(spm_path,'tpm','TPM.nii'), 2), 2, (False,False), (False, False))
    tissue3 = ((os.path.join(spm_path,'tpm','TPM.nii'), 3), 2, (False,False), (False, False))
    tissue4 = ((os.path.join(spm_path,'tpm','TPM.nii'), 4), 3, (False,False), (False, False))
    tissue5 = ((os.path.join(spm_path,'tpm','TPM.nii'), 5), 4, (False,False), (False, False))
    tissue6 = ((os.path.join(spm_path,'tpm','TPM.nii'), 6), 2, (False,False), (False, False))
    seg.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5, tissue6]    
    seg.inputs.affine_regularization = 'mni'
    seg.inputs.sampling_distance = 3
    seg.inputs.warping_regularization = [0, 0.001, 0.5, 0.05, 0.02]
    seg.inputs.write_deformation_fields = [False, False]
    seg_results = seg.run()
    
    # normalize bias corrected INV2
    norm_inv2_niimg = normalize(seg_results.outputs.bias_corrected_images)
                              
    # multiply normalized bias corrected INV2 with UNI
    uni_mprageized_niimg = multiply(norm_inv2_niimg,uni_file)
    return uni_mprageized_niimg

def cat12_seg(in_file):
    cat = cat12.CAT12Segment(in_files=in_file)
    cat_results = cat.run()
    return cat_results

def my_fs_recon_all(brainmask):
    pass


