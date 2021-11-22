from nipype.interfaces import spm
from nipype.interfaces import matlab
from nipype.interfaces import cat12
from nipype.interfaces.freesurfer import ApplyVolTransform, ApplyMask
import nipype.pipeline.engine as pe
import nibabel as nib
import numpy as np
import os
import shutil

spm_path = '/data/pt_02389/Software/spm12'
matlab.MatlabCommand.set_default_paths(spm_path)

def set_spm_path(new_spm_path):
    global spm_path
    spm_path = new_spm_path
    matlab.MatlabCommand.set_default_paths(spm_path)

def check_spm_path():
    print(spm_path)
    
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

def mprageize(inv2_file, uni_file, out_file=None):
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
    seg_results = seg.run(cwd=os.path.dirname(os.path.abspath(inv2_file)))
    
    # normalize bias corrected INV2
    norm_inv2_niimg = normalize(seg_results.outputs.bias_corrected_images)
                              
    # multiply normalized bias corrected INV2 with UNI
    uni_mprageized_niimg = multiply(norm_inv2_niimg,uni_file)

    if out_file:
        nib.save(uni_mprageized_niimg,out_file)
    return uni_mprageized_niimg

def cat12_seg(in_file):
    # currently runs external MATLAB script, TODO: change to nipype interface

    module_file_path = os.path.abspath(__file__)
    mfile_path = os.path.dirname(module_file_path)
    
    os.system(f"matlab" +
              f" -nodisplay -nodesktop -r " +
              f"\"addpath('{mfile_path}');" +
              f"cat12_seg('{in_file}','{spm_path}'); exit;\"")

    cwd = os.path.dirname(os.path.abspath(in_file))
    in_file_basename = os.path.basename(os.path.abspath(in_file))
    if os.path.isfile(os.path.join(cwd,'mri','p1' + in_file_basename)):
        return os.path.join(cwd,'mri','p1' + in_file_basename), os.path.join(cwd,'mri','p2' + in_file_basename)
    else:
        return os.path.join(cwd,'p1' + in_file_basename), os.path.join(cwd,'p2' + in_file_basename)
              
def mp2rage_recon_all(inv2_file,uni_file,output_fs_dir=None):

    # mprageize
    cwd = os.path.dirname(os.path.abspath(uni_file))
    uni_basename = os.path.basename(os.path.abspath(uni_file))
    
    uni_mprageized_file = os.path.join(cwd,uni_basename.replace('UNIT1','T1w'))
    uni_mprageized_brain_file =  os.path.join(cwd,uni_basename.replace('UNIT1','T1w_brain'))
    brainmask_file =  os.path.join(cwd,uni_basename.replace('UNIT1','brainmask'))
    if os.path.basename(uni_mprageized_file)==uni_basename:
        uni_mprageized_file = os.path.join(cwd,'T1w.nii')
        uni_mprageized_brain_file = os.path.join(cwd,'T1w_brain.nii')
        brainmask_file =  os.path.join(cwd,'brainmask.nii')
    
    mprageize(inv2_file,uni_file,uni_mprageized_file)

    # obtain brainmask using cat12
    wm_file, gm_file = cat12_seg(uni_mprageized_file)

    wm_nii = nib.load(wm_file)
    gm_nii = nib.load(gm_file)
    uni_mprageized_nii = nib.load(uni_mprageized_file)

    wm_data = wm_nii.get_fdata()
    gm_data = gm_nii.get_fdata()
    uni_mprageized_data = uni_mprageized_nii.get_fdata()

    brainmask_data = np.array(((wm_data > 0) | (gm_data > 0)),dtype=int)
    brainmask_nii = nib.Nifti1Image(brainmask_data,
                                    uni_mprageized_nii.affine,
                                    uni_mprageized_nii.header)
    nib.save(brainmask_nii,brainmask_file)
    
    uni_mprageized_brain_data = brainmask_data * uni_mprageized_data
    uni_mprageized_brain_nii = nib.Nifti1Image(uni_mprageized_brain_data,
                                               uni_mprageized_nii.affine,
                                               uni_mprageized_nii.header)
    nib.save(uni_mprageized_brain_nii,uni_mprageized_brain_file)

    # run recon-all
    if output_fs_dir:
         fs_dir = os.path.dirname(os.path.abspath(in_file))
         sub = os.path.basename(os.path.abspath(in_file))
    else:
         fs_dir = cwd
         sub = 'freesurfer'
        
    # autorecon1 without skullstrip removal (~11 mins)
    os.system("recon-all" + \
          " -i " + uni_mprageized_file + \
          " -hires" + \
          " -autorecon1" + \
          " -noskullstrip" + \
          " -sd " + fs_dir + \
          " -s " + sub + \
          " -parallel")

    # apply brain mask from CAT12
    transmask = ApplyVolTransform()
    transmask.inputs.source_file = brainmask_file
    transmask.inputs.target_file = os.path.join(fs_dir, sub, 'mri', 'orig.mgz')
    transmask.inputs.reg_header = True
    transmask.inputs.interp = "nearest"
    transmask.inputs.transformed_file = os.path.join(fs_dir, sub, 'mri', 'brainmask_mask.mgz')
    transmask.inputs.args = "--no-save-reg"
    transmask.run(cwd=cwd)

    applymask = ApplyMask()
    applymask.inputs.in_file = os.path.join(fs_dir, sub,'mri','T1.mgz')
    applymask.inputs.mask_file = os.path.join(fs_dir, sub, 'mri', 'brainmask_mask.mgz')
    applymask.inputs.out_file =  os.path.join(fs_dir, sub, 'mri', 'brainmask.mgz')
    applymask.run(cwd=cwd)

    shutil.copy2(os.path.join(fs_dir, sub, 'mri', 'brainmask.mgz'),
                 os.path.join(fs_dir, sub, 'mri','brainmask.auto.mgz'))

    # continue recon-all
    with open(os.path.join(cwd,'expert.opts'), 'w') as text_file:
        text_file.write('mris_inflate -n 100\n')
    # autorecon2 and 3
    os.system("recon-all" + \
              " -hires" + \
              " -autorecon2" + " -autorecon3"\
              " -sd " + fs_dir + \
              " -s " + sub + \
              " -expert " + os.path.join(cwd,'expert.opts') + \
              " -xopts-overwrite" + \
              " -parallel")
