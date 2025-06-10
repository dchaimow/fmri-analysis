#!/usr/bin/env python3
"""
Group analysis on fs_LR surface afer sampling volume data from individual subjects.

Steps are:
1. sample individual subject data to surfaces in native functional space. By default sample
    from entire cortical depth, optionally sample from one more depth range (layers).
    -> the result are two or more metric files per subject, one for each layer for each hemisphere
2. transform individual subject surface data to fs_LR space.
    -> the result are transformed metric files
3. generate one cifti file per subject, containing the two hemispheres
    -> the result is one or more cifti file per subject, containing the two hemispheres
4. process individual subject data in fs_LR space, e.g. smooth, possibly calculate layer contrasts
    -> the results is one processed cifti file per subject
5. calculate group contrasts in fs_LR space, e.g. average across subjects, calculate t-stats, etc.
6. optionally extract metric files from the cifti files, e.g. for further processing

Alternatively, we could generate the cifti after processing the individual subject data (step 4), or even
stick to metric files and not generate cifti files at all. 


Needed data:
- single subject volume files to sample (1st level estimates, fstats ...)
- single subject surfaces in functional space (pial and white)
- single subject transform to fs_LR space (ciftify)

Needed specifications:
- list of subjects to process
- single subjects laminar contrasts to compute
- group analysis contrasts to compute


If we want to calculated a depth contrast we need to specify the depth ranges of one or more layers to sample,
as well as the contrast(s) to compute from these layers.
    }
"""

import os
from joblib import Parallel, delayed
import layer_analysis as analysis
import numpy as np
from tempfile import TemporaryDirectory
import subprocess
import nibabel as nib

def sample_firstlevel_layer_contrast_to_fsLR(subject,firstlevel_analysis_dir,ciftify_dir,surf_dir,
                                             firstlevel_subpath,
                                             depth_ranges=None,
                                             layer_contrasts=None,
                                             group_contrasts=None,
                                             group_analysis_dir=None,
                                             smooth_sigma=None):
    
    with TemporaryDirectory() as tmpdirname:
        ciftify_dir = os.path.join(ciftify_dir, subject)
        surf_dir = os.path.join(surf_dir, subject)
        firstlevel_path = os.path.join(firstlevel_analysis_dir, subject, firstlevel_subpath)

        for hemi in ('L','R'):
            # set surface file names
            white_surf = os.path.join(surf_dir,f"{hemi}.white.func.surf.gii")
            pial_surf = os.path.join(surf_dir,f"{hemi}.pial.func.surf.gii")
            mid_surf_fsLR = os.path.join(ciftify_dir,"MNINonLinear",
                                            f"{subject}.{hemi}.midthickness.164k_fs_LR.surf.gii")

            firstlevel_fsLR_layer_data = np.zeros((len(depth_ranges), 163842))  # assuming 164k fs_LR vertices
            for layer_idx, depth_range in enumerate(depth_ranges):
    
                # set output file names
                firstlevel_fsLR_path = os.path.join(tmpdirname,f'{hemi}.{layer_idx}.164k_fs_LR.func.gii')                
                analysis.sample_layer_to_fs_LR(firstlevel_path, 
                                            firstlevel_fsLR_path,
                                            white_surf, pial_surf,
                                            ciftify_dir, 
                                            hemi,
                                            depth_range)
                
                # load sampled fsLR surface data
                firstlevel_fsLR_layer_data[layer_idx,:] =  nib.load(firstlevel_fsLR_path).darrays[0].data

            # calculate layer contrasts
            firstlevel_fsLR_layer_contrast_data = layer_contrasts @ firstlevel_fsLR_layer_data




                



def group_fslr_analysis(firstlevel_analysis_dir,ciftify_dir,surf_dir,firstlevel_subpath,subjects,
                        depth_ranges=None,layer_contrasts=None,group_contrasts=None,group_analysis_dir=None,
                        smooth_sigma=None):
    """
    Perform group analysis on fs_LR surface afer sampling volume data from individual subjects.

    Parameters
    ----------
    first_level_analysis_dir: str
        Path to the directory containing individual subjects' first level analysis subdirectories with first level analysis results.
    ciftify_dir: str
        Path to the directory containing individual subjects' ciftify subdirectories.
    firstlevel_subpath: 
        Relative path from each subjects first level analysis directory to the file to sample from.
    depth_ranges: list of tuples (default=None)
        List of depth ranges to sample from, e.g. [(0, 0.2), (0.2, 0.4)].
        If None, the entire cortical depth is sampled.
    layer_contrasts: ndarray
        List of contrasts to compute from the sampled layers, e.g. [[1, -1], [0, 1]].
        If None, each layer is processed separately.
    - group_contrasts: List of group contrasts to compute.
    - group_analysis_dir: Directory to store group analysis results.
    - smooth_sigma: Sigma value for smoothing the data.

    Returns:
    None
    """
    if depth_ranges is None:
        # Default to entire cortical depth
        depth_ranges = [(0, 1)]  

    if layer_contrasts is None:
        # each layer is processed separately
        layer_contrasts = np.eye(depth_ranges)

    # step 1 - sample individual subject data to surfaces in native functional space

    subject = 'sub-01'  # Example subject, replace with actual subject loop
    sample_firstlevel_layer_contrast_to_fsLR(subject,firstlevel_analysis_dir,ciftify_dir,surf_dir,
                                             firstlevel_subpath,
                                             depth_ranges,
                                             layer_contrasts)


study_data_dir = '/Users/dchaimow/data/finn-et-al-2019_replication/'
firstlevel_analysis_dir = os.path.join(study_data_dir,'derivatives','analysis')
ciftify_dir = os.path.join(study_data_dir,'derivatives','ciftify')
surf_dir = os.path.join(study_data_dir,'derivatives','ref_anat')
firstlevel_data_subpath = 'trialavg5/trialavg5_notnulled_alpharem_fstat.nii'
subjects = ['sub-01', 'sub-02']
depth_ranges = [(0,0.5),(0.5,1)]
layer_contrasts = np.array([[0.5, 0.5],[-1, 1]])
group_contrasts = np.array([1, 1, 1])/2

group_fslr_analysis(
    firstlevel_analysis_dir,
    ciftify_dir,
    surf_dir,
    firstlevel_data_subpath,
    subjects,
    depth_ranges,
    layer_contrasts,
    group_contrasts)