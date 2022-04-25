# Registration and transformation


### fs_surface_to_func

### ciftify_surface_to_func


### surftransform_gii

### surftransform_fs

### register_fs_to_vasot1

### apply_ants_transforms

### import_fs_ribbon_to_func


# Data processing

### process_vaso

### bold_correct
performs bold correction on nulled on notnulled data, optionally corrects for temporal dealy between both using 
fsl slice timer.
external calls:
- fsl slicetimer (using nipype interface)
- fslmaths (using nipype interface)

# ROI functions

## Top level functions:

### get_stat_cluster_roi
computes an roi from clustered activation 
calls:
- sample_surf_func_stat
- smooth_surf
- cluster_surf
- get_fs_roi

### get_stat_cluster_atlas 
computes an "atlas" of activation clusters 
(alternative to get_stat_cluster_roi; implements get_fs_roi for all clusters)
calls:
- sample_surf_func_stat
- smooth_surf
- cluster_surf
- fs_overlay_to_fs_volume
- apply_ants_transforms

### get_funcloc_roi
computes an roi from functional localizer activation 
(feat; typically low, res fullbrain)
calls:
- reg_feat_to_fs
- sample_surf_feat_stat
- smooth_surf
- cluster_surf
- get_fs_roi

### get_md_roi
computes an roi from MD atlas in fsLR (HCP) space
calls:
- get_fs_LR_atlas_roi

### get_glasser_roi
computes an roi from HCP MMP 1.0 (Glasser) atlas in fsLR (HCP) space
calls:
- get_fs_LR_atlas_roi

### sample_layer_to_fs_LR
Samples volume to fs_LR surface via native surfaces using specified layer depths by first creating volumetric 
intermediate surfaces. Use [0,1] for entire cortical ribbon. (0 = wm surface, 1 = pial surface)
external calls:
- wb_command -surface-cortex-layer
calls:
- sample_surf_hcp   
- transform_data_native_surf_to_fs_LR

## 2nd and 3rd level function

## reg_feat_to_fs
Registers a feat directory to freesurfer data set using bbr.
external calls:
- bbregister

### sample_surf_feat_stat
samples stat file from feat directory to fs surface (assume is has been registered to fs and feat2fs.lta is present)
external calls:
- mri_vol2surf

## sample_surf_func_stat
samples stat file to surface using a number of intermediate surfaces starting from white surface
calls:
external calls:
- mri_vol2surf
- mris_calc

## smooth_surf
smoothes freesurfer surface overlay (mgh)
external calls:
- mri_surf2surf

## cluster_surf
clusters freesurfer surface overlay and produces annotation file and converts it to labels mgh file
external calls:
- mri_surfcluster

## sample_surf_hcp
Samples volume to surface using arbitrary GIFTI surfaces using hcp tools (wb_command). 
Generates midthickness if file does not exists (and returns file name).
external calls:
- wb_command -volume-to-surface-mapping

## transform_data_native_surf_to_fs_LR
Transforms data from native space to fsLR space
external calls:
- wb_command -metric-resample



### get_fs_roi
Takes a fs roi defined as a surface overlay and transforms it to the functional volume.
calls:
- fs_overlay_to_fs_volume
- apply_ants_transforms
- index_roi

### fs_overlay_to_fs_volume
Tranforms freesurfer overlay file (where from?) to freesurfer volume
external calls:
- mri_surf2vol

### get_fs_LR_atlas_roi
Returns an ROI in functional space by transforming GIFTI label files in fs_LR space. 
calls:
- fs_LR_label_to_fs_volume
- apply_ants_transforms
- index_roi

### fs_LR_label_to_fs_volume
Transforms label .gii file in fs_LR space to Freesurfer volume.
external calls:
- wb_command -label-to-volume-mapping

## ROI from depth filling of activation maps

### get_funcact_roi_laynii
Depth filled ROI from activation using LN2_MASK (laynii).
external calls:
- LN2_COLUMNS
- LN2_MASK
- fslmaths

### get_funcact_roi_vfs
Depth filled ROI from activation using VDFS columns.

## ROI utils

### index_roi
Extracts ROI with a specific index from a multi-index label file.

### roi_and
Logical and combination of rois.

### roi_or
Logical or combination of rois.

### mask_image
(Check if it works on 4D and wheter it can be used to ignore affine on rois)
Mask data.

# Trial averaging

### calc_stim_times
Calculates stim onset times from a trial condition order list.

### write_stim_time_files
Writes files with stim times for multiple runs to be used by 3ddeconvolve.

### average_trials_3ddeconvolve
Performs trial averaging, potentially on multiple files (runs) using 3ddeconvolve.
Returns trial averages, baseline and fstat. Baseline is averaged from multiple runs.
calls:
- write_stim_time_files
external calls:
- 3ddeconvolve (using nipype interface)
- 3dTcat (using nipype interface)
- 3dTStat (using nipype interface)

### calc_percent_change_trialavg
Calculates percent change from trial average and baseline files.
external calls:
- 3dcalc (using nipype interface)


# sampling functions (layers and timecourses)

### calc_layers_laynii
calculates layers from gray mitter ribbon (run file) using laynii
external calls:
- LN2_LAYERS

### generate_two_layers
Generates superficial and deeper layer rois from cortical depths, possibly applying an roi and 
leaving an optional gap in depth between both

### average_roi
Averages data sampled from roi (assumes data and roi array correspond to each other, ignores affines).
Returns mean, std and number of timecourses.
calls:
- sample_roi

### sample_layer_profile
Samples data from roi with depths and then bins it into n layers, returns layer averaged data and 
depths of layer mid-points
calls:
- sample_depths

### sample_roi
samples data from roi using nilearn apply_mask (assumes data and roi array correspond to each other, ignores affines)

### sample_depths
Samples data and corresponding depths from roi using nilearn apply_mask (assumes data, depths and roi array correspond 
to each other, ignores affines).

### upsample
Uses 3dupsample to upsample a nifti file by some factor. Potentially to be used for laynii depth calculation.
Currently not used.

### get_labels
WIP to sample data according to labels (atlas rois) and save in a data frame.

### get_labels
WIP to sample data according to labels (atlas rois) and save in a data frame.

### get_labels_layers
WIP(still empty) to sample data according to labels (atlas rois) and layers and save in a data frame.

### get_labels_layers_masked
WIP(still empty) to sample data according to labels (atlas rois) and layers and using a mask (e.g activation)
and save in a data frame.

# Plotting time courses and layer profiles

### plot_cond_tcrs
(check if it is being used)
Plots time course for multiple conditions. 
All timecourses should have the same length.


### plot_profiles
Plots layer profiles computed from a list of images (one profile per image) using an roi and cortical depths image.


### plot_on_mmhcp_surface
WIP to generate glass brain visualization of values (provided as 1D vector) according to MMP-HCP (Glasser) parcellation.

# Util functions

### fsl_remove_ext

### add_prefix_to_nifti_basename

### add_postfix_to_nifti_basename

### reset_affine
sets affine (form and sform) of niimg to identity transform


# Finn replication (pilot) specific functions

plot_finn_panel

plot_finn_tcrses

finn_trial_averaging

finn_trial_averaging_with_boldcorrect

finn_trial_averaging_on_renzo_boldcorrect

paradigm

# Other

### feat_analysis

# potentially todo (currently empty)

### preprocess_funcloc

### fs_surf_to_fs_volume

### get_mni_coord_roi

### layer_extend_roi_laynii

### layer_extend_roi_vfs

### plot_timecourses

### normalize
Normalizes single voxel timecourses to change relative to a baselne

### avg_timcourse_nearest_volume

### get_funcact_roi_other_versions

### initialize_session

### motion_correction

### trial_averaging

### glm_analysis
