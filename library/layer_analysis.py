import subprocess
import os
import numpy as np
import nibabel as nib
from nipype.interfaces.afni import Deconvolve, TCatSubBrick, Refit, Calc, TStat
from nipype.interfaces.freesurfer import MRIConvert, MRIsConvert
from nipype.interfaces.fsl import SliceTimer
from nilearn.image import math_img
from nilearn.masking import apply_mask, intersect_masks
from collections import defaultdict
import nipype.interfaces.fsl as fsl
import nipype.pipeline.engine as pe
from niworkflows.interfaces.surf import CSVToGifti, GiftiToCSV
from nipype.interfaces.ants import ApplyTransformsToPoints
import pandas as pd
import seaborn as sns
from nilearn._utils import check_niimg
import copy
from itertools import zip_longest
import matplotlib.pyplot as plt
from shutil import move

def fsl_remove_ext(filename):
    result = subprocess.run(['remove_ext',filename],stdout=subprocess.PIPE)
    return result.stdout.strip().decode()


def surftransform_gii(gii_surf, transforms, invert_transform_flags,cwd=None):
    """Takes gifti surface and applies ants transforms
    """
    if cwd==None:
        cwd=os.path.dirname(os.path.normpath(gii_surf))
    out_file=os.path.basename(os.path.normpath(gii_surf))
    # convert gii to csv
    result_GiftiToCSV = GiftiToCSV(in_file=gii_surf,
                                   itk_lps=True).run(cwd=cwd)
    csv_surf = result_GiftiToCSV.outputs.out_file
    # apply transform
    result_ApplyTransformsToPoints = ApplyTransformsToPoints(
        dimension=3,
        input_file=csv_surf,
        transforms=transforms,
        invert_transform_flags=invert_transform_flags
    ).run(cwd=cwd)
    csv_surf_transformed = result_ApplyTransformsToPoints.outputs.output_file
    # convert csv to gii
    result_CSVToGifti = CSVToGifti(in_file=csv_surf_transformed,
                                   gii_file=gii_surf,
                                   itk_lps=True).run(cwd=cwd)
    gii_surf_transformed = result_CSVToGifti.outputs.out_file
    return gii_surf_transformed

def surftransform_fs(fs_surf, transforms, invert_transform_flags,out_file,cwd=None):
    if cwd==None:
        cwd=os.path.dirname(os.path.normpath(out_file))
    # convert fs to gii
    result_MRIsConvert = MRIsConvert(in_file = fs_surf,
                                     out_datatype = 'gii',
                                     to_scanner=True
                                     ).run(cwd=cwd)
    gii_surf = os.path.join(cwd,result_MRIsConvert.outputs.converted)
    gii_surf_transformed = surftransform_gii(gii_surf, transforms, invert_transform_flags)
    # convert gii to fs
    MRIsConvert(in_file=gii_surf_transformed,
                out_file=out_file,
                to_tkr=True
                ).run(cwd=cwd)
    return out_file

def fs_surface_to_func(fs_to_func_reg,fs_dir,analysis_dir=None):
    if analysis_dir == None:
        analysis_dir = os.path.join(fs_dir,'surf')
    transform_0_lin = fs_to_func_reg[2]
    transform_1_inversewarp = fs_to_func_reg[3]
    invert_transform_flags = [True, False]
    surf_trans_files = dict()
    for hemi in ['lh','rh']:
        for surf_type in ['white','pial']:        
            surf = os.path.join(fs_dir,'surf', hemi + '.' + surf_type)
            surf_trans = os.path.join(analysis_dir,hemi + '.' +surf_type + '_func')
            surf_trans_files[hemi,surf_type] = \
                surftransform_fs(surf,[transform_0_lin,transform_1_inversewarp],
                                 invert_transform_flags,out_file=surf_trans)
    return surf_trans_files

def ciftify_surface_to_func(fs_to_func_reg,ciftify_dir,analysis_dir=None):
    if analysis_dir == None:
        analysis_dir = os.path.join(ciftify_dir,'T1w','fsaverage_LR32k')
    ciftify_subject = os.path.basename(os.path.normpath(ciftify_dir))
    transform_0_lin = fs_to_func_reg[2]
    transform_1_inversewarp = fs_to_func_reg[3]
    invert_transform_flags = [True, False]
    surf_trans_files = dict()
    for hemi in ['L','R']:
        for surf_type in ['white','pial']:        
            surf = os.path.join(ciftify_dir,'T1w','fsaverage_LR32k',
                                ciftify_subject + '.' + hemi + '.' +
                                surf_type + '.32k_fs_LR.surf.gii')
            surf_trans = os.path.join(analysis_dir,
                                      ciftify_subject + '.' + hemi + '.' +
                                      surf_type + '.32k_fs_LR_func.surf.gii')
            out_file = surftransform_gii(surf,[transform_0_lin,transform_1_inversewarp],
                                  invert_transform_flags,cwd=analysis_dir)
            os.rename(out_file,surf_trans)
            surf_trans_files[hemi,surf_type] = surf_trans            
    return surf_trans_files
            

def process_vaso(session_dir,process_script,alpharem_runs,gonogo_runs,analysis_subdir='analysis'):
    analysis_dir = os.path.join(session_dir,analysis_subdir)
    if not os.path.isdir(os.path.join(analysis_dir)):
        os.mkdir(analysis_dir)
        with open(os.path.join(analysis_dir,'func_alpha-rem_task-runs.txt'),'w') as file:
            print(*alpharem_runs,sep='\n',file=file)
        with open(os.path.join(analysis_dir,'func_go-nogo_task-runs.txt'),'w') as file:
            print(*gonogo_runs,sep='\n',file=file)
        subprocess.run([process_script,session_dir,analysis_subdir])
    return analysis_dir


def register_fs_to_vasot1(fs_dir,analysis_dir, force=False):
    if not os.path.isfile(os.path.join(analysis_dir,'fs_to_func_0GenericAffine.mat')) \
       or force==True:
        subprocess.run(['register_fs-to-vasoT1.sh',
                        'func_all_T1.nii',
                        fs_dir,
                        'itksnap'],cwd=analysis_dir)

def apply_ants_transforms(vol_in, vol_out, ref_vol, affine, warp):
    """Applies ANTS non-linear registration transform, consisting of 1. warp and 2. affine to input volume.
    """
    subprocess.run(['antsApplyTransforms',
                    '--interpolation','BSpline''[5]''',
                    '-d','3',
                    '-i',vol_in,
                    '-r',ref_vol,
                    '-t',warp,
                    '-t',affine,
                    '-o',vol_out,
                    '-n','NearestNeighbor'])

def import_fs_ribbon_to_func(fs_dir,analysis_dir):
    """Calls shell script in fmri-analysis/library
    TODO: port to python function
    assume fs_to_func_reg files in analysis_dir
    """
    if subprocess.run(['/data/p_02389/code/fmri-analysis/library/import-fs-ribbon.sh',
                      fs_dir,
                      analysis_dir,
                      os.path.join(analysis_dir,'fs_t1_in-func.nii')]).returncode == 0:
        return os.path.join(analysis_dir,'rim.nii')
    

def index_roi(roi,idx):
    """Extracts ROI with a specific index from a multi-index label file.
    """
    return math_img(f'img=={idx}',img=roi)

def fs_LR_label_to_fs_volume(ciftify_dir,analysis_dir,labels,hemi,out_basename):
    """Transforms label .gii file in fs_LR space to Freesurfer volume.
    """
    ciftify_subject=os.path.basename(os.path.normpath(ciftify_dir))
    mid_surf=os.path.join(ciftify_dir,'T1w','fsaverage_LR32k',
                          ciftify_subject+'.'+hemi+'.midthickness.32k_fs_LR.surf.gii')
    white_surf=os.path.join(ciftify_dir,'T1w','fsaverage_LR32k',
                            ciftify_subject+'.'+hemi+'.white.32k_fs_LR.surf.gii')
    pial_surf=os.path.join(ciftify_dir,'T1w','fsaverage_LR32k',
                           ciftify_subject+'.'+hemi+'.pial.32k_fs_LR.surf.gii')
    volume=os.path.join(ciftify_dir,'T1w','T1w.nii.gz')
    volume_out=os.path.join(analysis_dir,out_basename+'_labels_'+hemi+'.nii')
    subprocess.run(['wb_command',
                    '-label-to-volume-mapping',
                    labels,
                    mid_surf,
                    volume,
                    volume_out,
                    '-ribbon-constrained',
                    white_surf,
                    pial_surf])
    return volume_out
    
def get_fs_LR_atlas_roi(parcel=None,atlas_labels=None,out_basename=None,analysis_dir=None,
                        ciftify_dir=None,fs_to_func_reg=None,force_import=False):
    """Returns an ROI in functional space by transforming GIFTI label files in fs_LR space. 
    ROI is specified using parcel=(hemi,idx).
    """
    hemi = parcel[0]
    parcel_idx = parcel[1]
    labels_in_func = os.path.join(analysis_dir,out_basename+'_labels_'+hemi+'_in-func.nii')
    if not os.path.isfile(labels_in_func) or force_import==True:
        labels_in_fs_individual = fs_LR_label_to_fs_volume(ciftify_dir,analysis_dir,
                                                           atlas_labels[hemi],hemi,out_basename)        
        apply_ants_transforms(vol_in=labels_in_fs_individual,
                              vol_out=labels_in_func,
                              ref_vol=fs_to_func_reg[0],
                              affine=fs_to_func_reg[1],
                              warp=fs_to_func_reg[2])
    roi = index_roi(labels_in_func,parcel_idx)
    return roi

def get_md_roi(parcel=None,analysis_dir=None,ciftify_dir=None,fs_to_func_reg=None,
               md_labels=None,force_import=False):
    """Returns a multiple-demand network ROI (Moataz et al. 2020) transformed to functional space
    """
    if md_labels==None:
        md_labels={'L':'/data/pt_02389/RL_analysis/ROIs/HCP_Glasser/Moataz/MD_L_0.2thresh.label.gii',
                   'R':'/data/pt_02389/RL_analysis/ROIs/HCP_Glasser/Moataz/MD_R_0.2thresh.label.gii'}
    out_basename='md'
    roi = get_fs_LR_atlas_roi(parcel,md_labels,out_basename,analysis_dir,ciftify_dir,
                              fs_to_func_reg,force_import)
    return roi
    
def get_glasser_roi(parcel=None,analysis_dir=None,ciftify_dir=None,fs_to_func_reg=None,
                    glasser_labels=None,force_import=False):
    """Returns a HCP MMP 1.0 atlas ROI (Glasser et al. 2016) transformed to functional space
    """
    if glasser_labels==None:
        glasser_labels={'L':'/data/pt_02389/FinnReplicationPilot/ROIs/GlasserAtlas.L.32k_fs_LR.label.gii',
                        'R':'/data/pt_02389/FinnReplicationPilot/ROIs/GlasserAtlas.R.32k_fs_LR.label.gii'}
    out_basename='glasser'
    roi = get_fs_LR_atlas_roi(parcel,glasser_labels,out_basename,analysis_dir,ciftify_dir,
                              fs_to_func_reg,force_import)
    return roi


def calc_stim_times(onset_delay, trial_duration, trial_order, condition_names=None):
    n = len(trial_order)
    t = np.arange(0,n) * trial_duration + onset_delay
    stim_times = dict()
    for condition in set(trial_order):
        stim_times[condition] = t[np.array(trial_order)==condition]
    return stim_times

def write_stim_time_files(stim_times_runs,cwd=None):
    if cwd==None:
        cwd=os.path.curdir()
    # find all conditions from all runs
    conditions = set.union(*[set(stim_times.keys()) for stim_times in stim_times_runs])
    # for each condition create a file and write a line of stim times for each run
    condition_stim_files=[]
    for condition in conditions:
        condition_stim_files.append([condition,
                                     os.path.join(cwd,'stim-times_' + str(condition) + '.txt')])
        with open(os.path.join(cwd,condition_stim_files[-1][1]),'w') as file:
            for stim_times in stim_times_runs:
                stim_times = defaultdict(list,stim_times)
                print(*stim_times[condition],file=file)
    return condition_stim_files

def average_trials_3ddeconvolve(in_files,stim_times_runs,trial_duration,
                                out_files_basename,polort=5,onset_shift=0,cwd=None):
    if cwd==None:
        cwd=os.path.dirname(os.path.normpath(in_files[0]))
    n_files = len(in_files)
    # returns estimated impules response components and baseline
    # set parameters
    tr = nib.load(in_files[0]).header.get_zooms()[3]
    n = np.ceil(trial_duration/tr)
    a = 0
    b = tr * (n - 1)
    # prepare stim times and model
    condition_stim_files = write_stim_time_files(stim_times_runs,cwd)
    n_conditions = len(condition_stim_files)
    stim_times=[]
    stim_label=[]
    i_condition = 0
    for condition, stim_file in condition_stim_files:
        i_condition = i_condition + 1
        stim_times.append((i_condition,stim_file,f'TENT({a},{b},{n})'))
        stim_label.append((i_condition,str(condition)))    

    #  cmdline_list = ['3dDeconvolve',
    #                 '-input', ' '.join(in_files),
    #                 '-overwrite',
    #                 '-cbucket', os.path.join(cwd,out_files_basename + '_cbucket.nii'),
    #                 '-fout',
    #                 '-local_times',
    #                 '-bucket',  os.path.join(cwd,out_files_basename + '_Deconvolve.nii'),
    #                 '-polort', str(polort),
    #                 '-stim_times_subtract', str(onset_shift),
    #                 '-num_stimts',str(n_conditions)]

    # for i in range(n_conditions):
    #     cmdline_list.extend(['-stim_times',str(stim_times[i][0]),stim_times[i][1],f"'{stim_times[i][2]}'"])
        
    # for i in range(n_conditions):
    #     cmdline_list.extend(['-stim_label',str(stim_label[i][0]),stim_label[i][1]])
    
    # # run deconvolve
    # subprocess.run(cmdline_list)b

    deconvolve = Deconvolve()
    deconvolve.inputs.in_files = in_files
    deconvolve.inputs.stim_times = stim_times
    deconvolve.inputs.stim_label = stim_label
    deconvolve.inputs.polort = polort
    deconvolve.inputs.local_times = True
    deconvolve.inputs.fout = True
    deconvolve.inputs.cbucket = os.path.join(cwd,out_files_basename + '_cbucket.nii.gz')
#    deconvolve.inputs.out_file = os.path.join(cwd,out_files_basename + '_Deconvolve.nii.gz')
    deconvolve.inputs.args ='-overwrite'
    deconvolve.inputs.stim_times_subtract = onset_shift
    #print(deconvolve.cmdline)
    result = deconvolve.run(cwd=cwd)
    #return
    # extract fstat
    result_fstat = TCatSubBrick(in_files=[(result.outputs.out_file,f"'[0]'")],
                                out_file = os.path.join(cwd,out_files_basename + '_fstat.nii'),
                                args='-overwrite').run()    
    # add back baseline  
    baseline_idcs = 0+(polort+1)*np.arange(n_files)
    baseline_idcs_str = ','.join([str(i) for i in baseline_idcs])
    result_baseline_vols = TCatSubBrick(in_files=[(result.outputs.cbucket,
                                              f"'[{baseline_idcs_str}]'")],
                                        out_file=os.path.join(cwd,out_files_basename + '_baseline_runs.nii'),
                                        args='-overwrite').run()
    result_baseline = TStat(in_file=os.path.join(cwd,out_files_basename + '_baseline_runs.nii'),
                            args='-mean -overwrite',
                            out_file=os.path.join(cwd,out_files_basename + '_baseline.nii')).run()
    trialavg_files=[]
    for i in range(n_conditions):
        condition = condition_stim_files[i][0]
        result_condition_diffresponse_timecourse = TCatSubBrick(
            in_files=[(result.outputs.cbucket,
                       f"'[{int((polort+1)*n_files+i*n)}..{int((polort+1)*n_files+(i+1)*n-1)}]'")],
            out_file=os.path.join(cwd,out_files_basename + f"_diffresponse_condition_{condition}.nii"),
            args='-overwrite').run()
        result_condition_response_timecourse = Calc(
            in_file_a=os.path.join(cwd,out_files_basename + f"_diffresponse_condition_{condition}.nii"),
            in_file_b=os.path.join(cwd,out_files_basename + '_baseline.nii'),
            out_file=os.path.join(cwd,out_files_basename + f"_response_condition_{condition}.nii"),
            expr='a+b',
            args='-overwrite').run()
        trialavg_files.append(result_condition_response_timecourse.outputs.out_file)
        
    baseline_file = os.path.join(cwd,out_files_basename + '_baseline.nii')
    fstat_file = os.path.join(cwd,out_files_basename + '_fstat.nii')
    return trialavg_files, baseline_file, fstat_file

def calc_percent_change_trialavg(trialavg_files,baseline_file,inv_change=False):    
    prc_change = []
    for trialavg_file in trialavg_files:
        if inv_change:
            expr='100-(100*a/b)'
        else:
            expr='100*a/b-100'
        trialavg_file_split = os.path.splitext(trialavg_file)
        out_file = trialavg_file_split[0] + '_prcchg' + trialavg_file_split[1]
        result_prcchg = Calc(
            in_file_a=trialavg_file,
            in_file_b=baseline_file,
            out_file=out_file,
            expr=expr,
            args='-overwrite').run()
        prc_change.append(result_prcchg.outputs.out_file)
    return prc_change


def plot_roi_tcrs(file_list,roi,xlabel='volume',ylabel='signal'):
    df_dict = dict()
    condition_idx = 0
    for file in file_list:
        condition_idx = condition_idx + 1
        condition_data = sample_timecourse(file,roi)
        df = pd.DataFrame(condition_data,index=np.arange(condition_data.shape[0])).stack()
        df_dict['cond_'+str(condition_idx)] = df
    data = pd.concat(df_dict).reset_index()
    data.columns=['condition','volume','voxel','signal']
    ax = sns.lineplot(x='volume',y='signal',hue='condition',data=data,ci=68)
    ax.legend(loc='best')

def plot_cond_tcrs(condition_data_list,t=None,TR=1,labels=None,colors=None,ax=None,periods=None,events=None):
    """ Plots time course for multiple conditions. All timecourses should have the same length.
    """
    if ax==None:
        ax = plt.axes()
    if t==None:
        N = len(condition_data_list[0])
        t = TR * np.arange(0,N)
    for period in (periods or []):
        ax.axvspan(period[0],period[1],alpha=0.1,color='gray')
    ax.axhline(0,color='gray',lw=0.5)
    line_handles = []
    for condition_data,color in zip_longest(condition_data_list,(colors or []),fillvalue=None):
        l, = ax.plot(t,condition_data,color=color)
        line_handles.append(l)
    ax.set_xticks(t)
    ax.axis()
    ax.set_xlim(min(t),max(t))
    y0,y1 = ax.get_ylim()
    if labels:
        ax.legend(line_handles,labels,loc='best')
    for event in (events or []):
        ax.axvline(event[1],color='gray',lw=0.5)
        ax.annotate(event[0],(event[1],y0),ha='center',va='bottom')
    return ax


def preprocess_funcloc(data):
    # motion correction
    
    # spatial smoothing

    # high pass filtering

    # register to freesurfer

    pass

def feat_analysis(feat_dir,fsf_template):
    pass



def reg_feat_to_fs(feat_dir,fs_dir):
    subject=os.path.basename(os.path.normpath(fs_dir))
    subjects_dir=os.path.dirname(os.path.normpath(fs_dir))
    my_env = os.environ.copy()
    my_env['SUBJECTS_DIR'] = subjects_dir
    if subprocess.run(['bbregister',
                       '--mov',os.path.join(feat_dir,'example_func.nii.gz'),
                       '--bold',
                       '--s',subject,
                       '--lta',os.path.join(feat_dir,'feat2fs.lta')],
                      env=my_env).returncode == 0:
        return os.path.join(feat_dir,'feat2fs.lta')

def sample_surf_feat_stat(feat_dir,stat_file,fs_dir,hemi):
    subjects_dir=os.path.dirname(os.path.normpath(fs_dir))
    my_env = os.environ.copy()
    my_env['SUBJECTS_DIR'] = subjects_dir
    surf_suffix = fsl_remove_ext(os.path.basename(os.path.normpath(stat_file))) + '.mgh'
    out_file = os.path.join(feat_dir,'stats',hemi + '.' + surf_suffix)
    if subprocess.run(['mri_vol2surf',
                       '--mov',os.path.join(feat_dir,'stats',stat_file),
                       '--reg',os.path.join(feat_dir,'feat2fs.lta'),
                       '--projfrac','0.5',
                       '--interp','nearest',
                       '--hemi',hemi,
                       '--o',out_file],
                      env=my_env).returncode == 0:
        return out_file

def smooth_surf(in_file, out_file=None, fs_dir=None, hemi=None, fwhm=0):
    if out_file==None:
        out_file=os.path.splitext(os.path.normpath(in_file))[0] + '_smooth.mgh'
    subject=os.path.basename(os.path.normpath(fs_dir))
    subjects_dir=os.path.dirname(os.path.normpath(fs_dir))
    my_env = os.environ.copy()
    my_env['SUBJECTS_DIR'] = subjects_dir
    if subprocess.run(['mri_surf2surf',
                       '--hemi',hemi,
                       '--s','freesurfer',
                       '--fwhm',str(fwhm),
                       '--cortex',
                       '--sval', in_file,
                       '--tval', out_file],
                      env=my_env).returncode == 0:
        return out_file


def surf_activation_clusters(surf,cluster_surf=None,threshold=3,fs_dir=None,hemi=None):
    TODO: WIP
    if cluster_surf==None:
        out_file=os.path.splitext(os.path.normpath(cluster_surf))[0] + '_smooth.mgh'
    subject=os.path.basename(os.path.normpath(fs_dir))
    subjects_dir=os.path.dirname(os.path.normpath(fs_dir))
    my_env = os.environ.copy()
    my_env['SUBJECTS_DIR'] = subjects_dir
  

def get_funcloc_roi(parcel=None,analysis_dir=None,fs_dir=None,fs_to_func_reg=None,feat_dir=None):
    stat_name='zstat1'
    fwhm=5
    if feat_dir == None:
        feat_dir = os.path.join(analysis_dir,'funcloc.feat')
    stat_file = os.path.join(feat_dir,'stats',stat_name+'.nii.gz')
    # follow previous code, but think about better organization of ROIs, possibly in a single file?
    # assume we have functional activation maps and a functional to FS registration
    if not os.path.isfile(os.path.join(feat_dir,'feat2fs.lta')):
        reg_feat_to_fs(feat_dir,fs_dir)
    for hemi in ['lh','rh']:
        # 1. take activation map and project to surface
        stat_surf = sample_surf_feat_stat(feat_dir,stat_file,fs_dir,hemi)
        # 2. smooth on surface
        stat_surf_smooth = smooth_surf(stat_surf, fs_dir=fs_dir, hemi=hemi, fwhm=fwhm)
        # 3. generate activation clusters
        
        # 4. transform cluster label files to volume
        # 5. non-linearily transform to func

    return None

def cluster_surf():
    subject=os.path.basename(os.path.normpath(feat_dir))
    subjects_dir=os.path.dirname(os.path.normpath(feat_dir))
    my_env = os.environ.copy()
    my_env['SUBJECTS_DIR'] = subjects_dir
    subprocess.run(['mri_surfcluster',
                    '--in',os.path.join(feat_dir,'stats',hemi + '.zstat1_smooth.mgh'),
                    '--thmin','3',
                    '--sign','pos',
                    '--hemi',hemi,
                    '--subject', subjects,
                    '--o',os.path.join(feat_dir,'stats',hemi + '.zstat1_smooth_clusters.mgh'),
                    '--olab',os.path.join(feat_dir,'stats',hemi+'.zstat1_smooth_clusters.label')],
                   env=my_env)


def fs_surf_to_fs_volume():
    pass

def get_mni_coord_roi():
    # generate an ROI based on mni coordinates and a radius around
    # (
    pass

def layer_extend_roi_laynii():
    pass
                    
def layer_extend_roi_vfs(roi):
    pass

def add_prefix_to_nifti_basename(path,prefix):
    norm_path = os.path.normpath(path)
    dir_name = os.path.dirname(norm_path)
    base_name = os.path.basename(norm_path)
    return(os.path.join(dir_name,prefix+base_name))    

def add_postfix_to_nifti_basename(path,postfix):
    norm_path = os.path.normpath(path)
    s = os.path.splitext(norm_path)
    if s[1] == '.gz':
        s = os.path.splitext(s[0])
        basename = s[0]
        extension = s[1] + '.gz'
    else:
        basename = s[0]
        extension = s[1]
    return basename + postfix + extension

    

def get_funcact_roi_laynii(act_file,rim_file,roi_out_file,n_columns=10000,threshold=1):
    columns_file=add_postfix_to_nifti_basename(rim_file,'_columns'+str(n_columns))
    if not os.path.isfile(columns_file):
        mid_gm_file=add_postfix_to_nifti_basename(rim_file,'_midGM_equidist')
        subprocess.run(['LN2_COLUMNS',
                        '-rim',rim_file,
                        '-midgm',mid_gm_file,
                        '-nr_columns',str(n_columns)])
    subprocess.run(['LN2_MASK',
                    '-scores', act_file,
                    '-columns', columns_file,
                    '-min_thr', str(threshold),
                    '-output',roi_out_file])
    subprocess.run(['fslmaths',roi_out_file,
                    '-bin',roi_out_file])
    return roi_out_file

def roi_and(roi1,roi2):
    
    return intersect_masks((roi1,roi2), threshold=1, connected=False)

def get_funcact_roi_other_versions():
    # not clear yet what to do here, possibly manual deliniation needed
    # then fill out entire GM
    # alternatively go to surface, smooth and back?
    # choose cluster within region/ close to coordinates?
    # apply activation mask to one of the above ROI definitions?
    pass


# functional processing
def initialize_session():
    pass

def motion_correction():
    pass

def trial_averaging():
    pass

def glm_analysis():
    pass

def bold_correct(nulled_file,notnulled_file,out_file,notnulled_shift=None):
    """ notnulled_shift should equal (positive) difference between readout blocks
    """
    if notnulled_shift is not None:
        slicetimer_result = SliceTimer(in_file=notnulled_file,
                                       global_shift=-notnulled_shift)
        notnulled_file = slicetimer_result.outputs.out_file

    my_env=os.environ.copy()
        
    if os.path.normpath(out_file)[-3:]=='nii':
        my_env['FSLOUTPUTTYPE'] = 'NIFTI'
    elif  os.path.normpath(out_file)[-6:]=='nii.gz':
        my_env['FSLOUTPUTTYPE'] = 'NIFTI_GZ'
            
    subprocess.run(['fslmaths',
                    nulled_file,'-div', notnulled_file,
                    '-max','0',
                    '-min','5',
                    out_file],env=my_env)
    return out_file

def normalize():
    # Normalizes single voxel timecourses to change relative to a baselne
    pass

def avg_timcourse_nearest_volume():
                    pass
                    
def sample_timecourse(func_filename,roi):
    masked_data=apply_mask(func_filename,roi)
    return masked_data

def calc_layers_laynii(rim_file,out_file_base=None,method='equidist',n_layers=3):
    # include upsampling methods?
    if out_file_base is None:
        out_file_base = fsl_remove_ext(rim_file)
        
    run_string_list = ['LN2_LAYERS',
                       '-rim',rim_file,
                       '-output',out_file_base,
                       '-nr_layers',str(n_layers)]
    if method == 'equivol':
        run_string_list.append('-equivol')

    subprocess.run(run_string_list)
    if method== 'equivol':
        return out_file_base + "_metric_equivol.nii"
    else:
        return out_file_base + "_metric_equidist.nii"

def generate_two_layers(analysis_dir,depths,delta=0,roi=None):
    test = intersect_masks([depths,roi])
    superficial = math_img(f'img<{0.5-delta/2}',img=depths)
    deeper = math_img(f'img>{0.5-delta/2}',img=depths)
    if roi is not None:
        #superficial.affine = roi.affine
        superficial = intersect_masks([superficial,roi],threshold=1)
        deeper = intersect_masks([deeper,roi],threshold=1)
    nib.save(superficial, os.path.join(analysis_dir,'superficial.nii'))
    nib.save(deeper, os.path.join(analysis_dir,'deeper.nii'))

def sample_depths(data,roi,depths):
    # assume voxel matrix of data, roi and depths corresponds to each other
    data = copy.deepcopy(check_niimg(data,dtype='auto'))
    roi = copy.deepcopy(check_niimg(roi,dtype='auto'))
    depths = copy.deepcopy(check_niimg(depths,dtype='auto'))

    data.set_qform(np.eye(4))
    data.set_sform(np.eye(4))
    
    roi.set_qform(np.eye(4))
    roi.set_sform(np.eye(4))

    depths.set_qform(np.eye(4))
    depths.set_sform(np.eye(4))
    
    masked_data=apply_mask(data,roi)
    masked_depths=apply_mask(depths,roi)
    return masked_data, masked_depths

def sample_layer_profile(data,roi,depths,n_layers):
    data, depths = sample_depths(data,roi,depths)
    depths = np.floor(depths*n_layers)
    y = np.zeros(n_layers)
    for i in np.arange(n_layers):
        y[i] = np.mean(data[depths==i])
    return y, (np.arange(n_layers)+0.5)/n_layers

def plot_timecourses():
    pass

def plot_profile(data,roi,depths,n_layers):
    ax = plt.axes()
    voxel_responses, voxel_depths = sample_depths(data,roi,depths)
    ax.plot(voxel_depths,voxel_responses,'.',alpha=0.1)
    layer_responses,layer_depths = sample_layer_profile(data,roi,depths,n_layers)
    ax.plot(layer_depths,layer_responses)
    return ax

def upsample(in_file, out_file, factor, method):
    voxel_widths = np.array(nib.load(in_file).header.get_zooms())
    scaled_voxel_widths = voxel_widths/factor
    subprocess.run(['3dresample',
                    '-dxyz',
                    str(scaled_voxel_widths[0]),
                    str(scaled_voxel_widths[1]),
                    str(scaled_voxel_widths[2]),
                    '-rmode',method,
                    '-overwrite',
                    '-prefix',out_file,
                    '-input',in_file])
    

### FinnReplicationPilot specifc functions
def plot_finn_panel(depths,roi,trialavg_data,run_type,layers,ax,d=0):
    condition_data = []
    for file in trialavg_data:
        sampled_data, sampled_depths = sample_depths(file,roi,depths)
        if layers=='deep':
            condition_data.append(np.mean(sampled_data[:,sampled_depths<(0.5-d)],axis=1))
        elif layers=='superficial':
            condition_data.append(np.mean(sampled_data[:,sampled_depths>(0.5-d)],axis=1))
            
    if run_type=='alpha-rem':
        labels=['rem','alpha']
        colors=['tab:green','tab:blue']
    elif run_type=='go-nogo':
        labels=['nogo','go']
        colors=['tab:orange','tab:red']
        
    plot_cond_tcrs(condition_data,TR=3.702,
                   labels=labels,
                   colors=colors,
                   periods=[[4,14],[14,20]],
                   events=[['Stim',0],['Cue',4],['Probe',14]],ax=ax)
    plt.title(layers)
    
def plot_finn_tcrses(depths,roi,trialavg_alpharem,trialavg_gonogo,d=0):
    if type(depths)==str:
        depths = nib.load(depths)
        
    fig=plt.figure(figsize=(10,8), dpi= 100, facecolor='w', edgecolor='k')
    ax = plt.subplot(2,2,1)
    plot_finn_panel(depths,roi,trialavg_alpharem,'alpha-rem','superficial',ax,d)
                    
    ax = plt.subplot(2,2,2)
    plot_finn_panel(depths,roi,trialavg_gonogo,'go-nogo','superficial',ax,d)
    
    ax = plt.subplot(2,2,3)
    plot_finn_panel(depths,roi,trialavg_alpharem,'alpha-rem','deep',ax,d)
    
    ax = plt.subplot(2,2,4)
    plot_finn_panel(depths,roi,trialavg_gonogo,'go-nogo','deep',ax,d)

def get_finn_tcrs_data(trial_averages,roi,):
    data = dict()
    for layer in layers:
        for condition in conditions:
            for modality in modalities:
                data[modality,layer,condition] = np.loadtxt(
                    fnamebase + modality + '_' + layer + '_'+  condition + '.1D')
    return data

def plot_finn_tcrs(fnamebase,modality,TR=3.702):
    data = get_finn_tcrs_data(fnamebase)
    fix, axes = plt.subplot(2,3,figsize=(15,5))
    periods=[[4,14],[14,20]]
    events=[['Stim',0],['Cue',4],['Probe',14]]
    for row, layer in enumerate(layers):
        plot_cond_tcrs([data[modality,layer,'alpha'],
                        data[modality,layer,'rem']],
                       colors=('tab:blue','tab:green'),
                       labels=('alpha','rem'),
                       periods=periods,
                       events=events,
                       TR=TR,
                       ax=axes[row,0])
        
        plot_cond_tcrs([data[modality,layer,'go'],
                        data[modality,layer,'nogo']],
                       colors=('tab:red','tab:orange'),
                       labels=('act','non-act'),
                       periods=periods,
                       events=events,
                       TR=TR,
                       ax=axes[row,0])

        plot_cond_tcrs([data[modality,layer,'alpha'] -  data[modality,layer,'rem'],
                        data[modality,layer,'go'] - data[modality,layer,'nogo']],                  
                       colors=('tab:purple','tab:cyan'),
                       labels=('alpha - rem','act - non-act'),
                       periods=periods,
                       events=events,
                       TR=TR,
                       ax=axes[row,0])
    axes[0,0].set_ylabel('signal change [%]')
    axes[1,0].set_ylabel('signal change [%]')
    axes[1,0].set_xlabel('trial time [s]')
    axes[1,1].set_xlabel('trial time [s]')
    axes[1,2].set_xlabel('trial time [s]')
        
    fig.text(0.08,0.45,'deeper',ha='right',weight='bold')
    fig.text(0.08,0.85,'superficial',ha='right',weight='bold')
    fig.suptitle(modality.upper(),weight='bold')

def finn_trial_averaging(run_type,analysis_dir):
    trial_duration = 32
    trial_order = paradigm(run_type)
    trialavg = dict()
    onset_delay = 8
    in_files_bold = [os.path.join(analysis_dir,f'func_{run_type}_notnulled_tshift.nii.gz')]
    in_files_vaso = [os.path.join(analysis_dir,f'func_{run_type}_vaso.nii.gz')]
    stim_times_runs = [calc_stim_times(onset_delay=8,trial_duration=trial_duration,
                                       trial_order=trial_order)]
    trialavg_files_bold, baseline_file_bold, fstat_file_bold = \
        average_trials_3ddeconvolve(in_files_bold,                                                                                                            stim_times_runs,
                                    trial_duration,
                                    out_files_basename='trialavg1_bold_' + run_type,
                                    polort=5)
    
    trialavg_files_vaso, baseline_file_vaso, fstat_file_vaso = \
        average_trials_3ddeconvolve(in_files_vaso,
                                    stim_times_runs,
                                    trial_duration,
                                    out_files_basename='trialavg1_vaso_' + run_type,
                                    polort=5)
    
    trialavg_bold_prcchg = calc_percent_change_trialavg(trialavg_files_bold,
                                                        baseline_file_bold,
                                                        inv_change=False)
    trialavg_vaso_prcchg = calc_percent_change_trialavg(trialavg_files_vaso,
                                                        baseline_file_vaso,
                                                        inv_change=True)

    return trialavg_bold_prcchg, trialavg_vaso_prcchg, fstat_file_bold, fstat_file_vaso

def finn_trial_averaging_with_boldcorrect(run_type,analysis_dir,TR1):
    trial_duration = 32
    trial_order = paradigm(run_type)
    trialavg = dict()
    onset_delay = 8
    in_files_nulled = [os.path.join(analysis_dir,f'func_{run_type}_nulled.nii')]
    in_files_notnulled = [os.path.join(analysis_dir,f'func_{run_type}_notnulled.nii')]
    stim_times_runs = [calc_stim_times(onset_delay=8,trial_duration=trial_duration,
                                       trial_order=trial_order)]

    trialavg_files_nulled, baseline_file_nulled, fstat_file_nulled = \
        average_trials_3ddeconvolve(in_files_nulled,
                                    stim_times_runs,
                                    trial_duration,
                                    out_files_basename='trialavg2_nulled',
                                    polort=5)
    
    trialavg_files_notnulled, baseline_file_notnulled, fstat_file_notnulled = \
        average_trials_3ddeconvolve(in_files_notnulled,
                                    stim_times_runs,
                                    trial_duration,
                                    out_files_basename='trialavg2_notnulled',
                                    polort=5,
                                    onset_shift=TR1)
  
    trialavg_files_vaso = [bold_correct(trialavg_files_nulled[0],trialavg_files_notnulled[0],
                                        trialavg_files_nulled[0].replace('nulled','vaso')),
                           bold_correct(trialavg_files_nulled[1],trialavg_files_notnulled[1],
                                        trialavg_files_nulled[1].replace('nulled','vaso'))]
                                                                                 
    baseline_file_vaso = bold_correct(baseline_file_nulled,baseline_file_notnulled,
                                      baseline_file_nulled.replace('nulled','vaso'))
        
    trialavg_bold_prcchg = calc_percent_change_trialavg(trialavg_files_notnulled,
                                                        baseline_file_notnulled,
                                                        inv_change=False)

    trialavg_vaso_prcchg = calc_percent_change_trialavg(trialavg_files_vaso,
                                                        baseline_file_vaso,
                                                        inv_change=True)

    fstat_file_bold = fstat_file_notnulled
    
    return trialavg_bold_prcchg, trialavg_vaso_prcchg, fstat_file_bold, fstat_file_nulled

# Define a function to obtain some paradigm related info. (For now trial order, TODO: trial period timings, GLM events, ...)
# we need:
# for trial averaging: exact trial onset times (what about VASO,GE-BOLD shifts?) with conditions
# + volume infor for period averagin?
# for glm analysis
# onsets and durations of all trial periods
def paradigm(run_type):
    conditionRem            = 2
    conditionAlpha          = 3
    conditionNogo           = 4
    conditionGo             = 5
    if run_type=='localizer':        
        letterStringDuration = 2.5
        fix1Duration         = 1.5
        cueDuration          = 1
        fixDelayDuration     = 9
        probeDuration        = 1
        interTrialDuration   = 5
        
        startBlankPeriod = 6
        
        trial_order=[3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]
    elif run_type in ['alpha-rem','go-nogo']:
        letterStringDuration = 2.5
        fix1Duration         = 1.5
        cueDuration          = 1
        fixDelayDuration     = 9
        probeDuration        = 2
        interTrialDuration   = 16    

        startBlankPeriod = 8
        
        if run_type=='alpha-rem':
            trial_order=[2,3,3,3,2,2,3,2,3,3,2,2,2,3,2,3,3,3,2,2]
        elif run_type=='go-nogo':
            trial_order=[4,5,5,5,4,4,5,4,5,5,4,4,4,5,4,5,5,5,4,4]
    
    return trial_order
