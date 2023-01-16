import copy
import os
import subprocess
import tempfile
from collections import defaultdict
from itertools import zip_longest
from shutil import rmtree, copy2

import numpy as np
import pandas as pd
import nibabel as nib
from nilearn._utils import check_niimg
from nilearn.image import math_img, mean_img, index_img
from nilearn.masking import apply_mask, intersect_masks
from nipype.interfaces.afni import Deconvolve, TCatSubBrick, Calc, TStat
from nipype.interfaces.ants import ApplyTransformsToPoints
from nipype.interfaces.freesurfer import MRIsConvert
from nipype.interfaces.fsl import SliceTimer
from niworkflows.interfaces.surf import CSVToGifti, GiftiToCSV
import glob


def fsl_remove_ext(filename):
    """
    Removes extension from filename using fsl's remove_ext
    :param filename:
    :return:
    """
    result = subprocess.run(["remove_ext", filename], stdout=subprocess.PIPE)
    return result.stdout.strip().decode()


### Registrations and Transform related functions
def surftransform_gii(gii_surf, transforms, invert_transform_flags, cwd=None):
    """
    Takes gifti surface and applies ants transforms
    :param gii_surf:
    :param transforms:
    :param invert_transform_flags:
    :param cwd:
    :return:
    """
    if cwd is None:
        cwd = os.path.dirname(os.path.normpath(gii_surf))
    out_file = os.path.basename(os.path.normpath(gii_surf))
    # convert gii to csv
    result_GiftiToCSV = GiftiToCSV(in_file=gii_surf, itk_lps=True).run(cwd=cwd)
    csv_surf = result_GiftiToCSV.outputs.out_file
    # apply transform
    result_ApplyTransformsToPoints = ApplyTransformsToPoints(
        dimension=3,
        input_file=csv_surf,
        transforms=transforms,
        invert_transform_flags=invert_transform_flags,
    ).run(cwd=cwd)
    csv_surf_transformed = result_ApplyTransformsToPoints.outputs.output_file
    # convert csv to gii
    result_CSVToGifti = CSVToGifti(
        in_file=csv_surf_transformed, gii_file=gii_surf, itk_lps=True
    ).run(cwd=cwd)
    gii_surf_transformed = result_CSVToGifti.outputs.out_file
    return gii_surf_transformed


def surftransform_fs(fs_surf, transforms, invert_transform_flags, out_file, cwd=None):
    """
    Takes freesurfer surface and applies ants transforms
    :param fs_surf:
    :param transforms:
    :param invert_transform_flags:
    :param out_file:
    :param cwd:
    :return:
    """
    if cwd == None:
        cwd = os.path.dirname(os.path.normpath(out_file))
    # convert fs to gii
    result_MRIsConvert = MRIsConvert(
        in_file=fs_surf, out_datatype="gii", to_scanner=True
    ).run(cwd=cwd)
    gii_surf = os.path.join(cwd, result_MRIsConvert.outputs.converted)
    gii_surf_transformed = surftransform_gii(
        gii_surf, transforms, invert_transform_flags
    )
    # convert gii to fs
    MRIsConvert(in_file=gii_surf_transformed, out_file=out_file, to_tkr=True).run(
        cwd=cwd
    )
    return out_file


def fs_surface_to_func(fs_to_func_reg, fs_dir, analysis_dir=None, force=False):
    """
    Transforms freesurfer surfaces to functional space using ANTs transfrom.
    :param fs_to_func_reg:
    :param fs_dir:
    :param analysis_dir:
    :param force:
    :return:
    """
    if analysis_dir == None:
        analysis_dir = os.path.join(fs_dir, "surf")
    transform_0_lin = fs_to_func_reg[1]
    transform_1_inversewarp = fs_to_func_reg[3]
    invert_transform_flags = [True, False]
    surf_trans_files = dict()
    for hemi in ["lh", "rh"]:
        for surf_type in ["white", "pial"]:
            surf = os.path.join(fs_dir, "surf", hemi + "." + surf_type)
            surf_trans = os.path.join(analysis_dir, hemi + "." + surf_type + "_func")
            if not os.path.isfile(surf_trans) or force == True:
                surf_trans_files[hemi, surf_type] = surftransform_fs(
                    surf,
                    [transform_0_lin, transform_1_inversewarp],
                    invert_transform_flags,
                    out_file=surf_trans,
                )
            else:
                surf_trans_files[hemi, surf_type] = surf_trans
    return surf_trans_files


def ciftify_surface_to_func(fs_to_func_reg, ciftify_dir, analysis_dir=None):
    """
    Transforms ciftify surfaces to functional space using ANTs transfrom.
    :param fs_to_func_reg:
    :param ciftify_dir:
    :param analysis_dir:
    :return:
    """
    # TODO: not working correctly -> check and fix
    if analysis_dir is None:
        analysis_dir = os.path.join(ciftify_dir, "T1w", "fsaverage_LR164k")
    ciftify_subject = os.path.basename(os.path.normpath(ciftify_dir))
    transform_0_lin = fs_to_func_reg[1]
    transform_1_inversewarp = fs_to_func_reg[3]
    invert_transform_flags = [True, False]
    surf_trans_files = dict()
    for hemi in ["L", "R"]:
        for surf_type in ["white", "pial"]:
            surf = os.path.join(
                ciftify_dir,
                "T1w",
                "fsaverage_LR164k",
                ciftify_subject + "." + hemi + "." + surf_type + ".164k_fs_LR.surf.gii",
            )
            surf_trans = os.path.join(
                analysis_dir,
                ciftify_subject
                + "."
                + hemi
                + "."
                + surf_type
                + ".164k_fs_LR_func.surf.gii",
            )
            out_file = surftransform_gii(
                surf,
                [transform_0_lin, transform_1_inversewarp],
                invert_transform_flags,
                cwd=analysis_dir,
            )
            os.rename(out_file, surf_trans)
            surf_trans_files[hemi, surf_type] = surf_trans
    return surf_trans_files


def process_vaso(
    session_dir,
    process_script,
    analysis_dir=None,
    alpharem_runs=None,
    gonogo_runs=None,
    analysis_subdir="analysis",
):
    """
    Wrapper function to run an external script that does VASO processing.
    :param session_dir:
    :param process_script:
    :param analysis_dir:
    :param alpharem_runs:
    :param gonogo_runs:
    :param analysis_subdir:
    :return:
    """
    if analysis_dir is None:
        analysis_dir = os.path.join(session_dir, analysis_subdir)
    if not os.path.isdir(analysis_dir):
        os.mkdir(analysis_dir)
        if alpharem_runs is not None:
            with open(
                os.path.join(analysis_dir, "func_alpha-rem_task-runs.txt"), "w"
            ) as file:
                print(*alpharem_runs, sep="\n", file=file)
        if gonogo_runs is not None:
            with open(
                os.path.join(analysis_dir, "func_go-nogo_task-runs.txt"), "w"
            ) as file:
                print(*gonogo_runs, sep="\n", file=file)
        subprocess.run([process_script, session_dir, analysis_dir])
    return analysis_dir


def register_fs_to_vasot1(fs_dir, analysis_dir, use_brain=False, force=False):
    """
    Wrapper function to run an external script that does registration of freesurfer dataset to VASO T1
    (functional space).
    :param fs_dir:
    :param analysis_dir:
    :param use_brain:
    :param force:
    """
    if (
        not os.path.isfile(os.path.join(analysis_dir, "fs_to_func_0GenericAffine.mat"))
        or force == True
    ):
        if use_brain == True:
            target = "func_all_T1_brain.nii"
        else:
            target = "func_all_T1.nii"

        subprocess.run(
            ["register_fs-to-vasoT1.sh", target, fs_dir, "itksnap"], cwd=analysis_dir
        )


def apply_ants_transforms(vol_in, vol_out, ref_vol, affine, warp):
    """
    Wrapper function to run an external script that applies a ANTS non-linear registration transform,
    consisting of 1. warp and 2. affine to input volume.
    :param vol_in:
    :param vol_out:
    :param ref_vol:
    :param affine:
    :param warp:
    """
    subprocess.run(
        [
            "antsApplyTransforms",
            "--interpolation",
            "BSpline" "[5]" "",
            "-d",
            "3",
            "-i",
            vol_in,
            "-r",
            ref_vol,
            "-t",
            warp,
            "-t",
            affine,
            "-o",
            vol_out,
            "-n",
            "NearestNeighbor",
        ]
    )
    # TODO: Check if the following comment can be removed.
    # NOTE: It seemed necessary, because the resulting affine from ANTS was not exactly the same as the ref volume
    # (they differ at the 8th decimal after the dot), but why are they not exatly the same?
    # nii_vol_out = nib.load(vol_out)
    # nii_ref_vol = nib.load(ref_vol)
    # nii_vol_out.header.set_qform(nii_ref_vol.header.get_qform())
    # nii_vol_out.header.set_sform(nii_ref_vol.header.get_sform())
    # nib.save(nii_vol_out,vol_out)


def import_fs_ribbon_to_func(fs_dir, analysis_dir, force=False):
    """
    Wrapper function to run external script that imports gray matter ribbon from freesurfer dataset and
    transforms it to functional space. Assumes there are fs-to-func registration files in analysis_dir
    :param fs_dir:
    :param analysis_dir:
    :param force:
    :return:
    """
    rim_file = os.path.join(analysis_dir, "rim.nii")
    if not os.path.isfile(rim_file) or force == True:
        if (
            subprocess.run(
                [
                    "/data/p_02389/code/fmri-analysis/library/import-fs-ribbon.sh",
                    fs_dir,
                    analysis_dir,
                    os.path.join(analysis_dir, "fs_t1_in-func.nii"),
                ]
            ).returncode
            != 0
        ):
            return None
    return rim_file


### ROI related functions

def generate_atlas_region_hcp(atlas_file, out_file, label_list):
    """
    generate a surface roi from a list of atlas labels
    """
    atlas = nib.load(atlas_file)
    roi_data = np.isin(atlas.darrays[0].data.copy(),label_list).astype(np.int32)
    roi_gii = nib.GiftiImage(header=atlas.header,
                             darrays = [nib.gifti.GiftiDataArray(data=roi_data,
                                                                 intent=1011)],
                             meta=atlas.meta)
    roi_gii.to_filename(out_file)
    return out_file

def index_roi_surflabel_hcp(label_file,roi_out,label):
    if type(label)==str:
        subprocess.run(['wb_command',
                        '-gifti-label-to-roi',
                        label_file,
                        roi_out,
                        '-name',
                        label],check=True)
    else:
        subprocess.run(['wb_command',
                        '-gifti-label-to-roi',
                        label_file,
                        roi_out,
                        '-key',
                        str(label)],check=True)
        return roi_out
                    

def index_roi(roi, idx):
    """
    Extracts ROI with a specific index from a multi-index label file.
    :param roi:
    :param idx:
    :return:
    """
    return math_img(f"img=={idx}", img=roi)


def fs_LR_label_to_fs_volume(ciftify_dir, analysis_dir, labels, hemi, out_basename):
    """Transforms label .gii file in fs_LR space to Freesurfer volume."""
    ciftify_subject = os.path.basename(os.path.normpath(ciftify_dir))
    mid_surf = os.path.join(
        ciftify_dir,
        "T1w",
        "fsaverage_LR164k",
        ciftify_subject + "." + hemi + ".midthickness.164k_fs_LR.surf.gii",
    )
    white_surf = os.path.join(
        ciftify_dir,
        "T1w",
        "fsaverage_LR164k",
        ciftify_subject + "." + hemi + ".white.164k_fs_LR.surf.gii",
    )
    pial_surf = os.path.join(
        ciftify_dir,
        "T1w",
        "fsaverage_LR164k",
        ciftify_subject + "." + hemi + ".pial.164k_fs_LR.surf.gii",
    )
    volume = os.path.join(ciftify_dir, "T1w", "T1w.nii.gz")
    volume_out = os.path.join(analysis_dir, out_basename + "_labels_" + hemi + ".nii")
    subprocess.run(
        [
            "wb_command",
            "-label-to-volume-mapping",
            labels,
            mid_surf,
            volume,
            volume_out,
            "-ribbon-constrained",
            white_surf,
            pial_surf,
            "-greedy",
        ]
    )
    return volume_out


def get_fs_LR_atlas_roi(
    parcel=None,
    atlas_labels=None,
    out_basename=None,
    analysis_dir=None,
    ciftify_dir=None,
    fs_to_func_reg=None,
    force=False,
):
    """Returns an ROI in functional space by transforming GIFTI label files in fs_LR space.
    ROI is specified using parcel=(hemi,idx).
    """
    hemi = parcel[0]
    parcel_idx = parcel[1]
    labels_in_func = os.path.join(
        analysis_dir, out_basename + "_labels_" + hemi + "_in-func.nii"
    )
    if not os.path.isfile(labels_in_func) or force == True:
        labels_in_fs_individual = fs_LR_label_to_fs_volume(
            ciftify_dir, analysis_dir, atlas_labels[hemi], hemi, out_basename
        )
        apply_ants_transforms(
            vol_in=labels_in_fs_individual,
            vol_out=labels_in_func,
            ref_vol=fs_to_func_reg[0],
            affine=fs_to_func_reg[1],
            warp=fs_to_func_reg[2],
        )
    roi = index_roi(labels_in_func, parcel_idx)
    return roi


def fs_overlay_to_fs_volume(
    overlay_file, fs_dir, analysis_dir, hemi, out_basename=None, force=False
):
    # assumes overlay belongs to {hemi}.white
    if out_basename == None:
        out_file = os.path.splitext(os.path.normpath(overlay_file))[0] + ".nii"
    else:
        out_file = os.path.join(analysis_dir, out_basename + "_labels_" + hemi + ".nii")

    subject = os.path.basename(os.path.abspath(fs_dir))
    subjects_dir = os.path.dirname(os.path.abspath(fs_dir))
    my_env = os.environ.copy()
    my_env["SUBJECTS_DIR"] = subjects_dir

    if not os.path.isfile(out_file) or force == True:
        if (
            subprocess.run(
                [
                    "mri_surf2vol",
                    "--so",
                    os.path.join(fs_dir, "surf", hemi + ".white"),
                    overlay_file,
                    "--subject",
                    subject,
                    "--o",
                    out_file,
                    "--ribbon",
                    os.path.join(fs_dir, "mri", "ribbon.mgz"),
                ],
                env=my_env,
            ).returncode
            != 0
        ):
            return None
        # remove -1 values
        out_file_nii = nib.load(out_file)
        out_file_data = out_file_nii.get_fdata()
        out_file_data[out_file_data == -1] = 0
        out_file_changed_nii = nib.Nifti1Image(
            out_file_data, out_file_nii.affine, out_file_nii.header
        )
        nib.save(out_file_changed_nii, out_file)
    return out_file


def get_fs_roi(
    parcel=None,
    overlay_files=None,
    out_basename=None,
    analysis_dir=None,
    fs_dir=None,
    fs_to_func_reg=None,
    force=False,
):
    hemi = parcel[0]
    parcel_idx = parcel[1]
    labels_in_func = os.path.join(
        analysis_dir, out_basename + "_labels_" + hemi + "_in-func.nii"
    )
    if not os.path.isfile(labels_in_func) or force == True:
        labels_in_fs_individual = fs_overlay_to_fs_volume(
            overlay_files[hemi], fs_dir, analysis_dir, hemi, out_basename, force
        )
        apply_ants_transforms(
            vol_in=labels_in_fs_individual,
            vol_out=labels_in_func,
            ref_vol=fs_to_func_reg[0],
            affine=fs_to_func_reg[1],
            warp=fs_to_func_reg[2],
        )
    roi = index_roi(labels_in_func, parcel_idx)
    return roi


def reg_feat_to_fs(feat_dir, fs_dir, force=False):
    subject = os.path.basename(os.path.normpath(fs_dir))
    subjects_dir = os.path.dirname(os.path.normpath(fs_dir))
    my_env = os.environ.copy()
    my_env["SUBJECTS_DIR"] = subjects_dir

    reg_file = os.path.join(feat_dir, "feat2fs.lta")

    if not os.path.isfile(reg_file) or force == True:
        if (
            subprocess.run(
                [
                    "bbregister",
                    "--mov",
                    os.path.join(feat_dir, "example_func.nii.gz"),
                    "--bold",
                    "--s",
                    subject,
                    "--lta",
                    os.path.join(feat_dir, "feat2fs.lta"),
                ],
                env=my_env,
            ).returncode
            != 0
        ):
            return None
    return reg_file


def sample_surf_feat_stat(feat_dir, stat_file, fs_dir, hemi, force=False):
    subjects_dir = os.path.dirname(os.path.normpath(fs_dir))
    my_env = os.environ.copy()
    my_env["SUBJECTS_DIR"] = subjects_dir
    surf_suffix = fsl_remove_ext(os.path.basename(os.path.normpath(stat_file))) + ".mgh"
    out_file = os.path.join(feat_dir, "stats", hemi + "." + surf_suffix)
    if not os.path.isfile(out_file) or force == True:
        if (
            subprocess.run(
                [
                    "mri_vol2surf",
                    "--mov",
                    os.path.join(feat_dir, "stats", stat_file),
                    "--reg",
                    os.path.join(feat_dir, "feat2fs.lta"),
                    "--projfrac",
                    "0.5",
                    "--interp",
                    "nearest",
                    "--hemi",
                    hemi,
                    "--o",
                    out_file,
                ],
                env=my_env,
            ).returncode
            != 0
        ):
            return None
    return out_file


def smooth_surf(in_file, out_file=None, fs_dir=None, hemi=None, fwhm=0, force=False):
    if out_file == None:
        out_file = os.path.splitext(os.path.normpath(in_file))[0] + "_smooth.mgh"
    subject = os.path.basename(os.path.abspath(fs_dir))
    subjects_dir = os.path.dirname(os.path.abspath(fs_dir))
    my_env = os.environ.copy()
    my_env["SUBJECTS_DIR"] = subjects_dir
    if not os.path.isfile(out_file) or force == True:
        if (
            subprocess.run(
                [
                    "mri_surf2surf",
                    "--hemi",
                    hemi,
                    "--s",
                    subject,
                    "--fwhm",
                    str(fwhm),
                    "--cortex",
                    "--sval",
                    in_file,
                    "--tval",
                    out_file,
                ],
                env=my_env,
            ).returncode
            != 0
        ):
            return None
    return out_file


def cluster_surf(
    in_file,
    out_file=None,
    fs_dir=None,
    hemi=None,
    threshold=10,
    force=False,
    sign="pos",
):
    if out_file == None:
        out_file = os.path.splitext(os.path.normpath(in_file))[0] + "_clusters.mgh"
    annot_file = os.path.splitext(os.path.normpath(out_file))[0] + ".annot"

    subject = os.path.basename(os.path.normpath(fs_dir))
    subjects_dir = os.path.dirname(os.path.normpath(fs_dir))
    my_env = os.environ.copy()
    my_env["SUBJECTS_DIR"] = subjects_dir

    if not os.path.isfile(out_file) or force == True:
        if (
            subprocess.run(
                [
                    "mri_surfcluster",
                    "--in",
                    in_file,
                    "--thmin",
                    str(threshold),
                    "--sign",
                    sign,
                    "--hemi",
                    hemi,
                    "--subject",
                    subject,
                    "--oannot",
                    annot_file,
                ],
                env=my_env,
            ).returncode
            != 0
        ):
            return None
    # convert to mgh format
    labels, _, _ = nib.freesurfer.io.read_annot(annot_file)
    in_mgh = nib.load(in_file)
    labels_mgh = nib.freesurfer.MGHImage(labels, in_mgh.affine, in_mgh.header)
    nib.save(labels_mgh, out_file)

    return out_file


def math_cifti(expr, cifti_out, **ciftis):
    """
    runs wb_command -cifti-math
    **ciftis can be : cifti1 = cifti1_file.nii, cifti2 = cifti2_file.nii, ...
    """
    cmd = ["wb_command", "-cifti-math", expr, cifti_out] + sum(
        [["-var", name, ciftis[name]] for name in ciftis], []
    )
    subprocess.run(cmd, check=True)
    return cifti_out


def math_metric(expr, metric_out, **metrics):
    """
    runs wb_command -metric math
    **metrics can be : metric1 = metric1_file.gii, metric2 = metric2_file.gii, ...
    """
    cmd = ["wb_command", "-metric-math", expr, metric_out] + sum(
        [["-var", name, metrics[name]] for name in metrics], []
    )
    subprocess.run(cmd, check=True,stdout=subprocess.DEVNULL)
    return metric_out

def stats_metric_hcp(metric_in, op, roi=None):
    cmd = ["wb_command",
           "-metric-stats",
           metric_in,
           "-reduce",op]
    if roi is not None:
        cmd += ["-roi",roi]
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE)
    return float(result.stdout)

    
def mask_metric_hcp(metric_in, metric_out, mask):
    subprocess.run(
        ["wb_command", "-metric-mask", metric_in, mask, metric_out], check=True
    )
    return metric_out


def find_clusters_hcp(metric_in, metric_out, mid_surf, threshold, min_area=0, roi=None):
    cmd = [
        "wb_command",
        "-metric-find-clusters",
        mid_surf,
        metric_in,
        str(threshold),
        str(min_area),
        metric_out,
    ]
    if roi:
        cmd += ["-roi", roi]
    subprocess.run(cmd, check=True)

def surf_to_vol_hcp(metric_in, volume_out, volume,
                    white_surf, pial_surf, mid_surf,greedy=False,is_roi=False):
    with tempfile.TemporaryDirectory() as tmpdirname:
        label_file = os.path.join(tmpdirname,'roi.label.gii')
        subprocess.run(["wb_command",
                        "-metric-label-import",
                        metric_in,
                        "",
                        label_file,
                        ])
        
        cmd = ["wb_command"]
        if is_roi:
            metric_in = label_file
            cmd += ["-label-to-volume-mapping"]
        else:
            cmd += ["-metric-to-volume-mapping"]
        
        cmd += [metric_in,
                mid_surf,
                volume,
                volume_out,
                "-ribbon-constrained",
                white_surf,
                pial_surf]
        if greedy==True:
            cmd += ['-greedy']
        subprocess.run(cmd,check=True)
    return volume_out
    
def sample_surf_hcp(
        volume_file, white_surf, pial_surf, mid_surf, outfile, mask_file=None, roi_out=None
):
    """
    Samples volume to surface using arbitrary GIFTI surfaces using hcp tools (wb_command).
    - generates midthickness if file does not exists
    :return:
    """
    # create midthickness
    if not os.path.isfile(mid_surf):
        subprocess.run(
            [
                "wb_command",
                "-surface-average",
                mid_surf,
                "-surf",
                pial_surf,
                "-surf",
                white_surf,
            ]
        )

    cmd_volume_to_surface = [
        "wb_command",
        "-volume-to-surface-mapping",
        volume_file,
        mid_surf,
        outfile,
        "-ribbon-constrained",
        white_surf,
        pial_surf,
    ]
    if roi_out is not None:
        cmd_volume_to_surface += ["-bad-vertices-out", roi_out]

    if mask_file is None:
        subprocess.run(cmd_volume_to_surface, check=True)
        return outfile, mid_surf
    else:
        cmd_volume_to_surface += ["-volume-roi", mask_file]
        cmd_fill_in_holes = [
            "wb_command",
            "-metric-dilate",
            outfile,
            mid_surf,
            str(0),
            outfile,
            "-nearest",
        ]

        subprocess.run(cmd_volume_to_surface, check=True)
        subprocess.run(cmd_fill_in_holes, check=True)
        return outfile, mid_surf


def smooth_surfmetric_hcp(metric_in, metric_out, mid_surf, fwhm):
    subprocess.run(
        [
            "wb_command",
            "-metric-smoothing",
            mid_surf,
            metric_in,
            str(fwhm),
            metric_out,
            "-fwhm",
            "-fix-zeros",
        ],
        check=True,
    )
    return metric_out


def transform_data_native_surf_to_fs_LR(
    data_native_surf, data_fs_LR_surf, native_mid_surf, hemi, ciftify_dir
):
    # find names (subject) of native anf fs_LR spheres
    native_sphere = glob.glob(
        os.path.join(
            ciftify_dir,
            "MNINonLinear",
            "Native",
            f"*.{hemi}.sphere.MSMSulc.native.surf.gii",
        )
    )[0]
    fs_LR_sphere = glob.glob(
        os.path.join(
            ciftify_dir, "MNINonLinear", f"*.{hemi}.sphere.164k_fs_LR.surf.gii"
        )
    )[0]

    # find name of fs_LR midthickness (for area correction)
    fs_LR_mid_surf = glob.glob(
        os.path.join(
            ciftify_dir, "MNINonLinear", f"*.{hemi}.midthickness.164k_fs_LR.surf.gii"
        )
    )[0]

    # find names of surface rois to exclude medial wall
    native_surf_roi = glob.glob(
        os.path.join(
            ciftify_dir, "MNINonLinear", "Native", f"*.{hemi}.roi.native.shape.gii"
        )
    )[0]
    fs_LR_surf_roi = glob.glob(
        os.path.join(
            ciftify_dir, "MNINonLinear", f"*.{hemi}.atlasroi.164k_fs_LR.shape.gii"
        )
    )[0]

    cmd1 = [
        "wb_command",
        "-metric-resample",
        data_native_surf,
        native_sphere,
        fs_LR_sphere,
        "ADAP_BARY_AREA",
        data_fs_LR_surf,
        "-area-surfs",
        native_mid_surf,
        fs_LR_mid_surf,
        "-current-roi",
        native_surf_roi,
    ]

    cmd2 = [
        "wb_command",
        "-metric-mask",
        data_fs_LR_surf,
        fs_LR_surf_roi,
        data_fs_LR_surf,
    ]

    subprocess.run(cmd1, check=True)
    subprocess.run(cmd2, check=True)
    return data_fs_LR_surf


def calc_area_hcp(roi, mid_surf):
    result = subprocess.run(
        [
            "wb_command",
            "-metric-weighted-stats",
            roi,
            "-sum",
            "-area-surface",
            mid_surf,
        ],
        check=True,
        stdout=subprocess.PIPE,
    )
    return float(result.stdout)


def sample_layer_to_fs_LR(
    volume_file,
    output_file,
    white_surf,
    pial_surf,
    ciftify_dir,
    hemi,
    depth_range=[0, 1],
    mask=None,
    depth_file=None,
):
    # if depth_file provided, use it to calculate mask, otherwise generate intermediate layer surfaces
    # depth_range: 0 = wm boundary, 1 = pial surface

    with tempfile.TemporaryDirectory() as tmpdirname:
        mid_surf = os.path.join(tmpdirname, "mid.surf.gii")
        data_native_surf = os.path.join(tmpdirname, "data_native_surf.func.gii")
        mask_file = os.path.join(tmpdirname, "mask.nii")

        # 1. generate boundary surfaces or compute layer mask
        if depth_file:
            if type(depth_file) == str:
                depth_file = nib.load(depth_file)
            # print(depth_file.affine)
            depth_surfs = [white_surf, pial_surf]
            layer_roi = math_img(
                f"(img>={depth_range[0]})& (img<={depth_range[1]})", img=depth_file
            )
            if mask is None:
                mask = layer_roi
            else:
                mask = roi_and((mask, layer_roi))
        else:
            depth_surfs = [
                os.path.join(tmpdirname, f"depth{i}.surf.gii") for i in [0, 1]
            ]
            cmds_depth_surf = [
                [
                    "wb_command",
                    "-surface-cortex-layer",
                    white_surf,
                    pial_surf,
                    str(depth_range[i]),
                    depth_surfs[i],
                ]
                for i in [0, 1]
            ]

            subprocess.run(cmds_depth_surf[0], check=True)
            subprocess.run(cmds_depth_surf[1], check=True)

        # 2. sample
        if isinstance(mask, nib.nifti1.Nifti1Image):
            nib.save(mask, mask_file)
            mask = mask_file
        # print(nib.load(mask_file).affine)
        # print(nib.load(volume_file).affine)
        # print(depth_surfs)
        data_native_surf, mid_surf = sample_surf_hcp(
            volume_file,
            depth_surfs[0],
            depth_surfs[1],
            mid_surf,
            outfile=data_native_surf,
            mask_file=mask,
        )
        # 3. resample to fs_LR
        output_file = transform_data_native_surf_to_fs_LR(
            data_native_surf, output_file, mid_surf, hemi, ciftify_dir
        )

    return output_file


def sample_surf_func_stat(
    stat_file,
    white_surf_file,
    thickness_file,
    out_file=None,
    n_depths=12,
    hemi=None,
    force=False,
    depths=None,
):
    # sample stat to a number of intermediate surfaces
    stat_file_dir = os.path.dirname(os.path.abspath(stat_file))
    stat_file_base = fsl_remove_ext(os.path.basename(os.path.abspath(stat_file)))

    if out_file is None:
        if hemi is not None:
            out_file = os.path.join(stat_file_dir, stat_file_base + "_" + hemi + ".mgh")
        else:
            out_file = os.path.join(stat_file_dir, stat_file_base + ".mgh")

    if depths is None:
        depths = np.linspace(0, 1, n_depths)

    if not os.path.isfile(out_file) or force == True:
        with tempfile.TemporaryDirectory() as tmpdirname:
            sample_file = os.path.join(tmpdirname, "sampled_depth.mgh")
            for i, depth in enumerate(depths):
                if (
                    subprocess.run(
                        [
                            "mri_vol2surf",
                            "--vol2surf",
                            stat_file,
                            white_surf_file,
                            "0",
                            str(depth),
                            thickness_file,
                            "regheader",
                            "novsm",
                            "5",
                            sample_file,
                        ]
                    ).returncode
                    != 0
                ):
                    return None
                if i == 0:
                    copy2(sample_file, out_file)
                else:
                    if (
                        subprocess.run(
                            [
                                "mris_calc",
                                "--output",
                                out_file,
                                out_file,
                                "add",
                                sample_file,
                            ]
                        ).returncode
                        != 0
                    ):
                        return None
        if (
            subprocess.run(
                ["mris_calc", "--output", out_file, out_file, "div", str(n_depths)]
            ).returncode
            != 0
        ):
            return None
    return out_file


def get_stat_cluster_roi(
    parcel=None,
    stat_file=None,
    analysis_dir=None,
    fs_dir=None,
    fs_to_func_reg=None,
    white_surf_files=None,
    thickness_files=None,
    fwhm=5,
    threshold=2,
    force=False,
):
    stat_cluster_labels = dict()
    for hemi in ["lh", "rh"]:
        # 1. take activation map and project to surface
        stat_surf = sample_surf_func_stat(
            stat_file,
            white_surf_files[hemi],
            thickness_files[hemi],
            hemi=hemi,
            force=force,
        )
        # 2. smooth on surface
        stat_surf_smooth = smooth_surf(
            stat_surf, fs_dir=fs_dir, hemi=hemi, fwhm=fwhm, force=force
        )
        # 3. generate activation clusters
        stat_cluster_labels[hemi] = cluster_surf(
            stat_surf_smooth, fs_dir=fs_dir, hemi=hemi, threshold=threshold, force=force
        )
        # 4. transform cluster label files to volume
        out_basename = fsl_remove_ext(stat_file)
    roi = get_fs_roi(
        parcel,
        stat_cluster_labels,
        out_basename,
        analysis_dir,
        fs_dir,
        fs_to_func_reg,
        force,
    )
    return roi
    # sample stat map to func surface (sample at multiple depths and average)

    # continue like in funcloc_roi:
    # 2. smooth on surface
    # 3. generate cluster (consider addiyional parameters, e.g.  minarea)
    # 4. transform cluster label files to fs volume
    # 5. transform to func space


def get_stat_cluster_atlas(
    hemi,
    stat_file=None,
    analysis_dir=None,
    fs_dir=None,
    fs_to_func_reg=None,
    white_surf_files=None,
    thickness_files=None,
    fwhm=5,
    threshold=2,
    force=False,
    dont_repeat_sample_and_smooth=False,
):

    stat_file_dir = os.path.dirname(os.path.abspath(stat_file))
    stat_file_base = fsl_remove_ext(os.path.basename(os.path.abspath(stat_file)))
    stat_surf_smooth = os.path.join(
        stat_file_dir, stat_file_base + "_" + hemi + "_smooth.mgh"
    )

    if (not dont_repeat_sample_and_smooth) or (not os.path.isfile(stat_surf_smooth)):
        # 1. take activation map and project to surface
        stat_surf = sample_surf_func_stat(
            stat_file,
            white_surf_files[hemi],
            thickness_files[hemi],
            hemi=hemi,
            force=force,
        )
        # 2. smooth on surface
        stat_surf_smooth = smooth_surf(
            stat_surf, fs_dir=fs_dir, hemi=hemi, fwhm=fwhm, force=force
        )

    # 3. generate activation clusters
    stat_cluster_labels = cluster_surf(
        stat_surf_smooth, fs_dir=fs_dir, hemi=hemi, threshold=threshold, force=force
    )
    # 4. transform cluster label files to volume
    out_basename = fsl_remove_ext(stat_file)
    labels_in_fs_individual = fs_overlay_to_fs_volume(
        stat_cluster_labels,
        fs_dir=fs_dir,
        analysis_dir=analysis_dir,
        hemi=hemi,
        out_basename=out_basename,
        force=force,
    )
    labels_in_func = os.path.join(
        analysis_dir, out_basename + "_labels_" + hemi + "_in-func.nii"
    )
    apply_ants_transforms(
        vol_in=labels_in_fs_individual,
        vol_out=labels_in_func,
        ref_vol=fs_to_func_reg[0],
        affine=fs_to_func_reg[1],
        warp=fs_to_func_reg[2],
    )

    return labels_in_func


def get_funcloc_roi(
    parcel=None,
    analysis_dir=None,
    fs_dir=None,
    fs_to_func_reg=None,
    feat_dir=None,
    stat_name="zstat1",
    fwhm=5,
    funcloc_labels=None,
    force=False,
):
    if feat_dir == None:
        feat_dir = os.path.join(analysis_dir, "funcloc.feat")
    stat_file = os.path.join(feat_dir, "stats", stat_name + ".nii.gz")

    if not os.path.isfile(os.path.join(feat_dir, "feat2fs.lta")) or force == True:
        reg_feat_to_fs(feat_dir, fs_dir)

    funcloc_labels = dict()
    for hemi in ["lh", "rh"]:
        # 1. take activation map and project to surface
        stat_surf = sample_surf_feat_stat(
            feat_dir, stat_file, fs_dir, hemi, force=force
        )
        # 2. smooth on surface
        stat_surf_smooth = smooth_surf(
            stat_surf, fs_dir=fs_dir, hemi=hemi, fwhm=fwhm, force=force
        )
        # 3. generate activation clusters
        funcloc_labels[hemi] = cluster_surf(
            stat_surf_smooth, fs_dir=fs_dir, hemi=hemi, force=force
        )
        # 4. transform cluster label files to volume
    out_basename = "funcloc"
    roi = get_fs_roi(
        parcel,
        funcloc_labels,
        out_basename,
        analysis_dir,
        fs_dir,
        fs_to_func_reg,
        force,
    )
    return roi


def get_md_roi(
    parcel=None,
    analysis_dir=None,
    ciftify_dir=None,
    fs_to_func_reg=None,
    md_labels=None,
    force=False,
):
    """Returns a multiple-demand network ROI (Moataz et al. 2020) transformed to functional space"""
    if md_labels == None:
        md_labels = {
            "L": "/data/pt_02389/FinnReplicationPilot/ROIs/MD_L_0.2thresh.label.gii",
            "R": "/data/pt_02389/FinnReplicationPilot/ROIs/MD_R_0.2thresh.label.gii",
        }
    out_basename = "md"
    roi = get_fs_LR_atlas_roi(
        parcel,
        md_labels,
        out_basename,
        analysis_dir,
        ciftify_dir,
        fs_to_func_reg,
        force,
    )
    return roi


def get_glasser_roi(
    parcel=None,
    analysis_dir=None,
    ciftify_dir=None,
    fs_to_func_reg=None,
    glasser_labels=None,
    force=False,
):
    """Returns a HCP MMP 1.0 atlas ROI (Glasser et al. 2016) transformed to functional space"""
    if glasser_labels == None:
        glasser_labels = {
            "L": "/data/pt_02389/FinnReplicationPilot/ROIs/GlasserAtlas.L.164k_fs_LR.label.gii",
            "R": "/data/pt_02389/FinnReplicationPilot/ROIs/GlasserAtlas.R.164k_fs_LR.label.gii",
        }
    out_basename = "glasser"
    roi = get_fs_LR_atlas_roi(
        parcel,
        glasser_labels,
        out_basename,
        analysis_dir,
        ciftify_dir,
        fs_to_func_reg,
        force,
    )
    return roi


def calc_stim_times(onset_delay, trial_duration, trial_order, condition_names=None):
    n = len(trial_order)
    t = np.arange(0, n) * trial_duration + onset_delay
    stim_times = dict()
    for condition in set(trial_order):
        stim_times[condition] = t[np.array(trial_order) == condition]
    return stim_times


def write_stim_time_files(stim_times_runs, cwd=None):
    if cwd == None:
        cwd = os.getcwd()
    # find all conditions from all runs
    conditions = set.union(*[set(stim_times.keys()) for stim_times in stim_times_runs])
    # for each condition create a file and write a line of stim times for each run
    condition_stim_files = []
    for condition in conditions:
        condition_stim_files.append(
            [condition, os.path.join(cwd, "stim-times_" + str(condition) + ".txt")]
        )
        with open(os.path.join(cwd, condition_stim_files[-1][1]), "w") as file:
            for stim_times in stim_times_runs:
                stim_times = defaultdict(list, stim_times)
                print(*stim_times[condition], file=file)
    return condition_stim_files


def average_trials_3ddeconvolve(
    in_files,
    stim_times_runs,
    trial_duration,
    out_files_basename,
    polort=5,
    onset_shift=0,
    cwd=None,
    force=None,
):
    if cwd == None:
        cwd = os.path.dirname(os.path.normpath(in_files[0]))
    n_files = len(in_files)
    # returns estimated impules response components and baseline
    # set parameters
    tr = nib.load(in_files[0]).header.get_zooms()[3]
    n = np.ceil(trial_duration / tr)
    a = 0
    b = tr * (n - 1)
    # prepare stim times and model
    condition_stim_files = write_stim_time_files(stim_times_runs, cwd)
    n_conditions = len(condition_stim_files)

    trialavg_files = []
    for i in range(n_conditions):
        condition = condition_stim_files[i][0]
        trialavg_files.append(
            os.path.join(
                cwd, out_files_basename + f"_response_condition_{condition}.nii"
            )
        )
    baseline_file = os.path.join(cwd, out_files_basename + "_baseline.nii")
    fstat_file = os.path.join(cwd, out_files_basename + "_fstat.nii")

    if (
        all([os.path.isfile(file) for file in trialavg_files])
        and os.path.isfile(baseline_file)
        and os.path.isfile(fstat_file)
        and force == False
    ):
        return trialavg_files, baseline_file, fstat_file

    stim_times = []
    stim_label = []
    i_condition = 0
    for condition, stim_file in condition_stim_files:
        i_condition = i_condition + 1
        stim_times.append((i_condition, stim_file, f"TENT({a},{b},{n})"))
        stim_label.append((i_condition, str(condition)))

    deconvolve = Deconvolve()
    deconvolve.inputs.in_files = in_files
    deconvolve.inputs.stim_times = stim_times
    deconvolve.inputs.stim_label = stim_label
    deconvolve.inputs.polort = polort
    deconvolve.inputs.local_times = True
    deconvolve.inputs.fout = True
    deconvolve.inputs.cbucket = os.path.join(
        cwd, out_files_basename + "_cbucket.nii.gz"
    )
    deconvolve.inputs.args = "-overwrite"
    deconvolve.inputs.stim_times_subtract = onset_shift
    result = deconvolve.run(cwd=cwd)
    # extract fstat
    result_fstat = TCatSubBrick(
        in_files=[(result.outputs.out_file, f"'[0]'")],
        out_file=os.path.join(cwd, out_files_basename + "_fstat.nii"),
        args="-overwrite",
    ).run()
    # add back baseline
    baseline_idcs = 0 + (polort + 1) * np.arange(n_files)
    baseline_idcs_str = ",".join([str(i) for i in baseline_idcs])
    result_baseline_vols = TCatSubBrick(
        in_files=[(result.outputs.cbucket, f"'[{baseline_idcs_str}]'")],
        out_file=os.path.join(cwd, out_files_basename + "_baseline_runs.nii"),
        args="-overwrite",
    ).run()
    result_baseline = TStat(
        in_file=os.path.join(cwd, out_files_basename + "_baseline_runs.nii"),
        args="-mean -overwrite",
        out_file=os.path.join(cwd, out_files_basename + "_baseline.nii"),
    ).run()
    for i in range(n_conditions):
        condition = condition_stim_files[i][0]
        result_condition_diffresponse_timecourse = TCatSubBrick(
            in_files=[
                (
                    result.outputs.cbucket,
                    f"'[{int((polort + 1) * n_files + i * n)}..{int((polort + 1) * n_files + (i + 1) * n - 1)}]'",
                )
            ],
            out_file=os.path.join(
                cwd, out_files_basename + f"_diffresponse_condition_{condition}.nii"
            ),
            args="-overwrite",
        ).run()
        result_condition_response_timecourse = Calc(
            in_file_a=os.path.join(
                cwd, out_files_basename + f"_diffresponse_condition_{condition}.nii"
            ),
            in_file_b=os.path.join(cwd, out_files_basename + "_baseline.nii"),
            out_file=os.path.join(
                cwd, out_files_basename + f"_response_condition_{condition}.nii"
            ),
            expr="a+b",
            args="-overwrite",
        ).run()

    baseline_file = os.path.join(cwd, out_files_basename + "_baseline.nii")
    fstat_file = os.path.join(cwd, out_files_basename + "_fstat.nii")
    return trialavg_files, baseline_file, fstat_file


def calc_percent_change_trialavg(
    trialavg_files,
    baseline_file,
    inv_change=False,
    force=False,
):
    if inv_change:
        expr = "100-(100*a/b)"
    else:
        expr = "100*a/b-100"
    prc_change = []
    for trialavg_file in trialavg_files:
        trialavg_file_split = os.path.splitext(trialavg_file)
        out_file = trialavg_file_split[0] + "_prcchg" + trialavg_file_split[1]
        if not os.path.isfile(out_file) or force == True:
            result_prcchg = Calc(
                in_file_a=trialavg_file,
                in_file_b=baseline_file,
                out_file=out_file,
                expr=expr,
                args="-overwrite",
            ).run()
            prc_change.append(result_prcchg.outputs.out_file)
        else:
            prc_change.append(out_file)
    return prc_change


def calc_normalized_trialavg(
    trialavg_files,
    baseline_vols=None,
    force=False,
):
    """
    Normalizes trial average by adding back baseline and dividing by mean from defined baseline volumes
    """
    expr = "(img1/img2[:,:,:,None])"
    normalized_imgs = []
    for trialavg_file in trialavg_files:
        trialavg_file_split = os.path.splitext(trialavg_file)
        out_file = trialavg_file_split[0] + "_norm" + trialavg_file_split[1]
        if not os.path.isfile(out_file) or force == True:
            baseline = mean_img(index_img(trialavg_file, baseline_vols))
            math_img(
                expr,
                img1=trialavg_file,
                img2=baseline,
            ).to_filename(out_file)
        normalized_imgs.append(out_file)
    return normalized_imgs


def plot_cond_tcrs(
    condition_data_list,
    t=None,
    TR=1,
    labels=None,
    colors=None,
    ax=None,
    periods=None,
    events=None,
):
    """Plots time course for multiple conditions. All timecourses should have the same length."""
    if ax == None:
        ax = plt.axes()
    if t == None:
        N = len(condition_data_list[0])
        t = TR * np.arange(0, N)
    for period in periods or []:
        ax.axvspan(period[0], period[1], alpha=0.1, color="gray")
    ax.axhline(0, color="gray", lw=0.5)
    line_handles = []
    for condition_data, color in zip_longest(
        condition_data_list, (colors or []), fillvalue=None
    ):
        (l,) = ax.plot(t, condition_data, color=color)
        line_handles.append(l)
    ax.set_xticks(t)
    ax.axis()
    ax.set_xlim(min(t), max(t))
    y0, y1 = ax.get_ylim()
    if labels:
        ax.legend(line_handles, labels, loc="best")
    for event in events or []:
        ax.axvline(event[1], color="gray", lw=0.5)
        ax.annotate(event[0], (event[1], y0), ha="center", va="bottom")

    # REMOVE ME LATER:
    # ax.axis([0,30,-0.4,1.2])
    return ax


def add_prefix_to_nifti_basename(path, prefix):
    norm_path = os.path.normpath(path)
    dir_name = os.path.dirname(norm_path)
    base_name = os.path.basename(norm_path)
    return os.path.join(dir_name, prefix + base_name)


def add_postfix_to_nifti_basename(path, postfix):
    norm_path = os.path.normpath(path)
    s = os.path.splitext(norm_path)
    if s[1] == ".gz":
        s = os.path.splitext(s[0])
        basename = s[0]
        extension = s[1] + ".gz"
    else:
        basename = s[0]
        extension = s[1]
    return basename + postfix + extension


def get_funcact_roi_laynii(
    act_file, rim_file, roi_out_file, n_columns=10000, threshold=1
):
    columns_file = add_postfix_to_nifti_basename(rim_file, "_columns" + str(n_columns))
    if not os.path.isfile(columns_file):
        mid_gm_file = add_postfix_to_nifti_basename(rim_file, "_midGM_equidist")
        subprocess.run(
            [
                "LN2_COLUMNS",
                "-rim",
                rim_file,
                "-midgm",
                mid_gm_file,
                "-nr_columns",
                str(n_columns),
            ]
        )
    subprocess.run(
        [
            "LN2_MASK",
            "-scores",
            act_file,
            "-columns",
            columns_file,
            "-min_thr",
            str(threshold),
            "-output",
            roi_out_file,
        ]
    )
    subprocess.run(["fslmaths", roi_out_file, "-bin", roi_out_file])
    return roi_out_file


def get_funcact_roi_vfs(act_file, columns_file, roi_out_file, threshold=1):
    scores = nib.load(act_file).get_fdata()

    nii_columns = nib.load(columns_file)
    columns = nii_columns.get_fdata()
    mask = np.zeros(nii_columns.shape)
    scores_thr = scores >= threshold

    act_columns = np.unique(columns[scores_thr])
    for column_idx in act_columns:
        if column_idx != 0:
            mask[columns == column_idx] = 1

    xform = nii_columns.affine
    mask_nii = nib.nifti2.Nifti2Image(mask, xform)
    nib.save(mask_nii, roi_out_file)
    return mask_nii


def roi_and(rois):
    # TODO: check out if using mask_image on rois is correct for ignoring affine.
    rois = [mask_image(rois[0], roi) for roi in rois]
    return intersect_masks(rois, threshold=1, connected=False)


def roi_or(rois):
    # TODO: check whether to use mask_image on rois for ignoring affine (see roi_and)
    return intersect_masks(rois, threshold=0, connected=False)


def mask_image(img, mask):
    if type(img) is str:
        img = nib.load(img)
    if type(mask) is str:
        mask = nib.load(mask)

    img_data = img.get_fdata()
    mask_data = mask.get_fdata()

    img_masked_data = img_data * (mask_data != 0)

    return nib.Nifti1Image(img_masked_data, img.affine, img.header)


def bold_correct(
    nulled_file, notnulled_file, out_file, notnulled_shift=None, force=None
):
    """notnulled_shift should equal (positive) difference between readout blocks"""
    if not os.path.isfile(out_file) or force == True:
        if notnulled_shift is not None:
            slicetimer_result = SliceTimer(
                in_file=notnulled_file, global_shift=-notnulled_shift
            ).run()
            notnulled_file = slicetimer_result.outputs.out_file

        my_env = os.environ.copy()

        if os.path.normpath(out_file)[-3:] == "nii":
            my_env["FSLOUTPUTTYPE"] = "NIFTI"
        elif os.path.normpath(out_file)[-6:] == "nii.gz":
            my_env["FSLOUTPUTTYPE"] = "NIFTI_GZ"

        subprocess.run(
            [
                "fslmaths",
                nulled_file,
                "-div",
                notnulled_file,
                "-max",
                "0",
                "-min",
                "5",
                out_file,
            ],
            env=my_env,
        )
    return out_file


def sample_timecourse(func_filename, roi):
    # assume voxel matrix of data, roi correspond to each other
    func_filename_reset = reset_affine(func_filename)
    roi_reset = reset_affine(roi)

    masked_data = apply_mask(func_filename_reset, roi_reset)
    return masked_data


def calc_layers_laynii(
    rim_file, out_file_base=None, method="equidist", n_layers=3, force=False
):
    # include upsampling methods?
    if out_file_base is None:
        out_file_base = fsl_remove_ext(rim_file)
    if method == "equivol":
        out_file = out_file_base + "_metric_equivol.nii"
    else:
        out_file = out_file_base + "_metric_equidist.nii"

    if not os.path.isfile(out_file) or force == True:
        run_string_list = [
            "LN2_LAYERS",
            "-rim",
            rim_file,
            "-output",
            out_file_base,
            "-nr_layers",
            str(n_layers),
        ]
        if method == "equivol":
            run_string_list.append("-equivol")
        subprocess.run(run_string_list)

    return out_file


def generate_two_layers(analysis_dir, depths, delta=0, roi=None):
    test = intersect_masks([depths, roi])
    superficial = math_img(f"img<{0.5 - delta / 2}", img=depths)
    deeper = math_img(f"img>{0.5 - delta / 2}", img=depths)
    if roi is not None:
        # superficial.affine = roi.affine
        superficial = intersect_masks([superficial, roi], threshold=1)
        deeper = intersect_masks([deeper, roi], threshold=1)
    nib.save(superficial, os.path.join(analysis_dir, "superficial.nii"))
    nib.save(deeper, os.path.join(analysis_dir, "deeper.nii"))


def reset_affine(img):
    img_reset = copy.deepcopy(check_niimg(img, dtype="auto"))
    img_reset.set_qform(np.eye(4))
    img_reset.set_sform(np.eye(4))
    return img_reset


def sample_roi(data, roi):
    # assume voxel matrix of data, roi correspond to each other
    data_reset = reset_affine(data)
    roi_reset = reset_affine(roi)

    masked_data = apply_mask(data_reset, roi_reset)
    return masked_data


def average_roi(data, roi):
    # assume voxel matrix of data, roi correspond to each other
    data_reset = reset_affine(data)
    roi_reset = reset_affine(roi)

    masked_data = sample_roi(data_reset, roi_reset)
    return np.mean(masked_data), np.std(masked_data), np.size(masked_data)


def sample_depths(data, roi, depths):
    # assume voxel matrix of data, roi correspond to each other
    data_reset = reset_affine(data)
    roi_reset = reset_affine(roi)
    depths_reset = reset_affine(depths)

    masked_data = apply_mask(data_reset, roi_reset)
    masked_depths = apply_mask(depths_reset, roi_reset)
    return masked_data, masked_depths


def sample_layer_profile(data, roi, depths, n_layers):
    data, depths = sample_depths(data, roi, depths)
    depths = np.floor(depths * n_layers)
    y = np.zeros(n_layers)
    for i in np.arange(n_layers):
        y[i] = np.mean(data[depths == i])
    return y, (np.arange(n_layers) + 0.5) / n_layers


def plot_profiles(data_list, roi, depths, n_layers, colors=None, labels=None):
    ax = plt.axes()
    line_handles = []
    for i, data in enumerate(data_list):
        voxel_responses, voxel_depths = sample_depths(data, roi, depths)
        if colors is not None:
            color = colors[i]
        else:
            prop_cycle = plt.rcParams["axes.prop_cycle"]
            color = prop_cycle.by_key()["color"][i]
        layer_responses, layer_depths = sample_layer_profile(
            data, roi, depths, n_layers
        )
        ax.plot(1 - voxel_depths, voxel_responses, ".", alpha=0.2, color=color)
        (l,) = ax.plot(1 - layer_depths, layer_responses, color=color, lw=2)
        line_handles.append(l)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CSF|GM", "GM|WM"])
        if labels:
            ax.legend(line_handles, labels, loc="best")
    return ax


def upsample(in_file, out_file, factor, method):
    voxel_widths = np.array(nib.load(in_file).header.get_zooms())
    scaled_voxel_widths = voxel_widths / factor
    subprocess.run(
        [
            "3dresample",
            "-dxyz",
            str(scaled_voxel_widths[0]),
            str(scaled_voxel_widths[1]),
            str(scaled_voxel_widths[2]),
            "-rmode",
            method,
            "-overwrite",
            "-prefix",
            out_file,
            "-input",
            in_file,
        ]
    )


def get_labels_data(
    data_file,
    labels_file,
    print_results=False,
    label_names=None,
    mask=None,
    layers=None,
):
    labels_data = nib.load(labels_file).get_fdata()
    labels = np.unique(labels_data)
    df = pd.DataFrame()
    for label in labels:
        if label != 0:
            if label_names is None:
                label_name = str(int(label))
            else:
                label_name = label_names[int(label) - 1]
            roi = index_roi(labels_file, label)
            if mask is not None:
                roi = roi_and((roi, mask))
            if np.any(roi.get_fdata()):
                d = sample_roi(data_file, roi)
                if layers is not None:
                    l = sample_roi(layers, roi)
                    df = df.append(
                        pd.DataFrame({"value": d, "layer": l}).assign(label=label_name)
                    )
                else:
                    df = df.append(pd.DataFrame({"value": d}).assign(label=label_name))
                if print_results:
                    m, s, n = average_roi(data_file, roi)
                    print(
                        f"{label_name}: {m:2.2f} +- {s / np.sqrt(n):2.2f} (n={int(n)})"
                    )
    return df


def get_labels_data_layers(
    data_file, labels_file, print_results=False, label_names=None
):
    pass


def get_labels_data_layers_masked(
    data_file, labels_file, print_results=False, label_names=None
):
    pass


#import matplotlib.pyplot as plt
#import nibabel as nib
#import nilearn.plotting as plotting
#import numpy as np
#import hcp_utils as hcp


def plot_on_mmhcp_surface(Xp):
    """Xp is a 1D Vector same size as hcp.mmp.labels."""
    mmp_labels = hcp.mmp.labels  # mmp = Glasser parcellation

    cm = "cold_hot"
    min_thresh = 0
    max_thresh = 0.1

    # 2D plot – I also detail here with an example how you can add subplots…
    fig = plt.figure(figsize=[20, 10])
    ax = fig.add_subplot(1, 4, 1, projection="3d")
    plotting.plot_surf_stat_map(
        hcp.mesh.inflated,
        hcp.cortex_data(hcp.unparcellate(Xp, hcp.mmp)),
        view="anterior",
        colorbar=True,
        threshold=min_thresh,
        vmax=max_thresh,
        bg_map=hcp.mesh.sulc,
        bg_on_data=True,
        darkness=0.3,
        axes=ax,
        figure=fig,
        cmap=cm,
        symmetric_cbar=True,
    )

    ax = fig.add_subplot(1, 4, 2, projection="3d")
    plotting.plot_surf_stat_map(
        hcp.mesh.inflated,
        hcp.cortex_data(hcp.unparcellate(Xp, hcp.mmp)),
        view="lateral",
        colorbar=True,
        threshold=min_thresh,
        vmax=max_thresh,
        bg_map=hcp.mesh.sulc,
        bg_on_data=True,
        darkness=0.3,
        axes=ax,
        figure=fig,
        cmap=cm,
        symmetric_cbar=True,
    )

    fig.suptitle("title", fontsize=16)
    plt.savefig("output.png", facecolor="white")
    plt.close()

    # or alternatively, you can get this in 3D and interact with it in html (on the cluster you need to use the Chrome browser)
    nn = plotting.view_surf(
        hcp.mesh.inflated,
        hcp.cortex_data(hcp.unparcellate(Xp, hcp.mmp)),
        bg_map=hcp.mesh.sulc,
        symmetric_cmap=False,
        vmax=max_thresh,
        vmin=min_thresh,
        title="title",
    )
    nn.save_as_html("output.html")


### FinnReplicationPilot specifc functions


def plot_finn_panel(depths, roi, trialavg_data, run_type, layers, ax, d=0, TR=3.702):
    condition_data = []
    for file in trialavg_data:
        sampled_data, sampled_depths = sample_depths(file, roi, depths)
        if layers == "deep":
            condition_data.append(
                np.mean(sampled_data[:, sampled_depths < (0.5 - d)], axis=1)
            )
        elif layers == "superficial":
            condition_data.append(
                np.mean(sampled_data[:, sampled_depths > (0.5 - d)], axis=1)
            )

    if run_type == "alpharem":
        labels = ["rem", "alpha"]
        colors = ["tab:green", "tab:blue"]
    elif run_type == "gonogo":
        labels = ["nogo", "go"]
        colors = ["tab:orange", "tab:red"]
    else:
        return None

    plot_cond_tcrs(
        condition_data,
        TR=TR,
        labels=labels,
        colors=colors,
        periods=[[4, 14], [14, 20]],
        events=[["Stim", 0], ["Cue", 4], ["Probe", 14]],
        ax=ax,
    )
    plt.title(layers)
    return condition_data


def plot_finn_tcrses(
    depths, roi, trialavg_alpharem, trialavg_gonogo, d=0, TR=3.702, ymin=-0.5, ymax=1.5
):
    if type(depths) == str:
        depths = nib.load(depths)

    fig = plt.figure(figsize=(10, 8), dpi=100, facecolor="w", edgecolor="k")
    ax = plt.subplot(2, 2, 1)
    plot_finn_panel(
        depths, roi, trialavg_alpharem, "alpharem", "superficial", ax, d, TR=TR
    )
    ax.axis([0, 30, ymin, ymax])

    ax = plt.subplot(2, 2, 2)
    plot_finn_panel(depths, roi, trialavg_gonogo, "gonogo", "superficial", ax, d, TR=TR)
    ax.axis([0, 30, ymin, ymax])

    ax = plt.subplot(2, 2, 3)
    plot_finn_panel(depths, roi, trialavg_alpharem, "alpharem", "deep", ax, d, TR=TR)
    ax.axis([0, 30, ymin, ymax])

    ax = plt.subplot(2, 2, 4)
    plot_finn_panel(depths, roi, trialavg_gonogo, "gonogo", "deep", ax, d, TR=TR)
    ax.axis([0, 30, ymin, ymax])


# def get_finn_tcrs_data(trial_averages, roi, layers=None, conditions=None, modalities=None):
#     data = dict()
#     for layer in layers:
#         for condition in conditions:
#             for modality in modalities:
#                 data[modality, layer, condition] = np.loadtxt(
#                     fnamebase + modality + '_' + layer + '_' + condition + '.1D')
#     return data


# def plot_finn_tcrs(fnamebase, modality, TR=3.702):
#     data = get_finn_tcrs_data(fnamebase)
#     fix, axes = plt.subplot(2, 3, figsize=(15, 5))
#     periods = [[4, 14], [14, 20]]
#     events = [['Stim', 0], ['Cue', 4], ['Probe', 14]]
#     for row, layer in enumerate(layers):
#         plot_cond_tcrs([data[modality, layer, 'alpha'],
#                         data[modality, layer, 'rem']],
#                        colors=('tab:blue', 'tab:green'),
#                        labels=('alpha', 'rem'),
#                        periods=periods,
#                        events=events,
#                        TR=TR,
#                        ax=axes[row, 0])
#
#         plot_cond_tcrs([data[modality, layer, 'go'],
#                         data[modality, layer, 'nogo']],
#                        colors=('tab:red', 'tab:orange'),
#                        labels=('act', 'non-act'),
#                        periods=periods,
#                        events=events,
#                        TR=TR,
#                        ax=axes[row, 0])
#
#         plot_cond_tcrs([data[modality, layer, 'alpha'] - data[modality, layer, 'rem'],
#                         data[modality, layer, 'go'] - data[modality, layer, 'nogo']],
#                        colors=('tab:purple', 'tab:cyan'),
#                        labels=('alpha - rem', 'act - non-act'),
#                        periods=periods,
#                        events=events,
#                        TR=TR,
#                        ax=axes[row, 0])
#     axes[0, 0].set_ylabel('signal change [%]')
#     axes[1, 0].set_ylabel('signal change [%]')
#     axes[1, 0].set_xlabel('trial time [s]')
#     axes[1, 1].set_xlabel('trial time [s]')
#     axes[1, 2].set_xlabel('trial time [s]')
#
#     fig.text(0.08, 0.45, 'deeper', ha='right', weight='bold')
#     fig.text(0.08, 0.85, 'superficial', ha='right', weight='bold')
#     fig.suptitle(modality.upper(), weight='bold')


def finn_trial_averaging(run_type, analysis_dir, out_dir=None, force=False):
    trial_duration = 32
    trial_order = paradigm(run_type)
    trialavg = dict()
    onset_delay = 8
    in_files_bold = [
        os.path.join(analysis_dir, f"func_{run_type}_notnulled_tshift.nii")
    ]
    in_files_vaso = [os.path.join(analysis_dir, f"func_{run_type}_vaso.nii")]

    if out_dir is None:
        out_dir = analysis_dir

    stim_times_runs = [
        calc_stim_times(
            onset_delay=8, trial_duration=trial_duration, trial_order=trial_order
        )
    ]
    (
        trialavg_files_bold,
        baseline_file_bold,
        fstat_file_bold,
    ) = average_trials_3ddeconvolve(
        in_files_bold,
        stim_times_runs,
        trial_duration,
        out_files_basename="trialavg1_bold_" + run_type,
        polort=5,
        cwd=out_dir,
        force=force,
    )

    (
        trialavg_files_vaso,
        baseline_file_vaso,
        fstat_file_vaso,
    ) = average_trials_3ddeconvolve(
        in_files_vaso,
        stim_times_runs,
        trial_duration,
        out_files_basename="trialavg1_vaso_" + run_type,
        polort=5,
        cwd=out_dir,
        force=force,
    )

    trialavg_bold_prcchg = calc_percent_change_trialavg(
        trialavg_files_bold, baseline_file_bold, inv_change=False, force=force
    )
    trialavg_vaso_prcchg = calc_percent_change_trialavg(
        trialavg_files_vaso, baseline_file_vaso, inv_change=True, force=force
    )

    return trialavg_bold_prcchg, trialavg_vaso_prcchg, fstat_file_bold, fstat_file_vaso


def finn_trial_averaging_with_boldcorrect(
    run_type, analysis_dir, TR1, out_dir=None, force=False
):
    trial_duration = 32
    trial_order = paradigm(run_type)
    trialavg = dict()
    onset_delay = 8
    in_files_nulled = [os.path.join(analysis_dir, f"func_{run_type}_nulled.nii")]
    in_files_notnulled = [os.path.join(analysis_dir, f"func_{run_type}_notnulled.nii")]
    stim_times_runs = [
        calc_stim_times(
            onset_delay=8, trial_duration=trial_duration, trial_order=trial_order
        )
    ]

    (
        trialavg_files_nulled,
        baseline_file_nulled,
        fstat_file_nulled,
    ) = average_trials_3ddeconvolve(
        in_files_nulled,
        stim_times_runs,
        trial_duration,
        out_files_basename="trialavg2_nulled_" + run_type,
        polort=5,
        force=force,
        cwd=out_dir,
    )

    (
        trialavg_files_notnulled,
        baseline_file_notnulled,
        fstat_file_notnulled,
    ) = average_trials_3ddeconvolve(
        in_files_notnulled,
        stim_times_runs,
        trial_duration,
        out_files_basename="trialavg2_notnulled_" + run_type,
        polort=5,
        onset_shift=TR1,
        cwd=out_dir,
        force=force,
    )

    trialavg_files_vaso = [
        bold_correct(
            trialavg_files_nulled[0],
            trialavg_files_notnulled[0],
            trialavg_files_nulled[0].replace("nulled", "vaso"),
            force=force,
        ),
        bold_correct(
            trialavg_files_nulled[1],
            trialavg_files_notnulled[1],
            trialavg_files_nulled[1].replace("nulled", "vaso"),
            force=force,
        ),
    ]

    baseline_file_vaso = bold_correct(
        baseline_file_nulled,
        baseline_file_notnulled,
        baseline_file_nulled.replace("nulled", "vaso"),
        force=force,
    )

    trialavg_bold_prcchg = calc_percent_change_trialavg(
        trialavg_files_notnulled, baseline_file_notnulled, inv_change=False, force=force
    )

    trialavg_vaso_prcchg = calc_percent_change_trialavg(
        trialavg_files_vaso, baseline_file_vaso, inv_change=True, force=force
    )

    fstat_file_bold = fstat_file_notnulled

    return (
        trialavg_bold_prcchg,
        trialavg_vaso_prcchg,
        fstat_file_bold,
        fstat_file_nulled,
    )


def finn_trial_averaging_on_renzo_boldcorrect(
    run_type, analysis_dir, out_dir=None, force=False
):
    trial_duration = 32
    trial_order = paradigm(run_type)
    trialavg = dict()
    onset_delay = 8
    in_files_bold = [os.path.join(analysis_dir, f"func_{run_type}_notnulled.nii")]
    in_files_vaso = [os.path.join(analysis_dir, f"func_{run_type}_rvaso.nii")]
    stim_times_runs = [
        calc_stim_times(
            onset_delay=8, trial_duration=trial_duration, trial_order=trial_order
        )
    ]
    (
        trialavg_files_bold,
        baseline_file_bold,
        fstat_file_bold,
    ) = average_trials_3ddeconvolve(
        in_files_bold,
        stim_times_runs,
        trial_duration,
        out_files_basename="trialavg3_bold_" + run_type,
        polort=5,
        cwd=out_dir,
        force=force,
    )

    (
        trialavg_files_vaso,
        baseline_file_vaso,
        fstat_file_vaso,
    ) = average_trials_3ddeconvolve(
        in_files_vaso,
        stim_times_runs,
        trial_duration,
        out_files_basename="trialavg3_vaso_" + run_type,
        polort=5,
        cwd=out_dir,
        force=force,
    )

    trialavg_bold_prcchg = calc_percent_change_trialavg(
        trialavg_files_bold, baseline_file_bold, inv_change=False, force=force
    )
    trialavg_vaso_prcchg = calc_percent_change_trialavg(
        trialavg_files_vaso, baseline_file_vaso, inv_change=True, force=force
    )

    return trialavg_bold_prcchg, trialavg_vaso_prcchg, fstat_file_bold, fstat_file_vaso


def finn_trial_averaging_on_boldcorrect_finn_baselining(
    run_type, analysis_dir, out_dir=None, force=False
):
    trial_duration = 32
    trial_order = paradigm(run_type)
    trialavg = dict()
    onset_delay = 8
    in_files_bold = [
        os.path.join(analysis_dir, f"func_{run_type}_notnulled_tshift.nii")
    ]
    in_files_vaso = [os.path.join(analysis_dir, f"func_{run_type}_vaso.nii")]
    stim_times_runs = [
        calc_stim_times(
            onset_delay=8, trial_duration=trial_duration, trial_order=trial_order
        )
    ]

    tr = nib.load(in_files_bold[0]).header.get_zooms()[3]
    trial_average_duration = tr * 12  # average 12 vols = 44.3s > 32s (trial length)
    baseline_vols = [6, 7, 8, 9, 10, 11]
    (
        trialavg_files_bold,
        baseline_file_bold,
        fstat_file_bold,
    ) = average_trials_3ddeconvolve(
        in_files_bold,
        stim_times_runs,
        trial_average_duration,
        out_files_basename="trialavg4_bold_" + run_type,
        polort=0,  # no detrending, just averaging
        onset_shift=tr,  # start average 1TR before stim onset
        cwd=out_dir,
        force=force,
    )

    (
        trialavg_files_vaso,
        baseline_file_vaso,
        fstat_file_vaso,
    ) = average_trials_3ddeconvolve(
        in_files_vaso,
        stim_times_runs,
        trial_average_duration,
        out_files_basename="trialavg4_vaso_" + run_type,
        polort=0,  # no detrending, just averaging
        onset_shift=tr,  # start average 1TR before stim onset
        cwd=out_dir,
        force=force,
    )

    trialavg_bold_norm = calc_normalized_trialavg(
        trialavg_files_bold, baseline_vols=baseline_vols, force=force
    )
    trialavg_vaso_norm = calc_normalized_trialavg(
        trialavg_files_vaso, baseline_vols=baseline_vols, force=force
    )

    return trialavg_bold_norm, trialavg_vaso_norm, fstat_file_bold, fstat_file_vaso


# Define a function to obtain some paradigm related info. (For now trial order, TODO: trial period timings, GLM events, ...)
# we need:
# for trial averaging: exact trial onset times (what about VASO,GE-BOLD shifts?) with conditions
# + volume infor for period averagin?
# for glm analysis
# onsets and durations of all trial periods
def paradigm(run_type):
    conditionRem = 2
    conditionAlpha = 3
    conditionNogo = 4
    conditionGo = 5
    if run_type == "localizer":
        letterStringDuration = 2.5
        fix1Duration = 1.5
        cueDuration = 1
        fixDelayDuration = 9
        probeDuration = 1
        interTrialDuration = 5

        startBlankPeriod = 6

        trial_order = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
    elif run_type in ["alpharem", "gonogo"]:
        letterStringDuration = 2.5
        fix1Duration = 1.5
        cueDuration = 1
        fixDelayDuration = 9
        probeDuration = 2
        interTrialDuration = 16

        startBlankPeriod = 8

        if run_type == "alpharem":
            trial_order = [2, 3, 3, 3, 2, 2, 3, 2, 3, 3, 2, 2, 2, 3, 2, 3, 3, 3, 2, 2]
        elif run_type == "gonogo":
            trial_order = [4, 5, 5, 5, 4, 4, 5, 4, 5, 5, 4, 4, 4, 5, 4, 5, 5, 5, 4, 4]
    else:
        return None

    return trial_order


### TODO or obsolete:


def preprocess_funcloc(data):
    # motion correction

    # spatial smoothing

    # high pass filtering

    # register to freesurfer

    pass


def feat_analysis(
    feat_template,
    func_file,
    output_dir,
    stim_timings_dir,
    smoothing_fwhm=0,
    overwrite=False,
):
    cwd = os.path.dirname(os.path.normpath(output_dir))
    feat_template_base = os.path.basename(os.path.normpath(feat_template))

    if overwrite == True:
        rmtree(output_dir)

    # copy
    subprocess.run(["cp", feat_template, "."], cwd=cwd)
    # edit
    subprocess.run(
        [
            "sed",
            "-i",
            "-e",
            f"s|templateVar_FuncFile|{func_file}|g",
            feat_template_base,
        ],
        cwd=cwd,
    )
    subprocess.run(
        [
            "sed",
            "-i",
            "-e",
            f"s|templateVar_OutputDir|{output_dir}|g",
            feat_template_base,
        ],
        cwd=cwd,
    )
    subprocess.run(
        [
            "sed",
            "-i",
            "-e",
            f"s|templateVar_StimTimingsDir|{stim_timings_dir}|g",
            feat_template_base,
        ],
        cwd=cwd,
    )
    subprocess.run(
        [
            "sed",
            "-i",
            "-e",
            f"s|templateVar_SmoothingFWHM|{smoothing_fwhm}|g",
            feat_template_base,
        ],
        cwd=cwd,
    )
    # run feat
    subprocess.run(["feat", feat_template_base], cwd=cwd)
    return output_dir


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


def plot_timecourses():
    pass


def normalize():
    # Normalizes single voxel timecourses to change relative to a baselne
    pass


def avg_timcourse_nearest_volume():
    pass


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
