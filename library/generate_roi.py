#! /usr/bin/env python3
""" 
Implements function find_roi and corresponding nipype interface FindRoi and command line interface. 
"""
import sys
import os
from scipy.optimize import bisect
from nipype.interfaces.workbench import MetricResample
from tempfile import TemporaryDirectory
import subprocess
sys.path.append("/data/p_02389/code/fmri-analysis/library/")
import layer_analysis as analysis

def find_roi(
    stat_file,
    anat_region,
    out_file,
    target,
    target_type,
    white_surf,
    pial_surf,
    hemi,
    anat_require_region=None,
    cluster=True,
    fwhm=3,
    cwd=None,
):
    """
    Finds an roi that most closely meets criteria based on predefined regions (e.g. atlas
    based and on activation maps.

    required inputs:

    stat_file           - volume file (TODO: also accept surface files?)
    anat_region         - surface mask
    out_file            - name of generated roi volume file
    target              - quantity of target measure for activation cluster/vertices/voxels
    target_type         - target measure, can be one of 'area','percentile','threshold',('voxels'?)
    white_surf          - white surface as .gii in native volume space
    pial_surf           - pial surface as .gii in native volume space
    hemi                - hemisphere 'L' or 'R'

    optional:

    anat_require_region - surface mask or None
    cluster             - if True use extract largest consecutive cluster of voxels
    fwhm                - smoothing kernel width (only used together with cluster)
    cwd                 - current working directory (derived from output file if None)
    """

    # set cwd to directory of output file, if not set otherwise
    if cwd == None:
        cwd = os.path.dirname(os.path.normpath(out_file))

    with TemporaryDirectory(dir=cwd) as tmpdirname:
        # 1. sample stat_file to surface, use wb_command (surface neeeded as input)
        mid_surf = os.path.join(tmpdirname, "mid.surf.gii")
        surf_stat_file = os.path.join(tmpdirname, "stat.func.gii")
        analysis.sample_surf_hcp(
            volume_file=stat_file,
            white_surf=white_surf,
            pial_surf=pial_surf,
            mid_surf=mid_surf,
            outfile=out_file,
            mask_file=None)

        # 2. if cluster: smooth surface
        surf_stat_smooth_file = os.path.join(tmpdirname, "stat.smooth.func.gii")
        if cluster:
            analysis.smooth_surfmetric_hcp(surf_stat_file,surf_stat_smooth_file,mid_surf,fwhm)

        # 3. calculate total area; if percentile, calculate area target
        A_total = analysis.calc_area_hcp(anat_region, mid_surf)
        if target_type == 'percentile':
            target = (percentile/100)*A_total
            target_type = 'area'
        
        # 4. use optimization function, that thresholds surface
        if cluster:
            # either by extracting largest above threshold cluster
            def calc_roi(threshold):
                roi_surf_file = os.path.join(tmpdirname,f'roi.{hash(threshold)}.shape.gii}')
                if not os.path.exists(roi_surf_file):
                    min_area = A_total/2

                    # find all cluster (use anat_region as ROI!!!)
                    clusters_file = os.path.join(tmpdirname, "stat.clusters.func.gii")
                    subprocess.run(
                        [
                            "wb_command",
                            "-metric-find-clusters",
                            mid_surf,
                            surf_stat_smooth_file,
                            str(threshold),
                            str(min_area),
                            clusters_file,
                        ],
                        check=True,
                    )
                    # - apply require_region as mask to cluster_file
                    # - load cluster file and get list of cluster numbers
                    # - for all cluster numbers extract cluster roi and calculate area
                    # - select largest cluster and return

                    
                    subprocess.run(
                        [
                            "wb_command",
                            "-metric-remove-islands",
                            mid_surf,
                            cluster_file,
                            roi_surf_file,
                        ],
                        check=True,
                    )
                return roi_file

        else:
            # or directly
            def calc_roi(threshold):
                roi_surf_file = os.path.join(tmpdirname,f'roi.{hash(threshold)}.shape.gii}')
                if not os.path.exists(roi_surf_file):
                    analysis.math_metric(
                        f"x > {threshold}",
                        metric_out=roi_surf_file,
                        x=surf_stat_file
                    )
                return roi_surf_file

        if target=='threshold':
            return calc_roi(threshold)
        else:
            if target=='area':
                max_area = 
                def f(threshold):
                    roi_file = calc_roi(threshold)
                    return target - area
            elif target=='percentile':
                def f(threshold):
                    roi_file = calc_roi(threshold)
                    return target - percentile
            elif target=='voxels':
                def f(threshold):
                    roi_file = calc_roi(threshold)
            b = MAXSTAT_IN_REGION
            threshold = bisect(f, 0, max_stat, xtol=0.01)

    # 4. return roi
    roi = None
    return roi
 

def test_find_roi():
    find_roi(
        stat_file="data/fstat.nii",
        anat_region="data/dlpfc_mask.nii",
        require_anat_region=None,
        out_file="data/roi.nii",
        target="100",
        target_type="voxels",
        cluster=True,
    )


class FindROIInputSpec(BaseInterfaceInputSpec):
    pass


class FinROIOutputSpec(TraitedSpec):
    pass


class FindROI:
    pass


if __name__ == "__main__":
    # process command line arguments and either call interface or implemented function
    pass
