#! /usr/bin/env python3
""" 
Implements function find_layer_contrast_roi and corresponding nipype interface FindLayerContrastRoi and command line interface. 
"""
import sys
import os
from scipy.optimize import bisect, brentq
from nipype.interfaces.workbench import MetricResample
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    TraitedSpec,
    traits,
    File,
    BaseInterface,
)
from tempfile import TemporaryDirectory
import subprocess
import nibabel as nib
import numpy as np
import shutil
from nilearn.image import math_img

# import dsargparse
import argparse
import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

sys.path.append("/data/p_02389/code/fmri-analysis/library/")
import layer_analysis as analysis

debug_print = True


def find_layer_contrast_roi(
    stat_file,
    anat_region,
    out_file,
    target_type,
    target,
    white_surf,
    pial_surf,
    hemi,
    depth_file,
    depth_ranges,
    depth_contrast,
    anat_require_region=None,
    cluster_roi=True,
    fwhm=3,
    cwd=None,
    keep_tmp=False,
):
    """
    Finds an roi that most closely meets criteria based on predefined regions (e.g. using atlas parcels and
    laminar contrasts of activation maps.)

    Args:
      stat_file: nifti volume file of activation statistic (only positive is used)
      anat_region: gifti surface mask of region (e.g. from atlas)
      out_file: name of generated roi nifti volume fil
      target_type: target measure, can be one of 'area','percentile','threshold','voxels' (not implemented)
      target: quantity of target measure for activation cluster/vertices/voxels
      white_surf: white surface as .gii in native volume space
      pial_surf: pial surface as .gii in native volume space
      hemi: hemisphere 'L' or 'R'
      depth_range: dictionary with depth labels (str) as keys and lists [min_dept, max_depth] as values
      depth_contrast: string specifying calculation to be used for math_metric

      cluster_roi: if True, use extract largest consecutive cluster of voxels
      anat_require_region: part of anat_region that cluster is strictly required to ovelap with (surface mask or None)
      fwhm: smoothing kernel width (only used together with cluster)
      cwd: current working directory (derived from output file if None)
    """

    # define function to calculate rois as a function of threshold
    def calc_roi(threshold):
        roi_surf_file = os.path.join(tmpdirname, f"roi.{hash(threshold)}.shape.gii")
        if not os.path.exists(roi_surf_file):
            if cluster_roi:
                # find all cluster (use anat_region as ROI)
                clusters_file = os.path.join(tmpdirname, "stat.clusters.func.gii")
                analysis.find_clusters_hcp(
                    surf_stat_file,
                    clusters_file,
                    mid_surf,
                    threshold,
                    min_area=0,
                    roi=anat_region,
                )

                clusters = nib.load(clusters_file).darrays[0].data
                if anat_require_region:
                    anat_require = (
                        nib.load(anat_require_region).darrays[0].data.astype(bool)
                    )
                    cluster_idcs = np.unique(clusters[anat_require])
                else:
                    cluster_idcs = np.unique(clusters)

                A_cluster = dict()
                for cluster_idx in cluster_idcs:
                    if cluster_idx > 0:
                        cluster = os.path.join(tmpdirname, "stat.cluster.func.gii")
                        analysis.math_metric(
                            f"clusters == {cluster_idx}",
                            cluster,
                            clusters=clusters_file,
                        )
                        A_cluster[cluster_idx] = analysis.calc_area_hcp(
                            cluster, mid_surf
                        )
                if A_cluster:
                    max_cluster_idx = max(A_cluster, key=A_cluster.get)
                    roi_surf_file = analysis.math_metric(
                        f"clusters == {max_cluster_idx}",
                        roi_surf_file,
                        clusters=clusters_file,
                    )
                else:
                    roi_surf_file = analysis.math_metric(
                        "0",
                        roi_surf_file,
                        clusters=clusters_file,
                    )
            else:
                analysis.math_metric(
                    f"(x > {threshold} ) && roi",
                    metric_out=roi_surf_file,
                    x=surf_stat_file,
                    roi=anat_region,
                )
        return roi_surf_file

    # set cwd to directory of output file, if not set otherwise
    if cwd == None:
        cwd = os.path.dirname(os.path.normpath(out_file))

    # process activation contrast (goal: laminar contrast on surface, maske by atlas requirement)
    with TemporaryDirectory(dir=cwd) as tmpdirname:
        # 1. sample stat_file to all surfaces, use wb_command (surface neeeded as input)
        mid_surf = os.path.join(tmpdirname, "mid.surf.gii")
        surf_stat_file = os.path.join(tmpdirname, "stat.func.gii")
        surf_stat_file_nosm = os.path.join(tmpdirname, "stat_nosm.func.gii")
        surf_stat_file_sm = os.path.join(tmpdirname, "stat_sm.func.gii")
        no_coverage_file = os.path.join(tmpdirname, "no_coverage.shape.gii")

        # for each depth sample to surface, keep surface activations, and no_coverage surface mask
        surf_stat_files_layers = dict()
        no_coverage_files_layers = dict()
        for layer in depth_ranges.keys():
            surf_stat_files_layers[layer] = os.path.join(tmpdirname, f"stat.{layer}.func.gii")
            no_coverage_files_layers[layer] = os.path.join(tmpdirname, f"no_coverage.{layer}.shape.gii")

            # compute layer masks
            depth_range = depth_ranges[layer]
            layer_roi = math_img(
                f"(img>={depth_range[0]}) & (img<={depth_range[1]})", img=depth_file
            )
            layer_roi.to_filename(os.path.join(tmpdirname, f"layer_{layer}.nii"))

            # sample surface from layer
            analysis.sample_surf_hcp(
                volume_file=stat_file,
                white_surf=white_surf,
                pial_surf=pial_surf,
                mid_surf=mid_surf,
                outfile=surf_stat_files_layers[layer],
                mask_file=os.path.join(tmpdirname, f"layer_{layer}.nii"),
                roi_out=no_coverage_files_layers[layer],
            )

        # combine no coverage files
        analysis.math_metric(' || '.join(no_coverage_files_layers.keys()), no_coverage_file, **no_coverage_files_layers)

        # calculate layer contrast and mask by no_coverage
        analysis.math_metric(f"({depth_contrast})*!no_coverage", surf_stat_file, no_coverage=no_coverage_file,**surf_stat_files_layers)

        # 2. if cluster: smooth surface
        if cluster_roi:
            shutil.copyfile(surf_stat_file, surf_stat_file_nosm)
            analysis.smooth_surfmetric_hcp(
                surf_stat_file, surf_stat_file, mid_surf, fwhm
            )
            shutil.copyfile(surf_stat_file, surf_stat_file_sm)

        # 4. Define optimization function
        if target_type == "threshold":
            tmp_roi = calc_roi(target)
        else:
            if target_type == "percentile":
                # calculate area of "covered" parcel
                anat_region_covered = analysis.math_metric(
                    "anat_region && !(not_covered)",
                    metric_out=os.path.join(
                        tmpdirname, "anat_region_covered.shape.gii"
                    ),
                    anat_region=anat_region,
                    not_covered=no_coverage_file,
                )
                A_total = analysis.calc_area_hcp(anat_region_covered, mid_surf)
                target = (target / 100) * A_total
                target_type = "area"
            if target_type == "area":

                def f(threshold):
                    roi_file = calc_roi(threshold)
                    area = analysis.calc_area_hcp(roi_file, mid_surf)
                    if debug_print:
                        print(f"A({threshold}) = {area}")
                    return target - area

            elif target_type == "voxels":
                raise NotImplementedError("voxel target not implemented yet!")

            max_stat = analysis.stats_metric_hcp(surf_stat_file, "MAX", anat_region)
            min_stat = analysis.stats_metric_hcp(surf_stat_file, "MIN", anat_region)
            if np.sign(f(min_stat)) == np.sign(f(max_stat)):
                tmp_roi = os.path.join(tmpdirname, f"roi.empty.shape.gii")
                analysis.math_metric(
                    "0",
                    tmp_roi,
                    img=surf_stat_file,
                )
            else:
                threshold = brentq(f, min_stat, max_stat)

                # 4. get final roi (on surface)
                tmp_roi = calc_roi(threshold)

        shutil.copyfile(tmp_roi, os.path.join(tmpdirname, "roi.final.shape.gii"))

        # sample surface roi to volume:
        analysis.surf_to_vol_hcp(
            tmp_roi,
            out_file,
            stat_file,
            white_surf,
            pial_surf,
            mid_surf,
            greedy=False,
            is_roi=True,
        )
        if keep_tmp:
            subprocess.run(
                "cp -r " + os.path.join(tmpdirname, "*") + " " + cwd,
                shell=True,
                check=True,
            )
    return out_file


def test_find_layer_contrast_roi(testdata_dir):
    # extract regions from atlas
    analysis.generate_atlas_region_hcp(
        os.path.join(testdata_dir, "sub-08.GlasserAtlas.native.label.gii"),
        os.path.join(testdata_dir, "dlpfc.shape.gii"),
        [83, 84, 73, 81, 82],
    )

    analysis.generate_atlas_region_hcp(
        os.path.join(testdata_dir, "sub-08.GlasserAtlas.native.label.gii"),
        os.path.join(testdata_dir, "p9_46v.shape.gii"),
        [83],
    )

    find_layer_contrast_roi(
        stat_file=os.path.join(testdata_dir, "fstat.nii"),
        anat_region=os.path.join(testdata_dir, "dlpfc.shape.gii"),
        out_file=os.path.join(testdata_dir, "roi.nii"),
        target_type="area",
        target=500,
        white_surf=os.path.join(testdata_dir, "lh.white_converted.transformed.gii"),
        pial_surf=os.path.join(testdata_dir, "lh.pial_converted.transformed.gii"),
        hemi="L",
        depth_file=os.path.join(testdata_dir, "vdfs_depths_equivol.nii"),
        depth_ranges={'deep':[0,0.5],'superficial':[0.5,1]},
        depth_contrast="deep - superficial",
        cluster_roi=True,
        anat_require_region=os.path.join(testdata_dir, "p9_46v.shape.gii"),
        fwhm=3,
        keep_tmp=True,
    )
    

class FindLayerContrastROIInputSpec(BaseInterfaceInputSpec):
    stat_file = File(
        exists=True,
        desc="nifti volume file of activation statistic (only positive is used)",
        mandatory=True,
    )
    anat_region = File(
        exists=True,
        desc="gifti surface mask of region (e.g. from atlas)",
        mandatory=True,
    )
    out_file = File(desc="name of generated roi nifti volume file", mandatory=True)
    target_type = traits.Enum(
        "area",
        "percentile",
        "thredhold",
        "voxels",
        desc="target measure",
        mandatory=True,
    )
    target = traits.Float(
        desc="quantity of target measure for activation cluster/vertices/voxels",
        mandatory=True,
    )
    white_surf = File(
        exists=True, desc="white surface as .gii in native volume space", mandatory=True
    )
    pial_surf = File(
        exists=True, desc="pial surface as .gii in native volume space", mandatory=True
    )
    hemi = traits.Enum("L", "R", desc="hemisphere", mandatory=True)
    depth_file = File(
        exists=True,desc="nifti volume file with depth information", mandatory=True
    )
    depth_ranges = traits.Dict(
        desc="dictionary with depth labels (str) as keys and lists [min_dept, max_depth] as values",
        mandatory=True,
    )
    depth_contrast = traits.Str(
        desc="string specifying calculation to be used for math_metric", mandatory=True
    )

    # optional
    cluster_roi = traits.Bool(
        mandatory=False,
        desc="if True, use extract largest consecutive cluster of voxels",
    )
    anat_require_region = File(
        desc="part of anat_region that cluster is strictly required to ovelap with (surface mask or None)"
    )
    fwhm = traits.Float(desc="smoothing kernel width (only used together with cluster)")


class FindLayerContrastROIOutputSpec(TraitedSpec):
    out_file = File(desc="the generated roi")


class FindLayerContrastROI(BaseInterface):
    input_spec = FindLayerContrastROIInputSpec
    output_spec = FindLayerContrastROIOutputSpec

    def _run_interface(self, runtime):
        find_layer_contrast_roi(
            stat_file=self.inputs.stat_file,
            anat_region=self.inputs.anat_region,
            out_file=self.inputs.out_file,
            target_type=self.inputs.target_type,
            target=self.inputs.target,
            white_surf=self.inputs.white_surf,
            pial_surf=self.inputs.pial_surf,
            hemi=self.inputs.hemi,
            depth_file=self.inputs.depth_file,
            depth_ranges=self.inputs.depth_ranges,
            depth_contrast=self.inputs.depth_contrast,
            cluster_roi=self.inputs.cluster_roi,
            anat_require_region=self.inputs.anat_require_region,
            fwhm=self.inputs.fwhm,
        )
        return runtime

    def _list_outputs(self):
        return {"out_file": self.inputs.out_file}


if __name__ == "__main__":
    # call test function with first command line argument as testdata_dir
    test_find_layer_contrast_roi('testdata')
    # process command line arguments and either call interface or implemented function
    # parser_description = "find roi"
    # parser = dsargparse.ArgumentParser(main=find_roi)
    # parser.add_argument("stat_file")
    # parser.add_argument("target")

    # parser.parse_and_run()

    # FindLayerContrastRoi(stat_file=args.stat_file).run()
