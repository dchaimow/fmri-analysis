#! /usr/bin/env python3
"""
Define interface and run in if in command line mode. Use a workflow?
"""
from nipype import Node, Workflow
from nipype.interfaces.freesurfer import MRIConvert
from nipype.interfaces.ants import Registration
from nipype.interfaces.io import FreeSurferSource

# TODO:
# - implement use of mask
# - implement use of initial_matrix
# - correct call using command line arguments
# = which files to output/save where (e.g. delete fs_T1.nii, fs_to_func_warped.nii,
# or specify files to output and where, working/output dir?)

def init_fs_to_epiT1_reg_wf(
    epi_t1_file, fs_subject, fs_subjects_dir, init_transform=None, mask=None
):
    # 0. get fs T1
    fs = FreeSurferSource()
    fs.inputs.subject_id = fs_subject
    fs.inputs.subjects_dir = fs_subjects_dir
    fs_node = Node(fs, name="fs_source")

    # 1. mri_convert
    conv = MRIConvert()
    conv.inputs.out_file = "fs_t1.nii"
    conv.inputs.out_type = "nii"
    conv_node = Node(conv, name="convert_to_nifti")

    # 2. antsRegistration
    reg = Registration()
    reg.inputs.fixed_image = epi_t1_file
    reg.inputs.float = False
    reg.inputs.interpolation = "BSpline"
    reg.inputs.interpolation_parameters = (5,)
    reg.inputs.use_histogram_matching = False
    reg.inputs.winsorize_lower_quantile = 0.005
    reg.inputs.winsorize_upper_quantile = 0.995

    reg.inputs.metric = ["MI", "MI", "CC"]
    reg.inputs.metric_weight = [1, 1, 1]
    reg.inputs.radius_or_number_of_bins = [32, 32, 4]
    reg.inputs.sampling_strategy = ["Regular", "Regular", None]
    reg.inputs.sampling_percentage = [0.25, 0.25, None]
    reg.inputs.transforms = ["Rigid", "Affine", "SyN"]
    reg.inputs.transform_parameters = [[0.1], [0.1], [0.1, 3, 0]]
    reg.inputs.smoothing_sigmas = [[4, 3, 2, 1], [4, 3, 2, 1], [5, 3, 2, 1, 0]]
    reg.inputs.shrink_factors = [[12, 8, 4, 2], [12, 8, 4, 2], [10, 6, 4, 2, 1]]
    reg.inputs.number_of_iterations = [
        [1000, 500, 250, 100],
        [1000, 500, 250, 100],
        [50, 50, 70, 50, 20],
    ]

    reg.inputs.output_transform_prefix = "fs_to_func_warped_"
    reg.inputs.output_warped_image = "fs_to_func_Warped.nii"
    reg.inputs.output_inverse_warped_image = "fs_to_func_InverseWarped.nii"
    reg_node = Node(reg, name="register_fs_to_vasot1")

    wf = Workflow(name="register_fs_to_vasot1", base_dir=None)
    wf.connect(
        [
            (fs_node, conv_node, [("T1", "in_file")]),
            (conv_node, reg_node, [("out_file", "moving_image")]),
        ]
    )

    return wf


if __name__ == "__main__":
    fs_to_epiT1_reg =  init_fs_to_epiT1_reg_wf(
        epi_t1_file, fs_subject, fs_subjects_dir, init_transform=None, mask=None)
):
