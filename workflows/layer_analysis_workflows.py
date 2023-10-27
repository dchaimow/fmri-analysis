from nipype.pipeline import engine as pe

# Registration and transformation

# Registration workflows
def init_surftransform_wf():
    """
    Worflow that transforms gifti or freesurfer surfaces using ants transform

    Parameters
    ----------

    Inputs
    ------
    - gifti or freesurfer surface
    - ants transform


    Outputs
    -------
    - transformed surface


    """
    # Nodes:

    # (MRIsConvert)

    # GiftiToCSV

    # ApplyTransformsToPoints

    # CSVToGifti

    # (MRIsConvert)

    pass


def init_process_vaso_wf():
    """
    Worflow to process functional vaso files

    Parameters
    ----------
    - first volume not-nulled
    - time between readouts
    - number of volumes to discard/overwrite
    - don't interpolate not-nulled


    Inputs
    ------
    - list (or single instance) of nulled/non-nulled pairs of filed
      - alternatively: list of files with interleaved nulled/not-nulled images each


    Outputs
    -------
    - processed vaso and bold files
     - optionally aligned by applying temporal shift to bold part
    - vaso t1
    - (vaso t1 skull stripped)
    - mean images
    - snr images
    - (skew, ...?)


    """
    # Nodes

    # something to overwrite/remove first volumes

    # something to possibly split volumes

    # motion correction wf

    # bold correction wf

    # vaso t1 workflow
    pass


def init_boldcorrect_wf():
    pass


def init_register_fs_to_vaso_wf():
    pass


def init_motion_correction_wf():
    # find_best_volume_wf (or choose alternative reference volume?)

    # 3dvolreg on all files (or choose alternative methods?)

    # generate plots (interface?)
    pass


def init_find_best_volume_wf():
    """
    Workflow to find least outlying volume from a list of (pairs of) runs.

    Parameters
    ----------

    Inputs
    ------
    - list of runs

    Outputs
    -------
    - number of run
    - number of volume

    """
    from nipype.interfaces.afni import OutlierCount, Eval

    outlier_count = pe.Node(
        OutlierCount(automask=True, fraction=True, polort=5, legendre=True),
        name="outlier_count",
    )

    add_counts = pe.Node(Eval(expr="a+b"), name="add_counts")

    pass


# ROI workflows
def init_fs_LR_volume_atlas_wf():
    pass


def init_stat_cluster_volume_atlas_wf():
    pass


# Surface sampling and transform workflows


# Analysis workflows (trial averaging and glm analysis)


# Layering wokflows
