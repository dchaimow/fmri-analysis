"""
Motion correction workflows

to consider
-----------
- single contrast and double contrast data (vaso), interleaved?
- afni, fsl and spm, (ants mc??)
- check fmriprep for ref volume generation
- for double contrast first mc each contrast tham register both contrasts and combine transformations?

"""

from nipype.pipeline import engine as pe
from nipype.interfaces.utility import IdentityInterface


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

    inputnode = pe.Node(IdentityInterface(fields=["in_files"], name="inputnode"))
    outputnode = pe.Node(
        IdentityInterface(fields=["run_idx, volume_idx"], name="outputnode")
    )
    outlier_count = pe.MapNode(
        OutlierCount(automask=True, fraction=True, polort=5, legendre=True),
        name="outlier_count",
        iterfield=['in_file'])
    )
    add_counts = pe.Node(Eval(expr="a+b"), name="add_counts")
    

    
    
    pass
