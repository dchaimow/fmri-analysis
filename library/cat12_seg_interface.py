from nipype.interfaces.base import traits
from nipype.interfaces.base import TraitedSpec, InputMultiPath
from nipype.interfaces.matlab import MatlabCommand, MatlabInputSpec
from nipype.interfaces.spm.base import ImageFileSPM

class cat12_seg_InputSpec(MatlabInputSpec):
    in_file = InputMultiPath(
        ImageFileSPM(exists=True),
        field="data",
        desc="file to segment",
        mandatory=True,
        copyfile=False,
    )

class cat12_seg_OutputSpec(TraitedSpec):
    matlab_output = traits.Str()


class cat12_seg(MatlabCommand):
    """
    CAT 12 Interface for GM and WM segmentation of 7T MPRAGE data using fixed parameters
    """
    
    input_spec = cat12_seg_InputSpec
    output_spec = cat12_seg_OutputSpec

    def _my_script(self):
        script = """
        matlabbatch{1}.spm.tools.cat.estwrite.data = {'%s,1'};
        matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
        matlabbatch{1}.spm.tools.cat.estwrite.nproc = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
        matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {'%s/tpm/TPM.nii'};
        matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
        matlabbatch{1}.spm.tools.cat.estwrite.opts.biasacc = 0.5;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.native = [];
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.setCOM = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.affmod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.spm_kamap = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASmyostr = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.WMHC = 2;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'%s/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'};
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.bb = 12;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.SRP = 22;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};
        matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0];
        matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;
        """ % (self.inputs.name, self.inputs.spm_path, self.inputs.spm_path)
        return script

    def run(self, **inputs):
        # Inject your script
        self.inputs.script = self._my_script()
        results = super(MatlabCommand, self).run(**inputs)
        stdout = results.runtime.stdout
        # Attach stdout to outputs to access matlab results
        results.outputs.matlab_output = stdout
        return results

    def _list_outputs(self):
        outputs = self._outputs().get()
        return outputs
