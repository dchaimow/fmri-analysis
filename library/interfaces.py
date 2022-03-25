# command line
from nipype.interfaces.base import CommandLine
from nipype.interfaces.base import CommandLineInputSpec, File
from nipype.interfaces.base import TraitedSpec, traits
from nipype.interfaces.utility import Function
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec
import re
import nibabel as nb
from scipy.io import savemat

from nipype.interfaces.base import CommandLineInputSpec, InputMultiObject

# interface for motion correction:
# - one interface for BOLD and VASO?
# - how to specifiy image pairs?
# - outputs:
#  - motion corrected data
#  - motion parameters
#  - plots?
# - check FSL, AFNI, ... interfaces
class MotionCorrectionInputSpec(CommandLineInputSpec):
    in_files = InputMultiObject(
        File(exists=True), mandatory=True, desc="list of input nifti files"))

class MotionCorrectionOutputSpec():
    pass

    
class MotionCorrection(CommandLine):
    pass





# 1. use CommandLine Interface class
nipype_ls = CommandLine('ls', args='-lh', terminal_output='allatonce')

# 2. create new class inheriting from CommandLine
class TransformInfoInputSpec(CommandLineInputSpec):
    in_file = File(exists=True, mandatory=True, argstr='%s',
                   position=0, desc='the input transform file')

class TransformInfoOutputSpec(TraitedSpec):
    translation = traits.List(traits.Float,
                              desc='the translation component of the input transform')

class TransformInfo(CommandLine):
    _cmd = 'antsTransformInfo'
    input_spec = TransformInfoInputSpec
    output_spec = TransformInfoOutputSpec

    # the next part is needed when providing additional outputs based on std output of run
    def _run_interface(self, runtime):
        import re

        # Run the command line as a natural CommandLine interface
        runtime = super(TransformInfo, self)._run_interface(runtime)

        # Search transform in the standard output
        expr_tra = re.compile('Translation:\s+\[(?P<translation>[0-9\.-]+,\s[0-9\.-]+,\s[0-9\.-]+)\]')
        trans = [float(v) for v in expr_tra.search(runtime.stdout).group('translation').split(', ')]

        # Save it for later use in _list_outputs
        setattr(self, '_result', trans)

        # Good to go
        return runtime

    def _list_outputs(self):
        # Get the attribute saved during _run_interface
        return {'translation': getattr(self, '_result')}

# 3. simpler case without extending run method
class CustomBETInputSpec(CommandLineInputSpec):
    in_file = File(exists=True, mandatory=True, argstr='%s', position=0, desc='the input image')
    mask = traits.Bool(mandatory=False, argstr='-m', position=2, desc='create binary mask image')

    # Do not set exists=True for output files!
    out_file = File(mandatory=True, argstr='%s', position=1, desc='the output image')

class CustomBETOutputSpec(TraitedSpec):
    out_file = File(desc='the output image')
    mask_file = File(desc="path/name of binary brain mask (if generated)")

class CustomBET(CommandLine):
    _cmd = 'bet'
    input_spec = CustomBETInputSpec
    output_spec = CustomBETOutputSpec

    def _list_outputs(self):
        # Get the attribute saved during _run_interface
        return {'out_file': self.inputs.out_file,
                'mask_file': self.inputs.out_file.replace('brain', 'brain_mask')}

# 4. wrap python code - simple approach
my_python_interface = Function(
    input_names=['img', 'translation', 'out_file'],
    output_names=['out_file'],
    function=translate_image
)

# 5. wrap -python code - complete approach
class TranslateImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='the input image')
    out_file = File(mandatory=True, desc='the output image') # Do not set exists=True !!
    translation = traits.List([50.0, 0.0, 0.0], traits.Float, usedefault=True,
                              desc='the translation component of the input transform')

class TranslateImageOutputSpec(TraitedSpec):
    out_file = File(desc='the output image')

class TranslateImage(BaseInterface):
    input_spec = TranslateImageInputSpec
    output_spec = TranslateImageOutputSpec

    def _run_interface(self, runtime):

        # Call our python code here:
        translate_image(
            self.inputs.in_file,
            self.inputs.translation,
            self.inputs.out_file
        )
        # And we are done
        return runtime

    def _list_outputs(self):
        return {'out_file': self.inputs.out_file}

# 6. matlab interface
class BrainVolumeMATLABInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True)
    script_file = File(exists=True, mandatory=True)

class BrainVolumeMATLABOutputSpec(TraitedSpec):
    volume = traits.Int(desc='brain volume')

class BrainVolumeMATLAB(BaseInterface):
    input_spec = BrainVolumeMATLABInputSpec
    output_spec = BrainVolumeMATLABOutputSpec

    def _run_interface(self, runtime):
        # Save the image in matlab format as tmp_image.mat
        tmp_image = 'tmp_image'
        data = nb.load(self.inputs.in_file).get_data()
        savemat(tmp_image, {b'data': data}, do_compression=False)

        # Load script
        with open(self.inputs.script_file) as script_file:
            script_content = script_file.read()

        # Replace the input_image.mat file for the actual input of this interface
        with open('newscript.m', 'w') as script_file:
            script_file.write(script_content.replace('input_image.mat', 'tmp_image.mat'))

        # Run a matlab command
        mlab = CommandLine('octave', args='newscript.m', terminal_output='stream')
        result = mlab.run()

        expr_tra = re.compile('total\ =\s+(?P<total>[0-9]+)')
        volume = int(expr_tra.search(result.runtime.stdout).groupdict()['total'])
        setattr(self, '_result', volume)
        return result.runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['volume'] = getattr(self, '_result')
        return outputs
