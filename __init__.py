'''BatchBrain
a Python batch processing interface for neuroimaging analysis. Basically a
fancy wrapper around functionality from FSL / mricron.
2012 J Carlin'''

__version__ = 'alpha'
# TODO:
# move code out of __init__.py 
# implement parallel support
import sys
import os
import platform

# Could also write outputs to a file...
output = sys.stdout.write
output('Welcome to BatchBrain.\nInitialising...\n')

# Initialise FSL
arch = platform.architecture()[0]
cp = os.environ['PATH']
# Check if fsl is already added
fslver = 'fsl-4.1.8'


# Note that the 64bit switching here will only work if your Python
# install is also 64bit
if fslver in cp:
	sys.stdout.write('FSL already on path\n')
else:
	sys.stdout.write('Adding FSL to system path...\n')
	if arch == '32bit':
		sys.stdout.write('Adding 32bit FSL...\n')
		fsldir = '/imaging/local/software/fsl/fsl32/fsl-4.1.8/fsl'
	elif arch == '64bit':
		sys.stdout.write('Adding 64bit FSL...\n')
		fsldir = '/imaging/local/software/fsl/fsl64/fsl-4.1.8/fsl'
	else:
		raise Exception('Unknown system architecture: %s' % arch)
	os.environ['PATH'] = os.path.join(fsldir,'bin:') + cp
	os.environ['FSLDIR'] = fsldir
	os.environ['FSLOUTPUTTYPE'] = 'NIFTI'

# Add dti_motion script to path
bbdir = '/home/jc01/code/UtilityScripts/Python/batchbrain:'
os.environ['PATH'] = bbdir + os.environ['PATH']

import batchbrain.base
import batchbrain.data
import batchbrain.processes
import batchbrain.runners

sys.stdout.write('Initialisation done.\n')

# TODO - parallelisation initialisation and settings here
