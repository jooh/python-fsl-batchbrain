'''Objects for processing - the intermediate step where the files in a 
StudyData-derived object are joined with relevant analyses and added to the
Tasks field in StudyData.'''

import os
import glob
import batchbrain.base
import pdb
import numpy
import csv

class DataProcess(batchbrain.base.ProcessingObject):
	"""Master object for things to do with processing imported data."""
	postrunfun = None

	def addPrefix(self,vol,pref):
		filepath,file = os.path.split(vol)
		if pref is None:
			return vol
		return os.path.join(filepath,'%s%s' % (pref,file))

	def processCommand(self,cmd):
		"""send a command line string to the shell. Easily parallelised."""
		os.system(cmd)
		return

	def ensureListofLen(self,inp,wantedlen):
		objinp = getattr(self,inp)
		if type(objinp) is list:
			if len(objinp) == wantedlen:
				return 
			if len(objinp) != 1:
				raise Exception('Lists must be len %d or 1' % wantedlen)
		else:
			objinp = [objinp]
		setattr(self,inp,objinp * wantedlen)
		return 

	def multiProcess(self,dt):
		'''Method for generating multi-in 1-out command line processes.
		Note that this _should_ work also with 1-in 1-out data.'''
		newjobs = 0
		# First, figure out how many 'multi' means
		# Get the len of any lists
		listcands = ['inprefix','targetsequence',
				'inext','acqrule']
		lens = [len(getattr(self,x)) for x in listcands \
				if type(getattr(self,x)) is list]
		if lens:
			maxlen = max(lens)
		else:
			maxlen = 1
		# Now make each listcands item a list if it isn't already
		for x in listcands:
			self.ensureListofLen(x,maxlen)

		for sub in dt.subjects.keys():
			# list of lists
			inputlist = []
			for ind in xrange(maxlen):
				invols = dt.getPaths(sub,self.targetsequence[ind],
						self.acqrule[ind],self.inprefix[ind],
						self.inext[ind])
				inputlist.append(invols)

			maxlistlen = len(max(inputlist,key=len))
			# We may have to check that we always get the same number of
			# volumes here (e.g., if different sequences), but for now,
			# let's assume all is well...
			for ind in xrange(maxlistlen):
				# Pull out invol list
				inputs = [x[ind] for x in inputlist]
				# Output checking is necessary to avoid re-running on every
				# iteration. Checks should be either based on checking dir
				# contents (current solution) OR checking attributes in
				# StudyData object.
				# Figure out what to call outvol
				output = self.getOutput(inputs,dt)
				alreadydone = self.checkOutput(output)
				# TODO actually make the outvol maybe. Also need to store
				# some kind of reference to outfield to set pars.
				if alreadydone and not self.overwrite:
					self.output('%s already exists, skipping.' % output)
					continue
				cmd = self.process(inputs,output)
				dt.addJobs(self.__class__.__name__,self.runner,cmd,
						postrunfun=self.postrunfun)
				newjobs+=1
		self.output('Done. Added %d new jobs to tasks.' % newjobs)
		return

class eddy_correct(DataProcess):

	def __init__(self,targetsequence,inprefix=None,outprefix='ec_', \
			refno=0,acqrule='all',**kwargs):
		self.processname = 'FSL Eddy Correct'
		self.targetsequence = targetsequence
		self.inprefix = inprefix
		self.outprefix = outprefix
		self.refno = refno
		self.runner = self.processCommand
		self.acqrule = acqrule
		self.logfn = self.outprefix + self.__class__.__name__ + '.log'
		self.checkInputs(kwargs)
		return

	def __call__(self,dt):
		self.multiProcess(dt)
		return

	def process(self,invol,outvol):
		path,filename = os.path.split(invol[0])
		logfile = os.path.join(path,self.logfn)
		self.initialiseLog(logfile)

		return 'eddy_correct %s %s %d >> %s' % (invol[0],outvol,self.refno,
				logfile)

class bet(DataProcess):
	'''Brain extraction with BET from FSL. This variant generates a brain
	mask only.'''

	def __init__(self,targetsequence,inprefix=None,outprefix='bet_', \
			f=0.5,acqrule='all',**kwargs):
		self.processname = 'FSL BET Brain Extraction'
		self.targetsequence = targetsequence
		self.inprefix = inprefix
		self.outprefix = outprefix
		self.f = f
		self.runner = self.processCommand
		self.acqrule = acqrule
		# the log file name ought to be the class name for transparency
		self.logfn = self.outprefix + self.__class__.__name__ + '.log'
		self.checkInputs(kwargs)
		return

	def __call__(self,dt):
		self.multiProcess(dt)
		return

	def process(self,invol,output):
		# Note that many other arguments are possible but not yet
		# implemented
		# Get base name
		path,filename = os.path.split(invol[0])
		logfile = os.path.join(path,self.logfn)
		self.initialiseLog(logfile)
		maincmd = 'bet %s %s -m -n -f %f >> %s' % (invol[0],output,
				self.f,logfile)
		# BET has an unfortunate tendency to enforce a certain filename
		# (output_mask.nii) - so rename output
		ind = output.index(self.outext)
		mvcmd = 'mv %s_mask%s %s' % (output[:ind],self.outext,output)
		return maincmd + ';' + mvcmd

class slicer(DataProcess):
	'''Interface to FSL slicer utility. Write out diagnostic axials.'''

	def __init__(self,targetsequence,inprefix=None,outprefix='slicer_',\
			nskipaxial=3,imwidth=400,acqrule='all',**kwargs):
		'''Identify 2 images (first is base, second overlay) by differences
		in targetsequence AND/OR inprefix, inex, acqrule.'''
		self.processname = 'FSL Slicer Image Output'
		self.targetsequence = targetsequence
		self.inprefix = inprefix
		self.outprefix = outprefix
		self.nskipaxial = nskipaxial
		self.imwidth = imwidth
		# NB slicer does NOT support other formats - don't change this!
		self.outext = '.ppm'
		self.acqrule = acqrule
		self.runner = self.processCommand
		# the log file name ought to be the class name for transparency
		self.logfn = self.outprefix + self.__class__.__name__ + '.log'
		self.checkInputs(kwargs)
		return
	
	def __call__(self,dt):
		self.multiProcess(dt)
		return

	def process(self,invol,output):
		'''FSL slicer method. input should be a list (first base, second
		overlay), output string.'''
		# Get base name
		path,filename = os.path.split(invol[0])
		logfile = os.path.join(path,self.logfn)
		self.initialiseLog(logfile)
		return 'slicer %s %s -S %d %d %s >> %s' % (invol[0],invol[1],
				self.nskipaxial,self.imwidth,output,logfile)

class dtifit(DataProcess):
	'''Interface to FSL dtifit.'''

	def __init__(self,targetsequence,inprefix=None,outprefix='dtifit_',\
			outext='*',acqrule='all',**kwargs):
		'''Identify 2 images (first is base, second overlay) by differences
		in targetsequence AND/OR inprefix, inex, acqrule.'''
		self.processname = 'FSL DTIFit'
		self.targetsequence = targetsequence
		self.inprefix = inprefix
		self.outprefix = outprefix
		self.acqrule = acqrule
		self.runner = self.processCommand
		self.outext = outext
		# the log file name ought to be the class name for transparency
		self.logfn = self.outprefix + self.__class__.__name__ + '.log'
		self.checkInputs(kwargs)

		return
	
	def __call__(self,dt):
		self.multiProcess(dt)
		return

	def process(self,invol,output):
		'''FSL DTIFit method. Input should be a list (first DTI, second
		mask), output string basename.'''
		# Get base name
		path,filename = os.path.split(invol[0])
		logfile = os.path.join(path,self.logfn)
		self.initialiseLog(logfile)
		# Find bvecs/bvals
		bvec = glob.glob(os.path.join(path,'*.bvec'))
		if len(bvec) != 1:
			raise Exception('Expected 1 bvec only, got: \n%s' % bvec)
		bvec = bvec[0]
		bval = glob.glob(os.path.join(path,'*.bval'))
		if len(bval) != 1:
			raise Exception('Expected 1 bval only, got: \n%s' % bval)
		bval = bval[0]
		return 'dtifit -k %s -m %s -r %s -b %s -o %s >> %s' % (invol[0],
				invol[1],bvec,bval,output.replace('*',''),logfile)

class dti_motion(DataProcess):
	'''Interface to FSL motion diagnostic bash script by Mark Jenkinson. NB
	dti_motion bash script MUST be on your path!'''

	def __init__(self,targetsequence,inprefix='ec_',outprefix='dtimotion_',\
			outext='*',acqrule='all',**kwargs):
		self.processname = 'FSL dti_motion'
		self.targetsequence = targetsequence
		self.inprefix = inprefix
		self.outprefix = outprefix
		self.acqrule = acqrule
		self.runner = self.processCommand
		self.outext = outext
		# Generated files that we rename
		self.outfilelist = ['ec_disp.png','ec_disp.txt','ec_rot.png',
				'ec_rot.txt','ec_trans.png','ec_trans.txt','grot_ts.txt']
		# the log file name ought to be the class name for transparency
		self.logfn = self.outprefix + self.__class__.__name__ + '.log'
		self.checkInputs(kwargs)
		return
	
	def __call__(self,dt):
		self.multiProcess(dt)
		return

	def process(self,invol,output):
		'''FSL motion estimate (dti_motion).'''
		# Get base name
		path,filename = os.path.split(invol[0])
		logfile = os.path.join(path,self.logfn)
		self.initialiseLog(logfile)
		# Find ecclog
		ecclog = glob.glob(os.path.join(path,'*.ecclog'))
		if len(ecclog) != 1:
			raise Exception('Expected 1 ecclog only, got: \n%s' % ecclog)
		ecclog = ecclog[0]
		# This is the command - but NB must be executed from subject dir to
		# avoid disastrous consequences
		cdcmd = 'cd %s;' % path
		inpath,infile = os.path.split(ecclog)
		infn,inext = os.path.splitext(infile)
		# Now with piping of log output
		basecmd = 'dti_motion %s >> %s;' % (infile,logfile)
		# Now one irritating thing about this script is that we can't
		# control output file naming, so need to rename after calling
		mvcmd = ''
		for file in self.outfilelist:
			fn,ext = os.path.splitext(file)
			mvcmd += 'mv %s %s;' % (file,self.outprefix+infn+'_'+fn+ext)
		return cdcmd+basecmd+mvcmd

class dti_motion_toparams(DataProcess):
	'''Place a stored list of motion parameters in the parameters field for
	each acquisition.'''

	def __init__(self,targetsequence,inprefix='dti_motion_', \
			inext=['*ec_trans.txt','*ec_rot.txt'],acqrule='all',**kwargs):
		self.processname = 'BB dti_motion to parameters'
		self.targetsequence = targetsequence
		self.inprefix = inprefix
		self.inext = inext
		# No output generated.
		self.outprefix = ''
		self.outext = ''
		self.acqrule = acqrule
		# NB things can go badly awry here if you mix up your inputs
		self.parstoset = {'translation': self.inext[0], 'rotation':
				self.inext[1]}
		# No need to define self.runner or self.postrunfun since we have
		# function handles that (presumably) overwrite

		# Generated files that we rename
		# the log file name ought to be the class name for transparency
		self.logfn = self.outprefix + self.__class__.__name__ + '.log'
		self.checkInputs(kwargs)
		return

	def __call__(self,dt):
		self.multiProcess(dt)
		return

	def getOutput(self,invol,sd):
		'''Use input volume to deduce acq details. Overwrites default
		DataProcess method that operates on files in the acq dir.'''
		# Figure out which acq object is relevant
		# Will be VERY interesting to see if changes to this object
		# actually carry over to main sd object
		return sd.pathToAcquisition(invol[0])

	def checkOutput(self,output):
		'''Check if the output (an acq object) contains the desired
		parameters. Overwrite default DataProcess method that operates on
		files in the acq dir.'''
		allouts = [output.parameters.has_key(x) for x in self.parstoset]
		# Only score as done if all pars were set
		return all(allouts)

	def process(self,infiles,output):
		'''Make a dict of parameter file names and their paths. Prepare a
		single variable to be called together with runner (below)'''
		pathdict = {}
		# Find each parameter's associated file
		for f,v in self.parstoset.iteritems():
			# Maybe a bit ambitious
			hits = self.findStrInList(v.replace('*',''),
					infiles)
			if len(hits) != 1:
				raise Exception('parstoset/inext mismatch for %s, got:\n' \
						% (f,hits))
			pathdict[f] = hits[0]
		return [pathdict,output]

	def runner(self,inputs):
		'''Get data out of inputs,store in parameters. Also summary stat
		(SSQ).'''
		# Unpack (deal) input list
		pathdict,acq = inputs
		resdict = {}
		for measure,filepath in pathdict.iteritems():
			# Load the data
			F = open(filepath,'r')
			R = csv.reader(F,delimiter=' ')
			# An annoying aspect of this is that a 4th garbage empty string
			# column is included for unknown reasons.
			resdict[measure] = numpy.array([i[0:3] for i in R],dtype='float')
		return [resdict,acq]

	def postrunfun(self,outputs):
		'''Add the list of resdicts to the SD object.'''
		# Unpack
		resdict,acq = outputs
		for k,v in resdict.iteritems():
			# Update acq object
			acq.parameters[k] = v
			# But is this in place?
			# Otherwise we can probably recover with pathToAcq
		return

class tbss_1_preproc(DataProcess):
	'''First step of TBSS preprocessing.'''

	def __init__(self,group,**kwargs):
		'''Initialise TBSS preprocessing. Unusual in that you only provide
		group object.'''
		self.processname = 'FSL TBSS Preprocessing 1'
		self.targetsequence = targetsequence
		self.inprefix = inprefix
		self.outprefix = outprefix
		self.f = f
		self.runner = self.processCommand
		self.acqrule = acqrule
		# the log file name ought to be the class name for transparency
		self.logfn = self.outprefix + self.__class__.__name__ + '.log'
		self.checkInputs(kwargs)
		return

	def __call__(self,dt):
		self.multiProcess(dt)
		return

	def process(self,invol,output):
		# Note that many other arguments are possible but not yet
		# implemented
		# Get base name
		path,filename = os.path.split(invol[0])
		logfile = os.path.join(path,self.logfn)
		self.initialiseLog(logfile)
		maincmd = 'bet %s %s -m -n -f %f >> %s' % (invol[0],output,
				self.f,logfile)
		# BET has an unfortunate tendency to enforce a certain filename
		# (output_mask.nii) - so rename output
		ind = output.index(self.outext)
		mvcmd = 'mv %s_mask%s %s' % (output[:ind],self.outext,output)
		return maincmd + ';' + mvcmd

