'''Objects for storing data.'''
import collections
import os
import glob
import numpy
import batchbrain.base 
import pdb
import shutil
import scipy.stats
import csv

class Task(batchbrain.base.ProcessingObject):
	"""Stores properties for a given task (ie, a function and a list of
	arguments to call it with).
	In the future, this will likely be useful for flow control
	(waitfornext etc)"""

	def __init__(self,taskname,procfunc,postrunfun=None,**kwargs):
		"""Initialise a new task object."""
		self.name = taskname
		self.func = procfunc
		self.postrunfun = postrunfun
		# Just return 
		if self.postrunfun is None:
			def postrunfun(out):
				return None
			self.postrunfun = postrunfun
		self.checkInputs(kwargs)
		self.jobs = []
		return

	def updateJobs(self,args):
		self.jobs.append(args)
		return

	def postRunUpdate(self,outputs):
		"""Deal with outputs from some processing."""
		for o in outputs:
			self.postrunfun(o)
		return

class AcquisitionData(batchbrain.base.ProcessingObject):
	"""Storage object for a single acquisition (smallest object in the
	hierarchy unless we go to VolumeData..."""

	def __init__(self,acqdir,parameters=None,**kwargs):
		# Check for DICOM import status
		self.checkInputs(kwargs)
		madedir = self.mkAnaDir(acqdir)
		self.imported = madedir==False
		self.path = acqdir
		self.name = os.path.split(acqdir)[1]
		# Store things like e.g. motion data
		self.parameters = parameters
		if self.parameters is None:
			self.parameters = {}
		return

	def findFiles(self,prefix,ext):
		"""Find some files in the acquisition dir with glob."""
		searchstr = os.path.join(self.path,
				'%s*%s' % (prefix,ext))
		return glob.glob(searchstr)

	def summariseTranslation(self):
		'''Calculate some summary statistic for how much a subject moved.
		Not sure how to do this ultimately. For now, use the sum of the
		absolute distance.'''
		# Get the scan-by-scan difference in xyz
		transdiff = numpy.diff(self.parameters['translation'],axis=0)
		# Compute euclidean distance of diffs
		transdist = numpy.sqrt(transdiff[:,0]**2 + transdiff[:,1]**2 + 
				transdiff[:,2]**2)
		self.parameters['translation_totdist'] = numpy.sum(transdist)
		return self.parameters['translation_totdist']

class SequenceData(batchbrain.base.ProcessingObject):
	"""Storage object for a sequence."""

	def __init__(self,**kwargs):
		'''Initialise a SequenceData object, setup data structures.'''
		self.acquisitions = collections.OrderedDict()
		self.checkInputs(kwargs)

		# Shortcut
		self.acq = self.acquisitions
		return

	def chooseLeastMover(self):
		score = []
		potacqs = self.acquisitions.keys()
		for acq in potacqs:
			score.append(self.acq[acq].summariseTranslation())
		# Pick the min
		ind = score.index(min(score))
		return self.acq[potacqs[ind]]

	def resolveAcquisition(self,acqrule):
		'''Return a list of acquisitions according to a rule.'''
		potacqs = self.acquisitions.keys()
		if acqrule is int:
			return [self.acq[potacqs[acqrule]]]
		if acqrule == 'all':
			return [self.acq[x] for x in self.acq]
		if acqrule == 'first':
			return [self.acq[potacqs[0]]]
		if acqrule == 'last':
			return [self.acq[potacqs[-1:]]]
		if acqrule == 'leasttranslation':
			return [self.chooseLeastMover()]
		if acqrule == 'leastrotation':
			raise Exception('No support yet for rotation!')
		raise Exception('Unknown acqrule: %s' % acqrule)
		return

class GroupData(batchbrain.base.ProcessingObject):
	"""Storage object for a group analysis."""

	def __init__(self,groupdir,**kwargs):
		"""Make a group analysis directory and set up internal data
		attributes."""
		self.path = groupdir
		self.name = os.path.split(groupdir)[1]
		self.subdata = collections.OrderedDict()
		self.checkInputs(kwargs)
		self.mkAnaDir(groupdir)
		return

	def addSubjectToGroup(self,sub,targetsequence,inprefix='', \
			inext='.nii',acqrule='last',outprefix='',outext='.nii'):
		"""Include a new subject in the group."""
		# Find the acquisition we want
		acq = sub.seq[targetsequence].resolveAcquisition(acqrule)
		if len(acq) > 1:
			raise Exception('No support for multivol per sub group analysis')
		acq = acq[0]
		# Find the volume
		invol = acq.findFiles(inprefix,inext)
		if len(invol) != 1:
			raise Exception('Expected 1 invol, got:\n%s' % invol)
		invol = invol[0]
		# Get outvol. getOutput method would be ideal here but can't
		# make it general enough
		fpath,ffile = os.path.split(invol)
		ffn,fext = os.path.splitext(ffile)
		outvol = os.path.join(self.path,
			outprefix+sub.name+'_'+ffn+outext)
		alreadydone = self.checkOutput(outvol)
		if alreadydone:
			if self.overwrite:
				pass
			elif not self.update:
				raise Exception('%s already exists.' % outvol)
		# Copy volume across
		shutil.copy(invol,outvol)
		return

class SubjectData(batchbrain.base.ProcessingObject):
	"""Storage object for individual subjects. Used for accessing paths and
	subject info (e.g. demographics)."""

	def __init__(self,subjdir,demographics=None,**kwargs):
		"""Make a subject directory and set up internal data attributes."""
		# Initialise other bits
		self.path = subjdir
		self.name = os.path.split(self.path)[1]
		self.sequences = collections.OrderedDict()
		self.demographics = demographics
		self.checkInputs(kwargs)
		self.mkAnaDir(subjdir)
		# TODO: update subject list in super. Not sure how to do this.

		# Shortcut
		self.seq = self.sequences
		return

class StudyData(batchbrain.base.ProcessingObject):
	"""Master object that keeps track of the data in the study."""

	def __init__(self,anadir,projectcodes=None,**kwargs):
		'''Initialise a new dataset, create main study directory.'''
		# Check if any interesting flags came in (e.g., overruling
		# overwrite setting)
		self.checkInputs(kwargs)
		self.output('Initialising StudyData object')
		# Set up main data directory
		self.mkAnaDir(anadir)
		self.anadir = anadir
		self.subjects = collections.OrderedDict()
		self.projectcodes = projectcodes
		self.groups = {}

		# Keep track of functions and inputs to map - by order of addition
		self.tasks = collections.OrderedDict()

		# shortcut
		self.sub = self.subjects
		return

	def addJobs(self,tn,procfunc,arguments,postrunfun=None):
		"""Add a given prepared process to the tasklist."""
		if not self.tasks.has_key(tn):
			# Initialise new task
			self.tasks[tn] = Task(tn,procfunc,update=self.update,
					overwrite=self.overwrite,verbose=self.verbose,
					postrunfun=postrunfun)
		else:
			# Ensure the user isn't doing something odd
			if not self.tasks[tn].func == procfunc:
				raise Exception('procfunc mismatch: old %s new %s' % (
					self.tasks[tn].func,procfunc))
		self.tasks[tn].updateJobs(arguments)
		return

	def printSummary(self):
		'''TODO: collection of n subjects by sequence etc.'''
		self.output('STUDY SUMMARY:')
		self.output('%d subjects' % len(self.subjects))
		self.output('%d sequences' % len(self.sequences))
		self.output('SUBJECT BY SEQUENCE BREAKDOWN:')
		# Longest sequence
		totlen = len(max(self.sequences+self.subjects,key=len))
		cn = len(self.sequences)+1
		rn = len(self.subjects)
		collabs = [s.center(totlen) for s in ['']+self.sequences]
		# First label row
		self.output('%s\t'*cn % tuple(collabs))
		for sub in self.subjects:
			# first the label
			row = sub.ljust(totlen)
			# Then each sequence
			for seq in self.sequences:
				if self.countdict[sub].has_key(seq):
					num = self.countdict[sub][seq]
				else:
					num = 0
				row += '\t' + str(num).center(totlen)
			self.output(row)
		return

	def addSubject(self,sub,demographics=None):
		'''Add a list of subjects to the object and initialise
		directories.'''
		subdir = os.path.join(self.anadir,sub)
		if demographics is None:
			demographics = {}
		self.subjects[sub] = SubjectData(subdir,demographics=demographics,
				overwrite=self.overwrite,
				update=self.update,
				verbose=self.verbose)
		return

	def getPaths(self,sub,seq,acq,prefix='',ext='.nii'):
		'''Retrieve a list of paths to volumes.'''
		# Get acq objects in a list
		acqs = self.sub[sub].seq[seq].resolveAcquisition(acq)
		vols = []
		for a in acqs:
			vols += a.findFiles(prefix,ext)
		return vols

	def pathToAcquisition(self,inpath):
		"""Deduce the likely acq object from a file path."""
		# Figure out a subject
		sub = self.findListItemInStr(inpath,self.subjects.keys())
		if len(sub) > 1:
			raise Exception('Multiple possible subjects:\n%s' % \
					sub)
		sub = sub[0]
		seq = self.findListItemInStr(inpath,self.sub[sub].seq.keys())
		if len(seq) > 1:
			raise Exception('Multiple possible sequences:\n%s' % \
					seq)
		seq = seq[0]
		acq = self.findListItemInStr(inpath,self.sub[sub].seq[seq].acq.keys())
		if len(acq) > 1:
			raise Exception('Multiple possible acquisitions:\n%s' % \
					acq)
		acq = acq[0]
		return self.sub[sub].seq[seq].acq[acq]

	def getSubjectsByDemographics(self,key,onlyvals=None):
		"""Return a dict of lists of SubjectData objects. The key will
		correspond to the unique entries in key."""
		outdict = {}
		for sub in self.subjects.values():
			if not sub.demographics.has_key(key):
				raise Exception('All subjects must have key: %s' % key)
			outk = sub.demographics[key]
			if onlyvals is not None and not outk in onlyvals:
				continue
			if outdict.has_key(outk):
				outdict[outk].append(sub)
			else:
				outdict[outk] = [sub]
		return outdict

	def compareSubjectsByDemographics(self,dvkey,ivkey,onlyivs=None):
		"""Use a between-samples T or F to compare subjects on scores in
		dvkey, sorted into groups by ivkey."""
		groups = self.getSubjectsByDemographics(ivkey)
		gdata = {}
		# Restrict to only certain iv levels. Useful for post hoc tests.
		if onlyivs is not None:
			for k in groups.keys():
				if not k in onlyivs:
					del groups[k]
		# Make a list for each group
		for g in groups:
			gdata[g] = []
			for sub in groups[g]:
				if sub.demographics[dvkey] == '':
					self.output('No demographic %s for %s (group %s), skipping' %(
						sub.name,dvkey,g))
					continue
				gdata[g].append(sub.demographics[dvkey])
		ngroups = len(gdata.keys())
		if ngroups < 2:
			raise Exception('Need more than 1 group!')
		elif ngroups == 2:
			compfun = scipy.stats.ttest_ind
			test = 'T'
		else:
			compfun = scipy.stats.f_oneway
			test = 'F'
		inference = dict(zip([test,'p'],compfun(*gdata.values())))
		for g in gdata.keys():
			inference[g+'raw'] = gdata[g]
			inference[g+'mean'] = numpy.mean(gdata[g])
			inference[g+'std'] = numpy.std(gdata[g])
			inference[g+'n'] = len(gdata[g])
			inference[g+'sterr'] = numpy.std(gdata[g]) / \
					numpy.sqrt(len(gdata[g]))
		return inference

	def getSubjectDemographics(self,key,subjectnames=None):
		"""Return a dict where each subject is a key, and the value is the
		entry in demographics[key]."""
		if subjectnames is None:
			subjectnames = self.subjects.keys()
		outdict = {}
		for subname in subjectnames:
			subobj = self.subjects[subname]
			if not subobj.demographics.has_key(key):
				raise Exception('All subjects must have key: %s' % key)
			outdict[subname] = subobj.demographics[key]
		return outdict

	def parametersToDemographics(self,seqname,acqrule='all'):
		"""Update each subject's demographics dict with the keys from the
		parameters dict in acq."""
		for s in self.sub.values():
			seq = s.seq[seqname]
			acqs = seq.resolveAcquisition(acqrule)
			for a in acqs:
				for k_in,v_in in a.parameters.iteritems():
					# Ugly conditional to avoid inconsistent column labels
					if acqrule == 'leasttranslation':
						k_out = '%s_leasttran_%s' % (seqname,k_in)
					else:
						k_out = '%s_%s' % (a.name,k_in)
					s.demographics[k_out] = v_in
		return

	def exportDemographics(self,outfile):
		"""Write a CSV table with demographics for all subjects."""
		F = open(outfile,'wb')
		W = csv.writer(F,dialect='excel')
		# First, make a master dict
		mdict = {}
		for subname,s in self.subjects.iteritems():
			for k,v in s.demographics.iteritems():
				# Skip data that isn't appropriate
				if type(v) is list or type(v) is numpy.ndarray:
					continue
				# Remove empty keys
				if k == '':
					continue
				if mdict.has_key(k):
					mdict[k][subname] = str(v)
				else:
					mdict[k] = {}
					mdict[k][subname] = str(v)
		# Now we know the labels
		collabs = mdict.keys()
		collabs.sort()
		# Begin writing out
		W.writerow(['bb_subname']+collabs)
		for subname in self.subjects.keys():
			row = [subname]
			for k in collabs:
				try:
					d = mdict[k][subname]
				except KeyError:
					d = ''
				row.append(d)
			W.writerow(row)
		del W
		F.close()
		return

	def addGroup(self,groupdir):
		'''Initialise a group analysis.'''
		if self.groups.has_key(groupdir):
			if self.overwrite:
				pass
			elif not self.update:
				raise Exception('%s already exists!' % groupdir)
			self.output('%s already exists, skipping...' % groupdir)
			return
		self.groups[groupdir] = GroupData(os.path.join(self.anadir,
			groupdir), overwrite=self.overwrite,  update=self.update,
			verbose=self.verbose)
		return

# An alternative function might be CopyData that simply grabs niftis from
# somewhere
class CBUMRI(StudyData):
	"""StudyData-derived object with functionality for CBU's
	mridata system."""

	def __init__(self,anadir,projectcodes=None,**kwargs):
		# Initialise with standard StudyData
		StudyData.__init__(self,anadir,projectcodes=projectcodes,\
				**kwargs)

		# sub-class specific information
		self.rootdir='/mridata/cbu'
		self.ncharignorefn=11
		self.indsseriesnfn=[7,10]
		self.indprojcodeneg=-7

		# Check if anything else came through
		self.checkInputs(kwargs)
		return

	def importCBUData(self,series,nslices=None):
		"""Iterate over subjects, converting DICOMs to NIFTI and copying to
		anadir. Also initialise SequenceData and AcquisitionData objects as
		necessary."""
		for sub in self.subjects.keys():
			# Find data
			# (need underscore because sometimes subject codes are LONGER than
			# meant to be)
			subdir = glob.glob(os.path.join(self.rootdir,sub+'_*'))

			# Catch bad data
			if not len(subdir):
				raise Exception('No subject dir found: %s' % sub)
			if len(subdir) > 1:
				if self.projectcodes is None:
					raise Exception('Duplicate subject directories:\n%s' % subdir)
				# Need to work out if remaining dirs belong to project
				# Attempt to disambiguate by project code
				check = numpy.array([s[self.indprojcodeneg:] in self.projectcodes for s in
					subdir],dtype='bool')
				subdir = list(numpy.array(subdir)[check])
				if not len(subdir):
					raise Exception('No subject dirs match project codes:\n%s'
							% subdir)

			# Now we have a list of subdir, usually len==1, sometimes more
			# Collect all series for this subject
			self.output('Identifying series in\n%s' % subdir)
			# Ascending project numbers should correspond to order of acq.
			subdir.sort()
			seriesdirs = []
			for sd in subdir:
				dirs = glob.glob(os.path.join(
					self.rootdir,sd,'*','*'))
				# Ensure in order of acquisition
				dirs.sort()
				seriesdirs += dirs
			# Now exact match to end of fn
			goodseries = [i for i in seriesdirs if i[-len(series):]==series]

			# Initialise empty list of acquisitions
			if self.sub[sub].seq.has_key(series):
				if self.overwrite:
					pass
				elif not self.update:
					raise Exception('%s and %s already exists!' % (sub,series))
			else:
				self.sub[sub].seq[series] = SequenceData(update=self.update,
						overwrite=self.overwrite,verbose=self.verbose)

			for ind,gs in enumerate(goodseries):
				# Check for completeness
				if nslices is not None:
					ndcm = os.listdir(gs)
					if len(ndcm) != nslices:
						self.output('%s is incomplete, skipping...' % gs)
						continue
				fn = '%s_acq%03d' % (series,ind+1)
				seriesoutdir = os.path.join(self.sub[sub].path,fn)
				self.sub[sub].seq[series].acq[fn] = AcquisitionData(
						seriesoutdir,update=self.update,overwrite=self.overwrite,
						verbose=self.verbose)
				# Don't reimport dicoms for existing series
				if self.sub[sub].seq[series].acq[fn].imported:
					continue
				self.output('converting %s' % gs)
				# TEMP HACK CODE
				cmd = 'dcm2nii -o %s %s %s' % (seriesoutdir,
						'-a N -c N -d N -e N -g N -i N -f N -p Y',
						gs)
				os.system(cmd)
				#mricron.dcm2nii(gs,seriesoutdir)
			self.output('finished conversion for %s' % sub)
