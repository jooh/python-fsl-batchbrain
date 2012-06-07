'''Base objects and functionality.'''

import sys
import shutil
import os
import time
import glob

class ProcessingObject():
	"""Base object with the properties shared by all processing
	objects."""
	overwrite = False
	waitfornext = True
	verbose = True
	channelout = sys.stdout.writelines
	update = False
	parallel = False
	inext = '.nii'
	outext = '.nii'

	def findStrInList(self,searchstr,targetlist):
		'''Use a list comp search method to find a target string in each of
		the strings in targetlist. Numerous assumptions: 1) error if no
		matches, 2) only one match should be found.'''
		hits = [x for x in targetlist if searchstr in x]
		if not hits: 
			raise Exception('No matches for %s' % searchstr)
		return hits

	def findListItemInStr(self,searchstr,targetlist):
		'''Use a list comp search method to find a target string in each of
		the strings in targetlist. Numerous assumptions: 1) error if no
		matches, 2) only one match should be found.'''
		hits = [x for x in targetlist if x in searchstr]
		if not hits: 
			raise Exception('No matches for %s' % searchstr)
		return hits

	def output(self,text):
		'''Write out text to e.g. stdout but only if verbose.'''
		if type(text) is not list:
			text = [text]
		if self.verbose:
			self.channelout(('%s\n' % s for s in text))
		return

	def checkInputs(self,argdict):
		for attribute in argdict:
			if hasattr(self,attribute):
				#self.output('Updating attribute: %s' % attribute)
				setattr(self,attribute,argdict[attribute])
			else:
				raise Exception('Unknown argument: %s' % attribute)
		return

	def getOutput(self,invol,sd):
		"""The default getOutput method converts an invol to an outvol
		based on various object settings. Inherited objects can overwrite
		with methods based on checking attributes internally instead."""
		outpath,outfile = os.path.split(invol[0])
		outfn,outext = os.path.splitext(outfile)
		return self.addPrefix(os.path.join(outpath,outfn + self.outext)
				,self.outprefix)

	def checkOutput(self,outvol):
		"""Use the output from getOutput to decide if the task has already
		been done (ie, there is output) or not. To be overwritten with
		getOutput by processes that don't use volumes as output."""
		# Support wild card matching for the adventurous
		if '*' in outvol:
			return len(glob.glob(outvol)) >= 1
		else:
			return os.path.exists(outvol)

	def setAllInputs(self,argdict):
		for attribute in argdict:
			setattr(self,attribute,argdict[attribute])
		return

	def mkAnaDir(self,anadir):
		'''Make a directory while remaining conscious of existing dirs and
		overwrite settings.'''
		# Check if we need to worry about overwriting
		if os.path.exists(anadir):
			if self.overwrite:
				self.output('Creating %s after removing old dir' % anadir)
				shutil.rmtree(anadir)
			elif self.update:
				self.output('%s already exists, continuing' % anadir)
				return False
			else:
				raise Exception('%s already exists!' % anadir)
		os.mkdir(anadir)
		return True

	def initialiseLog(self,logfile):
		"""Make a log file that stdout gets piped to."""
		timestr = time.strftime('%y%m%d_%H%M')
		if os.path.exists(logfile) and not self.overwrite:
			F = open(logfile,'a')
		else:
			F = open(logfile,'w')
		F.write('\n\nbatchbrain on %s...\n' % timestr)
		F.close()
		return


