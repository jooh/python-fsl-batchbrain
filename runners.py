'''Map functions/objects that take the particular input from StudyData.tasks and
runs with it.'''
import sys
import os
import random

def simpleRunner(sd):
	'''Simplest possible task runner: go through in serial order on current
	machine.'''
	nsteps = len(sd.tasks)
	if nsteps == 0:
		sd.output('No outstanding processing steps in tasks')
		return
	sd.output('Running %d processing steps...' % len(sd.tasks))
	for tn,to in sd.tasks.iteritems():
		sd.output('Running %s with %d jobs' % (tn,len(to.jobs)))
		outputs = map(to.func,to.jobs)
		to.postRunUpdate(outputs)
	# Clear out the task list 
	sd.output('Done.')
	sd.tasks.clear()
	return

class IPClusterRunner():
	'''CBU-specific code for starting a IPython cluster and running
	jobs.'''

	def __init__(self,nengines=None,targetmachines=None):
		'''Initialise the cluster, either by re-connecting to a current
		cluster or starting a new one.'''
		import IPython.kernel
		# Pull out the usual CBU machines
		if targetmachines is None:
			self.availablemachines = getCBUMachines()
		else:
			self.availablemachines = targetmachines
		# The machines get started without the appropriate path settings,
		# so need to bring the current machine's over
		self.localpath = os.environ['PATH']

		# In the absence of true load balance, we can at least randomise
		random.shuffle(self.availablemachines)
		self.mlist = self.availablemachines.copy()

		self.pophandles = []

		# See if we already have a controller
		try:
			rc = IPython.parallel.Client()
		except TimeOutError:
			# Otherwise give me one
			self.pophandles.append(self.startMachine(
					self.randomEngine(),iptype='ipcontroller'))
			rc = IPython.parallel.Client()
		nrunning = len(rc.ids)

		if nengines is None:
			if nrunning > 0:
				# Assume you already have the engines you want running
				return
			else:
				raise Exception('No engines available. Must \
						specify nengines>0')

		# Time to start some engines
		while len(self.pophandles)+1 < nengines:
			self.pophandles.append(self.startMachine(
					self.randomEngine(),iptype='ipengine'))
			# Need to slow down here so the ipcontroller can keep up
			time.sleep(1)
		return
	#### TODO TODO - map function from IP, make consistent with seq run
	
	def __call__(sd):
		'''Do a parallel map on each processing step in sequence.'''
		nsteps = len(sd.tasks)
		if nsteps == 0:
			sd.output('No outstanding processing steps in tasks')
			return
		sd.output('Running %d processing steps...' % len(sd.tasks))
		for func,args in sd.tasks.iteritems():
			self.map(func,args)
		sd.tasks.clear()
		return

	def randomEngine(self):
		try:
			return self.mlist.pop()
		except IndexError:
			self.mlist = self.availablemachines.copy()
			random.shuffle(self.mlist)
			return self.mlist.pop()
			
	def getCBUMachines(self):
		sys.path.append('/imaging/local/spm/loadshare')
		import loadsharesettings
		machines = loadsharesettings.machines()
		# Only machines 43 and above are 64bit
		return [i for i in machines if int(i[1:])>=43]

	def startMachine(self,id,iptype='ipengine'):
		import subprocess
		"""Start an ipengine or ipcontroller. Return Popen handle for xterm
		process."""
		xcmd = 'xterm -e ssh %s -x "' % id
		pathcmd = 'setenv PATH %s;' % self.localpath
		ipcmd = '%s"' % iptype
		fullcmd = xcmd + pathcmd + ipcmd
		return subprocess.Popen(fullcmd,shell=True)
