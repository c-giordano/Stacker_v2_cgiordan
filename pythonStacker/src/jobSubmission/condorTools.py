#
# functionality for making job description files for the condor submission system #
#
# documentation:
# - T2B wiki: https://t2bwiki.iihe.ac.be/HTCondor
# - HTCondor manual: https://htcondor.readthedocs.io/en/latest/users-manual/index.html

# general use:
# two ingredients are needed for a condor job:
# - a job description (.txt) file
# - an executable
# the functions in this tool allow creating an executable bash script (similar as for qsub)
# and its submission via a job description file

import os
import sys


def makeUnique(fname, scriptfolder=""):
    # make a file name unique by appending a number to it,
    # e.g. test.txt -> test1.txt (in case test.txt already exists)
    # Next part commented out, for consistency keep names always with numbers
    # if not os.path.exists(fname):
    #     return fname
    [name, ext] = os.path.splitext(fname)
    app = 0
    while app < 2500:
        tryname = name + str(app) + ext
        if not os.path.exists(os.path.join(scriptfolder, tryname)):
            return tryname
        app += 1
    print('# ERROR #: already 2500 files named {} exist.'.format(fname))
    print(' consider choosing more specific names, splitting in folders, etc.')
    sys.exit()


def initJobScript(name, scriptfolder="", cwd=None, hardbreak=False):
    # initialize an executable bash script by setting correct cms env
    # but copied here to be more standalone
    # parse argument
    name = os.path.splitext(name)[0]
    fname = os.path.join(scriptfolder, name + '.sh')
    if os.path.exists(fname):
        os.system('rm {}'.format(fname))
    if not os.path.exists(scriptfolder):
        os.makedirs(scriptfolder)
    if not cwd:
        cwd = os.path.abspath(os.getcwd())
    # write script
    with open(fname, 'w') as script:
        script.write('#!/bin/bash\n')
        if hardbreak:
            script.write('set -e\n')
        script.write(f'cd /user/nivanden/{os.getenv("CMSSW_VERSION")}/src\n')
        script.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        script.write('eval `scram runtime -sh`\n')
        script.write('export X509_USER_PROXY=/user/$USER/x509up_u$(id -u $USER)\n')
        script.write('cd {}\n'.format(cwd))
        script.write('echo "{}"\n\n'.format(fname))

    os.system('chmod +x ' + fname)
    print('initJobScript created {}'.format(fname))


def makeJobDescription(name, exe, argstring=None, stdout=None, stderr=None, log=None,
                       cpus=1, mem=1024, disk=10240, jdName=None, scriptfolder=""):
    # create a single job description txt file
    # note: exe can for example be a runnable bash script
    # note: argstring is a single string containing the arguments to exe (space-separated)
    # parse arguments:
    name = os.path.splitext(name)[0]

    if jdName:
        fname = jdName
    else:
        fname = name + '.sub'

    pathToFname = os.path.join(scriptfolder, fname)
    if os.path.exists(pathToFname):
        os.system('rm {}'.format(pathToFname))
    if stdout is None:
        stdout = name + '_$(ClusterId)_$(ProcId).out'
    if stderr is None:
        stderr = name + '_$(ClusterId)_$(ProcId).err'
    if log is None:
        log = name + '_$(ClusterId)_$(ProcId).log'
    # write file
    with open(pathToFname, 'w') as f:
        f.write('executable = {}\n'.format(os.path.join(scriptfolder, exe)))
        if argstring is not None:
            f.write('arguments = "{}"\n\n'.format(argstring))
        f.write(f"output = /user/{os.environ['USER']}/condor/output/{stdout}\n")
        f.write(f"error = /user/{os.environ['USER']}/condor/error/{stderr}\n")
        f.write(f"log = /user/{os.environ['USER']}/condor/logs/{log}\n\n")
        f.write("should_transfer_files = NO\n\n")
        f.write("job_machine_attrs = Machine\n")
        f.write("job_machine_attrs_history_length = 5\n")
        f.write("requirements = target.machine =!= MachineAttrMachine1 && target.machine =!= MachineAttrMachine2 && target.machine =!= MachineAttrMachine3  && target.machine =!= MachineAttrMachine4\n")

        f.write("max_retries = 5\n\n")
        #f.write("max_retries = 10\n\n")
        # f.write('request_cpus = {}\n'.format(cpus)) # Don't specify if not necessary
        # f.write('request_memory = {}\n'.format(mem)) # Don't specify if not necessary
        # f.write('request_disk = {}\n\n'.format(disk)) # Don't specify if not necessary
        # f.write('should_transfer_files = yes\n\n')
        # (not fully sure whether to put 'yes', 'no' or omit it completely)
        # f.write('queue\n\n')
    print('makeJobDescription created {}'.format(pathToFname))


def submitCondorJob(jobDescription, addArgs="", scriptfolder=""):
    # submit a job description file as a condor job
    fname = os.path.join(scriptfolder, os.path.splitext(jobDescription)[0] + '.sub')
    if not os.path.exists(fname):
        print('# ERROR #: job description file {} not found'.format(fname))
        sys.exit()
    # maybe later extend this part to account for failed submissions etc!
    os.system('condor_submit {} {}'.format(fname, addArgs))


def submitCommandAsCondorJob(name, command, stdout=None, stderr=None, log=None,
                        cpus=1, mem=1024, disk=10240, scriptfolder=""):
    # submit a single command as a single job
    # command is a string representing a single command (executable + args)
    submitCommandsAsCondorJobs(name, [[command]], stdout=stdout, stderr=stderr, log=log,
                               cpus=cpus, mem=mem, disk=disk, scriptfolder=scriptfolder)


def submitCommandsAsCondorCluster(name, commands, stdout=None, stderr=None, log=None,
                                  cpus=1, mem=1024, disk=10240, hardbreak=False):
    # run several similar commands within a single cluster of jobs
    # note: each command must have the same executable and number of args, only args can differ!
    # note: commands can be a list of commands (-> a job will be submitted for each command)
    # parse arguments
    name = os.path.splitext(name)[0]
    shname = makeUnique(name + '.sh')
    jdname = makeUnique(name + '.sub')
    [exe,argstring] = commands[0].split(' ',1) # exe must be the same for all commands
    nargs = len(argstring.split(' ')) # nargs must be the same for all commands
    # first make the executable
    initJobScript(shname, hardbreak=hardbreak)
    with open(shname,'a') as script:
        script.write(exe)
        # for i in range(nargs): script.write(' ${}'.format(i+1))
        script.write(' $@')
        script.write('\n')
    # then make the job description
    # first job:
    makeJobDescription(name, shname, argstring=argstring, stdout=stdout, stderr=stderr, log=log,
                       cpus=cpus, mem=mem, disk=disk, jdName=jdname)
    # add other jobs:
    with open(jdname, 'a') as script:
        for command in commands[1:]:
            [thisexe, thisargstring] = command.split(' ', 1)
            thisnargs = len(thisargstring.split(' '))
            if (thisexe != exe or thisnargs != nargs):
                print('# ERROR #: commands are not compatible to put in same cluster')
                return

            script.write('arguments = "{}"\n'.format(thisargstring))
            script.write('queue\n\n')
    # finally submit the job
    submitCondorJob(jdname)


def submitCommandsAsCondorJob(name, commands, stdout=None, stderr=None, log=None,
                              cpus=1, mem=1024, disk=10240):
    # submit a set of commands as a single job
    # commands is a list of strings, each string represents a single command (executable + args)
    # the commands can be anything and are not necessarily same executable or same number of args.
    submitCommandsAsCondorJobs(name, [commands], stdout=stdout, stderr=stderr, log=log,
                               cpus=cpus, mem=mem, disk=disk)


def submitCommandsAsCondorJobs(name, commands, stdout=None, stderr=None, log=None,
                               cpus=1, mem=1024, disk=10240, scriptfolder="", hardbreak=False):
    # submit multiple sets of commands as jobs (one job per set)
    # commands is a list of lists of strings, each string represents a single command
    # the commands can be anything and are not necessarily same executable or number of args.
    for commandset in commands:
        # parse arguments
        name = os.path.splitext(name)[0]
        shname = makeUnique(name + '.sh', scriptfolder=scriptfolder)
        jdname = makeUnique(name + '.sub', scriptfolder=scriptfolder)
        # first make the executable
        initJobScript(shname, scriptfolder=scriptfolder, hardbreak=hardbreak)
        with open(os.path.join(scriptfolder, shname), 'a') as script:
            for cmd in commandset:
                script.write(cmd + '\n')
        # then make the job description
        makeJobDescription(name,shname,stdout=stdout,stderr=stderr,log=log,
                           cpus=cpus,mem=mem,disk=disk,jdName=jdname, scriptfolder=scriptfolder)
        with open(os.path.join(scriptfolder, jdname), 'a') as jdScript:
            jdScript.write('queue\n\n')

        # finally submit the job
        submitCondorJob(jdname, addArgs='-batch-name "{}"'.format(name), scriptfolder=scriptfolder)


def submitCommandsetsAsCondorCluster(name, commands, stdout=None, stderr=None, log=None,
                                     cpus=1, mem=1024, disk=10240, scriptfolder="", cwd=None, hardbreak=False):
    # submit multiple sets of commands as one cluster (one job of the clusterper set)
    # commands is a list of lists of strings, each string represents a single command
    # the commands can be anything and are not necessarily same executable or number of args.
    # parse arguments
    name = os.path.splitext(name)[0]
    shname = makeUnique(name + '.sh', scriptfolder)
    jdname = makeUnique(name + '.sub', scriptfolder)
    # first make the executable
    print(shname)
    initJobScript(shname, scriptfolder, cwd=cwd, hardbreak=hardbreak)
    with open(os.path.join(scriptfolder, shname), 'a') as script:
        for i, commandset in enumerate(commands):
            script.write(f"if [ $1 -eq {i} ]; then\n")
            for cmd in commandset:
                script.write("  " + cmd + '\n')
            script.write("fi\n")

    # then make the job description
    # first job:
    makeJobDescription(name, shname, stdout=stdout, stderr=stderr, log=log,
                       cpus=cpus, mem=mem, disk=disk, jdName=jdname, scriptfolder=scriptfolder)

    with open(os.path.join(scriptfolder, jdname), 'a') as jdScript:
        jdScript.write("arguments = $(ProcId)\n")
        jdScript.write(f'queue {len(commands)}\n\n')

        # for commandset in commands[1:]:
        #     # parse arguments
        #     name = os.path.splitext(name)[0]
        #     shname = makeUnique(name + '.sh', scriptfolder=scriptfolder)
        #     # first make the executable
        #     initJobScript(shname, scriptfolder)
        #     with open(os.path.join(scriptfolder, shname), 'a') as script:
        #         for cmd in commandset:
        #             script.write(cmd + '\n')
        #     # then make the job description
        #     jdScript.write('executable = {}\n'.format(os.path.join(scriptfolder, shname)))
        #     jdScript.write('queue\n\n')
    # finally submit the job
    submitCondorJob(jdname, addArgs='-batch-name "{}"'.format(name), scriptfolder=scriptfolder)


def submitScriptAsCondorJob(scriptName):
    jdName = os.path.splitext(scriptName)[0]

    makeJobDescription(scriptName, scriptName)
    submitCondorJob(scriptName)
