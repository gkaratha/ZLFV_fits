# Create a condor script to process the skimming of the ntuples
import os
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("-n"     , dest="name"      , default="v01"                 , type=str, help="Ntupling version name")
parser.add_argument("--tree" , dest="tree"      , default="mutau"               , type=str, help="Selection tree to process (mutau or etau)")
parser.add_argument("--year" , dest="year"      , default="Run2"                , type=str, help="Data period to process (2016, 2017, 2018, or Run2)")
parser.add_argument("--tag"  , dest="tag"       , default=""                    , type=str, help="Dataset tag to process")
parser.add_argument("--input", dest="input_path", default="zprime_tau_rootfiles", type=str, help="Dataset directory to process")

args, unknown = parser.parse_known_args()
if len(unknown)>0: 
    print "Arguments not found:", unknown, "exiting"
    exit()
if args.tree not in ['mutau', 'etau', 'mutau_e', 'etau_mu']:
    print 'Unknown selection', args.tree
    exit()
additional_args = '--selection %s' % (args.tree)
if args.tree == 'mutau_e' or args.tree == 'etau_mu':
    args.tree = 'emu'

#-------------------------------------------------------------------------------------
# Setup the output area
#-------------------------------------------------------------------------------------

if args.name == '': args.name = 'test'
if not os.path.exists(args.name):
    os.mkdir(args.name)
os.chdir(args.name)


#-------------------------------------------------------------------------------------
# Get the input datasets
#-------------------------------------------------------------------------------------

input_file_path = '/eos/cms/store/group/phys_smp/ZLFV/' + args.input_path
if input_file_path[-1] != '/': input_file_path += '/'
print 'Using input file path', input_file_path

input_files = [f for f in os.listdir(input_file_path) if '.root' in f and (args.tag == '' or args.tag in f) and (args.year == 'Run2' or args.year in f)]

#-------------------------------------------------------------------------------------
# Setup the condor submission configuration
#-------------------------------------------------------------------------------------

name = args.name
sub_txt = ''
sub_txt += 'executable = condor_skim_%s.sh\n' % (name)
sub_txt += 'arguments = $(ProcId)\n'
sub_txt += 'output                = skim_%s.$(ClusterId).$(ProcId).out\n' % (name)
sub_txt += 'error                 = skim_%s.$(ClusterId).$(ProcId).err\n' % (name)
sub_txt += 'log                   = skim_%s.$(ClusterId).log\n\n' % (name)

# Send the job to Held state on failure.
sub_txt += 'on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n\n'

# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.
sub_txt += 'periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n\n'
sub_txt += '+DesiredOS = "SL7"\n'
sub_txt += '+JobFlavour = "longlunch"\n\n'

sub_txt += 'queue %i' % (len(input_files))
with open('condor_skim_%s.sub' % (name),'w') as fl:
    fl.write(sub_txt)
fl.close()

#-------------------------------------------------------------------------------------
# Setup the condor processing script
#-------------------------------------------------------------------------------------

script = ''
script += '#!/bin/sh\n'
script += 'ulimit -s unlimited\n'
script += 'set -e\n'
script += 'cd %s/src\n' % (os.getenv('CMSSW_BASE'))
script += 'export SCRAM_ARCH=%s\n' % (os.getenv('SCRAM_ARCH'))
script += 'source /cvmfs/cms.cern.ch/cmsset_default.sh\n'
script += 'eval `scramv1 runtime -sh`\n'
script += 'cd %s/..\n\n' % (os.getcwd())

for index, sample in enumerate(input_files):
    year = sample.split('_')[-1].split('.')[0]
    s_name = sample.split('_')[1]
    script += 'if [ $1 -eq %i ]; then\n' % (index)
    script += '  python skim_tau_ntuple.py -o %s --input-file %s --year %s --tree-name %s %s --input-directory %s --out-dir %s\n' % (s_name, sample, year, args.tree,
                                                                                                                                     additional_args,
                                                                                                                                     args.input_path, name)
    script += 'fi\n\n'

with open('condor_skim_%s.sh' % (name),'w') as fl:
    fl.write(script)
fl.close()

condor_dir = '%s/private/grid/zprime/%s/' % (os.getenv('HOME'), name)

# if not os.path.exists('~/private/grid'): os.mkdir('~/private/grid')
# if not os.path.exists('~/private/grid/zprime'): os.mkdir('~/private/grid/zprime')
if not os.path.exists(condor_dir): os.makedirs(condor_dir)
shutil.move("condor_skim_%s.sh" % (name), "%scondor_skim_%s.sh" % (condor_dir, name))
shutil.move("condor_skim_%s.sub" % (name), "%scondor_skim_%s.sub" % (condor_dir, name))

print 'Condor directory: %s' % (condor_dir)
print 'Submit using "condor_submit condor_skim_%s.sub"' % (name)
