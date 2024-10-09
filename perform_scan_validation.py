# Process validation tests on COMBINE datacards for the Z'->emu scan
import os
import argparse
import ROOT as rt
from array import array
from math import exp
from math import log10
from multiprocessing import Process

#----------------------------------------------------------------------------------------
# Define the sorting of the datacards
def file_sort(f):
   return int(f.split('_mp')[1].replace('.txt',''))

#----------------------------------------------------------------------------------------
# Process a single mass point
def process_datacard(card, directory, name, test = 'bias', toys = 100, signal_rate = 0., tag = '', verbose = 0):
   if not os.path.isfile(directory + card):
      print "Card %s not found" % (card)
   if verbose > -1: print 'Processing mass point', name, '(card =', card+')'

   base_dir = os.getenv('CMSSW_BASE') + '/src/ZLFV_fits'

   #----------------------------------------------------------------------------
   # Perform the validation test
   #----------------------------------------------------------------------------

   results = []
   if test == 'bias':
      command = base_dir + '/tests/bemu_bias.sh ' + card + ' -t ' + str(toys) + ' -g ' + str(toys) + ' --name ' + name + ' -r 30 --skipplots'
      if signal_rate > 0.: command += ' --genarg "--expectSignal ' + str(signal_rate) + '"'
      if tag != "":  command += ' --tag ' + tag
      output = 'fit_bias_%s%s.log' % (name, tag)
      if verbose > 1: print command
      os.system('cd %s; %s >| %s; cd ..' % (directory, command, output))

#----------------------------------------------------------------------------------------
# Retrieve fit information for a single mass point
def retrieve_info(card, directory, name, test = 'bias', signal_rate = 0., tag = '', figdir = '', verbose = 0):
   results = []
   if test == 'bias':
      fit_file = directory + 'fitDiagnostics.' + name + '_closure_test' + tag + '.root'
      f = rt.TFile.Open(fit_file, 'READ')
      t = f.Get('tree_fit_sb')
      t.Draw("r >> h")
      h = rt.gPad.GetPrimitive('h')
      mean = h.GetMean() - signal_rate
      c = rt.TCanvas()
      h.Draw('hist')
      h.SetLineWidth(2)
      c.SaveAs(figdir+name+'_r.png')

      h = rt.TH1F('hpull', 'Pull distribution', 30, -3, 3)
      t.Draw('(r - %.3f) / (r < 0. ? rHiErr : rLoErr) >> hpull' % (signal_rate), 'rHiErr > 0. && rLoErr > 0.')
      func = rt.TF1("func", "[2]*TMath::Gaus(x, [0], [1])", -3, 3);
      func.SetParameters(0., 1., h.Integral())
      h.Fit(func, 'R')
      pull  = func.GetParameter(0)
      width = func.GetParameter(1)
      # pull = h.GetMean()
      # width = h.GetStdDev()
      h.Draw('hist')
      h.SetLineWidth(2)
      func.Draw('same')
      c.SaveAs(figdir+name+'_bias.png')
      f.Close()
      results = [mean, pull, width]

   #----------------------------------------------------------------------------
   # Extract the mass
   #----------------------------------------------------------------------------

   mass = float(card.split('mass-')[1].split('_')[0])
      
   #----------------------------------------------------------------------------
   # Return the results
   #----------------------------------------------------------------------------

   return [results, mass]

#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="bdt_v01", type=str,help="datacard directory name")
parser.add_argument("-j", "--nthreads", dest="nthreads",default=8,type=int,help="Number of threads to process using")
parser.add_argument("--test", dest="test",default="bias", type=str,help="Validation test to perform: bias, GOF, impact")
parser.add_argument("-t", dest="toys",default=100, type=int,help="Toys for validation tests")
parser.add_argument("--signal-rate", dest="signal_rate",default=0., type=float, help="Signal rate to inject into the test")
parser.add_argument("--skip-fits", dest="skip_fits",default=False, action='store_true',help="Skip fits, assume already processed")
parser.add_argument("--max-steps", dest="max_steps",default=-1, type=int, help="Maximum steps to take in the scan")
parser.add_argument("--first-step", dest="first_step",default=0, type=int, help="First mass step to process")
parser.add_argument("--mass-point", dest="mass_point",default=-1, type=int, help="Mass point to process")
parser.add_argument("--card-tag", dest="card_tag",default="", type=str, help="Card name tag to process")
parser.add_argument("--tag", dest="tag",default="", type=str, help="Output directory tag")
parser.add_argument("-v","--verbose", dest="verbose",default=0, type=int,help="Add verbose printout")

args, unknown = parser.parse_known_args()


##### configuaration ########

### check input flags
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()

#----------------------------------------------
# Setup the tests
#----------------------------------------------

### default path
if args.tag != '': args.tag = '_' + args.tag
figdir = "./figures/val_%s_%s%s/" % (args.name, args.test, args.tag)
carddir = "./datacards/%s/" % (args.name)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
os.system("[ ! -d %s ] && mkdir -p %s" % (carddir, carddir))

rt.gROOT.SetBatch(True)

#----------------------------------------------
# Retrieve the cards to process
#----------------------------------------------

list_of_files = [f for f in os.listdir(carddir) if '.txt' in f]
# If not using a BDT score region tag, only take the merged files
if args.card_tag == "":  list_of_files = [f for f in list_of_files if '0d7' not in f]
else:                    list_of_files = [f for f in list_of_files if args.card_tag in f]

if len(list_of_files) == 0:
   print "No card files found!"
   exit()

# Sort the list by mass point
list_of_files.sort(key=file_sort)
if args.first_step > 0:
   list_of_files = list_of_files[args.first_step:]
if args.max_steps > 0 and args.max_steps < len(list_of_files):
   list_of_files = list_of_files[:args.max_steps]

#----------------------------------------------
# Perform the validation processing
#----------------------------------------------

# List of results
masses = array('d')
val_results = []

jobs = [] # For multithreaded processing

if not args.skip_fits:
   for f in list_of_files:
      mass_point = f.split('_mp')[1].split('.txt')[0]
      if args.mass_point >= 0 and mass_point != str(args.mass_point): continue
      job = Process(target = process_datacard, args=(f, carddir, args.name+ '_mp'+mass_point, args.test, args.toys, args.signal_rate, args.tag, args.verbose))
      jobs.append(job)
   if args.nthreads < 1: args.nthreads = 1
   print("Parallel processing using %i threads" % (args.nthreads))
   for ithread in range(0,len(jobs),args.nthreads):
      nthread = args.nthreads+ithread
      if (nthread>len(jobs)):
         nthread=len(jobs)
      print("Processing threads %i to %i out of %i jobs: %5.1f%% processed" % (ithread, nthread-1, len(jobs), (ithread*100./len(jobs))))
      for job in jobs[ithread:nthread]:
         job.start()
      for job in jobs[ithread:nthread]:    
         job.join()

for f in list_of_files:
   mass_point = f.split('_mp')[1].split('.txt')[0]
   if args.mass_point >= 0 and mass_point != str(args.mass_point): continue
   [result, mass] = retrieve_info(f, carddir, args.name+ '_mp'+mass_point, args.test, args.signal_rate, args.tag, figdir, args.verbose)

   # store the results
   masses.append(mass)
   print result
   val_results.append(result)

if args.mass_point >= 0: exit() # don't make a plot for the single point

#----------------------------------------------
# Combine and plot the results
#----------------------------------------------

rt.gStyle.SetOptStat(0)

if args.test == 'bias':
   h = rt.TH1D('hpull', 'Z prime mass fit pulls', 40, -1.5, 1.5)
   for bias in [res[1] for res in val_results]: h.Fill(bias)
   c = rt.TCanvas()
   h.Draw('hist')
   h.SetLineWidth(2)
   h.SetFillStyle(3005)
   h.SetFillColor(rt.kBlue)
   h.SetXTitle('Pull')
   c.SaveAs(figdir+'pulls.png')

   for index in range(len(masses)):
      if abs(val_results[index][1]) > 0.3:
         print 'Mass point', masses[index],'(index', index,'): Large bias, mean =', val_results[index][1], 'width =', val_results[index][2]

