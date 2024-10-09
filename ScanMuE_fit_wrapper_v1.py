# Perform a scan of the emu data, fitting for a potential signal resonance
import os
import argparse
import ROOT as rt
from array import array
from signal_model import *
from ScanMuE_wrapper_helper import *
from multiprocessing import Process
rt.gInterpreter.Declare('#include "SFBDT_weight_combined.h"')

#----------------------------------------------------------------------------------------------------
# Define a function to process the signal and background fitting calls
def proc_unit(sgn_argument, bkg_argument):
    if sgn_argument !='':
        os.system(sgn_argument)
    if bkg_argument != '':
        os.system(bkg_argument)

#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="test", type=str,help="output root name")
parser.add_argument("--data-file", dest="data_file",default="", type=str,help="data file")
parser.add_argument("--xgb-min", dest="xgb_min",default="0.70", type=str,help="BDT score minimum for the category")
parser.add_argument("--xgb-max", dest="xgb_max",default="1.00", type=str,help="BDT score maximum for the category")
parser.add_argument("--scan-min", dest="scan_min",default=110., type=float,help="Minimum mass hypothesis to scan")
parser.add_argument("--scan-max", dest="scan_max",default=500., type=float,help="Maximum mass hypothesis to scan")
parser.add_argument("--min-fit-mass", dest="min_fit_mass",default=95., type=float,help="Minimum mass used in the fits")
parser.add_argument("--max-fit-mass", dest="max_fit_mass",default=700., type=float,help="Maximum mass used in the fits")
parser.add_argument("--scan-step", dest="scan_step",default=1., type=float,help="Step size in the mass scan, in units of signal core sigma")
parser.add_argument("--year", dest="year",default="Run2", type=str,help="Data period to use (2016, 2017, 2018, or Run2)")
parser.add_argument("--component", dest="component",type=str,default="all",help="Only process the given fit component: all, sgn, bkg, or none")
parser.add_argument("--unblind", dest="unblind",default=False, action='store_true',help="Unblind the fits")
parser.add_argument("--skip-shape-dc", dest="skip_shape_dc",default=False, action='store_true',help="shape experiment")
parser.add_argument("--skip-correct", dest="skip_correct",default=False, action='store_true', help="Skip sample year efficiency corrections")
parser.add_argument("--fit-version", dest="ver",default="2", type=str,help="fit version")
parser.add_argument("--outvar", dest="outvar",default="mass_ll",type=str,help="Name out the output observable")
parser.add_argument("--param-name", dest="param_name",default="bin",type=str,help="Name of the COMBINE category")
parser.add_argument("-j", "--nthreads", dest="nthreads",default=8,type=int,help="Number of threads to process using")
parser.add_argument("--skip-fit", dest="skip_fit",default=False, action='store_true',help="fit skip")
parser.add_argument("--skip-sgn-syst", dest="skip_sgn_syst",default=False, action='store_true',help="shape experiment")
parser.add_argument("--skip-bkg-altfits", dest="skip_bkg_altfits",default=False, action='store_true',help="shape experiment")
parser.add_argument("--no-scales", dest="no_scales",default=False, action='store_true',help="Ignore signal scale factors")
parser.add_argument("--skip-dc", dest="skip_dc",default=False, action='store_true',help="Skip datacard creation")
parser.add_argument("--mass-point", dest="mass_point",default=-1,type=int,help="Single mass point to process")
parser.add_argument("--full-mass", dest="full_mass",default=False,action='store_true',help="Fit the entire mass distribution")
parser.add_argument("--use-gaus", dest="use_gaus",default=False,action='store_true',help="Model the signal with a Gaussian instead of a Crystal Ball")
parser.add_argument("--use-mc-signal", dest="use_mc_signal",default=False, action='store_true', help="Use the MC signal fit to the mass point instead of the interpolation")
parser.add_argument("--skip-mc-point", dest="skip_mc_point",default="", type=str, help="MC sample to ignore in interpolation")
parser.add_argument("--log-files", dest="log_files",default=False,action='store_true',help="Write individual mass point fits to log files")
parser.add_argument("--dry-run", dest="dry_run",default=False,action='store_true',help="Don't execute script calls")

args, unknown = parser.parse_known_args()


##### configuaration ########

### check input flags
if len(unknown)>0: 
    print "not found:",unknown,"exiting"
    exit()

MaxMasspoints=-1 #X for debug; -1 for run

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
altpath="/eos/user/m/mimacken/ZEMu/CMSSW_11_3_4/src/ZLFV_fits/trees/" #FIXME: Replace these with official versions
figdir = "./figures/%s/" % (args.name)
carddir = "./datacards/%s/" % (args.name)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
os.system("[ ! -d %s ] && mkdir -p %s" % (carddir, carddir))
os.system("[ ! -d log ] && mkdir log")
if not os.path.exists(carddir+"WorkspaceScanBKG"):
    os.symlink("../../WorkspaceScanBKG", carddir+"WorkspaceScanBKG")
if not os.path.exists(carddir+"WorkspaceScanSGN"):
    os.symlink("../../WorkspaceScanSGN", carddir+"WorkspaceScanSGN")
   

### MC signal mass points

# sgn_masspoints=["100","125","150","175","200","300","400","500","600","800","1000"]
sgn_masspoints=["100","125","150","175","200","300","400","500","600"]

# Define the signal samples by mass and period
signal_samples = {
    "100" : {
        "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM100_mcRun16.root", 94800, 2016, path),
        "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM100_mcRun17.root", 97800, 2017, path),
        "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM100_mcRun18.root", 99200, 2018, path),
    },
    "125" : {
        "2016" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM125_mcRun18.root", 98400, 2018, path),
        "2017" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM125_mcRun18.root", 98400, 2018, path),
        "2018" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM125_mcRun18.root", 98400, 2018, path),
    },
    "150" : {
        "2016" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM150_mcRun18.root", 90500, 2018, path),
        "2017" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM150_mcRun18.root", 90500, 2018, path),
        "2018" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM150_mcRun18.root", 90500, 2018, path),
    },
    "175" : {
        "2016" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM175_mcRun18.root", 89900, 2018, path),
        "2017" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM175_mcRun18.root", 89900, 2018, path),
        "2018" : sample("Meas_fullAndSFAndGenParts_bdt_v7_emu_scan_sgnM175_mcRun18.root", 89900, 2018, path),
    },
    "200" : {
        "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root", 96300, 2018, path),
        "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root", 96300, 2018, path),
        "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root", 96300, 2018, path),
    },
    "300" : {
        "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM300_mcRun18.root", 99700, 2018, path),
        "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM300_mcRun18.root", 99700, 2018, path),
        "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM300_mcRun18.root", 99700, 2018, path),
    },
    "400" : {
        "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root", 97600, 2018, path),
        "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root", 97600, 2018, path),
        "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root", 97600, 2018, path),
    },
    "500" : {
        "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM500_mcRun16.root", 81300, 2016, path),
        "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM500_mcRun17.root", 97800, 2017, path),
        "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM500_mcRun18.root", 98700, 2018, path),
    },
    "600" : {
        "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root", 97600, 2018, path),
        "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root", 97600, 2018, path),
        "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root", 97600, 2018, path),
    },
    "800" : {
        "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root", 97600, 2018, path),
        "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root", 97600, 2018, path),
        "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root", 97600, 2018, path),
    },
    "1000" : {
        "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root", 97600, 2018, path),
        "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root", 97600, 2018, path),
        "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root", 97600, 2018, path),
    },
}

### Data
data_files={"2016":"Meas_full_bdt_v7_emu_scan_data_Run16.root",
            "2017":"Meas_full_bdt_v7_emu_scan_data_Run17.root",
            "2018":"Meas_full_bdt_v7_emu_scan_data_Run18.root",
            "Run2":"Meas_full_bdt_v7_emu_scan_data_Run1*.root"}

### sf
# Scale factors
sf="Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt"
sf += "*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt"
sf += "*PU_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*Prefire_wt"
# sf += "*PtZ_wt*MixZ_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)" # LFV Z-specific weights

if args.no_scales:
    sf="1"

### arange flags
unblind="false"
if args.unblind: unblind="true"
shape_dc="true"
if args.skip_shape_dc: shape_dc="false"
do_sgn_syst="true"
if args.skip_sgn_syst: do_sgn_syst="false"
do_bkg_altfits="true"
if args.skip_bkg_altfits: do_bkg_altfits="false"
if args.data_file=="": args.data_file=path+data_files[args.year]

rt.gROOT.SetBatch(True)


#----------------------------------------------
# Get the signal parameterization vs. mass
#----------------------------------------------

# Check if testing a single MC mass point with the MC fit
tested_mass_point = ""
if args.use_mc_signal:
    tested_mass_point = "%.0f" % (args.scan_min)
    sgn_masspoints = [tested_mass_point]
if args.skip_mc_point != "":
    if args.skip_mc_point not in sgn_masspoints:
        print "Mass point %s not found in the signal MC list!" % (args.skip_mc_point)
        exit()
    sgn_masspoints = [mpoint for mpoint in sgn_masspoints if mpoint != args.skip_mc_point]


# Get the signal mass distribution for each Z prime MC sample
signal_distributions = []
cuts = sf+"*(Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<=xgb && xgb<"+args.xgb_max+")"
for mpoint in sgn_masspoints:
    h = rt.TH1F("hmass_"+mpoint,"Signal mass distribution",1100,0,1100)
    signal_distributions.append(signal_distribution(signal_samples[mpoint], h, "mass_ll", cuts, args.year, not args.skip_correct))

# Create a signal interpolation model
masses = array('d')
for mass in sgn_masspoints: masses.append(float(mass))
signal_model = create_signal_interpolation(masses, signal_distributions, args.use_gaus, figdir)
exit()

#----------------------------------------------
# Perform the signal scan at all mass points
#----------------------------------------------

###### main code ######
if MaxMasspoints>0:
    print "WILL PRODUCE",MaxMasspoints,"points"
NextPoint=True 

# set intitial scan mass
sr_center=args.scan_min

# List of jobs to process for multi-threaded processing
jobs = []

cnt=-1
script_head = 'echo ' if args.dry_run else ''

while (NextPoint):
    cnt+=1

    # calculate blind range and expected widths/yields for a mass
    sr_buffer = 1 #in signal width units
    region_buffer = 10 #in signal width units
    sr_approx_width = sr_center*(4./200.) # approximate as a linear function so it's stable between BDT categories where the width may vary

    # Get the signal model parameters
    sig_mass = sr_center
    sig_params = interpolate(signal_model, sig_mass)
    sig_yield  = sig_params[0]
    sig_mean   = sig_params[1]
    sig_width  = sig_params[2]
    sig_alpha1 = sig_params[3] if not args.use_gaus else 0.
    sig_alpha2 = sig_params[4] if not args.use_gaus else 0.
    sig_enne1  = sig_params[5] if not args.use_gaus else 0.
    sig_enne2  = sig_params[6] if not args.use_gaus else 0.

    # Define the search region using the approximate signal width
    sr_width = round(sig_width,2)
    sr_min = round(sr_center - sr_buffer*sr_approx_width,2)
    sr_max = round(sr_center + sr_buffer*sr_approx_width,2)
    sr_yld = round(sig_yield,2)
    cut_off_min = args.min_fit_mass # Don't use below 95 GeV due to Z->tautau contamination
    cut_off_max = args.max_fit_mass
    min_mass = round(sr_center - region_buffer*sr_approx_width,2) if not args.full_mass else cut_off_min
    max_mass = round(sr_center + region_buffer*sr_approx_width,2) if not args.full_mass else cut_off_max
    min_mass = max(cut_off_min, min_mass)
    max_mass = min(cut_off_max, max_mass)
    print "SR central",sr_center,"width",sr_width,"min",sr_min,"max",sr_max,"yield",sr_yld

    # Check if this mass point should be processed
    if(args.mass_point < 0 or cnt == args.mass_point):
        # create pdfs for the mass point
        if not args.skip_fit:
            sgn_argument = ''
            bkg_argument = ''
            if args.component == "sgn" or args.component == "all":
                tail = (' >| log/fit_sgn_%s_mp%i.log' % (args.name, cnt)) if args.log_files else ''

                if args.use_gaus: #Define the signal shape parameters (Gaussian or Crystal Ball)
                    sig_line = '{'+str(sr_center)+','+str(sr_width)+'}'
                else:
                    sig_line = '{%.4f, %.4f, %.4f, %.4f, %.4f, %.4f}' % (sig_mean, sig_width, sig_alpha1, sig_alpha2, sig_enne1, sig_enne2)
                sgn_argument = script_head + 'root -l -b -q ScanMuE_fit_sgn_v'+args.ver+'.C\'("' \
                    +args.name+"_mp"+str(cnt)+'",' \
                    +str(min_mass)+','+str(max_mass)+','+sig_line+',' \
                    +str(sr_yld)+','+shape_dc+',"'+args.outvar+'",'+do_sgn_syst+',"'+args.param_name+'")\'' + tail

            if args.component == "bkg" or args.component == "all":
                tail = (' >| log/fit_bkg_%s_mp%i.log' % (args.name, cnt)) if args.log_files else ''
                bkg_argument = script_head + 'root -l -b -q ScanMuE_fit_bkg_v'+args.ver+'.C\'("'+args.name+"_mp"+str(cnt)+'","'+args.data_file \
                    +'","'+args.xgb_min+'","'+args.xgb_max+'",'+str(min_mass)+','+str(max_mass)+','+str(sr_min)+','+str(sr_max) \
                    +','+unblind+','+shape_dc+',"'+args.outvar+'","'+args.param_name+'")\'' + tail
            job = Process(target = proc_unit, args=(sgn_argument,bkg_argument))
            jobs.append(job)


        # Create a corresponding datacard
        cardname = carddir + "datacard_zprime_" + args.name + ("_mass-%.1f" % (sr_center)) + "_mp" + str(cnt) + ".txt"
        sig_file = "WorkspaceScanSGN/workspace_scansgn_v" + args.ver + "_" + args.name + "_mp" + str(cnt) + ".root"
        bkg_file = "WorkspaceScanBKG/workspace_scanbkg_v" + args.ver + "_" + args.name + "_mp" + str(cnt) + ".root"
        if not args.dry_run and not args.skip_dc: print_datacard(cardname, sig_file, bkg_file, args.param_name, sr_center)


    # next iteration mass and exit conditions
    sr_center = round(sr_center +args.scan_step*sr_approx_width,2)
    if cnt>= MaxMasspoints and MaxMasspoints>0:
        print "Requested only",MaxMasspoints,"run"
        NextPoint=False
    if sr_center >= args.scan_max:
        NextPoint=False
  
  
if args.nthreads < 1: args.nthreads = 1
print("Parallel processing using %i threads" % (args.nthreads))
for ithread in range(0,len(jobs),args.nthreads):
    nthread = args.nthreads+ithread
    if nthread > len(jobs):
        nthread=len(jobs)
    print("Processing threads %i to %i out of %i jobs: %5.1f%% processed" % (ithread, nthread-1, len(jobs), (ithread*100./len(jobs))))
    for job in jobs[ithread:nthread]:
        job.start()
    for job in jobs[ithread:nthread]:    
        job.join()

print "Scannned",cnt+1,"points"


