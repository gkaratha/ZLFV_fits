# Apply the Z->e+mu selection and output a sparse TTree
import os
import time
import argparse
import ROOT as rt
from array import array
from signal_model import *
from math import sqrt,cos,sin,acos
rt.gROOT.SetBatch(True)

#----------------------------------------------
# Helper function to add TTree branches
def add_branch(tree, name, br_type = 'F'):
   if br_type == 'F':
      val = array('f', [0])
   tree.Branch(name, val, name + '/' + br_type)
   return val
      
#----------------------------------------------
# Delta phi calculation
def delta_phi(phi_1, phi_2):
   dphi = phi_1 - phi_2
   pi = 3.14159265358979323846
   if dphi <= -pi: dphi += 2.*pi
   if dphi >   pi: dphi -= 2.*pi
   return dphi
      
#----------------------------------------------
# Apply a correction to the MET
def met_correction(pt_change, phi, met, met_phi):
   # Apply the correction in the opposite direction to the MET
   x = met*cos(met_phi) - pt_change*cos(phi)
   y = met*sin(met_phi) - pt_change*sin(phi)
   met_corr = sqrt(x*x + y*y)
   met_corr_phi = acos(max(-1., min(1., x/met_corr)))
   if y < 0.: met_corr_phi *= -1
   return [met_corr, met_corr_phi]
      

#----------------------------------------------
# Transverse mass approximation
def mt(pt_1, phi_1, pt_2, phi_2):
   return sqrt(2.*pt_1*pt_2*(1.-cos(delta_phi(phi_1, phi_2))));
      
#----------------------------------------------
# Trigger matching
def trigger_matching(pt, eta, phi, flavor, year, tree):
   # Configure the trigger parameters
   trig_bit = 1 if flavor == 13 or year != 2017 else 10
   if flavor == 13:
      trig_pt = 27. if year == 2017 else 24.
   else:
      trig_pt = 27. if year == 2016 else 32.

   # Loop through the trigger objects
   for index in range(tree.nTrigObj):
      if abs(tree.TrigObj_id[index]) != flavor: continue # Correct physics object
      if tree.TrigObj_filterBits[index] & (1 << trig_bit) == 0: continue # Pass trigger ID
      if tree.TrigObj_pt[index] <= trig_pt: continue # Trigger pT threshold
      dphi = delta_phi(phi, tree.TrigObj_phi[index])
      deta = eta - tree.TrigObj_eta[index]
      if sqrt(dphi*dphi + deta*deta) > 0.1: continue # Fails Delta R matching
      # Accept the trigger match
      return True
   return False # No trigger matches found


#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",required=True,type=str,help="output root name")
parser.add_argument("--input-file", dest="input_file",required=True,type=str,help="Input MC/data file")
parser.add_argument("--input-directory", dest="input_directory",default="lfvanalysis_rootfiles",type=str,help="Input MC/data file directory")
parser.add_argument("--tree-name", dest="tree_name",default="Events",type=str,help="Input TTree name")
parser.add_argument("--year", dest="year",required=True,type=int,help="Data/MC processing year (2016, 2017, or 2018)")
parser.add_argument("--out-tree-name", dest="out_tree_name",default="mytree",type=str,help="Output TTree name")
parser.add_argument("--min-mass", dest="min_mass",default=95.,type=float,help="Minimum mass allowed")
parser.add_argument("--first-entry", dest="first_entry",default=0,type=int,help="First entry to process")
parser.add_argument("--max-entries", dest="max_entries",default=-1,type=int,help="Maximum entries to process")

args, unknown = parser.parse_known_args()

#----------------------------------------------
# Validate the input
#----------------------------------------------

# check input flags
if len(unknown)>0: 
   print "not found:",unknown,"exitting"
   exit()

if args.year not in [2016, 2017, 2018]:
   print "Unknown year %i" % (args.year)
   exit()

#----------------------------------------------
# Read in the input data
#----------------------------------------------

# default path
# path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
path="/eos/cms/store/group/phys_smp/ZLFV/%s/" % (args.input_directory)
figdir = "./figures/skims/%s_%i/" % (args.name, args.year)
outdir = "./trees/"
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir, figdir))
os.system("[ ! -d %s ] && mkdir -p %s" % (outdir, outdir))
os.system("[ ! -d log ] && mkdir log")

f_in = rt.TFile.Open(path + args.input_file, 'READ')
t_in = f_in.Get(args.tree_name)
nentries = t_in.GetEntries()

print "Input ntuple has", nentries, "entries"

first_entry = min(args.first_entry, nentries)
max_entries = nentries-first_entry if args.max_entries < 0 else min(nentries-first_entry, args.max_entries)
max_entry = first_entry + max_entries

t_in.SetBranchStatus('GenPart*', 0)
t_in.SetBranchStatus('Tau*', 0)
t_in.SetBranchStatus('Photon*', 0)
t_in.SetBranchStatus('GenMET*', 0)
t_in.SetBranchStatus('RawMET*', 0)
t_in.SetBranchStatus('RawPuppiMET*', 0)
t_in.SetBranchStatus('PV*', 0)
t_in.SetBranchStatus('Pileup_*', 0)
t_in.SetBranchStatus('LHE*', 0)


#----------------------------------------------
# Narrow to relevant events with TEventList
#----------------------------------------------

cuts = 'SelectionFilter_LepM > ' + str(args.min_mass-5.) + ' && (Muon_charge[0] * Electron_charge[0] < 0) && Muon_corrected_pt[0] > 20. && Electron_pt[0] > 20.'
print 'Apply event list selection', cuts
t_in.Draw(">>elist", cuts);
elist = rt.gDirectory.Get("elist");

#----------------------------------------------
# Create the output data structure
#----------------------------------------------

f_out = rt.TFile.Open(outdir + 'skim_' + args.name + '_' + str(args.year) + '.root', 'RECREATE')
t_out = rt.TTree(args.out_tree_name, 'Z prime events tree')


pt_ll          = add_branch(t_out, 'pt_ll'           )
eta_ll         = add_branch(t_out, 'eta_ll'          )
phi_ll         = add_branch(t_out, 'phi_ll'          )
mass_ll        = add_branch(t_out, 'mass_ll'         )
pt_l1          = add_branch(t_out, 'pt_l1'           )
eta_l1         = add_branch(t_out, 'eta_l1'          )
phi_l1         = add_branch(t_out, 'phi_l1'          )
mass_l1        = add_branch(t_out, 'mass_l1'         )
isMuon_l1      = add_branch(t_out, 'isMuon_l1'       )
pt_l2          = add_branch(t_out, 'pt_l2'           )
eta_l2         = add_branch(t_out, 'eta_l2'          )
phi_l2         = add_branch(t_out, 'phi_l2'          )
mass_l2        = add_branch(t_out, 'mass_l2'         )
isMuon_l2      = add_branch(t_out, 'isMuon_l2'       )
pt_ratio       = add_branch(t_out, 'ratio_ptl2_ptl1' )

met            = add_branch(t_out, 'met'             )
met_phi        = add_branch(t_out, 'met_phi'         )
mt_l1          = add_branch(t_out, 'mt_l1'           )
mt_l2          = add_branch(t_out, 'mt_l2'           )
met_sig        = add_branch(t_out, 'met_significance')
pt_j1          = add_branch(t_out, 'pt_j1'           )
ht             = add_branch(t_out, 'ht'              )
st             = add_branch(t_out, 'st'              )
njets          = add_branch(t_out, 'njets'           )
ratio_met_ptll = add_branch(t_out, 'ratio_met_ptll'  )
ratio_met_ht   = add_branch(t_out, 'ratio_met_ht'    )
dphi_met_ll    = add_branch(t_out, 'dphi_met_ll'     )

Flag_met       = add_branch(t_out, 'Flag_met'        )
Flag_muon      = add_branch(t_out, 'Flag_muon'       )
Flag_electron  = add_branch(t_out, 'Flag_electron'   )

evtnum         = add_branch(t_out, 'evtnum'          )

#----------------------------------------------
# Useful evaluated objects
#----------------------------------------------

lv1 = rt.TLorentzVector()
lv2 = rt.TLorentzVector()
ll  = rt.TLorentzVector()

accepted = 0
tot_wt = 0.

step = 0
prev_time = time.time()
start_time = prev_time
for entry in range(first_entry, max_entry):
   if step % 10000 == 0:
      curr_time = time.time()
      print "Processing event %7i (entry %7i): %5.1f%% complete, %5.1f%% accepted, %6.0f Hz" % (step, entry,
                                                                                                step*100./(max_entry-first_entry),
                                                                                                accepted*100./(step if step > 0 else 1.),
                                                                                                (10000./(curr_time-prev_time)) if step > 0 else 0.
                                                                                                )
      prev_time = curr_time
   step += 1

   if not elist.Contains(entry): continue

   t_in.GetEntry(entry)

   #-----------------------------------------
   # Evaluate input information
   #-----------------------------------------

   # Set the lepton kinematics
   ele_lead = t_in.Electron_pt[0] > t_in.Muon_corrected_pt[0]
   if ele_lead:
      pt_l1[0] = t_in.Electron_pt      [0]; eta_l1[0] = t_in.Electron_eta[0]; phi_l1[0] = t_in.Electron_phi[0]; mass_l1[0] = 0.511e-3;
      pt_l2[0] = t_in.Muon_corrected_pt[0]; eta_l2[0] = t_in.Muon_eta    [0]; phi_l2[0] = t_in.Muon_phi    [0]; mass_l2[0] = 0.10566;
   else:
      pt_l1[0] = t_in.Muon_corrected_pt[0]; eta_l1[0] = t_in.Muon_eta    [0]; phi_l1[0] = t_in.Muon_phi    [0]; mass_l1[0] = 0.10566;
      pt_l2[0] = t_in.Electron_pt      [0]; eta_l2[0] = t_in.Electron_eta[0]; phi_l2[0] = t_in.Electron_phi[0]; mass_l2[0] = 0.511e-3;
   lv1.SetPtEtaPhiM(pt_l1[0], eta_l1[0], phi_l1[0], mass_l1[0])
   lv2.SetPtEtaPhiM(pt_l2[0], eta_l2[0], phi_l2[0], mass_l2[0])
   isMuon_l1[0] = 0. if ele_lead else 1.
   isMuon_l2[0] = 1. if ele_lead else 0.
   ll = lv1+lv2
   mass_ll[0] = ll.M  ()
   pt_ll  [0] = ll.Pt ()
   eta_ll [0] = ll.Eta()
   eta_sc = t_in.Electron_deltaEtaSC[0] + t_in.Electron_eta[0]
   pt_ratio[0] = pt_l2[0] / pt_l1[0]

   # Set the MET kinematics
   met           [0] = t_in.PuppiMET_pt
   met_phi       [0] = t_in.PuppiMET_phi
   # Propagate the muon pT correction to the MET
   [met_corr, met_corr_phi] = met_correction(t_in.Muon_corrected_pt[0] - t_in.Muon_pt[0], t_in.Muon_phi[0], met[0], met_phi[0])
   met           [0] = met_corr
   met_phi       [0] = met_corr_phi
   met_sig       [0] = t_in.MET_significance
   mt_l1         [0] = mt(pt_l1[0], phi_l1[0], met[0], met_phi[0])
   mt_l2         [0] = mt(pt_l2[0], phi_l2[0], met[0], met_phi[0])
   dphi_met_ll   [0] = delta_phi(phi_ll[0], met_phi[0])
   ratio_met_ptll[0] = met[0] / pt_ll[0] if pt_ll[0] > 0. else 0.

   # Event variables FIXME: Likely need to clean jet collection first
   pt_j1        [0] = t_in.Jet_pt[0] if t_in.nJet > 0 else 0.
   ht           [0] = t_in.HT
   njets        [0] = t_in.nJet
   ratio_met_ht [0] = met[0] / ht[0] if ht[0] > 0. else 0.
   st           [0] = ht[0] + pt_l1[0] + pt_l2[0]
   Flag_met     [0] = t_in.Flag_METFilters
   Flag_muon    [0] = t_in.Flag_muonBadTrackFilter and t_in.Flag_BadPFMuonFilter
   Flag_electron[0] = t_in.Flag_EcalDeadCellTriggerPrimitiveFilter and t_in.Flag_eeBadScFilter
   evtnum       [0] = t_in.event


   #-----------------------------------------
   # Selection cuts
   #-----------------------------------------

   # Lepton kinematic cuts
   if mass_ll[0] < args.min_mass: continue
   if pt_l1[0] < 20. or pt_l2[0] < 20.: continue
   if abs(eta_l1[0]) >= 2.4 or abs(eta_l2[0]) >= 2.4: continue
   if abs(lv1.DeltaR(lv2)) < 0.3: continue
   if abs(eta_sc) >= 1.444 and abs(eta_sc) <= 1.566: continue

   # Lepton IDs isolation cuts
   if t_in.Electron_pfRelIso03_all[0] >= 0.15: continue
   if t_in.Muon_pfRelIso04_all[0] >= 0.15: continue
   if not t_in.Muon_mediumPromptId: continue
   if not t_in.Electron_mvaFall17V2noIso_WP90: continue

   # Lepton displacement cuts
   if t_in.Electron_dz[0] >= 0.5 or t_in.Electron_dxy[0] >= 0.2: continue
   if t_in.Electron_dzErr[0] <= 0. or t_in.Electron_dxyErr[0] <= 0.: continue
   if t_in.Muon_dz[0] >= 0.5 or t_in.Muon_dxy[0] >= 0.2: continue
   if t_in.Muon_dzErr[0] <= 0. or t_in.Muon_dxyErr[0] <= 0.: continue
   if t_in.Electron_dz[0]/t_in.Electron_dzErr[0] >= 4.7: continue
   if t_in.Muon_dz[0]/t_in.Muon_dzErr[0] >= 4.7: continue
   if t_in.Electron_dxy[0]/t_in.Electron_dxyErr[0] >= 3.0: continue
   if t_in.Muon_dxy[0]/t_in.Muon_dxyErr[0] >= 3.0: continue

   # Lepton sign selection
   if t_in.Muon_charge[0]*t_in.Electron_charge[0] > 0: continue

   # Trigger selection
   if args.year == 2016:
      elec_trig = t_in.HLT_Ele27_WPTight_Gsf and t_in.Electron_pt[0] > 29.
      muon_trig = t_in.HLT_IsoMu24 and t_in.Muon_pt[0] > 25.
   elif args.year == 2017:
      elec_trig = t_in.HLT_Ele32_WPTight_Gsf_L1DoubleEG and t_in.Electron_pt[0] > 35.
      muon_trig = t_in.HLT_IsoMu27 and t_in.Muon_pt[0] > 28.
   else: #2018
      elec_trig = t_in.HLT_Ele32_WPTight_Gsf and t_in.Electron_pt[0] > 34.
      muon_trig = t_in.HLT_IsoMu24 and t_in.Muon_pt[0] > 25.
   # Trigger matching
   muon_trig &= trigger_matching(t_in.Muon_corrected_pt[0], t_in.Muon_eta[0]    , t_in.Muon_phi[0]    , 13, args.year, t_in)
   elec_trig &= trigger_matching(t_in.Electron_pt[0]      , t_in.Electron_eta[0], t_in.Electron_phi[0], 11, args.year, t_in)
   if not elec_trig and not muon_trig: continue

   # Event flags:
   if not t_in.Flag_METFilters                         : continue
   if not t_in.Flag_BadChargedCandidateFilter          : continue
   if not t_in.Flag_goodVertices                       : continue
   if not t_in.Flag_HBHENoiseFilter                    : continue
   if not t_in.Flag_HBHENoiseIsoFilter                 : continue
   if not t_in.Flag_eeBadScFilter                      : continue
   if not t_in.Flag_muonBadTrackFilter                 : continue
   if not t_in.Flag_EcalDeadCellTriggerPrimitiveFilter : continue
   if not t_in.Flag_globalTightHalo2016Filter          : continue
   if not t_in.Flag_BadPFMuonFilter                    : continue


   # Event cuts
   if t_in.nBJet > 0: # Check for any b-tagged jet with the b-tagging region
      btagged = False
      for bjet in range(t_in.nBJet):
         btagged |= abs(t_in.BJet_eta[bjet]) < 2.4
      if btagged: continue
         

   #-----------------------------------------
   # Fill the outgoing data
   #-----------------------------------------

   accepted += 1
   t_out.Fill()

#-----------------------------------------
# Post-processing steps
#-----------------------------------------

total_time = time.time() - start_time
print "%i events accepted (efficiency = %.3f), processing took %.1f %s" % (accepted, accepted*1./max_entries,
                                                                           total_time/60. if total_time > 120. else total_time,
                                                                           "min" if total_time > 120. else "s")

t_out.Write()
f_out.Close()
