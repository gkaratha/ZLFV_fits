# Process COMBINE datacards for the Z'->emu scan
import os
import argparse
import ROOT as rt
from array import array
from math import exp
from math import log10

#----------------------------------------------------------------------------------------
# Define the sorting of the datacards
def file_sort(f):
   return int(f.split('_mp')[1].replace('.txt',''))

#----------------------------------------------------------------------------------------
# Process a single mass point
def process_datacard(card, directory, name, asimov = False, skip_fit = False, verbose = 0):
   if not os.path.isfile(directory + card):
      print "Card %s not found" % (card)
   if verbose > -1: print 'Processing mass point', name, '(card =', card+')'

   #----------------------------------------------------------------------------
   # Perform a signal rate fit
   #----------------------------------------------------------------------------

   if not skip_fit:
      command = 'combine -d %s -n .%s%s -M FitDiagnostics' % (card, name, "asimov" if asimov else "")
      if asimov: command += ' -t -1'
      #Allow negative measured signal rates
      command += ' --rMin -50 --rMax 50'
      #Additional commands to help fit converge properly
      command += ' --cminDefaultMinimizerStrategy 0'
      command += ' --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants'
      command += ' --cminRunAllDiscreteCombinations'
      output = 'fit_rate_%s.log' % (name)
      if verbose > 1: print command
      os.system('cd %s; %s >| %s; cd ..' % (directory, command, output))

   #----------------------------------------------------------------------------
   # Evaluate the signal rate significance
   #----------------------------------------------------------------------------

   if not skip_fit:
      command = 'combine -d %s -n .%s%s -M Significance --uncapped 1' % (card, name, "asimov" if asimov else "")
      if asimov: command += ' -t -1'
      #Allow negative measured signal rates
      command += ' --rMin -50 --rMax 50'
      #Additional commands to help fit converge properly
      command += ' --cminDefaultMinimizerStrategy 0'
      command += ' --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants'
      command += ' --cminRunAllDiscreteCombinations'
      output = 'fit_sig_%s.log' % (name)
      if verbose > 1: print command
      os.system('cd %s; %s >| %s; cd ..' % (directory, command, output))

   #----------------------------------------------------------------------------
   # Perform an upper limit evaluation
   #----------------------------------------------------------------------------

   if not skip_fit:
      command = 'combine -d %s -n .%s%s -M AsymptoticLimits' % (card, name, "asimov" if asimov else "")
      if asimov: command += ' -t -1'
      #Allow negative measured signal rates
      command += ' --rMin -200 --rMax 200'
      #Additional commands to help fit converge properly
      command += ' --cminDefaultMinimizerStrategy 0'
      command += ' --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --X-rtd MINIMIZER_multiMin_hideConstants'
      command += ' --cminRunAllDiscreteCombinations'
      output = 'fit_limit_%s.log' % (name)
      if verbose > 1: print command
      os.system('cd %s; %s >| %s; cd ..' % (directory, command, output))

   #----------------------------------------------------------------------------
   # Extract the mass FIXME: Make this more robust
   #----------------------------------------------------------------------------

   ws_name = name
   ws_name = name.replace('bdt_', 'bdt_0d7_1d0_')
   ws_file = 'WorkspaceScanSGN/workspace_scansgn_v2_%s.root' % (ws_name)
   if 'toy' in ws_file:
      ws_file = ws_file.split('_toy')[0] + '_mp' + ws_file.split('_mp')[1]
   f = rt.TFile.Open(ws_file, 'READ')
   ws = f.Get('workspace_signal')
   mass_var = ws.var('mean_bin2')
   mass = mass_var.getVal()
   f.Close()
   

   #----------------------------------------------------------------------------
   # Extract the results
   #----------------------------------------------------------------------------

   fit_file = '%sfitDiagnostics.%s%s.root' % (directory, name, "asimov" if asimov else "")
   f = rt.TFile.Open(fit_file, 'READ')
   t = f.Get('tree_fit_sb')
   t.GetEntry(0)
   r_fit = [t.r, t.rLoErr, t.rHiErr]
   if verbose > 0: print 'r fit results:', r_fit
   f.Close()

   sig_file = '%shiggsCombine.%s%s.Significance.mH120.root' % (directory, name, "asimov" if asimov else "")
   f = rt.TFile.Open(sig_file, 'READ')
   t = f.Get('limit')
   t.GetEntry(0)
   r_sig = t.limit
   if verbose > 0: print 'r significance results:', r_sig
   f.Close()
   
   lim_file = '%shiggsCombine.%s%s.AsymptoticLimits.mH120.root' % (directory, name, "asimov" if asimov else "")
   f = rt.TFile.Open(lim_file, 'READ')
   t = f.Get('limit')
   r_lim = []
   for index in range(t.GetEntries()):
      t.GetEntry(index)
      r_lim.append(t.limit)
   if verbose > 0: print 'r limit results:', r_lim
   f.Close()
   
   #----------------------------------------------------------------------------
   # Return the results
   #----------------------------------------------------------------------------

   return [r_fit, r_lim, r_sig, mass]


#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="bdt_v01", type=str,help="datacard directory name")
parser.add_argument("--skip-fits", dest="skip_fits",default=False, action='store_true',help="Skip fits, assume already processed")
parser.add_argument("--asimov", dest="asimov",default=False, action='store_true',help="Perform fits Asimov dataset")
parser.add_argument("--unblind", dest="unblind",default=False, action='store_true',help="Plot the observed limits")
parser.add_argument("--max-steps", dest="max_steps",default=-1, type=int, help="Maximum steps to take in the scan")
parser.add_argument("--first-step", dest="first_step",default=0, type=int, help="First mass step to process")
parser.add_argument("--tag", dest="tag",default="", type=str, help="Output directory tag")
parser.add_argument("-v", dest="verbose",default=0, type=int,help="Add verbose printout")

args, unknown = parser.parse_known_args()


##### configuaration ########

### check input flags
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()

#----------------------------------------------
# Setup the scan
#----------------------------------------------

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/scan_%s%s%s/" % (args.name, "_asimov" if args.asimov else "", args.tag)
carddir = "./datacards/%s/" % (args.name)
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))
os.system("[ ! -d %s ] && mkdir -p %s" % (carddir, carddir))
   
rt.gROOT.SetBatch(True)

#----------------------------------------------
# Perform the scan
#----------------------------------------------

list_of_files = [f for f in os.listdir(carddir) if '.txt' in f and '0d7' not in f]
# Sort the list by mass point
list_of_files.sort(key=file_sort)
if args.first_step > 0:
   list_of_files = list_of_files[args.first_step:]
if args.max_steps > 0:
   list_of_files = list_of_files[:args.max_steps]

asimov = args.asimov

# List of results
masses = array('d')
masses_errs = array('d')
r_fits = array('d')
r_lo_errs = array('d')
r_hi_errs = array('d')
r_lims = array('d')
r_exps = array('d')
r_exps_lo = array('d') #1 sigma band
r_exps_hi = array('d')
r_exps_lo_2 = array('d') #2 sigma band
r_exps_hi_2 = array('d')
r_sigs = array('d')


prev_mass = -1.
min_lim =  1.e10
max_lim = -1.e10
min_r =  1.e10
max_r = -1.e10
for f in list_of_files:
   if '.txt' not in f: continue
   # only process the merged fits
   if '0d7' in f: continue
   mass_point = f.split('_mp')[1].split('.txt')[0]
   [r_fit, r_lim, r_sig, mass] = process_datacard(f, carddir, args.name + '_mp'+mass_point, asimov, args.skip_fits, args.verbose)

   # store the results
   masses.append(mass)
   masses_errs.append((mass - prev_mass)/2. if prev_mass > 0. else 0.5) #rough mass bin estimate
   r_fits.append(r_fit[0])
   r_lo_errs.append(r_fit[1])
   r_hi_errs.append(r_fit[2])
   r_lims.append(r_lim[-1])
   r_exps.append(r_lim[2])  
   r_exps_lo.append(r_lim[2] - r_lim[1])
   r_exps_hi.append(r_lim[3] - r_lim[2])
   r_exps_lo_2.append(r_lim[2] - r_lim[0])
   r_exps_hi_2.append(r_lim[4]- r_lim[2])
   r_sigs.append(r_sig)

   prev_mass = mass
   min_lim = min(r_lim[0], r_lim[-1], min_lim)
   max_lim = max(r_lim[4], r_lim[-1], max_lim)
   min_r = min(r_fit[0] - r_fit[1], min_r)
   max_r = max(r_fit[0] + r_fit[2], max_r)

if len(masses) > 1: masses_errs[0] = (masses[1] - masses[0])/2.

#----------------------------------------------
# Plot the results
#----------------------------------------------

rt.gStyle.SetOptStat(0)

#----------------------------------------------
# Limit plot
#----------------------------------------------

g_exp   = rt.TGraph(len(masses), masses, r_exps)
g_exp_1 = rt.TGraphAsymmErrors(len(masses), masses, r_exps, masses_errs, masses_errs, r_exps_lo  , r_exps_hi  )
g_exp_2 = rt.TGraphAsymmErrors(len(masses), masses, r_exps, masses_errs, masses_errs, r_exps_lo_2, r_exps_hi_2)

c = rt.TCanvas('c_lim', 'c_lim', 800, 600)
g_exp_2.SetTitle("95% CL_{S} limit vs. Z' mass; Z' mass (GeV/c^{2}); #sigma(Z')*BR(Z'->e#mu) (fb)")
g_exp_2.SetFillColor(rt.kOrange)
g_exp_2.SetLineColor(rt.kOrange)
g_exp_1.SetFillColor(rt.kGreen+1)
g_exp_1.SetLineColor(rt.kGreen+1)
g_exp_2.Draw("AE3")
g_exp_1.Draw("E3")
g_exp.SetLineStyle(rt.kDashed)
g_exp.SetLineWidth(2)
g_exp.SetLineColor(rt.kBlack)
g_exp.Draw("XL")

if args.asimov or args.unblind:
   g_obs = rt.TGraph(len(masses), masses, r_lims)
   g_obs.SetMarkerStyle(20)
   g_obs.SetMarkerSize(0.8)
   g_obs.SetLineWidth(2)
   g_obs.SetMarkerColor(rt.kBlack)
   g_obs.SetLineColor(rt.kBlack)
   g_obs.Draw("PL")

g_exp_2.GetXaxis().SetRangeUser(masses[0], masses[-1])
g_exp_2.GetYaxis().SetRangeUser(0.8*min_lim, 1.4*max_lim)

leg = rt.TLegend(0.7, 0.7, 0.89, 0.89)
leg.SetLineWidth(0)
if args.asimov or args.unblind:
   leg.AddEntry(g_obs  , 'Observed', 'PL')
leg.AddEntry(g_exp  , 'Expected', 'L')
leg.AddEntry(g_exp_1, '#pm1#sigma' , 'F')
leg.AddEntry(g_exp_2, '#pm2#sigma' , 'F')
leg.Draw()


c.SaveAs(figdir+'limits.png')

#----------------------------------------------
# Signal rate plot
#----------------------------------------------

if not args.asimov and not args.unblind: exit()

g_r = rt.TGraphAsymmErrors(len(masses), masses, r_fits, masses_errs, masses_errs, r_lo_errs  , r_hi_errs)

c = rt.TCanvas('c_fit', 'c_fit', 800, 800)
pad1 = rt.TPad('pad1', 'pad1', 0., 0.3, 1.0, 1.0); pad1.Draw()
pad2 = rt.TPad('pad2', 'pad2', 0., 0.0, 1.0, 0.3); pad2.Draw()
pad1.SetBottomMargin(0.01)
pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.25)

# Add the fit distribution
pad1.cd()
g_r.SetTitle("Z' fit rate vs. Z' mass;; #sigma(Z')*BR(Z'->e#mu) (fb)")
g_r.SetMarkerStyle(20)
g_r.SetMarkerSize(0.8)
g_r.SetLineWidth(2)
g_r.SetMarkerColor(rt.kBlack)
g_r.SetLineColor(rt.kBlack)
g_r.Draw("APE1")

g_r.GetXaxis().SetRangeUser(masses[0], masses[-1])
g_r.GetYaxis().SetRangeUser(min_r - 0.05*(max_r - min_r), max_r + 0.1*(max_r - min_r))
g_r.GetXaxis().SetLabelSize(0.)

# Add an approximate significance distribution below the fit results
pad2.cd()
significances = array('d')
sig_half = array('d')
sig_err  = array('d')
m_errs = array('d') #create buffers between mass points
for index in range(len(r_fits)):
   r_err = (r_lo_errs[index] if r_fits[index] > 0. else r_hi_errs[index])
   if r_err <= 0.:
      print "Fit errors for mass point %.1f are 0" % (masses[index])
      r_err = max([0.1, r_lo_errs[index], r_hi_errs[index]])
   significances.append((r_fits[index]) / r_err)
   sig_half.append(significances[-1]/2.)
   sig_err.append(significances[-1]/2.)
   if index < len(masses) - 1:
      m_errs.append(0.8*(masses[index+1] - masses[index])/2.)
   elif len(masses) > 1:
      m_errs.append(0.8*(masses[index] - masses[index-1])/2.)
   else:
      m_errs.append(1.)

g_sig = rt.TGraphErrors(len(masses), masses, sig_half, m_errs, sig_err)
g_sig.SetTitle(";Z' mass (GeV/c^{2});#sigma_{r}")
g_sig.SetFillColor(rt.kAtlantic)
g_sig.Draw("AE2")
g_sig.GetXaxis().SetRangeUser(masses[0], masses[-1])
g_sig.GetYaxis().SetRangeUser(-4, 4)
g_sig.GetXaxis().SetLabelSize(0.08)
g_sig.GetYaxis().SetLabelSize(0.08)
g_sig.GetXaxis().SetTitleSize(0.1)
g_sig.GetXaxis().SetTitleOffset(0.9)
g_sig.GetYaxis().SetTitleSize(0.1)
g_sig.GetYaxis().SetTitleOffset(0.4)
line = rt.TLine(masses[0], 0., masses[-1], 0.)
line.SetLineWidth(2)
line.SetLineStyle(rt.kDashed)
line.SetLineColor(rt.kBlack)
line.Draw('same')

c.SaveAs(figdir+'fits.png')


#----------------------------------------------
# Significance plot
#----------------------------------------------

g_sig   = rt.TGraph(len(masses), masses, r_sigs)

c = rt.TCanvas('c_sig', 'c_sig', 800, 600)
g_sig.SetTitle("Measurement significance vs. Z' mass; Z' mass (GeV/c^{2}); #sigma(BR(Z'->e#mu))")
g_sig.SetMarkerStyle(20)
g_sig.SetMarkerSize(0.8)
g_sig.SetLineWidth(2)
g_sig.SetMarkerColor(rt.kBlack)
g_sig.SetLineColor(rt.kBlack)
g_sig.Draw("APL")

min_sig = min(r_sigs)
max_sig = max(r_sigs)
g_sig.GetXaxis().SetRangeUser(masses[0], masses[-1])
g_sig.GetYaxis().SetRangeUser(min_sig - 0.1*(max_sig-min_sig), max_sig + 0.1*(max_sig-min_sig))

c.SaveAs(figdir+'sig.png')


#----------------------------------------------
# p-value plot
#----------------------------------------------

pvals = array('d')
for sig in r_sigs: pvals.append(rt.RooStats.SignificanceToPValue(sig))

# Calculate the global significance by evaluating N(up crossings)
ref_u = 0.
n_ref_u = 0

for index in range(1,len(r_sigs)):
   if r_sigs[index-1] < ref_u and r_sigs[index] > ref_u: n_ref_u += 1

global_pvals = array('d')
for pval in pvals:
   sig = rt.RooStats.PValueToSignificance(pval)
   global_pval = min(1., pval + n_ref_u*exp(-(sig**2-ref_u**2)/2.))
   global_pvals.append(global_pval)
   print "p = %.4f, sig = %.2f, n_ref = %i, p_global = %.4f" % (pval, sig, n_ref_u, global_pvals[-1])

g_pval   = rt.TGraph(len(masses), masses, pvals)
g_global = rt.TGraph(len(masses), masses, global_pvals)

c = rt.TCanvas('c_pval', 'c_pval', 800, 600)
g_pval.SetTitle("Measurement p-value vs. Z' mass; Z' mass (GeV/c^{2}); p")
g_pval.SetMarkerStyle(20)
g_pval.SetMarkerSize(0.8)
g_pval.SetLineWidth(2)
g_pval.SetMarkerColor(rt.kRed)
g_pval.SetLineColor(rt.kRed)
g_pval.Draw("AL")
# g_pval.Draw("APL")

g_global.SetMarkerStyle(20)
g_global.SetMarkerSize(0.8)
g_global.SetLineWidth(2)
g_global.SetMarkerColor(rt.kRed)
g_global.SetLineColor(rt.kRed)
g_global.SetLineStyle(rt.kDashed)
g_global.Draw("L")

min_pval = min(pvals)
max_pval = max(pvals)
g_pval.GetXaxis().SetRangeUser(masses[0], masses[-1])
g_pval.GetYaxis().SetRangeUser(0.2*min_pval, 2.)
c.SetLogy()

# Add sigma lines to the p-value plot
sig_p_min = rt.RooStats.PValueToSignificance(0.2*min_pval)
lines = []
for sig in range(int(sig_p_min)):
   p_sig = rt.RooStats.SignificanceToPValue(sig)
   line = rt.TLine(masses[0], p_sig, masses[-1], p_sig)
   line.SetLineStyle(rt.kDashed)
   line.SetLineWidth(2)
   line.SetLineColor(rt.kBlack)
   line.Draw("same")
   lines.append(line)

c.SaveAs(figdir+'pval.png')

# Make just a histogram of log(p-value)
h = rt.TH1D('h_pval', '-log_{10}(p-value) distribution', 40, 0, 8)
for pval in pvals: h.Fill(-log10(abs(pval)))
h.SetLineColor(rt.kBlue)
h.SetFillColor(rt.kAtlantic)
h.SetLineWidth(2)
h.Draw("hist")
h.SetXTitle('-log_{10}(p)')
c.SaveAs(figdir+'pval_hist.png')
