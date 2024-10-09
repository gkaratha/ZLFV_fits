# Evaluate the signal rate effects by year
import os
import argparse
import ROOT as rt
from array import array
from signal_model import *
from math import sqrt
import ctypes

rt.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("--gaus", dest="gaus",default=False, action='store_true',help="Model the signal with a Gaussian instead of a Crystal Ball")
parser.add_argument("--xgb-min", dest="xgb_min",default=1. ,type=float ,help="BDT score minimum cut")
parser.add_argument("--xgb-max", dest="xgb_max",default=-1.,type=float ,help="BDT score maximum cut")
args, unknown = parser.parse_known_args()

gaus    = args.gaus
xgb_min = args.xgb_min
xgb_max = args.xgb_max

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/years_effect"
if gaus: figdir += "_gaus"
if xgb_min < xgb_max:
    figdir += "_%.0f_%.0f" % (xgb_min*100., xgb_max*100.)
figdir += "/"

print "Using figure directory:", figdir
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))

### MC signal mass points
sgn_masspoints=["100","500"]
# sgn_masspoints=["100","200","300","400","500","600","800","1000"]

# Define the signal samples by mass and period
signal_samples = {
   "100" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM100_mcRun16.root", 94800, 2016, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM100_mcRun17.root", 97800, 2017, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM100_mcRun18.root", 99200, 2018, path),
   },
   "125" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM125_mcRun18.root", 98400, 2018, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM125_mcRun18.root", 98400, 2018, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM125_mcRun18.root", 98400, 2018, path),
   },
   "150" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM150_mcRun18.root", 90500, 2018, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM150_mcRun18.root", 90500, 2018, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM150_mcRun18.root", 90500, 2018, path),
   },
   "175" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM175_mcRun18.root", 89900, 2018, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM175_mcRun18.root", 89900, 2018, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM175_mcRun18.root", 89900, 2018, path),
   },
   "200" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root", 96300, 2016, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root", 96300, 2017, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM200_mcRun18.root", 96300, 2018, path),
   },
   "300" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM300_mcRun18.root", 99700, 2016, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM300_mcRun18.root", 99700, 2017, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM300_mcRun18.root", 99700, 2018, path),
   },
   "400" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root", 97600, 2016, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root", 97600, 2017, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM400_mcRun18.root", 97600, 2018, path),
   },
   "500" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM500_mcRun16.root", 81300, 2016, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM500_mcRun17.root", 97800, 2017, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM500_mcRun18.root", 98700, 2018, path),
   },
   "600" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root", 97600, 2016, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root", 97600, 2017, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM600_mcRun18.root", 97600, 2018, path),
   },
   "800" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root", 97600, 2016, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root", 97600, 2017, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM800_mcRun18.root", 97600, 2018, path),
   },
   "1000" : {
      "2016" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root", 97600, 2016, path),
      "2017" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root", 97600, 2017, path),
      "2018" : sample("Meas_fullAndSF_bdt_v7_emu_scan_sgnM1000_mcRun18.root", 97600, 2018, path),
   },
}

# Loop through the signal files, retrieving the signal rates for each running period
yields = {"2016": array('d'), "2017": array('d'), "2018": array('d')}
yerrs  = {"2016": array('d'), "2017": array('d'), "2018": array('d')}
# Scale factors
sf="Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt"
sf += "*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt"
sf += "*PU_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*Prefire_wt"
# sf += "*PtZ_wt*MixZ_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)" # LFV Z-specific weights

cuts = sf + "*(Flag_met && Flag_muon && Flag_electron)"
if xgb_min < xgb_max: cuts += " && (xgb > %.2f && xgb <= %.2f)" % (xgb_min, xgb_max)
lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":137.6}
for mpoint in sgn_masspoints:
    for year in ["2016", "2017", "2018"]:
        # Retrieve the signal distribution
        hname = "hmass_"+mpoint
        h = rt.TH1F(hname,"Mass = %s GeV" % (mpoint),1500,0,1500)
        h = signal_distribution(signal_samples[mpoint], h, "mass_ll", cuts, year, False)
        # Remove the luminosity effect
        h.Scale(1./lumis[year])
        eff_err = ctypes.c_double(1.)
        yields[year].append(h.IntegralAndError(1, h.GetNbinsX(), eff_err))
        yerrs[year].append(eff_err.value)

# Add the x-axis information
masses = array('d')
xerrs  = array('d')
for mass in sgn_masspoints: masses.append(float(mass))
for index in range(len(masses)):
    mass_1 = masses[index]
    mass_2 = masses[index-1] if index > 0 else masses[index+1]
    xerrs.append(abs(mass_1-mass_2)/2.)
    

#--------------------------------------------------------------------
# Plot the efficiencies for each year
#--------------------------------------------------------------------

for year in ["2016", "2017", "2018"]:
    c = rt.TCanvas('c_yield', 'c_yield', 700, 500)
    g = rt.TGraphErrors(len(masses), masses, yields[year], xerrs, yerrs[year])
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.8)
    g.SetLineWidth(2)
    g.SetTitle("Signal efficiency vs. mass;Z prime mass (GeV/c^{2});efficiency")
    g.Draw("APE")
    c.SaveAs(figdir+'eff_' + year + '.png')

#--------------------------------------------------------------------
# Plot the efficiencies for each mass
#--------------------------------------------------------------------

for index in range(len(sgn_masspoints)):
    c = rt.TCanvas('c_yield', 'c_yield', 700, 500)
    years = array('d', [2016, 2017, 2018])
    years_errs = array('d', [0.5, 0.5, 0.5])
    effs_by_year = array('d', [yields[year][index] for year in ["2016","2017","2018"]])
    errs_by_year = array('d', [yerrs [year][index] for year in ["2016","2017","2018"]])
    g = rt.TGraphErrors(len(years), years, effs_by_year , years_errs, errs_by_year)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.8)
    g.SetLineWidth(2)
    g.SetTitle("Signal efficiency vs. year;;efficiency")
    g.Draw("APE")
    g.GetYaxis().SetRangeUser(0., 0.3)
    xaxis = g.GetXaxis()
    for ibin in [2016, 2017, 2018]:
        b = xaxis.FindBin(ibin)
        xaxis.SetBinLabel(b, str(ibin))
    xaxis.LabelsOption("h")
    xaxis.SetLabelSize(0.07)
    c.Modified()
    c.Update()
    c.SaveAs(figdir+'eff_M' + sgn_masspoints[index] + '.png')

#--------------------------------------------------------------------
# Plot the efficiencies relative to 2018
#--------------------------------------------------------------------

corrs = []
for year in ["2016", "2017"]:
    y    = array('d')
    yerr = array('d')
    for index in range(len(masses)):
        y.append(yields[year][index]/yields["2018"][index])
        # (x+s_x)/(y+s_y) ~ x/y*(1+s_x/x)*(1-s_y/y)
        yerr.append(y[-1]*sqrt((yerrs[year][index]/yields[year][index])**2 + (yerrs["2018"][index]/yields["2018"][index])**2))
    c = rt.TCanvas('c_yield', 'c_yield', 700, 500)
    g = rt.TGraphErrors(len(masses), masses, y, xerrs, yerr)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.8)
    g.SetLineWidth(2)
    g.SetTitle("Signal efficiency relative to 2018 vs. mass;Z prime mass (GeV/c^{2});eff / 2018 eff")
    g.Draw("APE")
    g.GetYaxis().SetRangeUser(0.8, 1.1)
    c.SaveAs(figdir+'rel_eff_' + year + '.png')
    corrs.append(y)
    print year, y

for index in range(len(corrs[0])):
    print "Average correction %i: %.3f" % (index, (corrs[0][index]+corrs[1][index])/2.)
