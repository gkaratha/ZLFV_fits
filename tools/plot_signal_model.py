# Plot the Z prime signal model interpolation
import os
import argparse
import ROOT as rt
from array import array
from signal_model import *

rt.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("--gaus", dest="gaus",default=False, action='store_true',help="Model the signal with a Gaussian instead of a Crystal Ball")
parser.add_argument("--xgb-min", dest="xgb_min",default=1. ,type=float ,help="BDT score minimum cut")
parser.add_argument("--xgb-max", dest="xgb_max",default=-1.,type=float ,help="BDT score maximum cut")
parser.add_argument("--year", dest="year",default="Run2", type=str,help="Data period to use (2016, 2017, 2018, or Run2)")
args, unknown = parser.parse_known_args()

gaus    = args.gaus
xgb_min = args.xgb_min
xgb_max = args.xgb_max

### default path
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
figdir = "./figures/signal_model"
if gaus: figdir += "_gaus"
if xgb_min < xgb_max:
    figdir += "_%.0f_%.0f" % (xgb_min*100., xgb_max*100.)
figdir += "/"

print "Using figure directory:", figdir
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))

### MC signal mass points
# sgn_masspoints=["200","400","600","800","1000"]
sgn_masspoints=["100","125","150","175","200","300","400","500","600","800","1000"]

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

# Loop through the signal files
hmasses = []
# Scale factors
sf="Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt"
sf += "*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt"
sf += "*PU_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*Prefire_wt"
# sf += "*PtZ_wt*MixZ_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)" # LFV Z-specific weights

cuts = sf + "*(Flag_met && Flag_muon && Flag_electron)"
if xgb_min < xgb_max: cuts += " && (xgb > %.2f && xgb <= %.2f)" % (xgb_min, xgb_max)

effs = array('d')
lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":137.6}
for mpoint in sgn_masspoints:
    # Retrieve the signal distribution
    hname = "hmass_"+mpoint
    h = rt.TH1F(hname,"Mass = %s GeV" % (mpoint),1500,0,1500)
    h = signal_distribution(signal_samples[mpoint], h, "mass_ll", cuts, args.year)
    hmasses.append(h)
    effs.append(h.Integral() / lumis[args.year])

# Create the signal interpolation model
masses = array('d')
for mass in sgn_masspoints: masses.append(float(mass))
signal_model = create_signal_interpolation(masses, hmasses, args.gaus, figdir)

#--------------------------------------------------------------------
# Plot the efficiency in addition to the yield
#--------------------------------------------------------------------

c = rt.TCanvas('c_eff', 'c_eff', 700, 500)
g = rt.TGraph(len(masses), masses, effs)
g.SetMarkerStyle(20)
g.SetMarkerSize(0.8)
g.SetLineWidth(2)
g.SetTitle("Signal efficiency vs. mass;Z prime mass (GeV/c^{2});efficiency")
g.Draw("APE")
c.SaveAs(figdir+'eff.png')

#--------------------------------------------------------------------
# Plot the mass distributions along with fits to them
#--------------------------------------------------------------------

c = rt.TCanvas('c_mass', 'c_mass', 700, 500)
haxis = hmasses[0]
colors = [rt.kGreen, rt.kRed, rt.kOrange, rt.kViolet+6, rt.kGreen-7, rt.kMagenta, rt.kAtlantic]
leg = rt.TLegend(0.1, 0.75, 0.9, 0.9)
leg.SetTextSize(0.04)
leg.SetNColumns(3)
for index in range(len(hmasses)):
    h = hmasses[index]
    if index == 0: h.Draw("hist")
    else         : h.Draw("hist same")
    h.SetLineColor(colors[index % len(colors)])
    h.SetLineWidth(2)
    leg.AddEntry(h)
leg.Draw()

haxis.SetTitle("Z prime mass distributions")
haxis.SetXTitle("Z prime mass (GeV/c^{2})")
max_val = max([h.GetMaximum() for h in hmasses])
haxis.GetYaxis().SetRangeUser(0., 1.3*max_val)
haxis.GetXaxis().SetRangeUser(0., 1100.)
rt.gStyle.SetOptStat(0)

c.SaveAs(figdir+'masses.png')


#--------------------------------------------------------------------
# Plot the mass distribution interpolations
#--------------------------------------------------------------------

obs = rt.RooRealVar("obs", "M_{e#mu}", 500., 0., 1100., "GeV/c^{2}")
obs.setBins(550)

frame = obs.frame(rt.RooFit.Title("Signal model interpolation"))

max_yield = 0.
for mpoint in sgn_masspoints:
    mass = float(mpoint)
    sig_params = interpolate(signal_model, mass)
    sig_yield  = sig_params[0]
    sig_mean   = sig_params[1]
    sig_width  = sig_params[2]
    sig_alpha1 = sig_params[3] if not args.gaus else 0.
    sig_alpha2 = sig_params[4] if not args.gaus else 0.
    sig_enne1  = sig_params[5] if not args.gaus else 0.
    sig_enne2  = sig_params[6] if not args.gaus else 0.

    mean   = rt.RooRealVar("mean"  +mpoint, "mean"  , sig_mean)
    sigma  = rt.RooRealVar("sigma" +mpoint, "sigma" , sig_width)
    if gaus:
        pdf    = rt.RooGaussian("pdf"+str(mass_point), "PDF", obs, mean, sigma)
    else:
        alpha1 = rt.RooRealVar("alpha1"+mpoint, "alpha1", sig_alpha1)
        alpha2 = rt.RooRealVar("alpha2"+mpoint, "alpha2", sig_alpha2)
        enne1  = rt.RooRealVar("enne1" +mpoint, "enne1" , sig_enne1)
        enne2  = rt.RooRealVar("enne2" +mpoint, "enne2" , sig_enne2)
        pdf    = rt.RooDoubleCrystalBall("pdf"+mpoint, "PDF", obs, mean, sigma, alpha1, enne1, alpha2, enne2)
    pdf.plotOn(frame, rt.RooFit.Normalization(sig_yield))
    max_yield = max(max_yield, sig_yield)

for h in hmasses:
    h.Rebin(2)
    dh = rt.RooDataHist("dh"+h.GetName(), "Mass distribution", obs, h)
    dh.plotOn(frame,rt.RooFit.MarkerSize(0.6), rt.RooFit.MarkerColor(h.GetLineColor()))

c = rt.TCanvas()
frame.Draw()
c.SaveAs(figdir+"interpolation.png")
frame.GetYaxis().SetRangeUser(1.e-5*max_yield,5.*max_yield)
c.SetLogy()
c.SaveAs(figdir+"interpolation_log.png")

