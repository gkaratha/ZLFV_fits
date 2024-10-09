# Test the flatness of the MC background sources
import os
import argparse
import ROOT as rt
from array import array


#---------------------------------------------------------------
# Create a mass histogram for the given process and cuts
def mass_distribution(name, file_path, cuts):
    # Create a TChain for the signal
    cc=rt.TChain("mytreefit")
    cc.Add(file_path)

    # Create a mass histogram
    h = rt.TH1F(name,"Mass distribution", 120, 0,600)
    cc.Draw("mass_ll>>"+name, cuts)
    if h.Integral() == 0.: return h
    scale = 1. / h.Integral()
    h.Scale(scale)
    return h

#---------------------------------------------------------------
# Plot the mass distributions
def plot_shapes(name, h_inc, h_low, h_high):
    first_val = h_inc.GetBinLowEdge(h_inc.FindFirstBinAbove(0.))
    xmin = first_val if first_val < 95. else 95.
    c = rt.TCanvas('c', 'c', 700, 500)
    leg = rt.TLegend(0.1, 0.75, 0.9, 0.9)
    leg.SetTextSize(0.04)
    leg.SetNColumns(3)
    h_inc.SetTitle("Inclusive")
    h_low.SetTitle("0.3-0.7")
    h_high.SetTitle("0.7-1.0")
    colors = [rt.kBlack, rt.kBlue, rt.kRed]
    hists = [h_inc, h_low, h_high]
    first = True
    for color,h in zip(colors, hists):
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.8)
        h.SetLineWidth(2)
        h.Draw("E1" if first else "E1 same")
        leg.AddEntry(h)
        first = False
    leg.Draw()
    max_val = max([h_inc.GetMaximum(), h_low.GetMaximum(), h_high.GetMaximum()])
    h_inc.SetTitle("Mass shape vs. BDT score")
    h_inc.SetXTitle("Mass (GeV/c^{2})")
    h_inc.GetXaxis().SetRangeUser(xmin, h_inc.GetXaxis().GetXmax())
    h_inc.GetYaxis().SetRangeUser(0., 1.3*max_val)
    c.SaveAs(name + "shape.png")
    h_inc.GetYaxis().SetRangeUser(1.e-5*max_val, 50*max_val)
    c.SetLogy()
    c.SaveAs(name + "shape_log.png")

    # Make a zoomed in on 95-200 GeV plot
    h_inc.GetXaxis().SetRangeUser(xmin, 150. if first_val < 95. else 200.)
    h_inc.GetYaxis().SetRangeUser(0., 1.3*max_val)
    c.SetLogy(0)
    c.SaveAs(name + "shape_zoom.png")
    h_inc.GetYaxis().SetRangeUser(1.e-2*max_val, 5*max_val)
    c.SetLogy()
    c.SaveAs(name + "shape_zoom_log.png")

#---------------------------------------------------------------
# Main processing

rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument("--tag", dest="tag",default="",type=str,help="Output figure directory tag")
parser.add_argument("--min-mass", dest="min_mass",default=95.,type=float,help="Minimum mass in the tests")
args, unknown = parser.parse_known_args()

### default path
# path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/"
path="/eos/user/m/mimacken/ZEMu/CMSSW_11_3_4/src/ZLFV_fits/trees/" #FIXME: Replace these with official versions
figdir = "./figures/flatness/"
if args.tag != "": figdir += args.tag + '/'
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir ))

### MC backgrounds
mcs=["dy","htautau"]
mc_files={
    "ww"     :"forMeas_bdt_v7_emu_scan_WW_*.root",\
    "ttbar"  :"forMeas_bdt_v7_emu_scan_ttbarlnu_mcRun201*.root",\
    "htautau":"forMeas_bdt_v7_emu_scan_ggFH-TauTau_mcRun201*.root",\
    "dy"     :"forMeas_bdt_v7_emu_scan_DY50*.root",\
}

# Loop through the MC processes
bdts = []
cdfs = []
base_cuts = "(mass_ll > " + str(args.min_mass) + " && Flag_met && Flag_muon && Flag_electron)"
for mc in mcs:
    # Get the mass distribution inclusive, 0.3-0.7, and 0.7-1.0
    h_inc  = mass_distribution(mc + "_inc" , path + mc_files[mc], base_cuts)
    h_low  = mass_distribution(mc + "_low" , path + mc_files[mc], base_cuts + " && (xgb > 0.3 && xgb <= 0.7)")
    h_high = mass_distribution(mc + "_high", path + mc_files[mc], base_cuts + " && (xgb > 0.7 && xgb <= 1.0)")
    if h_inc.Integral() == 0. or h_low.Integral() == 0. or h_high.Integral() == 0.: continue
    plot_shapes(figdir + mc + "_", h_inc, h_low, h_high)

