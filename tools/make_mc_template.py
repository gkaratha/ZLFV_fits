# Make PDF templates from an input histogram file
import os
import argparse
import ROOT as rt


#----------------------------------------------
# Fit an input histogram
def fit_hist(hist, figdir, min_mass = 95.):
    rt.gStyle.SetOptStat(0)
    name = hist.GetName()
    if name == 'top':
        func = rt.TF1('func', 'pol2(0)', min_mass, 150)
    elif name == 'higgs' or name == 'zmumu':
        func = rt.TF1('func', 'TMath::Exp([0] + (x/100.)*[1]) + [2] + [3]*(x/100.) + [4]*TMath::Power(x/100.,2)', min_mass, 150)
        func.SetParameters(10., -10., 0.1, -1., 1.)
    else:
        func = rt.TF1('func_' + name, 'expo(0)', min_mass, 150)
    hist.Fit(func, 'R 0')
    c = rt.TCanvas()
    hist.SetLineWidth(2)
    hist.SetLineColor(rt.kBlue)
    hist.SetMarkerColor(rt.kBlue)
    hist.SetMarkerSize(0.8)
    hist.SetMarkerStyle(20)
    hist.Draw('E1')
    func.Draw('same')
    hist.SetXTitle('M_{e#mu} (GeV/c^{2})')
    hist.SetYTitle('Events / %.1f GeV/c^{2}' % (hist.GetBinWidth(1)))
    hist.GetXaxis().SetRangeUser(min_mass, 150.)
    c.SaveAs(figdir+name+'_fit.png')
    c.SetLogy()
    c.SaveAs(figdir+name+'_fit_log.png')
    return func

rt.gROOT.SetBatch(True)

#----------------------------------------------
# Read in the input parameters
#----------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("--set", dest="set",default=11, type=int,help="Histogram file set")
parser.add_argument("--tag", dest="tag",default="v01", type=str,help="Output naming tag")
parser.add_argument("--mass", dest="mass",default=110., type=float,help="Signal mass for the template")
parser.add_argument("--min-mass", dest="min_mass",default=95., type=float,help="Minimum mass used in the data card")
parser.add_argument("--max-mass", dest="max_mass",default=700., type=float,help="Maximum mass used in the data card")

args, unknown = parser.parse_known_args()
if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()


#----------------------------------------------
# Output data paths
#----------------------------------------------

figdir = "./figures/templates/set_%i_%s/" % (args.set, args.tag)
outdir = "./templates/"
os.system("[ ! -d %s ] && mkdir -p %s" % (figdir , figdir))
os.system("[ ! -d %s ] && mkdir -p %s" % (outdir , outdir))

#----------------------------------------------
# Read in the input MC info
#----------------------------------------------

f_in = rt.TFile.Open('hists/zprime_lepm_%i_2016_2017_2018.hist' % (args.set), 'READ')
stack = f_in.Get('hstack')
higgs = stack.GetHists().At(0)
top = stack.GetHists().At(1)
ww = stack.GetHists().At(2)
zmumu = stack.GetHists().At(3) # Not used here
ztautau = stack.GetHists().At(4)
qcd = stack.GetHists().At(5)

higgs.SetName('higgs')
top.SetName('top')
ww.SetName('ww')
zmumu.SetName('zmumu')
ztautau.SetName('ztautau')
qcd.SetName('qcd')

higgs.SetTitle('H#rightarrow#tau#tau/WW')
top.SetTitle('Top')
ww.SetTitle('WW')
zmumu.SetTitle('Z#rightarrow#mu#mu')
ztautau.SetTitle('Z#rightarrow#tau#tau')
qcd.SetTitle('QCD')

#----------------------------------------------
# Process the histograms of interest
#----------------------------------------------

if args.min_mass < 95.: #only include Z->mumu below this threshold
    hists = [ww, top, qcd, zmumu, ztautau, higgs]
else:
    hists = [ww, top, qcd, ztautau, higgs]
fit_res = []

for hist in hists: fit_res.append(fit_hist(hist, figdir, args.min_mass))

#----------------------------------------------
# Make a background template from the fits
#----------------------------------------------

mass = args.mass
width = 1./50.*mass
min_mass = max(mass - 10.*width, args.min_mass)
max_mass = min(mass + 10.*width, args.max_mass)
nbins = int((max_mass - min_mass) / (0.5 * width))
hbkg = rt.TH1D('hbkg', 'Background template', nbins, min_mass, max_mass)
x_w_in = ww.GetBinWidth(1)
x_w_out = hbkg.GetBinWidth(1)
print x_w_in, x_w_out

for xbin in range(nbins):
    x = hbkg.GetBinCenter(xbin+1)
    bkg_val = 0.
    for fit in fit_res: bkg_val += fit.Eval(x)/x_w_in*x_w_out
    hbkg.SetBinContent(xbin+1, bkg_val)
    hbkg.SetBinError(xbin+1, 0.)

c = rt.TCanvas()
hbkg.Draw('hist')
hbkg.SetLineWidth(2)
hbkg.SetXTitle('M_{e#mu} (GeV/c^{2})')
hbkg.SetYTitle('Events / %.1f GeV/c^{2}' % (x_w_out))
c.SaveAs(figdir+'template.png')

#----------------------------------------------
# Make a data card for the MC template
#----------------------------------------------
cardname = 'datacard_zprime_mc_template_'  + args.tag + '_' + str(args.set) + ('_mass-%.2f' % (mass)) + '.txt'
wsname   = 'workspace_zprime_mc_template_' + args.tag + '_' + str(args.set) + ('_mass-%.2f' % (mass)) + '.root'

# Standard preamble
txt="# -*- mode: tcl -*-\n"
txt+="#Auto generated Z prime search COMBINE datacard\n"
txt+="#Using Z prime mass = %.3f\n\n" % (mass)
txt+="imax * number of channels\njmax * number of backgrounds\nkmax * number of nuisance parameters\n\n"

# Define the background and signal PDFs
txt+="#-----------------------------------------------------------------------------------------------------------\n"
txt+="shapes signal      bin1 %s ws:sig_pdf\n" % (wsname)
txt+="shapes background  bin1 %s ws:bkg_pdf\n" % (wsname)
txt+="shapes data_obs    bin1 %s ws:data_obs\n" % (wsname)

txt+="#-----------------------------------------------------------------------------------------------------------\n\n"
# Define the channel
txt+="bin         bin1\n"
txt+="observation  -1\n\n"
txt+="bin         bin1      bin1\n"
txt+="process    signal   background\n\n"
txt+="process       0         1\n\n"
txt+="rate          1         1\n" #Rate is taken from the _norm variables

# Define the uncertainties
txt+="#-----------------------------------------------------------------------------------------------------------\n"
txt+="SignalSys        lnN 1.10     -\n"

# Write the file
with open(outdir+cardname,'w') as fl:
    fl.write(txt)
fl.close()

#----------------------------------------------
# Make a workspace for the datacard
#----------------------------------------------

ws = rt.RooWorkspace('ws')
obs = rt.RooRealVar('mass_ll', 'Dilepton mass', (min_mass+max_mass)/2., min_mass, max_mass)
bkg_data = rt.RooDataHist('bkg_hist', 'Background template', obs, hbkg)
bkg_pdf = rt.RooHistPdf('bkg_pdf', 'Background PDF', obs, bkg_data)
bkg_norm = rt.RooRealVar('bkg_pdf_norm', 'Background normalization', hbkg.Integral(), 0., 1.e6)
ws.Import(bkg_pdf)
ws.Import(bkg_norm)
mc_data = bkg_pdf.generateBinned(obs, hbkg.Integral())
mc_data.SetName('data_obs')
ws.Import(mc_data)

# Create a signal PDF
sig_mean  = rt.RooRealVar('sig_mean', 'Signal mean', mass)
sig_width = rt.RooRealVar('sig_width', 'Signal width', mass/50.)
sig_mean.setConstant(1)
sig_width.setConstant(1)
sig_pdf = rt.RooGaussian('sig_pdf', 'Signal PDF', obs, sig_mean, sig_width)
sig_norm = rt.RooRealVar('sig_pdf_norm', 'Background normalization', 50.)
sig_norm.setConstant(1)
ws.Import(sig_pdf)
ws.Import(sig_norm)
ws.writeToFile(outdir+wsname)
