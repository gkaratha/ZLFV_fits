import os
import argparse
import ROOT as rt
rt.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("--signal", dest="signal", required=True, type=str, help="Signal file")
parser.add_argument("--background", dest="background", required=True, type=str, help="Background file")
parser.add_argument("--tree", dest="tree", default="mytreefit", type=str, help="Tree name")

args, unknown = parser.parse_known_args()
  
for arg in unknown:
  print "warning uknown parameter",arg
if len(unknown) > 0: exit()

f_s = rt.TFile.Open(args.signal    , 'READ')
f_b = rt.TFile.Open(args.background, 'READ')
if not f_s or not f_b: exit()

t_s = f_s.Get(args.tree)
t_b = f_b.Get(args.tree)

h_s = rt.TH1D('h_s', 'BDT score comparison', 30, 0., 1.)
h_b = rt.TH1D('h_b', 'BDT score comparison', 30, 0., 1.)
t_s.Draw('xgb >> h_s', '', 'goff')
t_b.Draw('xgb >> h_b', '', 'goff')

h_s.Scale(1./h_s.Integral())
h_b.Scale(1./h_b.Integral())

c = rt.TCanvas()
h_s.Draw('hist')
h_s.SetLineWidth(2)
h_s.SetLineColor(rt.kBlue)

h_b.Draw('hist same')
h_b.SetLineWidth(2)
h_b.SetLineColor(rt.kRed)

rt.gStyle.SetOptStat(0)
max_val = max(h_b.GetMaximum(), h_s.GetMaximum())
h_s.SetAxisRange(max_val/1.e4, 1.2*max_val, "Y")
c.SaveAs('compare_bdt.png')
h_s.SetAxisRange(max_val/1.e4, 3.*max_val, "Y")
c.SetLogy()
c.SaveAs('compare_bdt_log.png')

leg = rt.TLegend()
hists = []
c.SetLogy(False)
for index in range(4):
  bdt_min = index*0.25
  bdt_max = bdt_min + 0.25
  h = rt.TH1D('h_mcol_%i' % (index), 'Collinear mass vs. BDT', 30, 80, 300)
  t_b.Draw('mcol >> ' + h.GetName(), 'xgb > ' + str(bdt_min) + ' && xgb <= ' + str(bdt_max), 'goff')
  h.Scale(1./h.Integral())
  if index == 0: h.Draw('hist')
  else         : h.Draw('hist same')
  h.SetLineWidth(2)
  h.SetLineColor(index + 1)
  leg.AddEntry(h, str(bdt_min) + ' < BDT < ' + str(bdt_max), 'L')
  hists.append(h)
max_val = max([h.GetMaximum() for h in hists])
hists[0].GetYaxis().SetRangeUser(max_val/50., 1.1*max_val)
leg.Draw('same')
c.SaveAs('compare_bdt_mcol.png')
