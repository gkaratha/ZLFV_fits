import ROOT as rt


infiles =["template_zemu_embed_v3_bin1.root", "template_zemu_embed_v3_bin2.root", "template_zemu_embed_v3_bin3.root"]
outfiles =["bin1", "bin2", "bin3"]
Fzmm=1.00

rt.gROOT.SetBatch(True)


for idx in range(len(infiles)):

  tfile = rt.TFile(infiles[idx],"READ")  
  stack = tfile.Get("bkg_stack")
  th1s = stack.GetHists()
  
  c1 = rt.TCanvas("c"+str(idx),"c1",800,600)
  leg = rt.TLegend(0.7,0.7,1,1)
  htotal = rt.TH1F()
  for ith1,th1 in enumerate(th1s):
    print th1.GetName()
    leg.AddEntry(th1,th1.GetName())
    th1.SetFillColor(0)
    th1.SetLineWidth(2)
    if ith1==0:
      htotal=th1.Clone()
      htotal.SetName("hbkg_safe")
      htotal.SetTitle("hbkg_safe")
      th1.Draw("HIST")
      th1.GetYaxis().SetRangeUser(0,1000)
    else:
      if "Z->" in th1.GetName():
        th1.Scale(Fzmm)
      htotal.Add(th1)
      th1.Draw("HIST sames")
  
  leg.Draw("sames")
  
  htotal.SetLineColor(1)
  htotal.Draw("HIST sames")
  leg.AddEntry(htotal,"total")
  
  leg.Draw("sames")
  c1.SaveAs("ctmpl_embed_ZmmF_"+str(Fzmm)+"_"+outfiles[idx]+".png")
  
  fout = rt.TFile("template_zemu_v3_custom_ZmmF_"+str(Fzmm)+"_"+outfiles[idx]+".root","RECREATE")
  htotal.Write();

