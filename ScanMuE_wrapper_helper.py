from array import array
import ROOT as rt

#----------------------------------------------------------------------------------------
# Write a combine data card for a single bin
def print_datacard(name, sig_file, bkg_file, param_name, mass):
   # sig_file = "WorkspaceSGN/workspace_scansgn_v" + ver + "_" + name + "_mp" + str(cnt) + ".root"
   # bkg_file = "WorkspaceBKG/workspace_scanbkg_v" + ver + "_" + name + "_mp" + cnt + ".root"
      
   # Standard preamble
   txt="# -*- mode: tcl -*-\n"
   txt+="#Auto generated Z prime search COMBINE datacard\n"
   txt+="#Using Z prime mass = %.3f\n\n" % (mass)
   txt+="imax * number of channels\njmax * number of backgrounds\nkmax * number of nuisance parameters\n\n"

   # Define the background and signal PDFs
   txt+="### inputs-----------------------------------------------------------------------------------------------------------\n"
   txt+="shapes signal "+param_name+" %s ws_sgn:signal_pdf_%s\n" % (sig_file, param_name)
   txt+="shapes background "+param_name+" %s ws_bkg:multipdf_%s\n" % (bkg_file, param_name)
   txt+="shapes data_obs "+param_name+" %s ws_bkg:data_obs\n" % (bkg_file)

   txt+="### bins-----------------------------------------------------------------------------------------------------------\n\n"
   # Define the channel
   txt+="bin         "+param_name+"\n"
   txt+="observation  -1\n\n"
   txt+="bin         "+param_name+"    "+param_name+"\n"
   txt+="process    signal   background\n\n"
   txt+="process       0         1\n\n"
   txt+="rate          1         1\n" #Rate is taken from the _norm variables

   # Define the uncertainties
   txt+="### uncertainties-----------------------------------------------------------------------------------------------------------\n"
   #FIXME: Implement mass-dependent uncertainties
   txt+="ElectronID lnN 1.02     -\n"
   txt+="MuonID     lnN 1.02     -\n"
   txt+="Lumi       lnN 1.02     -\n"
   txt+="BTag       lnN 1.005    -\n"
   txt+="Theory     lnN 1.01     -\n"
   txt+="BDT        lnN 1.02     -\n"
   txt+="SigModel   lnN 1.08     -\n" #Uncertainty due to the interpolation

   txt+="#### scales-----------------------------------------------------------------------------------------------------------\n\n"
   # Scale uncertainties
   txt+="elec_ES_shift param 0 1 [-7, 7]\n"
   txt+="muon_ES_shift param 0 1 [-7, 7]\n"

   txt+="#### pdfs-----------------------------------------------------------------------------------------------------------\n\n"
   # Define the envelope discrete index to be scanned
   txt+="pdfindex_%s discrete\n" % (param_name)

   # Write the file
   with open(name,'w') as fl:
     fl.write(txt)
   fl.close()



#----------------------------------------------------------------------------------------
# Make a plot and save the figure, fit a 2nd order polynomial to it
def plot_graph(mpoint_array,signal_array,name,xaxis="",yaxis=""):
   signal_vs_mass = rt.TGraph(len(mpoint_array), mpoint_array, signal_array)

   fnc = rt.TF1("fnc", "pol2", 100, 500)
   par_yld = array('d',[0, 0, 0])
   cnv = rt.TCanvas("cnv_"+name,"",800,600)
   signal_vs_mass.Draw("A*")
   signal_vs_mass.SetTitle("#bf{CMS} #it{Preliminary};"+xaxis+";"+yaxis)
   signal_vs_mass.Fit(fnc,"R")
   cnv.SaveAs(name+".png")
   fnc.GetParameters(par_yld)
   return par_yld


#----------------------------------------------------------------------------------------
##### mass points analysis ######
### efficiency and yield vs mass
def mass_analysis(sgn_masspoints,path,sgn_masspoint_files, xgb_min, xgb_max, sf, ndens, lumi, figdir):
  nsgn_array = array('d',[])
  toteff_array=  array('d',[])
  bdteff_array=  array('d',[])
  mpoint_array = array('d',[])
  for mpoint in sgn_masspoints:
    cc=rt.TChain("mytreefit")
    cc.Add(path+"/"+sgn_masspoint_files[mpoint])
    htemp = rt.TH1F("htemp_"+mpoint,"",1,0,2)
    cc.Draw("1>>htemp_"+mpoint,sf+"*(Flag_met && Flag_muon && Flag_electron && "+xgb_min+"<=xgb && xgb<"+xgb_max+")")
    nnum = htemp.Integral()
    cross_section = 1. #in units femto-barns * BR(Z'->emu)
    nsgn_array.append(lumi*1*nnum/ndens[mpoint])
    mpoint_array.append(float(mpoint))
    toteff_array.append(nnum/(1.0*ndens[mpoint]))
    bdteff_array.append(nnum/(1.0*cc.GetEntries()))
    print "total pass ("+mpoint+") =",nnum,"eff",nnum/(1.0*ndens[mpoint]),"yield",nsgn_array[-1]

  par_yield = plot_graph(mpoint_array,nsgn_array,figdir+"yld_vs_mass","m(e,#mu)","Yield")
  par_bdt_eff = plot_graph(mpoint_array,bdteff_array,figdir+"bdteff_vs_mass","m(e,#mu)","BDT eff.")
  par_total_eff = plot_graph(mpoint_array,toteff_array,figdir+"totaleff_vs_mass","m(e,#mu)","Total eff.")

  ### get widths
  width_array = array('d',[])
  for mpoint in sgn_masspoints:
    cc=rt.TChain("mytreefit")
    cc.Add(path+"/"+sgn_masspoint_files[mpoint])
    htemp = rt.TH1F("htemp_"+mpoint,"",50,float(mpoint)*0.8,float(mpoint)*1.2)
    cc.Draw("mass_ll>>htemp_"+mpoint,sf+"*(Flag_met && Flag_muon && Flag_electron && "+xgb_min+"<=xgb && xgb<"+xgb_max+")")
    gauss = rt.TF1("gauss", "gaus", float(mpoint)*0.95,float(mpoint)*1.05)
    par_gs = array('d',[0, 0, 0])
    cnv = rt.TCanvas("cnv_"+mpoint,"",800,600)
    htemp.Draw("HIST")
    htemp.Fit(gauss,"LR")
    htemp.SetTitle("#bf{CMS} #it{Preliminary};m(e,#mu);")
    gauss.Draw("sames")
    cnv.SaveAs(figdir+"mass_fit_"+mpoint+".png")
    gauss.GetParameters(par_gs)
    width_array.append(par_gs[2])
  
  par_width = plot_graph(mpoint_array,width_array,figdir+"width_vs_mass","m(e,#mu)","Width")
  return  par_yield, par_width
