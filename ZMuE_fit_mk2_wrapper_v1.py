#Compitible with mk2 
import os
import argparse
import ROOT as rt
rt.gInterpreter.Declare('#include "SFBDT_weight_combined.h"')


parser = argparse.ArgumentParser()
parser.add_argument("-o", dest="name",default="test", type=str,help="output root name")
parser.add_argument("--data-file", dest="data_file",default="", type=str,help="data file")
parser.add_argument("--sgn-file", dest="sgn_file",default="", type=str,help="data file")
parser.add_argument("--zmm-file", dest="zmm_file",default="", type=str,help="Zmm file")
parser.add_argument("--bkg-file", dest="bkg_file",default="", type=str,help="Zmm file")
parser.add_argument("--xgb-min", dest="xgb_min",default="0.75", type=str,help="data file")
parser.add_argument("--xgb-max", dest="xgb_max",default="1.01", type=str,help="data file")
parser.add_argument("--add-pol-order", dest="add_pol_order",type=str, default="0",help="data file")
parser.add_argument("--add-exp-order", dest="add_exp_order",type=str, default="0",help="data file")
parser.add_argument("--add-plaw-order", dest="add_plaw_order",type=str, default="0",help="data file")
parser.add_argument("--create-shape-dc", dest="shape_dc",default=False, action='store_true',help="shape experiment")
parser.add_argument("--lumi", dest="lumi",default=0.0, type=float,help="lumi to compute signal")
parser.add_argument("--year", dest="year",default="Run2", type=str,help="year")
parser.add_argument("--fit-version", dest="ver",default="12", type=str,help="fit version")
parser.add_argument("--outvar", dest="outvar",default="mass_ll",type=str,help="data file")
parser.add_argument("--param-name", dest="param_name",default="bin",type=str,help="data file")
parser.add_argument("--skip-fit", dest="skip_fit",default=False, action='store_true',help="fit skip")
parser.add_argument("--component", dest="component",type=str,default="all",help="Options: all sgn zmm bkg data ")
parser.add_argument("--skip-sgn-syst", dest="skip_sgn_syst",default=False, action='store_true',help="shape experiment")
parser.add_argument("--run-pseudodata", dest="run_pseudodata",default=False, action='store_true',help="run pseudodata instead of data")
parser.add_argument("--run-histo", dest="run_histo",default=False, action='store_true',help="run pseudodata instead of data")
parser.add_argument("--histo-toy", dest="histo_toy",default=False, action='store_true')


args, unknown = parser.parse_known_args()


if len(unknown)>0: 
   print "not found:",unknown,"exiting"
   exit()


########################## cofiguration flags #################################
shape_dc = "false"
do_sgn_syst="true"
run_histo="false"
histo_toy="false"

if args.shape_dc:
   shape_dc="true"
if args.skip_sgn_syst:
   do_sgn_syst="false"
if args.run_histo:
   run_histo="true"
if args.histo_toy:
   histo_toy="true"

if args.run_histo:
   if "bin1" in args.data_file :
       args.xgb_min, args.xgb_max = "0.3","0.7"
   if "bin2" in args.data_file :
       args.xgb_min, args.xgb_max = "0.7","0.9"
   if "bin3" in args.data_file :
       args.xgb_min, args.xgb_max = "0.9","1.1" 

############################### file input ####################################
path="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/"
systematic_path="../BDT/Systematics_v2_no_met_significance/"

sgn_files={"2016":"Meas_fullAndSF_bdt_v7_signal_mcRun16*.root","2017":"Meas_fullAndSF_bdt_v7_signal_mcRun17*.root","2018":"Meas_fullAndSF_bdt_v7_signal_mcRun18*.root","Run2":"Meas_fullAndSF_bdt_v7_signal_mcRun1*.root"}
data_files={"2016":"Meas_full_bdt_v7_data_emu_Run16.root","2017":"Meas_full_bdt_v7_data_emu_Run17.root","2018":"Meas_full_bdt_v7_data_emu_Run18.root","Run2":"Meas_full_bdt_v7_data_emu_Run1*.root"}
ndens={"2016":180400.,"2017":187000.,"2018":194000.,"Run2":561400.}
lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":137.6}


pseudodata_r0_path = "./pseudo_data_from_MC_v2_r0.root"


if args.lumi==0:
   args.lumi=lumis[args.year]
   
if args.data_file=="":
   args.data_file=path+data_files[args.year]
print args.data_file
if args.sgn_file=="":
   args.sgn_file=path+sgn_files[args.year]



###################### expected signal yield ##################################
nden=ndens[args.year]
rt.gROOT.SetBatch(True)
sf="Muon_RecoID_wt*Muon_ID_wt*Muon_IsoID_wt*Muon_dxydz_wt*Electron_RecoID_wt*Electron_ID_wt*Electron_IsoID_wt*Electron_dxydz_wt*PU_wt*PtZ_wt*Trg_wt*SFbtag_wt*JetPUIDWeight*PtSignal_wt*MixZ_wt*Prefire_wt*(SFBDT_weight_Zmumu(xgb)/2.+SFBDT_weight_Zee(xgb)/2.)"
#sf="1"

nsgn=0.
if args.year =="Run2":
  for year in sgn_files.keys():
    if year=="Run2": continue
    cc=rt.TChain("mytreefit")
    cc.Add(path+"/"+sgn_files[year])
    htemp = rt.TH1F("htemp_"+year,"",1,0,2)
    cc.Draw("1>>htemp_"+year,sf+"*(70<mass_ll && mass_ll<110  && Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<xgb && xgb<="+args.xgb_max+")")
    nnum = htemp.Integral()
    print "total pass count ("+year+") =",nnum,"den",ndens[year],"ae",nnum/ndens[year]
    nsgn += lumis[year]*6077220.0/(3*0.0336)*2.62e-7*nnum/ndens[year]
else:
  cc=rt.TChain("mytreefit")
  cc.Add(args.sgn_file)
  htemp = rt.TH1F("htemp","",1,0,2)
  cc.Draw("1>>htemp",sf+"*(70<mass_ll && mass_ll<110  && Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<xgb && xgb<="+args.xgb_max+")")
  nnum = htemp.Integral()
  print "total pass =",nnum
  nsgn = args.lumi*6077220.0/(3*0.0336)*2.62e-7*nnum/nden

print "sgn= ",nsgn


########################## pseudodata normalization ##########################
run_on_pseudodata="false"
pseudodata_norm=-1.0
if args.run_pseudodata:
   run_on_pseudodata="true"
   cc_real=rt.TChain("mytreefit")
   cc_real.Add(path+data_files["Run2"])
   htemp_real = rt.TH1F("htemp_real","",1,0,2)
   cc_real.Draw("1>>htemp_real","(70<mass_ll && mass_ll<110  && Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<=xgb && xgb<"+args.xgb_max+")")
   cc_pseudo=rt.TChain("mytreefit")
   cc_pseudo.Add(pseudodata_r0_path)
   htemp_pseudo = rt.TH1F("htemp_pseudo","",1,0,2)
   cc_pseudo.Draw("1>>htemp_pseudo","NormGen_wt*( 70<mass_ll && mass_ll<110  && Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<xgb && xgb<="+args.xgb_max+")")
   pseudodata_norm= htemp_real.Integral()/htemp_pseudo.Integral()
   print "BKG normalization factor",pseudodata_norm,"pseudo",htemp_pseudo.GetBinContent(1),"real",htemp_real.GetBinContent(1),"entries real",cc_real.GetEntries(),"entries gen",cc_pseudo.GetEntries()
   hmll_real = rt.TH1F("hmll_real","",40,70,110)
   hmll_pseudo = rt.TH1F("hmll_pseudo","",40,70,110)
   cc_real.Draw("mass_ll>>hmll_real","(70<mass_ll && mass_ll<110  && Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<=xgb && xgb<"+args.xgb_max+")")
   cc_pseudo.Draw("mass_ll>>hmll_pseudo",str(pseudodata_norm)+"*NormGen_wt*( 70<mass_ll && mass_ll<110  && Flag_met && Flag_muon && Flag_electron && "+args.xgb_min+"<xgb && xgb<="+args.xgb_max+")")
   c1 = rt.TCanvas("c1","",700,700)
   hmll_real.Draw()
   hmll_pseudo.Draw("sames")
   c1.SaveAs("cfit_wrapper_norm_"+args.name+".png")
   


###############################################################################
############################## FIT ############################################
if args.skip_fit:
   exit()

os.system(". ./Load_cmssw.sh")
if args.component == "sgn" or args.component == "all":
  os.system('root -l -b -q ZMuE_fit_mk2_sgn_v'+args.ver+'.C\'("'+args.name+'","'+args.sgn_file+'","'+args.xgb_min+'","'+args.xgb_max+'",'+shape_dc+','+str(nsgn)+',"'+args.outvar+'",'+do_sgn_syst+',"'+args.param_name+'")\'')

if args.component == "zmm" or args.component == "all":
  os.system('root -l -b -q ZMuE_fit_mk2_Zmm_v'+args.ver+'.C\'("'+args.name+'","'+args.zmm_file+'","'+args.xgb_min+'","'+args.xgb_max+'",'+shape_dc+',"'+args.outvar+'","'+args.param_name+'",'+str(pseudodata_norm)+','+run_histo+','+histo_toy+')\'')   

if args.component == "bkg" or args.component == "all":
  os.system('root -l -b -q ZMuE_fit_mk2_bkg_v'+args.ver+'.C\'("'+args.name+'","'+args.bkg_file+'","'+args.xgb_min+'","'+args.xgb_max+'",'+shape_dc+',"'+args.outvar+'","'+args.param_name+'",'+run_on_pseudodata+','+str(pseudodata_norm)+','+run_histo+','+histo_toy+','+args.add_pol_order+','+args.add_exp_order+','+args.add_plaw_order+')\'')

if args.component == "data" or args.component == "all":
  os.system('root -l -b -q ZMuE_fit_mk2_datagen_v'+args.ver+'.C\'("'+args.name+'","'+args.data_file+'","'+args.bkg_file+'","'+args.xgb_min+'","'+args.xgb_max+'","'+args.outvar+'","'+args.param_name+'",'+str(pseudodata_norm)+','+run_histo+')\'')   

