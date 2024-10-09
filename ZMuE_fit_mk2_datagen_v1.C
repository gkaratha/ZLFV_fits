#include <fstream>
#include "fit_helper.h"
#include "fit_helper_core.h"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// mk2: seperating fit in components
///////////// MC PDF and data obs workspace v1
//////////////////// +add more options for data obs
//////////////////// +



int ZMuE_fit_mk2_datagen_v1(TString name="bin1_r2",
    TString data_file="pseudo_data_from_MC_v2_r0.root",
    TString mc_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSF_bdt_v7_signal_mcRun1*.root",
    TString xgbmin="0.3",TString xgbmax="0.7", TString outvar="mass_ll",
    TString varname="bin", float pseudodata_norm=-1.0,
    bool histo_mc_input=false){

   //////////////////////////////////// configuration /////////////////////////
   gROOT->SetBatch(true);

   TString cuts = xgbmin+"<xgb && xgb<"+xgbmax+" && Flag_met && Flag_muon && Flag_electron"; // not add mass_ll here when run systematics
   TString tree_name="mytreefit";

   int nbin_data=80;
   double min_fit_range=70;
   double max_fit_range=110;

   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   ///////////////////////////////////////////////////////////////////////////
   cout<<" \n ********************************************************** "<<endl;
   cout<<" ********************************************************** "<<endl;
   cout<<" ***ZMuE mk2 fit v1: data obs and mc pdf prep ***Output name:"+name<<endl;
   cout<<" ********************************************************** "<<endl;


   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   dilep_mass_out.setBins(nbin_data);

   // read real data
   TTree * data_tree = get_tree("mytreefit",data_file,cuts);
   TH1F * hdata = new TH1F("hdata","",nbin_data, min_fit_range , max_fit_range);
   cout << "*** Retrieving the data...\n";
   data_tree->Draw("mass_ll>>hdata",cuts);

   RooDataHist * dhist_data = new RooDataHist("real_data_"+varname,"real_data_"+varname,RooArgSet(dilep_mass_out),hdata);

   // read MC
   TH1F * hmc = new TH1F("hmc","",nbin_data,min_fit_range,max_fit_range);
   TH1F * hmc_nozmm = new TH1F("hmc_nozmm","",nbin_data,min_fit_range,max_fit_range);
   if (histo_mc_input) {
      TFile * fmc = new TFile(mc_file,"READ");
      THStack * stack = (THStack *) fmc->Get("bkg_stack");
      auto * hists = stack -> GetHists();
      for (auto * hist: *hists){
         TH1F* th1 = (TH1F*) hist;
         hmc->Add(th1);
         TString hist_name = th1->GetName();
         if (!hist_name.BeginsWith("Z"))
            hmc_nozmm->Add(th1);
      }
   }
   else{
       TChain * cc = new TChain(tree_name);
       cc->Add(mc_file);
       if (pseudodata_norm>0){
          cout << "*** Retrieving the total MC...\n";
          cc->Draw("mass_ll>>hmc",TString(std::to_string(pseudodata_norm))+"*Total_wt*( "+cuts+" )");
          cout << "*** Retrieving the non-Zmm MC...\n";
          cc->Draw("mass_ll>>hmc_nozmm",TString(std::to_string(pseudodata_norm))+"*Total_wt*( "+cuts+" && !IsZmm  && !IsZee )");
       } else{
          cout << "*** Retrieving the total MC...\n";
          cc->Draw("mass_ll>>hmc","Total_wt*( "+cuts+" )");
          cout << "*** Retrieving the non-Zmm MC...\n";
          cc->Draw("mass_ll>>hmc_nozmm","Total_wt*( "+cuts+"  && !IsZmm  && !IsZee )");
       }
   }

   cout << "*** Retrieved all inputs!\n";
   RooDataHist * dhist_mc = new RooDataHist("mc_dataset_"+varname,"mc_dataset_"+varname,RooArgSet(dilep_mass_out),hmc);
   RooDataHist * dhist_mc_nozmm = new RooDataHist("mc_nozmm_dataset_"+varname,"mc_nozmm_dataset_"+varname,RooArgSet(dilep_mass_out),hmc_nozmm);

   // create mc pdf
   RooHistPdf pdf_template("mc_pdf_"+varname,"",RooArgSet(dilep_mass_out),*dhist_mc);
   RooRealVar pdf_template_norm("mc_pdf_"+varname+"_norm","",dhist_data->sumEntries(),0,2*dhist_data->sumEntries());
   RooHistPdf pdf_nozmm_template("mc_nozmm_pdf_"+varname,"",RooArgSet(dilep_mass_out),*dhist_mc_nozmm);
   RooRealVar pdf_nozmm_template_norm("mc_nozmm_pdf_"+varname+"_norm","",dhist_data->sumEntries(),0,2*dhist_data->sumEntries());

   // write everything
   RooWorkspace *wspace = new RooWorkspace("workspace_data","workspace_data");
   cout << "*** Importing the data...\n";
   wspace->import(*dhist_data);

   // create single mc toy
   const bool add_toy = false;

   if(add_toy) {
     cout << "*** Generate toy dataset\n";
     RooDataSet * mc_toy = pdf_template.generate(RooArgSet(dilep_mass_out),dhist_data->sumEntries());
     RooDataHist dhist_mc_toy("mc_toy_"+varname,"mc_toy_"+varname,dilep_mass_out,*mc_toy);
     cout << "*** Generate no Zmm toy dataset\n";
     RooDataSet * mc_nozmm_toy = pdf_nozmm_template.generate(RooArgSet(dilep_mass_out),dhist_data->sumEntries());
     RooDataHist dhist_mc_nozmm_toy("mc_nozmm_toy_"+varname,"mc_nozmm_toy_"+varname,dilep_mass_out,*mc_nozmm_toy);

     cout << "*** Importing the MC...\n";
     wspace->import(*dhist_mc);
     cout << "*** Importing the non-Zmm MC...\n";
     wspace->import(*dhist_mc_nozmm);

     wspace->import(dhist_mc_toy);
     wspace->import(dhist_mc_nozmm_toy);
   }

   cout << "*** Importing the templates...\n";
   wspace->import(pdf_template);
   wspace->import(pdf_template_norm);
   wspace->import(pdf_nozmm_template);
   wspace->import(pdf_nozmm_template_norm);

   cout << "*** Writing to disk...\n";
   wspace->writeToFile("workspace_mk2dataAndMC_v1_"+name+".root");


 return 0;
}
