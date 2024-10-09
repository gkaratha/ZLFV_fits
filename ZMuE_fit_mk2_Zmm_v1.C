#include <fstream>
#include "fit_helper.h"
#include "fit_helper_core.h"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////// mk2: seperating fit in components
///////////// Z->mumu component v1
//////////////////// + Same as v401



int ZMuE_fit_mk2_Zmm_v1(TString name="bin1_r2",
    TString zmm_file="pseudo_data_from_MC_v2_r0.root",
    TString xgbmin="0.3",TString xgbmax="0.7", bool create_dc_input=false,
    TString outvar="mass_ll", TString varname="bin", float pseudodata_norm=-1.0,
    bool histo_input=false, bool histo_toy=false){

   //////////////////////////////////// configuration /////////////////////////
   gROOT->SetBatch(true);


   TString cuts = xgbmin+"<xgb && xgb<"+xgbmax+" && Flag_met && Flag_muon && Flag_electron"; // not add mass_ll here when run systematics

   TString tree_name="mytreefit";
   TString histo_name = "bkg_stack";
   double  min_fit_range = 70.;
   double  max_fit_range = 110.;
   double blind_min = 86.;
   double blind_max = 96.;
   int nbin_data = 80;
   bool Verbose=false; // RooFit verbosity
   TString constrain_type_mean = "fixed"; // fixed, linear, gaussian
   int printout_levels=1; // 0: Print only final fits parameters, 1: Print all fits from F test tests

   if (constrain_type_mean !="fixed" && constrain_type_mean !="linear" && constrain_type_mean !="gaussian"){
       cout<<" wrong constraint on mean"<<endl;
       return 0;
   }

   ///////////////////////////////////////////////////////////////////////////
   cout<<" \n ********************************************************** "<<endl;
   cout<<" ********************************************************** "<<endl;
   cout<<" ***ZMuE mk2 fit v1: Z->mumu part  ***Output name:"+name<<endl;
   cout<<" ********************************************************** "<<endl;

   if (!Verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   // read trees
   RooRealVar dilep_mass("mass_ll","m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");

   dilep_mass.setBins(nbin_data);
   dilep_mass_out.setBins(nbin_data);

   ///////////////////////////////// Data/bkg fit ////////////////////////////
   cout<<"\n *********************** Data fit ********************** "<<endl;

   ////// get Zmm and bkg from embeded
   TH1F * hzmm = new TH1F("hzmm","",nbin_data,min_fit_range,max_fit_range);

   if (histo_input) {
     TFile * fzmm = new TFile(zmm_file,"READ");
     THStack * stzmm = (THStack *) fzmm->Get("bkg_stack");
     auto * hists = stzmm -> GetHists();
     for (auto * hist: *hists){
       TString hist_name =hist->GetName();
       TH1F* th1 = (TH1F*) hist;
       if (hist_name.BeginsWith("Z"))
         hzmm->Add(th1);
     }
   } else {
     TChain * cc = new TChain(tree_name);
     cc->Add(zmm_file);
     if (pseudodata_norm>0)
       cc->Draw("mass_ll>>hzmm",TString(std::to_string(pseudodata_norm))+"*NormGen_wt*( (IsZmm || IsZee) && "+cuts+" )");
     else
       cc->Draw("mass_ll>>hzmm","NormGen_wt*( (IsZmm || IsZee) && "+cuts+" )");
   }

   TCanvas * c1 = new TCanvas("c1","",800,600);
   hzmm->Draw();
   hzmm->SetLineWidth(2);
   c1->SaveAs("mk2zmm_v1_templ_"+name+".png");

   RooDataHist * dhist_zmm;
   if (!histo_toy)
       dhist_zmm = new RooDataHist("dhist_zmm_"+varname,"dhist_zmm_"+varname,RooArgSet(dilep_mass),hzmm);
   else{
     RooDataHist * dhist_original_zmm = new RooDataHist("dhist_original_zmm","dhist_original_zmm",RooArgSet(dilep_mass),hzmm);
     RooHistPdf zmm_template("zmm_template","",RooArgSet(dilep_mass),*dhist_original_zmm);
     RooDataSet * dhist_zmm_toy = zmm_template.generate(RooArgSet(dilep_mass),dhist_original_zmm->sumEntries());
     dhist_zmm = new RooDataHist("dhist_zmm_"+varname,"dhist_zmm_"+varname,dilep_mass,* dhist_zmm_toy);
   }

   /// fit ds crystal ball on Z->mm fakes
   RooRealVar* zmm_dcb_mean   = new RooRealVar("zmm_dcb_mean_"+varname , "zmm_dcb_mean_"+varname  , 84.4, 70., 90.);
   RooRealVar* zmm_dcb_sigma  = new RooRealVar("zmm_dcb_sigma_"+varname, "zmm_dcb_sigma_"+varname , 4.48,  3., 10.);
   RooRealVar* zmm_dcb_a1 = new RooRealVar("zmm_dcb_a1_"+varname, "zmm_dcb_a1_"+varname, 1.22, 0.5,  5.);
   RooRealVar* zmm_dcb_a2 = new RooRealVar("zmm_dcb_a2_"+varname, "zmm_dcb_a2_"+varname, 1.78, 0.5,  5.);
   RooRealVar* zmm_dcb_n1  = new RooRealVar("zmm_dcb_n1_"+varname, "zmm_dcb_n1_"+varname , 1.5, 0., 10.);
   RooRealVar* zmm_dcb_n2  = new RooRealVar("zmm_dcb_n2_"+varname , "zmm_dcb_n2_"+varname , 9.14, 0., 10.);

   RooDoubleCrystalBall zmm_dcb_pdf("zmm_dcb_pdf", "zmm_dcb_pdf", dilep_mass, *zmm_dcb_mean, *zmm_dcb_sigma, *zmm_dcb_a1, *zmm_dcb_n1, *zmm_dcb_a2, *zmm_dcb_n2);
   RooRealVar zmm_dcb_pdf_norm("zmm_dcb_pdf_norm" , "zmm_dcb_pdf_norm" , 100, 0, 10000000.);
   RooAddPdf zmm_dcb_epdf("zmm_dcb_epdf","zmm_dcb_epdf",RooArgList(zmm_dcb_pdf),RooArgList(zmm_dcb_pdf_norm));

   RooFitResult * zmm_fit_result = zmm_dcb_epdf.fitTo(*dhist_zmm, RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range , max_fit_range),RooFit::PrintLevel(-1));

   print_details (zmm_fit_result);
   auto zmm_frame = dilep_mass.frame();
   dhist_zmm->plotOn(zmm_frame,RooFit::Binning(nbin_data),RooFit::Name("data"));
   zmm_dcb_pdf.plotOn(zmm_frame,RooFit::Name("zmm_dcb_pdf"));
   save_plot(zmm_frame,"m(e,#mu)","mk2zmm_v1_dcb_"+name,new TLegend(),NULL, false, false);
   save_plot_and_band(zmm_frame,dilep_mass,{"zmm_dcb_pdf"},"m(e,#mu)","mk2zmm_v1_dcb_band_"+name);

   /// fix parameters based on fit
   if (constrain_type_mean == "fixed")
      zmm_dcb_mean->setConstant(true);

   if (constrain_type_mean == "linear")
      zmm_dcb_mean->setRange(zmm_dcb_mean->getVal()-zmm_dcb_mean->getError(),zmm_dcb_mean->getVal()+zmm_dcb_mean->getError());

   // mean formula for gaussian
   RooRealVar *zmm_scale = new RooRealVar("nuisance_zmm_scale_"+varname,"",0,-1,1);

   RooFormulaVar *zmm_mean_formula=new RooFormulaVar("zmm_mean_formula_"+varname, "", "@0*(1+"+TString(std::to_string(zmm_dcb_mean->getError()))+"*@1)", RooArgList(*zmm_dcb_mean,*zmm_scale));

   zmm_dcb_sigma->setConstant(true);
   zmm_dcb_a1->setConstant(true);
   zmm_dcb_a2->setConstant(true);
   zmm_dcb_n1->setConstant(true);
   zmm_dcb_n2->setConstant(true);

   if (create_dc_input){
     RooWorkspace *wspace = new RooWorkspace("workspace_zmm","workspace_zmm");

     RooDoubleCrystalBall *zmm_dcb_pdf_out;
     if ( constrain_type_mean == "gaussian")
       zmm_dcb_pdf_out = new RooDoubleCrystalBall("zmm_dcb_pdf_"+varname, "zmm_dcb_pdf_"+varname, dilep_mass_out, *zmm_mean_formula, *zmm_dcb_sigma, *zmm_dcb_a1, *zmm_dcb_n1, *zmm_dcb_a2, *zmm_dcb_n2);
     else
       zmm_dcb_pdf_out = new RooDoubleCrystalBall("zmm_dcb_pdf_"+varname, "zmm_dcb_pdf_"+varname, dilep_mass_out, *zmm_dcb_mean, *zmm_dcb_sigma, *zmm_dcb_a1, *zmm_dcb_n1, *zmm_dcb_a2, *zmm_dcb_n2);

     RooRealVar zmm_dcb_pdf_out_norm("zmm_dcb_pdf_"+varname+"_norm" , "zmm_dcb_pdf_"+varname+"_norm" , zmm_dcb_pdf_norm.getVal(), 0, 10000000.);
     zmm_dcb_pdf_out_norm.setConstant(true);
     wspace->import(*zmm_dcb_pdf_out);
     wspace->import(zmm_dcb_pdf_out_norm);
     wspace->writeToFile("workspace_mk2zmm_v1_"+name+".root");

   }

 return 0;
}
