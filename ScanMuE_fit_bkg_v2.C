#include <fstream>
#include "fit_helper.h"
#include "fit_helper_core.h"

bool force_standard_env_ = false; //force a specific background envelope with a specific asimov template

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// Z prime v2:
//////////////////// +based on v15
//////////////////// +divided in signal and bkg part: here BKG



int ScanMuE_fit_bkg_v2(TString name="bin1_r2",
                       TString data_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Scan_Zprime/Meas_full_bdt_v7_emu_scan_data_Run1*.root",
                       TString xgbmin="0.7",TString xgbmax="1.01", double  min_fit_range=90,
                       double  max_fit_range = 135, double blind_min=106, double blind_max=114,
                       bool unblind=false, bool create_dc_input=false, TString outvar="mass_ll",
                       TString varname="bin",TString outdir="WorkspaceScanBKG/"){

  //////////////////////////////////// configuration /////////////////////////
  gROOT->SetBatch(true);
  if(outdir == "") outdir = "WorkspaceScanBKG/";
  if(!outdir.EndsWith("/")) outdir += "/";
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", outdir.Data(), outdir.Data()));
  const bool mmackenz = TString(gSystem->Getenv("USER")).Contains("macken"); //preferred structure for Michael
  if(mmackenz) {
    figdir_ = "figures/" + name + "/";
    figdir_.ReplaceAll("_mp", "/mp"); //put mass point fits within the same sub-directory of the base processing name
  }
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir_.Data(), figdir_.Data()));

  /// bkg functions
  bool Fit_cheb=true;
  int min_cheb_order=1,max_cheb_order=(max_fit_range > 500. && min_fit_range < 150.) ? 6 : 4;
  int add_orders_cheb=0;
  bool Fit_sumexp=true;
  int min_sumexp_order=1,max_sumexp_order=4;
  int add_orders_sumexp=0;
  bool sumexp_recurse_coef=true;
  bool Fit_sumplaw=true;
  int min_sumplaw_order=1,max_sumplaw_order=4;
  int add_orders_sumplaw=0;
  bool sumplaw_recurse_coef=true;

  /// generic
  TString cuts = xgbmin+"<xgb && xgb<"+xgbmax+" && Flag_met && Flag_muon && Flag_electron"; // not add mass_ll here when run systematics
  TString tree_name="mytreefit";
  //bin data such that a single bin is ~1/2 of the signal width --> +-1 sigma = 4 bins
  int nbin_data = (max_fit_range-min_fit_range)/(0.5*(blind_max-blind_min)/2.); //(blind_max - blind_min) = 2*width
  bool Verbose=false; // RooFit verbosity
  int printout_levels=1; // 0: Print only final fits parameters, 1: Print all fits from F test tests


  ///////////////////////////////////////////////////////////////////////////
  cout<<" \n ********************************************************** "<<endl;
  cout<<" ********************************************************** "<<endl;
  cout<<" ***ScanMuE fit v2: bkg part ***Output name:"+name<<endl;
  cout<<" ********************************************************** "<<endl;

  if (!Verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);


  int nbin_blind =0;
  if ( !unblind )
    nbin_blind = (blind_max-blind_min) /( max_fit_range - min_fit_range )*nbin_data;

  // read trees
  RooRealVar dilep_mass ("mass_ll","m(e,#mu)", min_fit_range, min_fit_range , max_fit_range, "GeV/c^{2}");
  RooRealVar dilep_mass_out(outvar,"m(e,#mu)", min_fit_range, min_fit_range , max_fit_range, "GeV/c^{2}");

  dilep_mass.setRange("left",min_fit_range, blind_min);
  dilep_mass.setRange("right",blind_max, max_fit_range);
  dilep_mass.setRange("sr",blind_min,blind_max);
  dilep_mass.setRange("full",min_fit_range, max_fit_range);


  ///////////////////////////////// Data/bkg fit ////////////////////////////
  cout<<"\n *********************** Data fit ********************** "<<endl;

  dilep_mass.setBins(nbin_data);

  TH1F * hbkg = new TH1F("hbkg","",nbin_data,min_fit_range,max_fit_range);

  TChain * cc = new TChain("mytreefit");
  cc->Add(data_file);
  cc->Draw("mass_ll>>hbkg",cuts);

  RooDataHist * dhist_bkg = new RooDataHist("dhist_bkg","dhist_bkg",RooArgSet(dilep_mass),hbkg);

  RooDataHist * dhist_data = new RooDataHist("dhist_data","dhist_data",RooArgSet(dilep_mass_out),hbkg);

  cout<<"dataset = "<<dhist_bkg->sumEntries()<<endl;


  /////// chebychev
  std::vector<RooAbsPdf*> bkg_cheb_pdfs;
  std::vector<RooRealVar*> bkg_cheb_ampl;
  std::vector<TString> bkg_cheb_names;
  std::vector<TString> bkg_cheb_legs;
  std::vector<int> bkg_cheb_orders;

  for (int iorder=min_cheb_order; iorder<max_cheb_order+1; iorder++){
    TString sorder(std::to_string(iorder));
    bkg_cheb_pdfs.push_back(CreateChebychev( "cheb"+sorder+"_bkgPDF", iorder, dilep_mass));
    bkg_cheb_ampl.push_back(new RooRealVar("cheb"+sorder+"_bkgPDF_norm","",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
    bkg_cheb_names.push_back("cheb"+sorder+"_"+name);
    bkg_cheb_legs.push_back("Chebychev "+sorder);
    bkg_cheb_orders.push_back(iorder);
  }

  ///// ftest
  TString cfg_tag = (mmackenz) ? "" : outdir+"Ftest_figs/cfit";
  FtestStruct cheb_Ftest;
  cheb_Ftest.success=false;
  const double min_p_value = 0.001;
  const double ftest_step = 0.05; //information gain requirement
  const bool   force_inclusion = true; //force at least one order of each family to be included in the envelope, independent of p-value
  if (Fit_cheb) {
    cout<<" ************************ Chebychev begin ************************ "<<endl;
    cheb_Ftest =  HistFtest(bkg_cheb_pdfs, bkg_cheb_ampl,  dhist_bkg, dilep_mass, bkg_cheb_orders,
                            bkg_cheb_names, bkg_cheb_legs, nbin_data,"scanbkg_v2_cheb_"+name,
                            ftest_step, min_p_value,printout_levels, force_inclusion, force_standard_env_, cfg_tag);
    cout<<" ************************ Chebychev end ************************ "<<endl;
  }



  /////// sum exp
  std::vector<RooAbsPdf*> bkg_sumexp_pdfs;
  std::vector<RooRealVar*> bkg_sumexp_ampl;
  std::vector<TString> bkg_sumexp_names;
  std::vector<TString> bkg_sumexp_legs;
  std::vector<int> bkg_sumexp_orders;
  for (int iorder=min_sumexp_order; iorder<max_sumexp_order+1; iorder++){
    TString sorder(std::to_string(iorder));
    bkg_sumexp_pdfs.push_back(CreateSumExpo( "sumexp"+sorder+"_bkgPDF", iorder, dilep_mass, sumexp_recurse_coef));
    bkg_sumexp_ampl.push_back(new RooRealVar("sumexp"+sorder+"_bkgPDF_norm","",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
    bkg_sumexp_names.push_back("sumexp"+sorder+"_"+name);
    bkg_sumexp_legs.push_back("#Sigma Expo "+sorder);
    bkg_sumexp_orders.push_back(iorder);
  }

  FtestStruct sumexp_Ftest;
  sumexp_Ftest.success=false;
  if (Fit_sumexp) {
    ////// ftest
    cout<<" ************************ Ftest SumExp begin ************************ "<<endl;
    sumexp_Ftest =  HistFtest(bkg_sumexp_pdfs, bkg_sumexp_ampl,  dhist_bkg, dilep_mass, bkg_sumexp_orders,
                              bkg_sumexp_names, bkg_sumexp_legs, nbin_data, "scanbkg_v2_exp_"+name,
                              ftest_step, min_p_value,printout_levels, force_inclusion, force_standard_env__, cfg_tag);
    cout<<" ************************ Ftest SumExp end ************************ "<<endl;
  }


  /////// sum power law
  std::vector<RooAbsPdf*> bkg_sumplaw_pdfs;
  std::vector<RooRealVar*> bkg_sumplaw_ampl;
  std::vector<TString> bkg_sumplaw_names;
  std::vector<TString> bkg_sumplaw_legs;
  std::vector<int> bkg_sumplaw_orders;
  for (int iorder=min_sumplaw_order; iorder<max_sumplaw_order+1; iorder++){
    TString sorder(std::to_string(iorder));
    bkg_sumplaw_pdfs.push_back(CreateSumPower( "sumplaw"+sorder+"_bkgPDF", iorder, dilep_mass,sumplaw_recurse_coef));
    bkg_sumplaw_ampl.push_back(new RooRealVar("sumplaw"+sorder+"_bkgPDF_norm","",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
    bkg_sumplaw_names.push_back("sumplaw"+sorder+"_"+name);
    bkg_sumplaw_legs.push_back("#Sigma power law "+sorder);
    bkg_sumplaw_orders.push_back(iorder);
  }

  ////// ftest
  FtestStruct sumplaw_Ftest;
  sumplaw_Ftest.success=false;
  if (Fit_sumplaw){
    cout<<" ************************ Ftest Sum Power Law begin ************************ "<<endl;
    sumplaw_Ftest =  HistFtest(bkg_sumplaw_pdfs, bkg_sumplaw_ampl,  dhist_bkg, dilep_mass, bkg_sumplaw_orders,
                               bkg_sumplaw_names, bkg_sumplaw_legs, nbin_data, "scanbkg_v2_sumplaw_"+name,
                               ftest_step, min_p_value,printout_levels, force_inclusion, force_standard_env__, cfg_tag);
    cout<<" ************************ Ftest Sum Power Law end ************************ "<<endl;
  }



  //////////////////// fit all candidate functions ////////////////

  // For each accepted function, create a new PDF and fit this to the data sidebands
  std::vector<RooAbsPdf*> bkg_pdfs; //background PDF
  std::vector<RooRealVar*> bkg_ampl; //normalization for each PDF
  std::vector<TString> bkg_names; //name to write to the workspace
  std::vector<TString> bkg_legs; //title for the legend

  if (cheb_Ftest.success) {
    for (int i=0; i<cheb_Ftest.getAllOrder.size(); i++){
      int iord =cheb_Ftest.getAllOrder[i];
      TString sord(std::to_string(iord));
      bkg_pdfs.push_back(CreateChebychev( varname+"_cheb"+sord+"_pdf", iord, dilep_mass));
      bkg_ampl.push_back(new RooRealVar(varname+"_cheb"+sord+"_pdf_norm", varname+"_cheb"+sord+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
      bkg_names.push_back(varname+"_cheb"+sord+"_pdf");
      bkg_legs.push_back("Chebychev "+sord);
    }
  }

  if (sumexp_Ftest.success) {
    for (int i=0; i<sumexp_Ftest.getAllOrder.size(); i++){
      int iord =sumexp_Ftest.getAllOrder[i];
      TString sord(std::to_string(iord));
      bkg_pdfs.push_back( CreateSumExpo( "sumexp"+sord+"_pdf", iord, dilep_mass, sumexp_recurse_coef) );
      bkg_ampl.push_back(new RooRealVar(varname+"_sumexp"+sord+"_pdf_norm", varname+"_sumexp"+sord+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
      bkg_names.push_back(varname+"_sumexp"+sord+"_pdf");
      bkg_legs.push_back("SumExp "+sord);
    }
  }

  if (sumplaw_Ftest.success) {
    for (int i=0; i<sumplaw_Ftest.getAllOrder.size(); i++){
      int iord =sumplaw_Ftest.getAllOrder[i];
      TString sord(std::to_string(iord));
      bkg_pdfs.push_back( CreateSumPower( "sumplaw"+sord+"_pdf", iord, dilep_mass,sumplaw_recurse_coef) );
      bkg_ampl.push_back(new RooRealVar(varname+"_sumplaw"+sord+"_pdf_norm", varname+"_sumplaw"+sord+"_pdf_norm",dhist_bkg->sumEntries(),0,2*dhist_bkg->sumEntries()));
      bkg_names.push_back(varname+"_sumplaw"+sord+"_pdf");
      bkg_legs.push_back("SumPower "+sord);
    }
  }


  //////////// fit here
  cout<<" Fit of best variables "<<endl;
  TString fig_tag = "cfit";
  std::vector<std::vector<float>> final_results =  FitHistBkgFunctions(bkg_pdfs, bkg_ampl, dhist_bkg, dilep_mass,
                                                                       bkg_names, bkg_legs, nbin_data, unblind,
                                                                       "scanbkg_v2_best_"+name, true,fig_tag,true,unblind);

  if (!create_dc_input)
    return 0;
  //////////////////// create bkg workspace ////////////////////////

  RooCategory cat("pdfindex_"+varname, "");
  RooArgList models_out; //container for multpdf
  RooWorkspace *wspace = new RooWorkspace("ws_bkg","ws_bkg"); //background & data workspace
  int best_index = 0; //track the best fit result, set this as the default
  float best_pvalue = 0.;

  // Create new PDFs initialized to the fit parameter results returned by FitHistBkgFunctions
  // FIXME: Why re-create the PDFs?
  int nstartFNC1=0;
  if (cheb_Ftest.success) {
    for (int i=0; i<cheb_Ftest.getAllOrder.size(); i++){
      std::vector<float> param;
      TString sord (std::to_string(cheb_Ftest.getAllOrder[i]));
      for(int j=2; j<final_results[i].size(); j++) //skip indices 0 and 1 which contain chi^2 and N(DOF), respectively
        param.push_back(final_results[i][j]);
      auto pdf = CreateChebychev( "bkg_cheb"+sord+"_pdf_"+varname, cheb_Ftest.getAllOrder[i], dilep_mass_out);
      RooArgList param_list(*(pdf->getVariables()));
      for(int index = 0; index < param_list.getSize()-1; ++index) {
        ((RooRealVar*) param_list.at(index))->setVal(param[index]);
      }
      models_out.add(*pdf);
      pdf->Print();
      const float pvalue = ROOT::Math::chisquared_cdf_c(final_results[nstartFNC1][0],final_results[nstartFNC1][1]);
      if(pvalue > best_pvalue) {
        best_index  = nstartFNC1;
        best_pvalue = pvalue;
      }
      nstartFNC1+=1;
    }
    for( int iord =cheb_Ftest.getBestOrder+1; iord<cheb_Ftest.getBestOrder+add_orders_cheb+1; iord++){
      TString sord (std::to_string(iord));
      models_out.add( *( CreateChebychev( "bkg_cheb"+sord+"_pdf_"+varname, iord, dilep_mass_out) ));
    }
  }

  int nstartFNC2=nstartFNC1;
  if (sumexp_Ftest.success) {
    for (int i=0; i<sumexp_Ftest.getAllOrder.size(); i++){
      std::vector<float> param;
      for(int j=2; j<final_results[i+nstartFNC1].size(); j++)
        param.push_back(final_results[i+nstartFNC1][j]);
      TString sord (std::to_string(sumexp_Ftest.getAllOrder[i]));
      auto pdf = CreateSumExpo( "bkg_sumexp"+sord+"_pdf_"+varname,  sumexp_Ftest.getAllOrder[i], dilep_mass_out, sumexp_recurse_coef);
      RooArgList param_list(*(pdf->getVariables()));
      for(int index = 0; index < param_list.getSize()-1; ++index) {
        ((RooRealVar*) param_list.at(index))->setVal(param[index]);
      }
      models_out.add(*pdf);
      const float pvalue = ROOT::Math::chisquared_cdf_c(final_results[nstartFNC2][0],final_results[nstartFNC2][1]);
      if(pvalue > best_pvalue) {
        best_index = nstartFNC2;
        best_pvalue = pvalue;
      }
      nstartFNC2+=1;
    }
    for( int iord =sumexp_Ftest.getBestOrder+1; iord<sumexp_Ftest.getBestOrder+add_orders_sumexp+1; iord++){
      TString sord (std::to_string(iord));
      models_out.add( *CreateSumExpo( "bkg_sumexp"+sord+"_pdf_"+varname, iord, dilep_mass_out, sumexp_recurse_coef) );
    }
  }

  int nstartFNC3=nstartFNC2;
  if (sumplaw_Ftest.success) {
    for (int i=0; i<sumplaw_Ftest.getAllOrder.size(); i++){
      std::vector<float> param;
      for(int j=2; j<final_results[i+nstartFNC2].size(); j++)
        param.push_back(final_results[i+nstartFNC2][j]);
      TString sord (std::to_string(sumplaw_Ftest.getAllOrder[i]));
      auto pdf = (CreateSumPower( "bkg_sumplaw"+sord+"_pdf_"+varname, sumplaw_Ftest.getAllOrder[i], dilep_mass_out,sumplaw_recurse_coef));
      RooArgList param_list(*(pdf->getVariables()));
      for(int index = 0; index < param_list.getSize()-1; ++index) {
        ((RooRealVar*) param_list.at(index))->setVal(param[index]);
      }
      models_out.add(*pdf);
      const float pvalue = ROOT::Math::chisquared_cdf_c(final_results[nstartFNC3][0],final_results[nstartFNC3][1]);
      if(pvalue > best_pvalue) {
        best_index = nstartFNC3;
        best_pvalue = pvalue;
      }
      nstartFNC3+=1;
    }
    for( int iord =sumplaw_Ftest.getBestOrder+1; iord<sumplaw_Ftest.getBestOrder+add_orders_sumplaw+1; iord++){
      TString sord (std::to_string(iord));
      models_out.add( *(CreateSumPower( "bkg_sumplaw"+sord+"_pdf_"+varname, iord, dilep_mass_out,sumplaw_recurse_coef)) );
    }
  }

  //Create a RooMultiPdf envelope with all of the accepted functions and the associated function index
  RooMultiPdf multipdf("multipdf_"+varname, "", cat, models_out);
  RooRealVar norm_out("multipdf_"+varname+"_norm","",dhist_bkg->sumEntries(),0.5*dhist_bkg->sumEntries(),10*dhist_bkg->sumEntries());
  if(force_standard_env_) cat.setIndex(1); //use the 2nd order Chebychev for the Asimov template if using fixed envelopes
  else {
    cat.setIndex(best_index); //initialize to the best fit result
    cout << "Using best index " << best_index << " with p-value " << best_pvalue << " (" << models_out.at(best_index)->GetName() << ")\n";
  }

  models_out.Print();
  wspace->import(cat);
  wspace->import(multipdf);
  wspace->import(norm_out);
  dhist_data->SetName("data_obs");
  wspace->import(*dhist_data); //import the data as well
  wspace->writeToFile(outdir+"workspace_scanbkg_v2_"+name+".root"); // write the output



  return 0;
}
