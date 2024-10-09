#include <fstream>
#include "fit_helper.h"
#include "fit_helper_core.h"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// mk2: seperating fit in components
///////////// Signal component v1
//////////////////// + Same thing as v15 etc



int ZMuE_fit_mk2_sgn_v1(TString name="bin1_r2",
    TString sgn_file="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSF_bdt_v7_signal_mcRun1*.root",
    TString xgbmin="0.3",TString xgbmax="0.7",  bool create_dc_input=false,
    float expected_Nsgn=0, TString outvar="mass_ll", bool syst_sgn=false,
    TString varname="bin"){

   //////////////////////////////////// configuration /////////////////////////
   gROOT->SetBatch(true);
   TString sgn_file_path = "/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/";
   TString sgn_mu_up_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_MU_SCALEUP_bdt_v7_signal_mcRun1*.root";
   TString sgn_mu_down_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_MU_SCALEDOWN_bdt_v7_signal_mcRun1*.root";
   TString sgn_ele_up_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_ELE_SCALEUP_bdt_v7_signal_mcRun1*.root";
   TString sgn_ele_down_file="../BDT/Systematics_v2_no_met_significance/Meas_fullAndSF_syst_ELE_SCALEDOWN_bdt_v7_signal_mcRun1*.root";
   TString dilep_var_name="mass_ll";
   TString dilep_muup_name="mass_ll_Muon_scale_up";
   TString dilep_mudown_name="mass_ll_Muon_scale_down";
   TString dilep_eleup_name="mass_ll_Electron_scale_up";
   TString dilep_eledown_name="mass_ll_Electron_scale_down";
   std::vector<TString> syst_names{"Muon_scale_up", "Muon_scale_down", "Electron_scale_up", "Electron_scale_down"};
   std::vector<TString> syst_files{sgn_mu_up_file, sgn_mu_down_file, sgn_ele_up_file, sgn_ele_down_file};

   ///////////////////////////////////////////////////////////////////////////
   cout<<" \n ********************************************************** "<<endl;
   cout<<" ********************************************************** "<<endl;
   cout<<" ***ZMuE mk2 fit v1: signal part ***Output name:"+name<<endl;
   cout<<" ********************************************************** "<<endl;


   TString cuts = xgbmin+"<xgb && xgb<"+xgbmax+" && Flag_met && Flag_muon && Flag_electron"; // not add mass_ll here when run systematics
   TString tree_name="mytreefit";

   int nbin_data=80;
   double min_fit_range=70;
   double max_fit_range=110;
   bool Verbose = false;

   if (!Verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);


   // read trees
   TTree * sgn_tree = get_tree("mytreefit",sgn_file,cuts);
   RooRealVar dilep_mass("mass_ll","m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");
   RooDataSet sgn_dataset("sgn_dataset","sgn_dataset",RooArgSet(dilep_mass),RooFit::Import(*sgn_tree));
   RooRealVar dilep_mass_out(outvar,"m(e,#mu)", 110., min_fit_range , max_fit_range, "GeV/c^{2}");



   /////////////////////////////// Signal Fit ////////////////////////////////
   cout<<"\n ************** Signal fit: ************** "<<endl;
   cout<<" Signal MC events: "<< sgn_tree->GetEntries()<<" norm to "<<expected_Nsgn<<endl;

   RooWorkspace ws ("ws");
   ws.import(dilep_mass);

   ///// parameters
   RooRealVar sgn_cb_mean("sgn_cb_mean_"+varname,"",91.0e+00, 84.0e+00, 100.0e+00);
   RooRealVar sgn_cb_width("sgn_cb_width_"+varname,"",5., 0.1, 10.);
   RooRealVar sgn_cb_a1("sgn_cb_a1_"+varname,"",1.0, 0.1, 100.0);
   RooRealVar sgn_cb_n1("sgn_cb_n1_"+varname,"",1.0, 0.00001, 100.0);
   RooRealVar sgn_cb_a2("sgn_cb_a2_"+varname,"",1.0, 0.1, 100.0);
   RooRealVar sgn_cb_n2("sgn_cb_n2_"+varname,"",1.0, 0.1, 100.0);

   ///// function
   RooDoubleCB sgn_PDF("sgn_PDF","cb",dilep_mass,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);

   ///// extend
   RooRealVar n_sgn("sgn_PDF_norm", "",sgn_dataset.sumEntries(),-10*sgn_dataset.sumEntries(),10*sgn_dataset.sumEntries());
   RooAddPdf esgn_PDF("sgn_epdf","esgnPDF",RooArgList(sgn_PDF),RooArgList(n_sgn));

   ///// fit
   RooFitResult * sgn_result = esgn_PDF.fitTo(sgn_dataset,RooFit::Extended(1),RooFit::Save(),RooFit::Range(min_fit_range, max_fit_range),RooFit::PrintLevel(-1));
   ///// fix shape
   sgn_cb_mean.setConstant(true);   sgn_cb_width.setConstant(true);
   sgn_cb_a1.setConstant(true);     sgn_cb_n1.setConstant(true);
   sgn_cb_a2.setConstant(true);     sgn_cb_n2.setConstant(true);
   RooDoubleCB sgn_PDF_out("sgn_PDF_"+varname,"cb",dilep_mass_out,sgn_cb_mean,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);

   auto sgn_frame = dilep_mass.frame();
   print_details(sgn_result);
   int n_param_sgn = sgn_result->floatParsFinal().getSize();
   float chi2_sgn = get_chi_squared(dilep_mass, &esgn_PDF, sgn_dataset, true, nbin_data, n_param_sgn);
   sgn_dataset.plotOn(sgn_frame,RooFit::Binning(nbin_data),RooFit::Name("data"));
   esgn_PDF.plotOn(sgn_frame,RooFit::LineColor(kBlue),RooFit::Name("sgn_epdf"));

   n_sgn.setVal(expected_Nsgn);

   save_plot(sgn_frame,"m(#mu,e)","mk2sgn_v1_dcb_"+name);
   save_plot_and_band(sgn_frame,dilep_mass,{"sgn_epdf"},"m(e,#mu)","mk2sgn_v1_dcb_band_"+name);

   cout<<" Result  whole range    |    85-95 only"<<endl;
   std::pair<double,double> nSgn = yield_calc( n_sgn.getVal(), dilep_mass, &sgn_PDF);
   cout<<" - Expected nSgn "<<nSgn.first<<"   |  "<<nSgn.second<<endl;
   cout<<" ************************************************* "<<endl;

   /////// signal systematics
   double mean_shift_mu_up=0,width_shift_mu_up=0,mean_shift_mu_down=0,width_shift_mu_down=0;
   double mean_shift_ele_up=0,width_shift_ele_up=0,mean_shift_ele_down=0,width_shift_ele_down=0;

   if (syst_sgn){
     for (int isyst=0; isyst<syst_files.size(); isyst++){
        // FIXME: Put the systematic ntuples in a publicly viewable place
        if(false) {
          std::pair<double,double> max_mean_width = SignalSystematicsMaxMeanWidth(sgn_file_path+syst_files[isyst], cuts, syst_names[isyst], min_fit_range, max_fit_range, nbin_data, name, "mk2sgn_v1") ;
          cout<<"signal systematic "+syst_names[isyst]<<" mean "<<max_mean_width.first<<" width "<<max_mean_width.second<<endl;
          if (isyst==0){ mean_shift_mu_up=max_mean_width.first; width_shift_mu_up=max_mean_width.second;}
          if (isyst==1){ mean_shift_mu_down=max_mean_width.first; width_shift_mu_down=max_mean_width.second;}
          if (isyst==2){ mean_shift_ele_up=max_mean_width.first; width_shift_ele_up=max_mean_width.second;}
          if (isyst==3){ mean_shift_ele_down=max_mean_width.first; width_shift_ele_down=max_mean_width.second;}
        } else {
          if (isyst==0){ mean_shift_mu_up=1.001*sgn_cb_mean.getVal(); width_shift_mu_up=1.001*sgn_cb_width.getVal();}
          if (isyst==1){ mean_shift_mu_down=0.999*sgn_cb_mean.getVal(); width_shift_mu_down=0.999*sgn_cb_width.getVal();}
          if (isyst==2){ mean_shift_ele_up=1.002*sgn_cb_mean.getVal(); width_shift_ele_up=1.001*sgn_cb_width.getVal();}
          if (isyst==3){ mean_shift_ele_down=0.992*sgn_cb_mean.getVal(); width_shift_ele_down=0.999*sgn_cb_width.getVal();}
        }
     }
   }
   /////// create signal workspace
   if (create_dc_input){
      RooRealVar n_sgn_out("sgn_PDF_"+varname+"_norm", "",expected_Nsgn,-10*expected_Nsgn,10*expected_Nsgn);
      n_sgn_out.setConstant(true);

      RooWorkspace *wspace_sgn = new RooWorkspace("workspace_signal","workspace_signal");
      wspace_sgn->import(sgn_PDF_out);
      wspace_sgn->import(n_sgn_out);

      double max_shift_muon=0;
      double max_shift_electron=0;
      if (syst_sgn){
        max_shift_muon= (fabs(sgn_cb_mean.getVal()-mean_shift_mu_down) + fabs(sgn_cb_mean.getVal()-mean_shift_mu_up) )/(2*sgn_cb_mean.getVal());
        max_shift_electron= fabs(mean_shift_ele_down-mean_shift_ele_up)/(2*sgn_cb_mean.getVal());
        cout<<"delta/M muon "<<max_shift_muon<<" electron "<<max_shift_electron<<endl;
      }
      RooRealVar *mu_scale = new RooRealVar("nuisance_mu_scale_"+varname,"",0,-1,1);
      RooRealVar *ele_scale = new RooRealVar("nuisance_ele_scale_"+varname,"",0,-1,1);
      mu_scale->setConstant(true);
      ele_scale->setConstant(true);
      RooFormulaVar *mean_formula=new RooFormulaVar("mean_formula_"+varname, "", "@0*(1+"+TString(std::to_string(max_shift_muon))+"*@1)*(1+"+TString(std::to_string(max_shift_muon))+"*3*@2)", RooArgList(sgn_cb_mean,*mu_scale,*ele_scale));
      RooDoubleCB syst_sgn_PDF_out ("syst_sgn_PDF_"+varname,"cb",dilep_mass_out,*mean_formula,sgn_cb_width,sgn_cb_a1,sgn_cb_n1,sgn_cb_a2,sgn_cb_n2);
      RooRealVar n_syst_sgn_out("syst_sgn_PDF_"+varname+"_norm", "",expected_Nsgn,-10*expected_Nsgn,10*expected_Nsgn);
      n_syst_sgn_out.setConstant(true);
      if (syst_sgn){
        wspace_sgn->import(syst_sgn_PDF_out);
        wspace_sgn->import(n_syst_sgn_out);
      }
      wspace_sgn->writeToFile("workspace_mk2sgn_v1_"+name+".root");
   }

 return 0;
}
