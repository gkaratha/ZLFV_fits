#include <fstream>
#include "fit_helper.h"
#include "fit_helper_core.h"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Main code ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////// Z prime v2:
//////////////////// +based on v15
////////////////////

//--------------------------------------------------------------------------------------------
// Model the Z prime signal with a Crystal Ball or Gaussian distribution
int ScanTau_fit_sig(TString name="bin1_r2",
                    float min_mass=100, float max_mass=500,
                    vector<float> signal_parameters = {}, //Signal shape parameters
                    float expected_Nsgn=65, bool create_dc_input=true, TString outvar="mcol",
                    bool syst_sgn=false, TString varname="bin", TString outdir="tau_workspaces/"){

  if(signal_parameters.size() != 6 && signal_parameters.size() != 2) {
    printf("%s: Error! Parameter vector does not have the correct size\n", __func__);
    return -1;
  }

  const bool mutau = name.Contains("mutau");

  // Setup output directories
  gROOT->SetBatch(true);
  if(outdir == "") outdir = "tau_workspaces/";
  if(!outdir.EndsWith("/")) outdir += "/";
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", outdir.Data(), outdir.Data()));
  figdir_ = "figures/" + name + "/";
  figdir_.ReplaceAll("_mp", "/mp"); //put mass point fits within the same sub-directory of the base processing name
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir_.Data(), figdir_.Data()));

  // Create the signal model
  cout<<" \n ********************************************************** "<<endl;
  cout<<" ********************************************************** "<<endl;
  cout<<" ***ScanTau fit: signal part ***Output name:"+name<<endl;
  cout<<" ********************************************************** "<<endl;

  RooRealVar dilep_mass_out(outvar,"m(l,#tau)", (max_mass+min_mass)/2., min_mass , max_mass, "GeV/c^{2}");

  const bool use_gaus = signal_parameters.size() == 2;
  if(use_gaus) { //push dummy values for the rest
    for(unsigned i = signal_parameters.size(); i < 6; ++i) signal_parameters.push_back(0.);
  }
  RooRealVar mean  ("mean_"  +varname, "mean"  , signal_parameters[0]); mean  .setConstant(true);
  RooRealVar sigma ("sigma_" +varname, "sigma" , signal_parameters[1]); sigma .setConstant(true);
  RooRealVar alpha1("alpha1_"+varname, "alpha1", signal_parameters[2]); alpha1.setConstant(true);
  RooRealVar alpha2("alpha2_"+varname, "alpha2", signal_parameters[3]); alpha2.setConstant(true);
  RooRealVar enne1 ("enne1_" +varname, "enne1" , signal_parameters[4]); enne1 .setConstant(true);
  RooRealVar enne2 ("enne2_" +varname, "enne2" , signal_parameters[5]); enne2 .setConstant(true);

  const bool use_energy_scale = true; //include energy scale uncertainties
  RooRealVar lep_ES_shift("lep_ES_shift", "lepton ES shift", 0., -5., 5.); lep_ES_shift.setConstant(true); //correlated between categories
  RooRealVar tau_ES_shift("tau_ES_shift", "Tau ES shift"    , 0., -5., 5.); tau_ES_shift.setConstant(true);
  //Define the electron scale uncertainty as 0.5%, muon as 0.3%, and tau as 1%, assume the effect ~sqrt(sys/nom)
  // --> mean shift = 0.25% for electron ES, 0.15% for muon ES, and 0.5% for tau ES
  RooFormulaVar mean_func("mean_func_"+varname, "mean with offset",
                          (mutau) ? "@0*(1 + @1*0.0015 + @2*0.0050)" : "@0*(1 + @1*0.0025 + @2*0.0050)",
                          RooArgList(mean, lep_ES_shift, tau_ES_shift));

  RooAbsPdf * signal_pdf(nullptr);
  if(use_gaus) signal_pdf = new RooGaussian("signal_pdf_"+varname,"Gaussian PDF",dilep_mass_out,mean_func,sigma);
  else signal_pdf = new RooDoubleCrystalBall("signal_pdf_"+varname, "Crystal Ball PDF", dilep_mass_out, mean_func, sigma, alpha1, enne1, alpha2, enne2);

  RooRealVar n_sgn("signal_pdf_"+varname+"_norm","",expected_Nsgn,0,100000);
  n_sgn.setConstant(true);

  auto sgn_frame = dilep_mass_out.frame(RooFit::Title(""));
  signal_pdf->plotOn(sgn_frame,RooFit::LineColor(kBlue),RooFit::Normalization(n_sgn.getVal(), RooAbsReal::NumEvent));
  save_plot(sgn_frame,"m(l,#tau)","sig_"+name);

  //free the energy scale uncertainties if being included
  if(use_energy_scale) {
    lep_ES_shift.setConstant(false);
    tau_ES_shift.setConstant(false);
  }

  if (create_dc_input){
    RooWorkspace *wspace_sgn = new RooWorkspace("ws_sgn","ws_sgn");
    wspace_sgn->import(*signal_pdf);
    wspace_sgn->import(n_sgn);
    wspace_sgn->writeToFile(outdir+"workspace_sig_"+name+".root");
  }

  return 0;
}
