#ifndef FIT_HELPER_H
#define FIT_HELPER_H

#include <fstream>
#include <random>

TString figdir_ = "./"; //Figure directory

/////// read ttree
TTree * get_tree(TString tree_name, TString path, TString cuts="",Long64_t nmax=TTree::kMaxEntries){
    TChain * tree = new TChain(tree_name);
    tree->Add(path);
    if (cuts=="")
       cuts="1>0";
    if(tree->GetEntries() == 0) {
      cout << __func__ << ": Tree " << tree_name.Data() << " (path = " << path.Data() << ") has no entries\n";
      return nullptr;
    }
    TTree * tree_cut = tree->CopyTree(cuts,"",nmax);
    return tree_cut;
}


TLatex * CMS_lumi(bool IsData){
    TLatex *mark = new TLatex();

    mark->SetNDC();
    TString lumistamp = "Run 2 (13 TeV)";
    float cmsTextSize = 0.042 * 1.25;
    float extraOverCmsTextSize  = 0.76;
    float extraTextSize = extraOverCmsTextSize*cmsTextSize;

    mark->SetTextAlign(11);
    mark->SetTextSize(cmsTextSize);
    mark->SetTextFont(61);
    mark->DrawLatex(gPad->GetLeftMargin()+0.09, 1 - (gPad->GetTopMargin() - 0.017), "CMS");
    mark->SetTextSize(0.042);
    mark->SetTextFont(52);
    if (IsData)
       mark->DrawLatex(gPad->GetLeftMargin() + 0.2, 1 - (gPad->GetTopMargin() - 0.017),  "Preliminary");
    else
       mark->DrawLatex(gPad->GetLeftMargin() + 0.2, 1 - (gPad->GetTopMargin() - 0.017),  "Simulation");
    mark->SetTextSize(extraTextSize);
    mark->SetTextFont(42);
    mark->SetTextAlign(31);
    mark->DrawLatex(1 - gPad->GetRightMargin(), 1 - ( gPad->GetTopMargin() - 0.017), lumistamp);
    return mark;
}


//-----------------------------------------------------------------------------------------------------------------------------------
// Count free PDF params
int count_pdf_params(RooAbsPdf* pdf) {
  int nfree = 0;
  auto vars = pdf->getVariables();
  auto itr = vars->createIterator();
  auto var = itr->Next();
  while(var) {
    if(!((RooAbsReal*) var)->isConstant()) ++nfree;
    var = itr->Next();
  }
  return max(nfree-1,0); //remove the observable from counting
}

//-----------------------------------------------------------------------------------------------------------------------------------
// List of PDF parameters, optionally sorted alphabetically
list<RooRealVar*> list_pdf_params(RooAbsPdf* pdf, RooRealVar& obs, bool sorted = false) {
  auto vars = pdf->getVariables();
  auto itr = vars->createIterator();
  auto var = itr->Next();
  list<RooRealVar*> results;
  while(var) {
    TString name = var->GetName();
    if(name != obs.GetName()) results.push_back((RooRealVar*) var);
    var = itr->Next();
  }
  if(sorted) {
    results.sort([](const RooRealVar* a, const RooRealVar* b) { return std::string(a->GetName()) < std::string(b->GetName()); });
  }
  return results;
}

double get_manual_subrange_chisquare(RooRealVar& obs, RooAbsPdf* pdf, RooDataSet& data,
				     const char* range = nullptr, const char* norm_range = nullptr,
				     bool norm_skip = true, int* nbins = nullptr, bool IsPseudodata=false) {
  TH1* htmp_pdf  = pdf->createHistogram("htmp_chisq_pdf" , obs);
  TH1* htmp_data = data.createHistogram("htmp_chisq_data", obs);
  if(htmp_pdf->GetNbinsX() != htmp_data->GetNbinsX()) {
    cout << __func__ << ": PDF and data don't have the same number of bins ("
         << htmp_pdf->GetNbinsX() << " vs " << htmp_data->GetNbinsX() << ")!\n";
    delete htmp_pdf;
    delete htmp_data;
    return -1.;
  }

  //Create the PDF normalization by matching it to the data, skipping the region requested
  const double xmin_norm = (norm_range) ? obs.getMin(norm_range) : (norm_skip) ?  1. : obs.getMin();
  const double xmax_norm = (norm_range) ? obs.getMax(norm_range) : (norm_skip) ? -1. : obs.getMax();
  const int bin_norm_lo = (xmin_norm < xmax_norm) ? htmp_data->GetXaxis()->FindBin(xmin_norm) : (norm_skip) ? -1 : 1;
  const int bin_norm_hi = (xmin_norm < xmax_norm) ? htmp_data->GetXaxis()->FindBin(xmax_norm) : (norm_skip) ? -1 : htmp_data->GetNbinsX();
  double pdf_norm = 0.;
  double data_norm = 0.;
  for(int ibin = 1; ibin <= htmp_data->GetNbinsX(); ++ibin) {
    bool in_norm = ibin >= bin_norm_lo && ibin <= bin_norm_hi;
    if(norm_skip && in_norm) continue; //if given a range to skip
    if(!norm_skip && !in_norm) continue; //if given a range for normalizing
    const double npdf = htmp_pdf->GetBinContent(ibin)*htmp_pdf->GetBinWidth(ibin); //expected N(events)
    const double ndata = htmp_data->GetBinContent(ibin); //observed N(events)
    pdf_norm  += npdf;
    data_norm += ndata;
   // cout << "Norm Bin " << ibin << " (" << htmp_data->GetBinLowEdge(ibin) << " - " << htmp_data->GetBinLowEdge(ibin) + htmp_data->GetBinWidth(ibin) << "): Data = " << ndata << " PDF = " << npdf  << endl;
  }
  if(pdf_norm <= 0.) {
    cout << __func__ << ": PDF normalization is non-positive!\n";
    delete htmp_pdf;
    delete htmp_data;
    return -1.;
  }
  htmp_pdf->Scale(data_norm / pdf_norm);
// cout << __func__ << ": Scaling PDF:\n" << "#### Data norm = " << data_norm << endl << "#### PDF norm = " << pdf_norm << endl << "#### Scaling the PDF histogram by " << data_norm / pdf_norm << endl;

  const double xmin = obs.getMin(range);
  const double xmax = obs.getMax(range);

  double chisq = 0.;
  const int bin_lo = max(1, htmp_data->GetXaxis()->FindBin(xmin));
  const int bin_hi = min(htmp_data->GetNbinsX(), htmp_data->GetXaxis()->FindBin(xmax));
  for(int ibin = bin_lo; ibin <= bin_hi; ++ibin) {
    const double x_data = htmp_data->GetBinCenter(ibin);
    const double x_pdf  = htmp_pdf ->GetBinCenter(ibin);
    if(x_data != x_pdf) {
      cout << __func__ << ": Warning! Data center = " << x_data << " but PDF center = " << x_pdf << endl;
    }
    const double npdf = htmp_pdf->GetBinContent(ibin)*htmp_pdf->GetBinWidth(ibin);
    const double error = htmp_data->GetBinError  (ibin);

    const double ndata = htmp_data->GetBinContent(ibin);
    const double val  = ndata - npdf;
    //cout<<"val "<<val<<" error "<<error<<" npdf "<<npdf<<endl;
    const double sigma = val*val / ((IsPseudodata) ? error*error : (npdf <= 0. ? 1.e-5 : npdf));
   // const double sigma = val*val/(error*error);
    chisq += sigma;

 //     cout << "Bin " << ibin << " (" << htmp_data->GetBinLowEdge(ibin) << " - "<< htmp_data->GetBinLowEdge(ibin) + htmp_data->GetBinWidth(ibin) << "): Data = " << ndata << " PDF = " << npdf<< " --> sigma = " << sigma << endl;
    }

//  cout << __func__ << ": Total chi^2 = " << chisq << " / " << bin_hi - bin_lo+1 << " bins\n";
  if(nbins) *nbins = bin_hi - bin_lo+1;
  delete htmp_pdf;
  delete htmp_data;
  return chisq;
}

double get_manual_subrange_chisquare(RooRealVar& obs, RooAbsPdf* pdf, RooDataHist& data,
				     const char* range = nullptr, const char* norm_range = nullptr,
				     bool norm_skip = true, int* nbins = nullptr, bool IsPseudodata=false) {
  TH1* htmp_pdf  = pdf->createHistogram("htmp_chisq_pdf" , obs);
  TH1* htmp_data = data.createHistogram("htmp_chisq_data", obs);
  if(htmp_pdf->GetNbinsX() != htmp_data->GetNbinsX()) {
    cout << __func__ << ": PDF and data don't have the same number of bins ("
         << htmp_pdf->GetNbinsX() << " vs " << htmp_data->GetNbinsX() << ")!\n";
    delete htmp_pdf;
    delete htmp_data;
    return -1.;
  }

  //Create the PDF normalization by matching it to the data, skipping the region requested
  const double xmin_norm = (norm_range) ? obs.getMin(norm_range) : (norm_skip) ?  1. : obs.getMin();
  const double xmax_norm = (norm_range) ? obs.getMax(norm_range) : (norm_skip) ? -1. : obs.getMax();
  const int bin_norm_lo = (xmin_norm < xmax_norm) ? htmp_data->GetXaxis()->FindBin(xmin_norm) : (norm_skip) ? -1 : 1;
  const int bin_norm_hi = (xmin_norm < xmax_norm) ? htmp_data->GetXaxis()->FindBin(xmax_norm) : (norm_skip) ? -1 : htmp_data->GetNbinsX();
  double pdf_norm = 0.;
  double data_norm = 0.;
  for(int ibin = 1; ibin <= htmp_data->GetNbinsX(); ++ibin) {
    bool in_norm = ibin >= bin_norm_lo && ibin <= bin_norm_hi;
    if(norm_skip && in_norm) continue; //if given a range to skip
    if(!norm_skip && !in_norm) continue; //if given a range for normalizing
    const double npdf = htmp_pdf->GetBinContent(ibin)*htmp_pdf->GetBinWidth(ibin); //expected N(events)
    const double ndata = htmp_data->GetBinContent(ibin); //observed N(events)
    pdf_norm  += npdf;
    data_norm += ndata;
   // cout << "Norm Bin " << ibin << " (" << htmp_data->GetBinLowEdge(ibin) << " - " << htmp_data->GetBinLowEdge(ibin) + htmp_data->GetBinWidth(ibin) << "): Data = " << ndata << " PDF = " << npdf  << endl;
  }
  if(pdf_norm <= 0.) {
    cout << __func__ << ": PDF normalization is non-positive!\n";
    delete htmp_pdf;
    delete htmp_data;
    return -1.;
  }
  htmp_pdf->Scale(data_norm / pdf_norm);
// cout << __func__ << ": Scaling PDF:\n" << "#### Data norm = " << data_norm << endl << "#### PDF norm = " << pdf_norm << endl << "#### Scaling the PDF histogram by " << data_norm / pdf_norm << endl;

  const double xmin = obs.getMin(range);
  const double xmax = obs.getMax(range);
  // cout << __func__ << ": Using observable range " << xmin << " - " << xmax << endl;

  double chisq = 0.;
  const int bin_lo = max(1, htmp_data->GetXaxis()->FindBin(xmin));
  const int bin_hi = min(htmp_data->GetNbinsX(), htmp_data->GetXaxis()->FindBin(xmax));
  // cout << __func__ << ": Hist binning corresponds to " << htmp_data->GetBinLowEdge(bin_lo)
  //      << " - " << htmp_data->GetXaxis()->GetBinUpEdge(bin_hi) << " (widths = " << htmp_data->GetBinWidth(1) << ")" << endl;
  for(int ibin = bin_lo; ibin <= bin_hi; ++ibin) {
    const double x_data = htmp_data->GetBinCenter(ibin);
    const double x_pdf  = htmp_pdf ->GetBinCenter(ibin);
    if(x_data != x_pdf) {
      cout << __func__ << ": Warning! Data center = " << x_data << " but PDF center = " << x_pdf << endl;
    }
    const double npdf = htmp_pdf->GetBinContent(ibin)*htmp_pdf->GetBinWidth(ibin);
    const double error = htmp_data->GetBinError  (ibin);

    const double ndata = htmp_data->GetBinContent(ibin);
    const double val  = ndata - npdf;
    //cout<<"val "<<val<<" error "<<error<<" npdf "<<npdf<<endl;
    const double sigma = val*val / ((IsPseudodata) ? error*error : (npdf <= 0. ? 1.e-5 : npdf));
   // const double sigma = val*val/(error*error);
    chisq += sigma;

 //     cout << "Bin " << ibin << " (" << htmp_data->GetBinLowEdge(ibin) << " - "<< htmp_data->GetBinLowEdge(ibin) + htmp_data->GetBinWidth(ibin) << "): Data = " << ndata << " PDF = " << npdf<< " --> sigma = " << sigma << endl;
    }

//  cout << __func__ << ": Total chi^2 = " << chisq << " / " << bin_hi - bin_lo+1 << " bins\n";
  if(nbins) *nbins = bin_hi - bin_lo+1;
  delete htmp_pdf;
  delete htmp_data;
  return chisq;
}


double get_chi_squared(RooRealVar& obs, RooAbsPdf* pdf, RooDataSet& data, bool unblind,
		       int nbins_data, int n_param, bool ReturnNorm=true, bool IsPseudodata=false){

    if(!unblind) {
      int nbins = 0; //count of total bins used
      int nbin_running = 0; //bins used in a single subrange
      double chi_sq = get_manual_subrange_chisquare(obs, pdf, data, "left", "full", false, &nbin_running, IsPseudodata);
      nbins += nbin_running;
 //     cout<<">>>> chi2 (1): "<<chi_sq<<" bins "<<nbins_data<<" n_param "<<n_param<<" chi2/ndof "<<chi_sq/(nbins_data - (n_param))<<endl;
      chi_sq += get_manual_subrange_chisquare(obs, pdf, data, "right","full", false, &nbin_running, IsPseudodata);
      nbins += nbin_running;
      cout<<__func__<<">>>> chi2: "<<chi_sq<<" bins "<<nbins<<" n_param "<<n_param<<" chi2/ndof "<<chi_sq/(nbins - (n_param))<<endl;
      if (ReturnNorm) chi_sq/=(nbins - n_param);
      return chi_sq;
    } else {
      double chi_sq= get_manual_subrange_chisquare(obs, pdf, data, "full",nullptr, true, nullptr, IsPseudodata);
      cout<<">>>> chi2: "<<chi_sq<<" bins "<<nbins_data<<" n_param "<<n_param<<" chi2/ndof "<<chi_sq/(nbins_data - (n_param))<<endl;
      if (ReturnNorm) chi_sq/=(nbins_data - n_param);
      return chi_sq;
    }

}

double get_chi_squared(RooRealVar& obs, RooAbsPdf* pdf, RooDataHist& data, bool unblind,
		       int& nbins, int n_param, bool ReturnNorm=true, bool IsPseudodata=false){

    if(!unblind) {
      nbins = 0; //count of total bins used
      int nbin_running = 0; //bins used in a single subrange
      double chi_sq = get_manual_subrange_chisquare(obs, pdf, data, "left", "sr", true, &nbin_running, IsPseudodata);
      nbins += nbin_running;
      chi_sq += get_manual_subrange_chisquare(obs, pdf, data, "right","sr", true, &nbin_running, IsPseudodata);
      nbins += nbin_running;
      cout<<">>>> chi2: "<<chi_sq<<" bins "<<nbins<<" n_param "<<n_param<<" chi2/ndof "<<chi_sq/(nbins - (n_param))<<endl;
      if (ReturnNorm) chi_sq/=(nbins - n_param);
      return chi_sq;
    } else {
      double chi_sq= get_manual_subrange_chisquare(obs, pdf, data, "full",nullptr, true, &nbins, IsPseudodata);
      cout<<">>>> chi2: "<<chi_sq<<" bins "<<nbins<<" n_param "<<n_param<<" chi2/ndof "<<chi_sq/(nbins - n_param)<<endl;
      if (ReturnNorm) chi_sq/=(nbins - n_param);
      return chi_sq;
    }

}

TCanvas * create_canvas(TString name){
    TCanvas * ctemp = new TCanvas(name,"",50,50,800,600);
    ctemp->SetFillColor(0);
    ctemp->SetBorderMode(0);
    ctemp->SetFrameFillStyle(0);
    ctemp->SetLeftMargin( 0.12 );
    ctemp->SetRightMargin( 0.04 );
    ctemp->SetTopMargin(0.08);
    ctemp->SetBottomMargin(0.12);
    ctemp->SetTickx(0);
    ctemp->SetTicky(0);
    return ctemp;
}

void save_plot( RooPlot * xframe, TString xaxis_title, TString name, TLegend* leg=new TLegend(),TPaveText * pt=NULL, bool IsData=true,bool Logy=false){

    TCanvas * ctemp = create_canvas("c"+name);
    xframe->SetMinimum(0);
    xframe->Draw();
    xframe->GetXaxis()->SetTitle(xaxis_title);
    xframe->GetYaxis()->SetLabelSize(0.05);
    xframe->GetYaxis()->SetLabelOffset(0.01);
    xframe->GetYaxis()->SetTitleSize(0.06);
    xframe->GetYaxis()->SetTitleOffset(0.90);
    xframe->GetXaxis()->SetLabelSize(0.05);
    xframe->GetXaxis()->SetLabelOffset(0.007);
    xframe->GetXaxis()->SetTitleSize(0.07);
    xframe->GetXaxis()->SetTitleOffset(0.80);
    TLatex * logo = CMS_lumi(IsData);
    logo->Draw("sames");
    if (Logy){
       ctemp->SetLogy(true);
       xframe->SetMinimum(1);
    }
    if (leg->GetNRows()>0)
       leg->Draw("sames");
    if (pt!=NULL)
       pt->Draw("sames");
    ctemp->SaveAs(figdir_+name+".png");
}

void save_plot_and_band( RooPlot * xframe,  RooRealVar var, std::vector<TString> functions, TString xaxis_title, TString name, TLegend* leg=new TLegend(),TPaveText * pt=NULL, bool IsData=true,bool Logy=false){

    TCanvas * ctemp = create_canvas("c"+name);
    TPad * pad1 = new TPad("pad1","pad1",0,0.24,1,1.0);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);
    pad1->SetBottomMargin(0.11);
    pad1 ->SetFillColor(0);
    pad1->Draw();
    pad1->cd();
    pad1->SetTickx(0);
    pad1->SetTicky(0);
    xframe->SetMinimum(0);
    xframe->SetTitle("");
    xframe->GetYaxis()->SetTitle(Form("Events / ( %.2f GeV/c^{2} )", xframe->GetXaxis()->GetBinWidth(1)));
    xframe->GetYaxis()->SetLabelSize(0.05);
    xframe->GetYaxis()->SetLabelOffset(0.01);
    xframe->GetYaxis()->SetTitleSize(0.06);
    xframe->GetYaxis()->SetTitleOffset(0.95);
    xframe->GetYaxis()->SetMaxDigits(3);
    xframe->GetXaxis()->SetLabelSize(0);
    xframe->GetXaxis()->SetLabelOffset(0.007);
    xframe->GetXaxis()->SetTitleSize(0.07);
    xframe->GetXaxis()->SetTitleOffset(0.80);
    xframe->Draw("PE1");
    TLatex * logo = CMS_lumi(IsData);
    logo->Draw("sames");
    if (Logy){
       ctemp->SetLogy(true);
       xframe->SetMinimum(1);
    }
    if (leg->GetNRows()>0) {
       xframe->SetMaximum(1.3*xframe->GetMaximum()); //create more room for the legend
       leg->Draw("sames");
    }
    if (pt!=NULL)
       pt->Draw("sames");
    ctemp->cd();
    TPad * pad2 = new TPad("pad2","pad2",0,0.0,1,0.32);
    pad2->SetTopMargin(0.03);
    pad2->SetLeftMargin(pad1->GetLeftMargin());
    pad2->SetRightMargin(pad1->GetRightMargin());
    pad2->SetBottomMargin(0.35);
    pad2 ->SetFillColor(0);
    pad2->SetTickx(0);
    pad2->SetTicky(0);
    pad2->Draw();
    pad2->cd();
    auto xframe3 = var.frame();

  //  for (auto fnc :functions)
//      cout<<"fnc "<<fnc<<endl;
    for (int i =0; i<functions.size(); i++){
      auto hpull = xframe->pullHist("data",functions[i]);
      hpull->SetName("ratio_fnc");
      hpull->SetLineColor(i+1 + (i >= 4)); //skip yellow due to the difficulty to see
      hpull->SetLineWidth(2);
      hpull->SetMarkerSize(0);
      xframe3->addPlotable(hpull,"PE1");
    }
    xframe3->GetYaxis()->SetRangeUser(-4,4);
    xframe3->GetYaxis()->SetNdivisions(5);
    xframe3->GetYaxis()->SetLabelSize(0.125);
    xframe3->GetYaxis()->SetLabelOffset(0.01);
    xframe3->GetYaxis()->SetTitleSize(0.145);
    xframe3->GetYaxis()->SetTitleOffset(0.35);
    xframe3->GetXaxis()->SetLabelSize(0.15);
    xframe3->GetXaxis()->SetLabelOffset(0.008);
    xframe3->GetXaxis()->SetTitleSize(0.2);
    xframe3->GetXaxis()->SetTitleOffset(0.72);

    xframe3->SetTitle("");
    xframe3->GetYaxis()->SetTitle("#frac{N_{data} - N_{fit}}{#sigma_{data}}");
    xframe3->GetXaxis()->SetTitle(xaxis_title);
    xframe3->Draw();
    TLine line(xframe3->GetXaxis()->GetXmin(), 0., xframe3->GetXaxis()->GetXmax(), 0.);
    line.SetLineWidth(2);
    line.SetLineStyle(kDashed);
    line.SetLineColor(kBlack);
    line.Draw("same");
//    TLine *line = new TLine(gPad->GetUxmin()+0.12, 0.665, gPad->GetUxmax()-0.04, 0.665);
  //  line->SetNDC(kTRUE);
   // line->Draw();
    ctemp->SaveAs(figdir_+name+".png");
}



void save_pull( RooPlot * xframe, RooRealVar var, TString xaxis_title, TString name){

    TCanvas * cresid = create_canvas("cresid_"+name);
    auto hresid = xframe->residHist();
    auto xframe2 = var.frame();
    xframe2->addPlotable(hresid,"P");
    xframe2->Draw();
    xframe2->SetTitle(xaxis_title);
    cresid->SaveAs(figdir_+"cresid_"+name+".png");
    TCanvas * cpull= create_canvas("cpull_"+name);
    auto hpull = xframe->pullHist();
    auto xframe3 = var.frame();
    hpull->GetYaxis()->SetRangeUser(-5,5);
    xframe3->addPlotable(hpull,"P");
    xframe3->GetYaxis()->SetRangeUser(-5,5);
    xframe3->Draw();
    xframe3->GetXaxis()->SetTitle(xaxis_title);
    cpull->SaveAs(figdir_+"cpull_"+name+".png");
}

void print_details(RooFitResult * result){
    result->Print();
    cout<<"edm "<<result->edm()<<" log "<<result->minNll()<<endl;
}


void create_pdf_from_mc(TString tree_name, TString path_to_pdf_mc, TString cuts, int nbin, float min_fit_range,float max_fit_range, TString dilep_var_name, TString pdf_outfile_name ){
   TTree * mc_tree = get_tree(tree_name,path_to_pdf_mc,cuts);
   TH1F* hmc_pdf = new TH1F("hpdf","",nbin,min_fit_range,max_fit_range);
   mc_tree->Draw(dilep_var_name+">>hpdf",cuts);
   TFile * fmc_out = new TFile(pdf_outfile_name,"RECREATE");
   hmc_pdf->Write();
   fmc_out->Close();
}


void save_th1(std::vector<TH1F*> histos, std::vector<TString> legs, TString xaxis, TString name, bool FitGauss=false){

   TCanvas * cnv = create_canvas("c"+name);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);

   TLegend * leg = new TLegend(0.7,0.7,0.9,0.9);
   for (int ih=0; ih<histos.size(); ih++){
     leg->AddEntry(histos[ih],legs[ih]);
     if (ih ==0)
       histos[ih]->Draw("PE");
     else
       histos[ih]->Draw("PE sames");
     histos[ih]->SetMarkerSize(2);
     histos[ih]->SetLineWidth(2);
     histos[ih]->SetMarkerColor(ih+1);
     histos[ih]->SetLineColor(ih+1);
     histos[ih]->GetXaxis()->SetTitle(xaxis);
     histos[ih]->GetYaxis()->SetLabelSize(0.05);
     histos[ih]->GetYaxis()->SetLabelOffset(0.01);
     histos[ih]->GetYaxis()->SetTitleSize(0.06);
     histos[ih]->GetYaxis()->SetTitleOffset(0.90);
     histos[ih]->GetXaxis()->SetLabelSize(0.05);
     histos[ih]->GetXaxis()->SetLabelOffset(0.007);
     histos[ih]->GetXaxis()->SetTitleSize(0.07);
     histos[ih]->GetXaxis()->SetTitleOffset(0.80);

     TLatex * logo = CMS_lumi(false);
     logo->Draw("sames");

     if (FitGauss){
       TF1 *gauss = new TF1("gauss","gaus",-5,5);
       gauss->SetLineColor(kRed);
       histos[ih]->Fit("gauss","RL");
       double param_pull[3];
       double param_errors_pull[] = {gauss->GetParError(0), gauss->GetParError(1), gauss->GetParError(2)};
       gauss->GetParameters(&param_pull[0]);
       TPaveText *  pt_pull = new TPaveText(0.2,0.7,0.5,0.9,"tlNDC");
       pt_pull->SetFillColor(0);
       pt_pull->SetFillStyle(0);
       pt_pull->SetBorderSize(0);
       pt_pull->AddText("Mean "+TString(std::to_string(param_pull[1]))+" #pm"+TString(std::to_string( param_errors_pull[1] )));
       pt_pull->AddText("Sigma "+TString(std::to_string(param_pull[2]))+" #pm"+TString(std::to_string( param_errors_pull[2] )));
       pt_pull->Draw("sames");
     }
  }

    if (histos.size()>1)
       leg->Draw("sames");

     cnv->SaveAs(figdir_+name+".png");
}

#endif
