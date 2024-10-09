//Draw workspace PDFs
TString get_pdf_title(TString pdf_name) {
  TString title = "";
  if(pdf_name.Contains("_gs"))   title += "Gaussian + ";
  if(pdf_name.Contains("pol1"))  title += "Poly (1)";
  if(pdf_name.Contains("pol2"))  title += "Poly (2)";
  if(pdf_name.Contains("pol3"))  title += "Poly (3)";
  if(pdf_name.Contains("pol4"))  title += "Poly (4)";
  if(pdf_name.Contains("cheb1")) title += "Chebychev (1)";
  if(pdf_name.Contains("cheb2")) title += "Chebychev (2)";
  if(pdf_name.Contains("cheb3")) title += "Chebychev (3)";
  if(pdf_name.Contains("cheb4")) title += "Chebychev (4)";
  if(pdf_name.Contains("cheb5")) title += "Chebychev (5)";
  if(pdf_name.Contains("plaw1")) title += "Power (1)";
  if(pdf_name.Contains("plaw2")) title += "Power (2)";
  if(pdf_name.Contains("plaw3")) title += "Power (3)";
  if(pdf_name.Contains("exp1"))  title += "Expo (1)";
  if(pdf_name.Contains("exp2"))  title += "Expo (2)";
  if(pdf_name.Contains("exp3"))  title += "Expo (3)";
  return title;
}

int debug_zemu_ws(const char* bin = "bin3") {

  // Get the name of the observable
  TString bin_s(bin);
  TString obs_name = "lepm_1";
  if(bin_s.BeginsWith("bin")) obs_name = obs_name + bin_s[bin_s.Sizeof()-2];
  else obs_name = "lepm_" + bin_s;

  ////////////////////////////////////////
  // Retrieve the input background PDFs
  ////////////////////////////////////////

  TFile* f_bkg = TFile::Open(Form("workspace_mk2bkg_v1_%s.root", bin), "READ");
  if(!f_bkg) return 1;
  auto ws = (RooWorkspace*) f_bkg->Get("workspace_background");
  if(!ws) return 2;
  auto multipdf = (RooMultiPdf*) ws->pdf(Form("multipdf_%s", bin));
  if(!multipdf) return 3;
  auto cat = (RooCategory*) ws->obj(Form("pdfindex_%s", bin));
  if(!cat) return 4;
  cout << "Using obs name " << obs_name.Data() << endl;
  auto obs = (RooRealVar*) ws->var(obs_name);
  if(!obs) return 5;
  const float xmin(70.f), xmax(110.f), blind_min(86.f), blind_max(96.f);
  obs->setRange("full", xmin, xmax);
  obs->setRange("LowSideband", xmin, blind_min-1.e-4);
  obs->setRange("HighSideband", blind_max+1.e-4, xmax);
  obs->setRange("BlindRegion", blind_min, blind_max);

  ////////////////////////////////////////
  // Retrieve the input data
  ////////////////////////////////////////

  // workspace_mk2dataAndMC_v1_bin1.root workspace_data:real_data_bin1
  TFile* f_data = TFile::Open(Form("workspace_mk2dataAndMC_v1_%s.root", bin), "READ");
  if(!f_data) return 11;
  auto ws_data = (RooWorkspace*) f_data->Get("workspace_data");
  if(!ws_data) return 12;
  auto data = ws_data->data(Form("real_data_%s", bin));
  if(!data) return 13;

  ////////////////////////////////////////
  // Retrieve the Z->mumu PDF
  ////////////////////////////////////////

  // workspace_mk2zmm_v1_bin1.root workspace_zmm:zmm_dcb_pdf_bin1
  TFile* f_zmm = TFile::Open(Form("workspace_mk2zmm_v1_%s.root", bin), "READ");
  if(!f_zmm) return 21;
  auto ws_zmm = (RooWorkspace*) f_zmm->Get("workspace_zmm");
  if(!ws_zmm) return 22;
  auto zmm = ws_zmm->pdf(Form("zmm_dcb_pdf_%s", bin));
  if(!zmm) return 23;
  auto zmm_norm = (RooRealVar*) ws_zmm->var(Form("zmm_dcb_pdf_%s_norm", bin));
  if(!zmm_norm) return 24;

  ////////////////////////////////////////
  // Make a plot of the inputs
  ////////////////////////////////////////

  //Create a legend
  TLegend* leg = new TLegend(0.11, 0.65, 0.94, 0.89);
  leg->SetNColumns(3);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  auto frame = obs->frame();
  data->plotOn(frame, RooFit::Name("data"));

  const int ncats = cat->numTypes();
  const int colors[] = {kRed-6, kGreen-6, kOrange-4, kViolet-4};
  // const int colors[] = {kRed-6, kGreen-6, kOrange-4, kViolet-4, kYellow-3, kAtlantic, kOrange+1, kSpring+5};
  const int ncolors = sizeof(colors)/sizeof(*colors);
  vector<TString> titles;
  for(int ipdf = 0; ipdf < ncats; ++ipdf) {
    auto pdf = multipdf->getPdf(ipdf);
    TString title = get_pdf_title(pdf->GetName());
    titles.push_back(title);
    pdf->plotOn(frame, RooFit::LineColor(colors[ipdf % ncolors]), RooFit::LineStyle((ipdf >= ncolors) ? kDashed : kSolid),
                RooFit::Title(title), RooFit::Name(Form("pdf_%i", ipdf)));
  }

  TCanvas c; c.SetRightMargin(0.05);
  frame->Draw();
  leg->AddEntry("data", "Data", "PL");
  for(int ipdf = 0; ipdf < ncats; ++ipdf) leg->AddEntry(Form("pdf_%i", ipdf), titles[ipdf], "L");
  leg->Draw("same");
  frame->GetYaxis()->SetRangeUser(0.1, 1.4*frame->GetMaximum());
  c.SaveAs(Form("inputs_%s.png", bin));

  // Draw the Z->mumu PDF
  frame = obs->frame();
  zmm->plotOn(frame);
  frame->Draw();
  c.SaveAs(Form("inputs_zmm_%s.png", bin));

  ////////////////////////////////////////
  // Re-fit the Z->mumu+PDFs to the data
  ////////////////////////////////////////

  // Blind the data for this plot
  TH1* h_data = data->createHistogram("hdata", *obs);
  for(int ibin = h_data->FindBin(blind_min+0.01); ibin <= h_data->FindBin(blind_max-0.1); ++ibin) {
    h_data->SetBinContent(ibin, 0.);
    h_data->SetBinError(ibin, 0.);
  }
  auto blind_data = new RooDataHist("blind_data", "Blinded data", *obs, h_data);

  TCanvas c2("c2","c2", 700, 800);
  TPad pad1("pad1","pad1",0.,0.4,1.,1.0); pad1.Draw();
  TPad pad2("pad2","pad2",0.,0.0,1.,0.4); pad2.Draw();
  pad1.SetTopMargin(0.1);
  pad1.SetRightMargin(0.05);
  pad1.SetBottomMargin(0.01);
  pad2.SetTopMargin(0.03);
  pad2.SetRightMargin(0.05);
  pad2.SetBottomMargin(0.20);

  pad1.cd();
  frame = obs->frame();
  auto frame2 = obs->frame();
  blind_data->plotOn(frame, RooFit::Name("data"));
  const double frac_zmm = zmm_norm->getVal() / data->sumEntries();
  RooRealVar* pdf_frac = new RooRealVar("pdf_frac", "PDF fraction", 1. - frac_zmm);
  RooRealVar* zmm_val = new RooRealVar("zmm_val", "Z->mumu events", zmm_norm->getVal());
  pdf_frac->setConstant(true);
  zmm_val->setConstant(true);
  for(int ipdf = 0; ipdf < ncats; ++ipdf) {
    auto pdf = multipdf->getPdf(ipdf);
    // RooAddPdf* tot_pdf = new RooAddPdf(Form("tot_pdf_%i", ipdf), titles[ipdf], *pdf, *zmm, *pdf_frac);
    RooRealVar* pdf_val = new RooRealVar(Form("pdf_val_%i", ipdf), "PDF events", data->sumEntries() - zmm_norm->getVal());
    RooAddPdf* tot_pdf = new RooAddPdf(Form("tot_pdf_%i", ipdf), titles[ipdf], RooArgList(*pdf, *zmm), RooArgList(*pdf_val, *zmm_val));
    tot_pdf->fitTo(*data, RooFit::PrintLevel(-1), RooFit::Warnings(0), RooFit::PrintEvalErrors(-1),
                   RooFit::Range("LowSideband,HighSideband"), RooFit::NormRange("LowSideband,HighSideband"));
    tot_pdf->plotOn(frame, RooFit::LineColor(colors[ipdf % ncolors]), RooFit::LineStyle((ipdf >= ncolors) ? kDashed : kSolid),
                    RooFit::Title(titles[ipdf]), RooFit::Name(Form("pdf_%i", ipdf)), RooFit::Range("full"), RooFit::NormRange("LowSideband,HighSideband"));
    if(ipdf == 0) {
      tot_pdf->plotOn(frame, RooFit::Components(zmm->GetName()), RooFit::LineColor(kGreen), RooFit::Name("zmm"),
                      RooFit::Range("full"), RooFit::NormRange("LowSideband,HighSideband"));
    }
    auto h = frame->pullHist("data", Form("pdf_%i", ipdf));
    h->SetLineWidth(2);
    h->SetLineColor(colors[ipdf % ncolors]);
    h->SetLineStyle((ipdf >= ncolors) ? kDashed : kSolid);
    h->SetMarkerSize(0.1);
    frame2->addPlotable(h, "PE1");
  }
  frame->Draw();
  leg->SetTextSize(0.036);
  leg->Draw("same");
  frame->GetYaxis()->SetRangeUser(0.1, 1.4*frame->GetMaximum());
  frame->GetXaxis()->SetLabelSize(0.);

  pad2.cd();
  frame2->Draw();

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("");
  frame2->SetTitle("");
  frame2->GetYaxis()->SetTitle("#frac{N_{data} - N_{fit}}{#sigma_{data}}");

  frame->GetYaxis()->SetLabelSize(0.04);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetTitleOffset(1.2);

  frame2->GetYaxis()->SetRangeUser(-4., 4.);
  frame2->GetYaxis()->SetNdivisions(5);
  frame2->GetYaxis()->SetLabelSize(0.07);
  frame2->GetYaxis()->SetLabelOffset(0.01);
  frame2->GetYaxis()->SetTitleSize(0.07);
  frame2->GetYaxis()->SetTitleOffset(0.58);
  frame2->GetXaxis()->SetLabelSize(0.07);
  frame2->GetXaxis()->SetLabelOffset(0.008);
  frame2->GetXaxis()->SetTitleSize(0.07);
  frame2->GetXaxis()->SetTitleOffset(1.0);

  TLine line(frame->GetXaxis()->GetXmin(), 0., frame->GetXaxis()->GetXmax(), 0.);
  line.SetLineWidth(2);
  line.SetLineColor(kBlack);
  line.SetLineStyle(kDashed);
  line.Draw("same");
  pad2.SetGridy();

  c2.SaveAs(Form("inputs_refit_%s.png", bin));
  return 0;
}
