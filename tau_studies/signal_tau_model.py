# Fit the signal MC and create an interpolation model
import ROOT as rt
from array import array
from math import erf
import ctypes

#----------------------------------------------------------------------------------------
# Signal sample structure
class sample:
    def __init__(self, file_path, n_gen, year, base_path, mass = 0.):
        self.file_path_ = file_path
        self.n_gen_ = n_gen
        self.year_ = year
        self.base_path_ = base_path
        self.full_path_ = base_path
        if base_path[-1] != "/": self.full_path_ += "/"
        self.full_path_ += file_path

        if mass <= 0.: # try to get the mass from the file name
            mass = float(file_path.split('_sgnM')[1].split('_mcRun')[0])
        self.mass_ = mass

        if n_gen <= 0:
            # Try to get the normalization
            cc=rt.TChain("norm")
            cc.Add(self.full_path_)
            n_gen = 0
            for entry in range(cc.GetEntries()):
                cc.GetEntry(entry)
                n_gen += cc.nevents - 2*cc.nnegative
            self.n_gen_ = n_gen


    def __repr__(self):
        return "Sample:\n file = %s\n n_gen = %i\n year = %i\n base_path = %s\n" % (self.file_path_, self.n_gen_, self.year_, self.base_path_)

#----------------------------------------------------------------------------------------
# Retrieve the signal distribution, properly weighing each year's component
def signal_distribution(sample_map, h, var, cuts, period = "Run2", correct_samples = True):
    lumis={"2016":36.33,"2017":41.48,"2018":59.83,"Run2":137.6}
    periods = ["2016", "2017", "2018"] if period == "Run2" else [period]
    for p in periods:
        cc=rt.TChain("mytree")
        cc.Add(sample_map[p].full_path_)
        if cc.GetEntries() == 0:
            print "Sample path", sample_map[p].full_path_, "has no entries!"
        htmp = h.Clone("htmp")
        htmp.Reset()
        cc.Draw(var + ">>htmp", cuts)
        cross_section = 1. #Use 1 fb cross section
        scale = lumis[p]/sample_map[p].n_gen_
        htmp.Scale(scale)
        h.Add(htmp)
    if h.Integral() <= 0.:
        print "Histogram", var, cuts, period, "has no integral!"
    return h


#----------------------------------------------------------------------------------------
# Fit the signal mass distributions and create an interpolation map
def interpolate(param_fits, mass):
    params = []
    mass = mass/1000.
    index = 0
    for param_fit in param_fits:
        if index == 0: #yield is modeled with a function that includes the turn on
            params.append(param_fit[0]*(1. + erf(param_fit[1]*mass + param_fit[2])) + param_fit[3]*mass)
        else: params.append(param_fit[0] + param_fit[1]*mass + param_fit[2]*mass*mass)
        index += 1
    return params

#----------------------------------------------------------------------------------------
# Fit the signal mass distributions and create an interpolation map
def create_signal_interpolation(masses, distributions, use_gaus = False, figdir = ''):
   # List of distributions includes the overall efficiency for the selection

   # Loop through each distribution, fitting the MC
   mass_errs = array('d')
   fit_results = [] # Store the fit parameters
   fit_errs = []
   for h in distributions:
      mass_errs.append(0.) # For plotting purposes
      # Create a RooFit dataset and corresponding PDF to fit the signal distribution
      h_mean = h.GetBinCenter(h.GetMaximumBin())
      h_width = h.GetStdDev()

      # Define the fit range (only the core for Gaussian modeling)
      if use_gaus:
         min_mass = h_mean - 0.75*h_width
         max_mass = h_mean + 0.75*h_width
      else:
         min_mass = h_mean - 5.*h_width if h_mean > 95. else 70.
         max_mass = h_mean + 5.*h_width if h_mean > 95. else 110.
      min_mass = max(70., min_mass)

      # Create a RooFit setup to perform the fit
      obs = rt.RooRealVar("obs", "obs", h_mean, min_mass, max_mass, "GeV/c^{2}")
      dh = rt.RooDataHist("dh", "Mass distribution", obs, h)

      mean   = rt.RooRealVar("mean", "mean", h_mean, 0.9*h_mean, 1.1*h_mean)
      sigma  = rt.RooRealVar("sigma", "sigma", h_mean/20., h_mean/100., h_mean);
      alpha1 = rt.RooRealVar("alpha1", "alpha1", 1., 0.5, 10.);
      alpha2 = rt.RooRealVar("alpha2", "alpha2", 1., 0.5, 10.);
      enne1  = rt.RooRealVar("enne1", "enne1", 5., 0.1, 10.);
      enne2  = rt.RooRealVar("enne2", "enne2", 5., 0.1, 10.);
      if use_gaus:
         pdf = rt.RooGaussian("pdf", "PDF", obs, mean, sigma)
      else:
         pdf = rt.RooDoubleCrystalBall("pdf", "PDF", obs, mean, sigma, alpha1, enne1, alpha2, enne2)

      # Fit the PDF to the signal MC
      pdf.fitTo(dh, rt.RooFit.PrintLevel(-1), rt.RooFit.Warnings(0), rt.RooFit.PrintEvalErrors(-1), rt.RooFit.SumW2Error(1))

      # Print the fit results
      frame = obs.frame(rt.RooFit.Title("Signal model fit"))
      dh.plotOn(frame)
      pdf.plotOn(frame)

      c = rt.TCanvas()
      frame.Draw()
      c.SaveAs(figdir+h.GetName()+"_fit.png")
      max_val = h.GetMaximum()
      frame.GetYaxis().SetRangeUser(1.e-2*max_val if use_gaus else 1.e-4*max_val, 5.*max_val)
      c.SetLogy()
      c.SaveAs(figdir+h.GetName()+"_fit_log.png")

      # Include the efficiency in the yield
      eff_err = ctypes.c_double(1.)
      eff = h.IntegralAndError(1, h.GetNbinsX(), eff_err)
      fit_results.append([eff          , mean.getVal()  , sigma.getVal()  , alpha1.getVal()  , alpha2.getVal()  , enne1.getVal()  , enne2.getVal()  ])
      fit_errs   .append([eff_err.value, mean.getError(), sigma.getError(), alpha1.getError(), alpha2.getError(), enne1.getError(), enne2.getError()])

   # Create the interpolation by fitting the parameters as a function of mass
   if use_gaus:
      param_names = ["yield", "#mu", "#sigma"]
   else:
      param_names = ["yield", "#mu", "#sigma", "#alpha_{1}", "#alpha_{2}", "n_{1}", "n_{2}"]

   signal_params = []

   # Fit each parameter as a function of mass
   for param in range(len(param_names)):
      name = param_names[param]
      vals = array('d')
      errs = array('d')
      for index in range(len(fit_results)):
         vals.append(fit_results[index][param])
         errs.append(fit_errs   [index][param])

      # Plot the parameters and the fit result
      c = rt.TCanvas('c_param', 'c_param', 700, 500)
      c.SetBottomMargin(0.13)
      c.SetRightMargin(0.05)
      g = rt.TGraphErrors(len(vals), masses, vals, mass_errs, errs)
      g.SetTitle("Model %s parameter vs. mass;Z prime mass (GeV/c^{2});%s" % (name, name))
      g.SetMarkerStyle(20)
      g.SetMarkerSize(0.8)
      g.SetLineWidth(2)
      g.Draw("APE")
      min_val = min([vals[index]-errs[index] for index in range(len(vals))])
      max_val = max([vals[index]+errs[index] for index in range(len(vals))])
      g.GetYaxis().SetRangeUser(min_val - 0.1*(max_val-min_val), max_val + 0.1*(max_val-min_val))
      g.GetYaxis().SetLabelSize(0.04)
      g.GetXaxis().SetLabelSize(0.04)
      g.GetYaxis().SetTitleSize(0.06)
      g.GetYaxis().SetTitleOffset(0.65)
      g.GetXaxis().SetTitleSize(0.06)
      g.GetXaxis().SetTitleOffset(0.75)

      # Divide mass by 1,000 for numerical stability
      if name == "yield": #Use a higher order function for the yield turn on
          func = rt.TF1("func", "[0]*(1. + erf([1]*x/1000. + [2])) + [3]*x/1000", 0., 1000.)
          func.SetParameters(0.5, 1., 0., 0.)
      else:
          func = rt.TF1("func", "[0] + [1]*(x/1000.) + [2]*pow(x/1000.,2)", 0., 1000.)
          func.SetParameters(0.5, 0., 0.)
    
      if len(vals) == 1: # just use this value if only 1 data point
          if name == "yield": #Use a higher order function for the yield turn on
              func.SetParameters(vals[0], 0., 0., 0.)
          elif name == "#mu" or name == "#sigma":
              func.SetParameters(0., vals[0]*1000./masses[0], 0.)
          else:
              func.SetParameters(vals[0], 0., 0.)
      else:
          g.Fit(func, "R")
      func.Draw("same")
      c.SaveAs("%sparam_%i.png" % (figdir, param))

      # Store the fit result
      if name == "yield":
          signal_params.append([func.GetParameter(0), func.GetParameter(1), func.GetParameter(2), func.GetParameter(3)])
          print signal_params[-1]
      else:
          signal_params.append([func.GetParameter(0), func.GetParameter(1), func.GetParameter(2)])
   

   # Return the list of interpolation fit results, one fit per parameter
   return signal_params
