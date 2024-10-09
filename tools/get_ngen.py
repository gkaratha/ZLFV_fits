# Get the N(gen) information for the Z prime MC samples for normalization purposes
import os
import ROOT as rt

rt.gROOT.SetBatch(True)

#-----------------------------------------------------------------
# Signal sample definitions
#-----------------------------------------------------------------

# Default ntuple path
path="/eos/cms/store/group/phys_smp/ZLFV/MC_generation/"

# Signal samples to process
sgn_masspoints=["125","150","175"]
years = ["2018"]

# Signal samples by year and mass
sgn_masspoint_files = {"2016": {}, "2017": {}, "2018": {}}

sgn_masspoint_files["2016"]={
    "100" :"ZEMu_NANO_M100_2016_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M100_2016_100k4/240906_125113/0000/*.root",
    "500" :"ZEMu_NANO_2016_M500_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M500_100k4/240830_140542/0000/*.root",
}
sgn_masspoint_files["2017"]={
    "100" :"ZEMu_NANO_2017_M100_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M100_100k4/240830_171613/0000/*.root",
    "500" :"ZEMu_NANO_2017_M500_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M500_100k4/240830_173217/0000/*.root",
}
sgn_masspoint_files["2018"]={
    "100" :"ZEMu_NANO_2018_M100_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M100_100k4/240830_135648/0000/*.root",
    "125" :"ZEMu_NANO_M125_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M125_2018_100k4/240917_113503/0000/*.root",
    "150" :"ZEMu_NANO_M150_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M150_2018_100k4/240917_113734/0000/*.root",
    "175" :"ZEMu_NANO_M175_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M175_2018_100k4/240917_113816/0000/*.root",
    "200" :"ZEMu_NANO_M200_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M200_2018_100k4/240710_144622/0000/*.root",
    "300" :"ZEMu_NANO_2018_M300_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M300_100k4/240830_170925/0000/*.root",
    "400" :"ZEMu_NANO_M400_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M400_2018_100k4/240710_144531/0000/*.root",
    "500" :"ZEMu_NANO_2018_M500_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M500_100k4/240830_124627/0000/*.root",
    "600" :"ZEMu_NANO_M600_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M600_2018_100k4/240710_144439/0000/*.root",
    "800" :"ZEMu_NANO_M800_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M800_2018_100k4/240710_143659/0000/*.root",
    "1000":"ZEMu_NANO_M1000_2018_100k4/CRAB_UserFiles/ZLFVAnalysis_NANO_M1000_2018_100k4/240710_143332/0000/*.root",
}

#-----------------------------------------------------------------
# Process each requested sample
#-----------------------------------------------------------------

for year in years:
    for mpoint in sgn_masspoints:
        if mpoint not in sgn_masspoint_files[year]: continue
        # Create a TChain for the signal
        cc=rt.TChain("Events")
        filepath = path+sgn_masspoint_files[year][mpoint]
        cc.Add(filepath)
        norm = cc.GetEntries()
        print "Norm for %s mass point %-4s: %i" % (year, mpoint, norm)
