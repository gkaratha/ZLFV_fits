from array import array
import ROOT as rt

#----------------------------------------------------------------------------------------
# Write a combine data card for a single bin
def print_datacard(name, sig_file, bkg_file, param_name, mass):
   # Standard preamble
   txt="# -*- mode: tcl -*-\n"
   txt+="#Auto generated Z prime search COMBINE datacard\n"
   txt+="#Using Z prime mass = %.3f\n\n" % (mass)
   txt+="imax * number of channels\njmax * number of backgrounds\nkmax * number of nuisance parameters\n\n"

   # Define the background and signal PDFs
   txt+="#-----------------------------------------------------------------------------------------------------------\n"
   txt+="shapes signal "+param_name+" %s ws_sgn:signal_pdf_%s\n" % (sig_file, param_name)
   txt+="shapes background "+param_name+" %s ws_bkg:multipdf_%s\n" % (bkg_file, param_name)
   txt+="shapes data_obs "+param_name+" %s ws_bkg:data_obs\n" % (bkg_file)

   txt+="#-----------------------------------------------------------------------------------------------------------\n\n"
   # Define the channel
   txt+="bin         "+param_name+"\n"
   txt+="observation  -1\n\n"
   txt+="bin         "+param_name+"    "+param_name+"\n"
   txt+="process    signal   background\n\n"
   txt+="process       0         1\n\n"
   txt+="rate          1         1\n" #Rate is taken from the _norm variables

   # Define the uncertainties
   txt+="#-----------------------------------------------------------------------------------------------------------\n"
   #FIXME: Implement mass-dependent uncertainties
   txt+="ElectronID lnN 1.02     -\n"
   txt+="MuonID     lnN 1.02     -\n"
   txt+="TauID      lnN 1.05     -\n"
   txt+="Lumi       lnN 1.02     -\n"
   txt+="BTag       lnN 1.005    -\n"
   txt+="Theory     lnN 1.01     -\n"
   txt+="BDT        lnN 1.02     -\n"

   txt+="#-----------------------------------------------------------------------------------------------------------\n\n"
   # Scale uncertainties
   txt+="lep_ES_shift param 0 1 [-7, 7]\n"
   txt+="tau_ES_shift param 0 1 [-7, 7]\n"

   txt+="#-----------------------------------------------------------------------------------------------------------\n\n"
   # Define the envelope discrete index to be scanned
   txt+="pdfindex_%s discrete\n" % (param_name)

   # Write the file
   # with open(carddir + "/datacard_zprime_" + name +"_mass-"+str(mass)+ "_mp" + cnt + ".txt",'w') as fl:
   with open(name,'w') as fl:
     fl.write(txt)
   fl.close()

