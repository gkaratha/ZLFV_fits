# LFV Z resonance fits

Code base to perform LFV Z decay resonance fits.

## Z to e+mu decay resonance search

The Z to e+mu decay signature is a resonance at the Z mass in the e-mu data. This framework uses a BDT to separate the signal from the background,
without disturbing this resonance signature, and then perform resonance fits simultaneously in categories of the BDT score.

### Signal model

The signal is modeled with a double-sided Crystal Ball, fit to the signal MC separately for each BDT score region.
The signal MC fits are peformed using [ZMuE_fit_mk2_sgn_v1.C](ZMuE_fit_mk2_sgn_v1.C).

### Z->mumu model

The Z->mumu background is modeled with a double-sided Crystal Ball fit to the MC.
The fits are performed using [ZMuE_fit_mk2_Zmm_v1.C](ZMuE_fit_mk2_Zmm_v1.C)

### Parametric background model

The continuum of e-mu background events are modeled with data-driven parametric function fits.
The background is modeled using a broad Gaussian to describe Z->tau_e+tau_mu and a second function
chosen from polynomial, exponential, and power law families using a fit quality selection and an F-test.

The fits are performed using [ZMuE_fit_mk2_bkg_v1.C](ZMuE_fit_mk2_bkg_v1.C)

### Retrieve MC templates

Smoothed MC templates are used for bias tests and studying background function choices.
These can be retrieved using [CreatePseudoData_from_MC_v2.py](CreatePseudoData_from_MC_v2.py).

```
python CreatePseudoData_from_MC_v2.py
ls -lt pseudo_data_*.root | head -n 1

# Get the smoothed templates
cp /afs/cern.ch/user/m/mimacken/public/forGeorge/zemu_smoothed_embed_v03_1?.root ./
rename "zemu_smoothed_embed_v03_1" "template_zemu_embed_v3_bin" zemu_smoothed_embed_v03_1?.root
```

### Build the model

The total model is built using [ZMuE_fit_mk2_wrapper_v1.py](ZMuE_fit_mk2_wrapper_v1.py).
The Z->e+mu search uses three exclusive BDT score-defined categories: 0.3-0.7, 0.7-0.9, and 0.9-1.0.

Running with fits to the data:
```
OUTDIR="zemu_data/"
[ ! -d ${OUTDIR} ] && mkdir ${OUTDIR}

MCFILE="pseudo_data_from_MC_v2_r0_ZmmR1.25_updateID.root"
DATAFILE="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_full_bdt_v7_data_emu_Run1*.root"
# FIXME: ZMuE_fit_mk2_datagen_v1.C attempts to use the background fit file for an MC dataset, which fails on data
python ZMuE_fit_mk2_wrapper_v1.py -o bin1 --fit-version 1 --skip-sgn-syst --zmm-file ${MCFILE} --bkg-file "${DATAFILE}" --xgb-min 0.3 --xgb-max 0.7 --param-name bin1 --outvar lepm_11 --create-shape-dc
python ZMuE_fit_mk2_wrapper_v1.py -o bin2 --fit-version 1 --skip-sgn-syst --zmm-file ${MCFILE} --bkg-file "${DATAFILE}" --xgb-min 0.7 --xgb-max 0.9 --param-name bin2 --outvar lepm_12 --create-shape-dc
python ZMuE_fit_mk2_wrapper_v1.py -o bin3 --fit-version 1 --skip-sgn-syst --zmm-file ${MCFILE} --bkg-file "${DATAFILE}" --xgb-min 0.9 --xgb-max 1.1 --param-name bin3 --outvar lepm_13 --create-shape-dc

mv *mk2*_bin?.png ${OUTDIR}
mv *mk2*_bin?.root ${OUTDIR}
# Copy the corresponding COMBINE cards
cp /eos/cms/store/group/phys_smp/ZLFV/datacards/zemu/*.txt ${OUTDIR}
# Copy the signal workspaces with the correct systematics
cp /eos/cms/store/group/phys_smp/ZLFV/datacards/zemu/workspace_mk2sgn_v1_*.root ${OUTDIR}
```

Running with fits to the MC:
```
OUTDIR="zemu_mc/"
[ ! -d ${OUTDIR} ] && mkdir ${OUTDIR}

MCFILE="pseudo_data_from_MC_v2_r0_ZmmR1.25_updateID.root"
python ZMuE_fit_mk2_wrapper_v1.py -o bin1 --fit-version 1 --skip-sgn-syst --zmm-file ${MCFILE} --bkg-file ${MCFILE} --xgb-min 0.3 --xgb-max 0.7 --param-name bin1 --outvar lepm_11 --create-shape-dc
python ZMuE_fit_mk2_wrapper_v1.py -o bin2 --fit-version 1 --skip-sgn-syst --zmm-file ${MCFILE} --bkg-file ${MCFILE} --xgb-min 0.7 --xgb-max 0.9 --param-name bin2 --outvar lepm_12 --create-shape-dc
python ZMuE_fit_mk2_wrapper_v1.py -o bin3 --fit-version 1 --skip-sgn-syst --zmm-file ${MCFILE} --bkg-file ${MCFILE} --xgb-min 0.9 --xgb-max 1.1 --param-name bin3 --outvar lepm_13 --create-shape-dc

mv *mk2*_bin?.png ${OUTDIR}
mv *mk2*_bin?.root ${OUTDIR}
# Copy the corresponding COMBINE cards
cp /eos/cms/store/group/phys_smp/ZLFV/datacards/zemu/*.txt ${OUTDIR}
# Copy the signal workspaces with the correct systematics
cp /eos/cms/store/group/phys_smp/ZLFV/datacards/zemu/workspace_mk2sgn_v1_*.root ${OUTDIR}
```

Running with the smoothed templates:
```
OUTDIR="zemu_embed/"
[ ! -d ${OUTDIR} ] && mkdir ${OUTDIR}

DATAFILE="/eos/cms/store/cmst3/user/gkaratha/ZmuE_forBDT_v7_tuples/BDT_outputs_v7/Meas_fullAndSF_bdt_v7_signal_mcRun1*.root"
FITARGS="--run-histo --histo-toy --add-pol-order 1 --add-exp-order 1 --add-plaw-order 1"
MCFILE="template_zemu_embed_v3_bin1.root"
python ZMuE_fit_mk2_wrapper_v1.py -o bin1 --fit-version 1 --skip-sgn-syst ${FITARGS} --zmm-file ${MCFILE} --bkg-file ${MCFILE} --xgb-min 0.3 --xgb-max 0.7 --param-name bin1 --outvar lepm_11 --create-shape-dc
MCFILE="template_zemu_embed_v3_bin2.root"
python ZMuE_fit_mk2_wrapper_v1.py -o bin2 --fit-version 1 --skip-sgn-syst ${FITARGS} --zmm-file ${MCFILE} --bkg-file ${MCFILE} --xgb-min 0.7 --xgb-max 0.9 --param-name bin2 --outvar lepm_12 --create-shape-dc
MCFILE="template_zemu_embed_v3_bin3.root"
python ZMuE_fit_mk2_wrapper_v1.py -o bin3 --fit-version 1 --skip-sgn-syst ${FITARGS} --zmm-file ${MCFILE} --bkg-file ${MCFILE} --xgb-min 0.9 --xgb-max 1.1 --param-name bin3 --outvar lepm_13 --create-shape-dc
mv *mk2*_bin?.png ${OUTDIR}
mv *mk2*_bin?.root ${OUTDIR}
# Copy the corresponding COMBINE cards
cp /eos/cms/store/group/phys_smp/ZLFV/datacards/zemu/*.txt ${OUTDIR}
# Copy the signal workspaces with the correct systematics
cp /eos/cms/store/group/phys_smp/ZLFV/datacards/zemu/workspace_mk2sgn_v1_*.root ${OUTDIR}


python CreatePseudoData_from_EmbedTH1.py
mv ctmpl_embed*.png ${OUTDIR}/
```

## Z prime scan

The Z prime scan searches for a narrow resonance in the e-mu data using the Z->e+mu analysis framework and BDT.

### To do list

- Account for the effects of different running periods: see [evaluate_years_effect.py](tools/evaluate_years_effects.py) for a starting point, added linear correction to [signal_model.py](signal_model.py)
- Add uncertainty for efficiency turn-on (perhaps maximum difference from the fit)
- Validate systematic uncertainty values: see [eval_zprime_unc.py](tools/eval_zprime_unc.py) for a starting point
- Code updates: Merge ScanMuE_fit_wrapper_v2.py changes, add ntupling code to repo, add BDT score code to repo
- Update documentation: AN and paper have a preliminary discussion added
- Find a relevant theory model to compare to (e.g. EXO-19-014 Z' model with B = 1e-4 or something)
- Prepare for mini-preapproval presentation on Oct. 1st

### Create signal and background PDFs for a range of mass points for a single BDT category
```
MINMASS=110
MAXMASS=500
MINBDT=0.70
MAXBDT=1.00
NAME=bdt_0d7_1d0_v01
time python ScanMuE_fit_wrapper_v2.py -o ${NAME} --scan-min ${MINMASS} --scan-max ${MAXMASS} --xgb-min ${MINBDT} --xgb-max ${MAXBDT} --log-files --scan-step 1 [-j N]
ls -l datacards/${NAME}/datacard_zprime_${NAME}_mass-*_mp*.txt | head -n 2
ls -l WorkspaceScanSGN/workspace_scansgn_v2_${NAME}_mp*.root | head -n 2
ls -l WorkspaceScanBKG/workspace_scanbkg_v2_${NAME}_mp*.root | head -n 2
#figures are printed to: figures/${NAME}/ (signal) and figures/${NAME}_mp*/ (background/data)
```

### Create standard BDT categories
```
MINMASS=110
MAXMASS=500
NAME=v01
time ./make_scan_cards.sh --min-mass ${MINMASS} --max-mass ${MAXMASS} --tag ${NAME} --scan-arg "--scan-step 1 --log-files"
ls -l datacards/bdt_${NAME}/datacard_zprime_${NAME}_mass-*_mp*.txt | head -n 2
```

### Scan the mass points, evaluating signal rates and upper limits
```
NAME=v01
time python perform_scan.py -o bdt_${NAME} [--asimov] [--unblind] [--smooth-expected]
ls -l figures/scan_bdt_${NAME}[_asimov]/*.png
```

### Generate a toy dataset and run a scan over the toy dataset
This assumes the nominal scan is already processed on data with corresponding COMBINE cards available.

```
# Fit the data in the entire mass range for toy generation (only needed once, all toys can be generated from this initial fit)
python ScanMuE_fit_wrapper_v2.py -o bdt_0d3_0d7_LEE --full-mass --scan-min 300 --scan-max 300.1 --scan-step 1 --xgb-min 0.30 --xgb-max 0.70 --param-name bin1 --component bkg
python ScanMuE_fit_wrapper_v2.py -o bdt_0d7_1d0_LEE --full-mass --scan-min 300 --scan-max 300.1 --scan-step 1 --xgb-min 0.70 --xgb-max 1.01 --param-name bin2 --component bkg

# Generate a single toy dataset for each BDT region
python create_toy.py --fit-file WorkspaceScanBKG/workspace_scanbkg_v2_bdt_0d3_0d7_LEE_mp0.root -o toy_0d3_0d7 --toy 2 --param bin1 --seed 90
python create_toy.py --fit-file WorkspaceScanBKG/workspace_scanbkg_v2_bdt_0d7_1d0_LEE_mp0.root -o toy_0d7_1d0 --toy 2 --param bin2 --seed 90

# Create toy scan COMBINE cards from an existing scan dataset
TOY="2" #Toy number
MAIN="bdt_v01" #Main directory used for cloning
./clone_cards_for_toy.sh datacards/${MAIN}/ datacards/${MAIN}_toy_${TOY}/ WorkspaceScanTOY/toy_file_toy_0d3_0d7_${TOY}.root WorkspaceScanTOY/toy_file_toy_0d7_1d0_${TOY}.root

# Run the scan over the toy datacards
time python perform_scan.py -o ${MAIN}_toy_${TOY} --unblind --smooth-expected
```

### Validation studies

Validation studies for the entire scan are performed using the [perform_scan_validation.py](perform_scan_validation.py) tool.

```
# Perform a bias test for each mass point in the scan:
NTOYS=1000
NAME="bdt_v01"
time python perform_scan_validation.py -o ${NAME} -t ${NTOYS} --test bias
ls -l figures/val_${NAME}_bias/pulls.png
ls -l figures/val_${NAME}_bias/${NAME}_mp*_bias.png | head -n 5
```

Underlying tools to perform a single validation check:
- [bemu_bias.sh](tests/bemu_bias.sh): Perform a self-closure bias test on a single data card.
- [bemu_gen_fit_test.sh](tests/bemu_gen_fit_test.sh): Perform a closure test using (potentially) different generation and fit data cards.
- [mc_template_bias.sh](tests/mc_template_bias.sh): Perform a bias test generating with a fit MC template and fitting with an envelope data card
- [do_goodness_of_fit.sh](tests/do_goodness_of_fit.sh): Perform a goodness-of-fit test.
- [impacts.sh](tests/impacts.sh): Evaluate impacts (including condor configurations if requested).

#### Individual validation tests

Test the interpolation by comparing a mass point fit to the interpolation omitting this point:
```
[ ! -d biases ] && mkdir biases
GENSIGNAL=5
for MPOINT in "125" "300" "500"
do
  ./make_scan_cards.sh --min-mass ${MPOINT} --max-mass ${MPOINT}.1 --tag mass_${MPOINT}_mc --scan-arg "--use-mc-signal"
  ./make_scan_cards.sh --min-mass ${MPOINT} --max-mass ${MPOINT}.1 --tag mass_${MPOINT}_int --scan-arg "--skip-mc-point ${MPOINT}"
  CARD1="datacards/bdt_mass_${MPOINT}_mc/datacard_zprime_mass_${MPOINT}_mc_mass-${MPOINT}.0_mp0.txt"
  CARD2="datacards/bdt_mass_${MPOINT}_int/datacard_zprime_mass_${MPOINT}_int_mass-${MPOINT}.0_mp0.txt"
  time tests/bemu_gen_fit_test.sh --card_1 ${CARD1} --card_2 ${CARD2} -t 1000 -g 1000 -r 20 --gensignal "${GENSIGNAL}" --name mc_vs_int_${MPOINT} --skipfullplots
  mv bias_mc_vs_int_${MPOINT}_test.png biases/
  GENSIGNAL=$((${GENSIGNAL}-2))
done
```

### Additional useful tools/studies

- [plot_signal_model.py](tools/plot_signal_model.py): Plot the signal model interpolation, including the line shape and the overall efficiency.
- [plot_zprime_bdt.py](tools/plot_zprime_bdt.py): Plot the Z prime BDT score distribution and CDF as a function of mass.
- [eval_zprime_unc.py](tools/eval_zprime_unc.py): Evaluate standard uncertainties on the Z prime signal efficiency.
- [evaluate_years_effect.py](tools/evaluate_years_effect.py): Evaluate the difference in signal yield by year vs. mass
- [get_ngen.py](tools/get_ngen.py): Retrieve the signal ntuple normalization info.
- [test_mc_flatness.py](tools/test_mc_flatness.py): Test how flat the MC backgrounds are expected to be
- [make_mc_templates.py](tools/make_mc_templates.py): Make an MC template data card/workspace using fits to the MC process histograms
- [compare_scans.py](tools/compare_scans.py): Plot the ratio of the expected limit from two (identical mass point) scans as a function of mass to compare their sensitivity.
- [optimize_binning.py](tools/optimize_binning.py): Rough optimization of a single BDT bin using cut-and-count limits without systematic uncertainties.
- [skim_ntuple.py](skim_ntuple.py): Create a sparse TTree from the skimmed NANOAOD ntuples (not the official version, not fully validated).
- [evaluate_bdt.py](evaluate_bdt.py): Evaluate the signal BDT and output a TTree with this score included.
