#! /bin/bash
# Testing Z->emu fit results
# Usage: ./test_zemu_limits.sh <card>

Help() {
    echo "Testing Z->emu fit results"
    echo "Usage: ./test_zemu_limits.sh <card>"
}

CARD=$1
RMIN=-5
RMAX=5

if [[ "${CARD}" == "" ]]; then
    Help
    exit
fi
if [ ! -f ${CARD} ]; then
    echo "No input card ${CARD} found"
    Help
    exit
fi

ARGUMENTS="--cminDefaultMinimizerStrategy 0 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
ARGUMENTS="${ARGUMENTS} --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstant --X-rtd MINIMIZER_multiMin_maskConstraints"
ARGUMENTS="${ARGUMENTS} --rMin ${RMIN} --rMax ${RMAX} --X-rtd MINIMIZER_multiMin_maskChannels=2"

# Debug: setting initial parameters (and some ranges):
PARAMS="--setParameters "
ARGRANGES=""

# Bin 1
if [[ "${CARD}" == *"bin1"* ]] || [[ "${CARD}" == *"total"* ]]; then
    PARAMS="${PARAMS},bkg_gspol1_pdf_bin1_b0=-0.00639766,bkg_gspol1_pdf_bin1_mu=55.3031,bkg_gspol1_pdf_bin1_wd=12.1883"
    PARAMS="${PARAMS},bkg_gspol2_pdf_bin1_b0=0.0696676,bkg_gspol2_pdf_bin1_b1=-0.000496288,bkg_gspol2_pdf_bin1_mu=57.9333,bkg_gspol2_pdf_bin1_wd=12.2978"
    PARAMS="${PARAMS},bkg_gsexp1_pdf_bin1_mu=69,bkg_gsexp1_pdf_bin1_wd=150,bkg_gsexp1_pdf_bin1_x0=-0.105078"
    PARAMS="${PARAMS},bkg_gsexp2_pdf_bin1_c0=0.77646,bkg_gsexp2_pdf_bin1_mu=69,bkg_gsexp2_pdf_bin1_wd=5.47864,bkg_gsexp2_pdf_bin1_x0=-0.047225,bkg_gsexp2_pdf_bin1_x1=0.0242386"
    PARAMS="${PARAMS},bkg_gsplaw1_pdf_bin1_a0=-1.20994,bkg_gsplaw1_pdf_bin1_mu=60.4097,bkg_gsplaw1_pdf_bin1_wd=10.993"
    PARAMS="${PARAMS},bkg_gsplaw2_pdf_bin1_a0=4.91931,bkg_gsplaw2_pdf_bin1_a1=-3.33746,bkg_gsplaw2_pdf_bin1_c0=0.0832076,bkg_gsplaw2_pdf_bin1_mu=69,bkg_gsplaw2_pdf_bin1_wd=5.57072"
    PARAMS="${PARAMS},bkg_cheb3_pdf_bin1_0=-0.744645,bkg_cheb3_pdf_bin1_1=0.31931,bkg_cheb3_pdf_bin1_2=-0.0970132"
    PARAMS="${PARAMS},bkg_cheb4_pdf_bin1_0=-0.744743,bkg_cheb4_pdf_bin1_1=0.320949,bkg_cheb4_pdf_bin1_2=-0.0979635,bkg_cheb4_pdf_bin1_3=0.00261585"

    ARGRANGES="${ARGRANGES},bkg_cheb5_pdf_bin1_0=-3,3,bkg_cheb5_pdf_bin1_1=-2,2,bkg_cheb5_pdf_bin1_2=-1,1,bkg_cheb5_pdf_bin1_3=-0.5,0.5,bkg_cheb5_pdf_bin1_4=-0.2,0.2"
fi

# Bin 2
if [[ "${CARD}" == *"bin2"* ]] || [[ "${CARD}" == *"total"* ]]; then
    PARAMS="${PARAMS},bkg_gspol1_pdf_bin2_b0=-0.0061644,bkg_gspol1_pdf_bin2_mu=64.5875,bkg_gspol1_pdf_bin2_wd=9.03786"
    PARAMS="${PARAMS},bkg_gspol2_pdf_bin2_b0=0.0671617,bkg_gspol2_pdf_bin2_b1=-0.000521754,bkg_gspol2_pdf_bin2_mu=64.1955,bkg_gspol2_pdf_bin2_wd=9.26246"
    PARAMS="${PARAMS},bkg_gsexp1_pdf_bin2_mu=64.8665,bkg_gsexp1_pdf_bin2_wd=8.8785,bkg_gsexp1_pdf_bin2_x0=-0.0168155"
    PARAMS="${PARAMS},bkg_gsexp2_pdf_bin2_c0=0.704532,bkg_gsexp2_pdf_bin2_mu=68.9973,bkg_gsexp2_pdf_bin2_wd=7.52105,bkg_gsexp2_pdf_bin2_x0=-0.0162419,bkg_gsexp2_pdf_bin2_x1=-0.181696"
    PARAMS="${PARAMS},bkg_gsplaw2_pdf_bin2_a0=-12.3945,bkg_gsplaw2_pdf_bin2_a1=-1.41593,bkg_gsplaw2_pdf_bin2_c0=0.329129,bkg_gsplaw2_pdf_bin2_mu=68.9977,bkg_gsplaw2_pdf_bin2_wd=7.43669"
    PARAMS="${PARAMS},bkg_gsplaw3_pdf_bin2_a0=-14.7373,bkg_gsplaw3_pdf_bin2_a1=-6.1344,bkg_gsplaw3_pdf_bin2_a2=-1.33806,bkg_gsplaw3_pdf_bin2_c0=0.222879,bkg_gsplaw3_pdf_bin2_c1=0.16461,bkg_gsplaw3_pdf_bin2_mu=68.9999,bkg_gsplaw3_pdf_bin2_wd=7.44786"
    PARAMS="${PARAMS},bkg_cheb5_pdf_bin2_0=-1.22725,bkg_cheb5_pdf_bin2_1=0.639567,bkg_cheb5_pdf_bin2_2=-0.195609,bkg_cheb5_pdf_bin2_3=0.00271248,bkg_cheb5_pdf_bin2_4=0.0307185"
    ARGRANGES="${ARGRANGES},bkg_cheb5_pdf_bin2_0=-3,3,bkg_cheb5_pdf_bin2_1=-2,2,bkg_cheb5_pdf_bin2_2=-1,1,bkg_cheb5_pdf_bin2_3=-0.5,0.5,bkg_cheb5_pdf_bin2_2=-0.2,0.2"
fi

# Bin 3
if [[ "${CARD}" == *"bin3"* ]] || [[ "${CARD}" == *"total"* ]]; then
    PARAMS="${PARAMS},bkg_gspol1_pdf_bin3_b0=-0.008133,bkg_gspol1_pdf_bin3_mu=65.5644,bkg_gspol1_pdf_bin3_wd=8.84455"
    PARAMS="${PARAMS},bkg_gspol2_pdf_bin3_b0=0.059798,bkg_gspol2_pdf_bin3_b1=-0.000567103,bkg_gspol2_pdf_bin3_mu=65.3745,bkg_gspol2_pdf_bin3_wd=8.96678"
    PARAMS="${PARAMS},bkg_gsexp1_pdf_bin3_mu=66.2777,bkg_gsexp1_pdf_bin3_wd=8.43507,bkg_gsexp1_pdf_bin3_x0=-0.0492584"
    PARAMS="${PARAMS},bkg_gsexp2_pdf_bin3_c0=0.11467,bkg_gsexp2_pdf_bin3_mu=66.3163,bkg_gsexp2_pdf_bin3_wd=8.41403,bkg_gsexp2_pdf_bin3_x0=-0.0508909,bkg_gsexp2_pdf_bin3_x1=-0.0494316"
    PARAMS="${PARAMS},bkg_gsplaw1_pdf_bin3_a0=-5.03222,bkg_gsplaw1_pdf_bin3_mu=66.7112,bkg_gsplaw1_pdf_bin3_wd=8.25766"
    PARAMS="${PARAMS},bkg_gsplaw2_pdf_bin3_a0=-14.2918,bkg_gsplaw2_pdf_bin3_a1=-5.08336,bkg_gsplaw2_pdf_bin3_c0=0.00223443,bkg_gsplaw2_pdf_bin3_mu=66.7952,bkg_gsplaw2_pdf_bin3_wd=8.21747"
    PARAMS="${PARAMS},bkg_cheb5_pdf_bin3_0=-1.49644,bkg_cheb5_pdf_bin3_1=0.772098,bkg_cheb5_pdf_bin3_2=-0.222063,bkg_cheb5_pdf_bin3_3=-0.0116984,bkg_cheb5_pdf_bin3_4=0.0391811"

    ARGRANGES="${ARGRANGES},bkg_cheb5_pdf_bin3_0=-3,3,bkg_cheb5_pdf_bin3_1=-2,2,bkg_cheb5_pdf_bin3_2=-1,1,bkg_cheb5_pdf_bin3_3=-0.5,0.5,bkg_cheb5_pdf_bin3_4=-0.2,0.2"
fi

NAME="test"
if [[ "${CARD}" == *"bin1"* ]]; then
    NAME="bin1"
elif [[ "${CARD}" == *"bin2"* ]]; then
    NAME="bin2"
elif [[ "${CARD}" == *"bin3"* ]]; then
    NAME="bin3"
elif [[ "${CARD}" == *"total"* ]]; then
    NAME="total"
fi

# Flag to do parameter scans, not really working well right now
DOSCAN=""
# Flag to add discrete index cycling flag: Doesn't work well with the combined fit
DOCYCLE="d"
# Flag to add pre-fit arguments to all fits
DOPREFIT="d"

if [[ "${PREFIT}" != "" ]]; then
    ARGUMENTS="${ARGUMENTS} --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1"
fi
if [[ "${DOCYCLE}" != "" ]]; then
    ARGUMENTS="${ARGUMENTS} --cminRunAllDiscreteCombinations"
fi

echo "Performing original fits"
combine -d ${CARD} ${ARGUMENTS}
combine -d ${CARD} -M FitDiagnostics ${ARGUMENTS}
if [[ "${DOSCAN}" != "" ]]; then
    combine -d ${CARD} -M MultiDimFit ${ARGUMENTS} --saveNLL --setParameterRanges r=-1,1 -n .original_${NAME} --algo grid --points 10
fi

exit

echo "Performing updated fits"
combine -d ${CARD} ${ARGUMENTS} ${PARAMS} --setParameterRanges ${ARGRANGES}  --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1
combine -d ${CARD} -M FitDiagnostics ${ARGUMENTS} ${PARAMS} --setParameterRanges ${ARGRANGES}  --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1
if [[ "${DOSCAN}" != "" ]]; then
    combine -d ${CARD} -M MultiDimFit ${ARGUMENTS} ${PARAMS} --saveNLL --setParameterRanges r=-1,1 -n .updated_${NAME} --algo grid --points 10
fi

# Compare the parameter scan results
if [[ "${DOSCAN}" != "" ]]; then
    python ../../../../CombineHarvester/CombineTools/scripts/plot1DScan.py higgsCombine.original_${NAME}.MultiDimFit.mH120.root --others higgsCombine.updated_${NAME}.MultiDimFit.mH120.root:Updated:2 -o compare_scan_${NAME}
    python ../../../../CombineHarvester/CombineTools/scripts/plot1DScan.py higgsCombine.updated_${NAME}.MultiDimFit.mH120.root --others higgsCombine.original_${NAME}.MultiDimFit.mH120.root:Original:2 -o compare_scan_${NAME}
fi
