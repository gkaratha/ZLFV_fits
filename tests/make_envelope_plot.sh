#! /bin/bash

# Run combine scans with each individual background function, generating a plot of the envelope method

Help() {
    echo "Make an envelope plot for a single Z->emu fit category"
    echo " Usage: make_envelope_plot.sh <card> [options]"
    echo " --cat      : Category number for PDF index variable (default = 13)"
    echo " --maxpdfs  : Maximum PDF index (default = 2)"
    echo " --rrange   : r parameter range (default = 10)"
    echo " --rgen     : r for generated toy (default = 0)"
    echo " --index_gen: PDF index for generated toy"
    echo " --npoints  : N(scan points) (default = 100)"
    echo " --obs      : Use the observed dataset"
    echo " --plotonly : Only print plots"
    echo " --dryrun   : Only setup commands"
    echo " --tag      : Tag for output plot"
    echo " --seed     : Seed for generation"
}

CARD=""
CAT=13
MAXPDFS=2
RRANGE=10
RGEN=0
TAG=""
GENINDEX=""
NPOINTS=100
SEED=90
OBS=""
PLOTONLY=""
DRYRUN=""

iarg=1
while [ ${iarg} -le $# ]
do
    eval "var=\${${iarg}}"
    if [[ "${var}" == "--help" ]] || [[ "${var}" == "-h" ]]
    then
        Help
        exit
    elif [[ "${var}" == "--rrange" ]] || [[ "${var}" == "-r" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        RRANGE=${var}
    elif [[ "${var}" == "--rgen" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
    elif [[ "${var}" == "--npoints" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        NPOINTS=${var}
    elif [[ "${var}" == "--tag" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        TAG=${var}
    elif [[ "${var}" == "--maxpdfs" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        MAXPDFS=${var}
    elif [[ "${var}" == "--index_gen" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        GENINDEX=${var}
    elif [[ "${var}" == "--genarg" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        GENARG=${var}
    elif [[ "${var}" == "--cat" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        CAT=${var}
    elif [[ "${var}" == "--seed" ]] || [[ "${var}" == "-s" ]]
    then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        SEED=${var}
    elif [[ "${var}" == "--obs" ]]; then
        OBS="d"
    elif [[ "${var}" == "--plotonly" ]]; then
        PLOTONLY="d"
    elif [[ "${var}" == "--dryrun" ]]; then
        DRYRUN="d"
    elif [[ "${CARD}" != "" ]]
    then
        echo "Arguments aren't configured correctly!"
        Help
        exit
    else
        CARD=${var}
    fi
    iarg=$((iarg + 1))
done

if [[ "${CARD}" == "" ]]; then
    echo "No input Combine card/workspace provided!"
    exit 1
fi

if [[ "${CAT}" == "" ]]; then
    CAT=13
fi

if [[ "${MAXPDFS}" == "" ]]; then
    MAXPDFS=10
fi

RMIN=-${RRANGE}
RMAX=${RRANGE}
ARGRANGES=" --setParameterRanges r=${RMIN},${RMAX}"
ARGUMENTS="--algo grid --cminDefaultMinimizerStrategy 0 --saveNLL -m 125 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1"
ARGUMENTS="${ARGUMENTS} --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstant --X-rtd MINIMIZER_multiMin_maskConstraints"
ARGUMENTS="${ARGUMENTS} --cminApproxPreFitTolerance 0.1 --cminPreScan --cminPreFit 1 --X-rtd MINIMIZER_multiMin_maskChannels=2"
if [[ "${CARD}" != *"total"* ]]; then
    ARGUMENTS="${ARGUMENTS} --cminRunAllDiscreteCombinations"
fi
ARGUMENTS="${ARGUMENTS} --rMin ${RMIN} --rMax ${RMAX}"

#Generate a toy dataset if needed
if [[ "${PLOTONLY}" == "" ]]; then
    if [[ "${OBS}" == "" ]]; then
        if [[ "${DRYRUN}" == "" ]]; then
            combine -M GenerateOnly -d ${CARD} -t 1 -m 91 -s ${SEED} -n _env_${CAT} --saveToys
        fi
        TOYDATA=higgsCombine_env_${CAT}.GenerateOnly.mH91.${SEED}.root
        ls -l ${TOYDATA}
        TOYARG="--toysFile=${TOYDATA} -t -1"
        TOYTAG=".${SEED}"
    else
        TOYARG=""
        TOYTAG=""
    fi
fi

# Perform the fit for each envelope function
OUTPUTLIST="{"
for (( PDF=0; PDF<${MAXPDFS}; PDF++ )); do
    CATVAR="pdfindex_bin${CAT}"

    if [[ "${PLOTONLY}" == "" ]]; then
        echo combine -M MultiDimFit -d ${CARD} ${ARGUMENTS} ${ARGRANGES} --points ${NPOINTS} ${TOYARG} --setParameters ${CATVAR}=${PDF} --freezeParameters "${CATVAR}" -n "_env_${CAT}_cat_${PDF}"
        if [[ "${DRYRUN}" == "" ]] ; then
            combine -M MultiDimFit -d ${CARD} ${ARGUMENTS} ${ARGRANGES} --points ${NPOINTS} ${TOYARG} --setParameters ${CATVAR}=${PDF} --freezeParameters "${CATVAR}" -n "_env_${CAT}_cat_${PDF}"
        fi
    fi
    OUTFILE="higgsCombine_env_${CAT}_cat_${PDF}.MultiDimFit.mH125${TOYTAG}.root"
    if [ ! -f ${OUTFILE} ]; then
        echo "File ${OUTFILE} not found!"
    fi
    if [[ "${OUTPUTLIST}" == "{" ]]; then
        OUTPUTLIST="{\"${OUTFILE}\""
    else
        OUTPUTLIST="${OUTPUTLIST},\"${OUTFILE}\""
    fi
done

# Perform the total envelope fit
if [[ "${PLOTONLY}" == "" ]]; then
    echo combine -M MultiDimFit -d ${CARD} ${ARGUMENTS} ${ARGRANGES} --points ${NPOINTS} -n "_env_${CAT}_tot" ${TOYARG}
    if [[ "${DRYRUN}" == "" ]]; then
        combine -M MultiDimFit -d ${CARD} ${ARGUMENTS} ${ARGRANGES} --points ${NPOINTS} -n "_env_${CAT}_tot" ${TOYARG}
    fi
fi
OUTFILE="higgsCombine_env_${CAT}_tot.MultiDimFit.mH125${TOYTAG}.root"
if [ ! -f ${OUTFILE} ]; then
    echo "File ${OUTFILE} not found!"
    break
fi
if [[ "${OUTPUTLIST}" == "{" ]]; then
    OUTPUTLIST="{\"${OUTFILE}\""
else
    OUTPUTLIST="${OUTPUTLIST},\"${OUTFILE}\""
fi
OUTPUTLIST="${OUTPUTLIST}}"
echo ${OUTPUTLIST}

# Plot the results
if [[ "${DRYRUN}" == "" ]]; then
    if [[ "${OBS}" == "" ]]; then
        root.exe -q -b "${CMSSW_BASE}/src/ZLFV_fits/tools/plot_envelope.C(${CAT}, ${OUTPUTLIST}, \"${TAG}\", true)"
    else
        root.exe -q -b "${CMSSW_BASE}/src/ZLFV_fits/tools/plot_envelope.C(${CAT}, ${OUTPUTLIST}, \"${TAG}\", true, true)"
    fi
fi
