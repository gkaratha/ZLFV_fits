#! /bin/bash
# Create A' prime scan COMBINE datacards

Help() {
    echo "Create Z' prime scan COMBINE datacards"
    echo "Options:"
    echo " --min-mass       : Minimum mass to scan (default = 110)"
    echo " --max-mass       : Maximum mass to scan (default = 500)"
    echo " --scan-arg       : Arguments to pass to scan python script"
    echo " --tag            : Tag for output results (default = v01)"
    echo " --skip-fits      : Skip fits and initial datacard creation"
    echo " --dry-run        : Don't execute commands"
    echo " --help      (-h ): Print this information"
}

TAG="v01"
SIGNAL=""
ARG=""
MINMASS="110"
MAXMASS="500"
SKIPFITS=""
DRYRUN=""

iarg=1
while [ ${iarg} -le $# ]
do
    eval "var=\${${iarg}}"
    if [[ "${var}" == "--help" ]] || [[ "${var}" == "-h" ]]; then
        Help
        exit
    elif [[ "${var}" == "--scan-arg" ]]; then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        ARG=${var}
    elif [[ "${var}" == "--tag" ]]; then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        TAG=${var}
    elif [[ "${var}" == "--min-mass" ]]; then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        MINMASS=${var}
    elif [[ "${var}" == "--max-mass" ]]; then
        iarg=$((iarg + 1))
        eval "var=\${${iarg}}"
        MAXMASS=${var}
    elif [[ "${var}" == "--skip-fits" ]]; then
        SKIPFITS="d"
    elif [[ "${var}" == "--dry-run" ]]; then
        DRYRUN="d"
    else
        echo "Arguments aren't configured correctly! (at ${var})"
        Help
        exit
    fi
    iarg=$((iarg + 1))
done

if [[ "${DRYRUN}" != "" ]]; then
    ARG="${ARG} --dry-run"
fi
# Create the standard BDT score-defined regions
if [[ "${SKIPFITS}" == "" ]]; then
    [ ! -d log ] && mkdir log
    echo "Running fits and creating datacards..."
    python ScanMuE_fit_wrapper_v2.py -o bdt_0d3_0d7_${TAG} --scan-min ${MINMASS} --scan-max ${MAXMASS} --xgb-min 0.30 --xgb-max 0.70 --param-name bin1 ${ARG}
    python ScanMuE_fit_wrapper_v2.py -o bdt_0d7_1d0_${TAG} --scan-min ${MINMASS} --scan-max ${MAXMASS} --xgb-min 0.70 --xgb-max 1.00 --param-name bin2 ${ARG}
fi

# make a combined directory
DIR="datacards/bdt_${TAG}/"
[ ! -f ${DIR} ] && mkdir -p ${DIR}

if [[ "${DRYRUN}" != "" ]]; then
    exit
fi

# Create combined datacards
for CARD in `ls -d datacards/bdt_0d3_0d7_${TAG}/*.txt`
do
    if [[ "${CARD}" != *"datacard_"* ]]; then
        continue
    fi
    MASSPOINT=`echo ${CARD} | sed "s|datacards/bdt_0d3_0d7_${TAG}/datacard_zprime_bdt_0d3_0d7_${TAG}_||" | sed 's/.txt//'`
    CARD2=datacards/bdt_0d7_1d0_${TAG}/datacard_zprime_bdt_0d7_1d0_${TAG}_${MASSPOINT}.txt

    cp ${CARD}  ${DIR}
    cp ${CARD2} ${DIR}
    cd ${DIR}
    [ ! -e WorkspaceScanBKG ] && ln -s ../../WorkspaceScanBKG WorkspaceScanBKG
    [ ! -e WorkspaceScanSGN ] && ln -s ../../WorkspaceScanSGN WorkspaceScanSGN

    CARD1=datacard_zprime_bdt_0d3_0d7_${TAG}_${MASSPOINT}.txt
    CARD2=datacard_zprime_bdt_0d7_1d0_${TAG}_${MASSPOINT}.txt
    CARDOUT=datacard_zprime_${TAG}_${MASSPOINT}.txt
    echo "Creating merged datacard ${CARDOUT}..."
    echo '# -*- mode: tcl -*-' >| ${CARDOUT}
    echo '#Merged Z prime search COMBINE datacard' >> ${CARDOUT}
    cat ${CARD1} | grep 'Z prime mass' >> ${CARDOUT}
    echo "" >> ${CARDOUT}
    combineCards.py bin_1=${CARD1} bin_2=${CARD2} >> ${CARDOUT}
    cd ../..
done
