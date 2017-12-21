#!/bin/sh

KRN=5
ESPMAP=3
CALIB=24
THRESH=0.002
NITER=15

echo "----    Arguments   ----"
echo " $# arguments : $@"

while getopts "r:k:m:w:i:d" opt; do
	case $opt in
	r)
		CALIB=$OPTARG
                echo "Calibration size     :${CALIB}"
	;;
	k)
		KRN=$OPTARG
                echo "ESPIRiT kernel size  :${KRN}"
	;;
	m)
		ESPMAP=$OPTARG
                echo "ESPIRiT map number   :${ESPMAP}"
	;;
	w)
		THRESH=$OPTARG
                echo "PICS L1 regularization weight:${THRESH}"
        ;;
        i)
                NITER=$OPTARG
                echo "PICS iteration number:${NITER}"
	;;
	\?)
		echo "Invalid option       : -$OPTARG" >&2
	;;
	esac
done

shift $((OPTIND-1))

if [ $# -lt 1 ] ; then
        echo "Usage: $0 <kspace>"
        exit 1
fi

kspace=$(readlink -f "$1")

echo "---- Reconstruction ----"

echo "----Step 1: ESPIRiT Calibration               ----"
/home/amax/bart/bart ecalib -r${CALIB} -k${KRN} -m${ESPMAP} -S -t0.0005 -c0.9 ${kspace} maps
echo "----Step 2: L1-SENSE Reconstruciton on GPU    ----"
/home/amax/bart/bart pics -S -l1 -r${THRESH} -i${NITER} ${kspace} maps ims_soft_sense
echo "----Step 3: Fake kspace with Data consistency ----"
/home/amax/bart/bart fakeksp -r ims_soft_sense ${kspace} maps fakekspace
