#!/bin/bash

#source setup.sh

OUTDIR=output
mkdir -p ${OUTDIR}

GEOMETRY_FILE=nd280geometry.root
FLUX_FILE=/cygdrive/c/T2K/beam/BeamMC/data/11a/nd34/nu.ing_flukain.0.root
NEUTRATE_FILE=${OUTDIR}/nu.ing_flukain.0_neutrate.root

VOLUMES="+ingrid"
NDINDEX=3

FILEROOT_NAME=ingrid_test_neut
POT=0.0001 

NPROCS=8
iproc=0
for (( file=0; file<${NPROCS}; file++ ))
do
    
    ${NEUT_ROOT}/src/neutgeom/./genev -w 1 -g ${GEOMETRY_FILE} -v ${VOLUMES} -j ${FLUX_FILE} -i ${NEUTRATE_FILE} -o ${OUTDIR}/${FILEROOT_NAME}_${file}.root -e ${POT} -d ${NDINDEX} -r $RANDOM > /dev/null &

iproc=$(( $iproc + 1 ))

done
