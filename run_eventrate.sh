#!/bin/bash

#source setup.sh

OUTDIR=output
mkdir -p ${OUTDIR}

GEOMETRY_FILE=nd280geometry.root
FLUX_FILE=/cygdrive/c/T2K/beam/BeamMC/data/11a/nd34/nu.ing_flukain.0.root
NEUTRATE_FILE=${OUTDIR}/nu.ing_flukain.0_neutrate.root

VOLUMES="+ingrid"
NDINDEX=3
 
export RANFILE=../neutsmpl/random.tbl
   
./event_rate -g ${GEOMETRY_FILE} -v ${VOLUMES} -f ${FLUX_FILE} -o ${NEUTRATE_FILE} -d ${NDINDEX}
