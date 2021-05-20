#!/bin/bash
###############################################################################
# Shell script for doing initial processing of FLEXPART/FLEXDUST model output #
# This might take a while to run...                                           #
# Dependencies, DUST python package, process_flexpart.py script               #
#                                                                             #
# Author: Ove Haugvaldstad                                                    #
###############################################################################
shopt -s extglob
shopt -s globstar

export PATH=$PATH:/projects/NS2806K/ovewh/DUST/DUST_scripts/DUST

DATAPATH='/projects/NS2806K/ovewh/tracing_the_winds/flexpart/production_run/'
OUTPATH='/projects/NS2806K/ovewh/tracing_the_winds/flexpart/production_run/source_contribution'
PATHFLEXDUST='/projects/NS2806K/ovewh/tracing_the_winds/FLEXDUST_emission_flux/FLEXDUST1999_2019/'
SUBDIRS=('drydep' 'wetdep')

ps=('20micron' '2micron')


POSITIONAL=()
while [[ $# -gt 0 ]]
do
    arg="$1"
    case $arg in
        -op|--out_path)
        OUTPATH="$2"
        shift #past agrument
        shift #past value
        ;;
        -dp|--data_path)
        DATAPATH="$2"
        shift
        shift
        ;;
        -sd|--sub_dir)
        SUBDIRS="$2"
        shift
        shift
        ;;
    esac
done

echo "Processing data in ${DATAPATH}"
echo "Output stored at :${OUTPATH}"
echo "Sub dir : ${SUBDIRS}"


proc_output () {
    echo "processing: $1, creating output at $3"
    process_flexpart.py $1 $2 --op $3 --x0 73 --x1 115 --y0 30 --y1 50
}

N=2 #Number of processes to run the background at one time

[[ -d ${OUTPATH} ]] || mkdir ${OUTPATH}

echo ${SUBDIRS}

for dir in "${SUBDIRS[@]}"; do
echo ${dir}
[[ -d ${OUTPATH}/${dir}  ]] || mkdir ${OUTPATH}/${dir}  
   for size in "${ps[@]}"; do
       [[ -d ${OUTPATH}/${dir}/${size} ]] || mkdir ${OUTPATH}/${dir}/${size}
       
       cd  ${DATAPATH}${dir}/${size}
       for year in {1999..2019}; do
           [[ -d ${OUTPATH}/${dir}/${size}/${year} ]] || mkdir ${OUTPATH}/${dir}/${size}/${year}
           if [[ -d ./surface_sensitivity/${year} ]]; then
                    
               cd surface_sensitivity/
               for ncf in ./${year}/*.nc; do
                    ((i=i%N)); ((i++==0)) && wait
                    proc_output ${ncf} ${PATHFLEXDUST}${year}/ ${OUTPATH}/${dir}/${size}/${year} & 
                done
                cd ..
           else
                for ncf in ./${year}/**/output/*.nc; do
                    ((i=i%N)); ((i++==0)) && wait
                    proc_output ${ncf} ${PATHFLEXDUST}${year}/ ${OUTPATH}/${dir}/${size}/${year} &
                done
           fi
            
        done
    done
done



