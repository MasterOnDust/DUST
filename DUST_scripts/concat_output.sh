#!/bin/bash

##############################################################################
# This is for concatenating FLEXPART output from simulations setup according #
# to multiple receptor location in the RELEASE FILE, and new simulation      #
# every timestep.                                                            #
# Data which needs to concatenated into one file before subsequent analysis  #
# can be done.                                                               #
##############################################################################


shopt -s extglob
shopt -s globstar




concat_output () {
    N=4
    output_dir=$3
    [[ -d ${output_dir} ]] || mkdir ${output_dir} 

    year=$2
    temp_dir="$1$2"
    for month in {3..5}; do
        if [[ ${month} -eq 3 ]]; then
            for p in {0..6}; do
                ((i=i%N)); ((i++==0)) && wait
                python concat_output.py ${temp_dir}/${year}0${month}[06-09,10-31]*/output/*.nc --x0 73 \
                --x1 115 --y0 30 --y1 50 --op ${output_dir} --location $p --height 100& 
            done
        else 
            for p in {0..6}; do
                ((i=i%N)); ((i++==0)) && wait
                python concat_output.py ${temp_dir}/${year}0${month}[01-09,10-31]*/output/*.nc --x0 73 \
                --x1 115 --y0 30 --y1 50 --op ${output_dir} --location $p --height 100&
            done 
        fi
    done
}

TARGETDIR_CONC='/projects/NS2806K/ovewh/tracing_the_winds/flexpart/FLEXPART_spring/emission_sensitivities/Conc/2micron/'


OUTPUT_DIR_CONC='/projects/NS2806K/ovewh/tracing_the_winds/flexpart/FLEXPART_spring/emission_sensitivities/Conc/2micron/surface_sensitivities_test'
[[ -d ${OUTPUT_DIR_CONC} ]] || mkdir ${OUTPUT_DIR_CONC} 
# Years are in per timestep format
TARGETYEARS_CONC=(1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010
    2011 2012 2013 2014 2016 2017 2018 2019)



TARGETDIR_DRYDEP='/projects/NS2806K/ovewh/tracing_the_winds/flexpart/FLEXPART_spring/emission_sensitivities/DryDep/20micron/'
OUTPUT_DIR_DRYDEP='/projects/NS2806K/ovewh/tracing_the_winds/flexpart/FLEXPART_spring/emission_sensitivities/DryDep/20micron/surface_sensitivities_test'
TARGETYEARS_DRYDEP=2019

[[ -d ${OUTPUT_DIR_DRYDEP} ]] || mkdir ${OUTPUT_DIR_DRYDEP} 


#concat_output ${TARGETDIR_DRYDEP} ${TARGETYEARS_DRYDEP} \
#${OUTPUT_DIR_DRYDEP}/${TARGETYEARS_DRYDEP}

for year in ${TARGETYEARS_CONC[@]}; do
    concat_output ${TARGETDIR_CONC} ${year} ${OUTPUT_DIR_CONC}/${year}
    \ ${OUTPUT_DIR_CONC/${year}
}
done

