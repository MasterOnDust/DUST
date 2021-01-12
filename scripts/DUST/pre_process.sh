# Shell script for doing initial processing of FLEXPART/FLEXDUST model output
# This might take a while to run...
# Dependencies, DUST python package, process_flexpart.py script
#
# Author: Ove Haugvaldstad
shopt -s extglob
shopt -s globstar

DATAPATH='/projects/NS2806K/ovewh/tracing_the_winds/flexpart/FLEXPART_spring/emission_sensitivities/'

#SUBDIRS=('Conc' 'DryDep' 'WetDep')

SUBDIRS=('Conc')

ps=('20micron' '2micron')

locations=('SACOL' 'BADOE' 'YINCHUAN' 'LINGTAI' 'LUOCHUAN' 'LANTIAN' 'SHAPOTOU')


proc_output () {
    echo "process_flexpart.py $1 --op $2 --x0 73 --x1 115 --y0 30 --y1 50"
}

for dir in "${SUBDIRS[@]}"; do
    for l in "${locations[@]}"; do
        for size in "${ps[@]}"; do
            cd  ${DATAPATH}${dir}/${size}
            for year in {1999..2019}; do
                #ls 
                if [[ -d ./surface_sensitivity/gridtime_${year}0306-00_${year}0331-21 ]]; then
                    
                    cd surface_sensitivity/
                    for ncf in ./gridtime_${year}*/*.nc; do
                        echo $ncf
                    done
                    #ls ./gridtime_${year}*/*.nc
                    cd ..
                else
                    echo "HAHA ${year} ${size}"
                fi
            done
        done
    done
done



