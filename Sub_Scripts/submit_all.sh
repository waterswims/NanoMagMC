source load_compilers.sh

for SIZE in 20 25 30
do
    sed "s/SIZE .*/SIZE ${SIZE}/g" < INPUT.dat > INPUT_${SIZE}_X_X
    for PER in 0 1
    do
        sed "s/ISPERIO .*/ISPERIO ${PER}/g" < INPUT_${SIZE}_X_X > INPUT_${SIZE}_${PER}_X
        for PROT in 1 2 4
        do
            sed "s/PROTOCOL .*/PROTOCOL ${PROT}/g" < INPUT_${SIZE}_${PER}_X > INPUT_${SIZE}_${PER}_${PROT}
            sed "s/run INPUT.dat/run INPUT_${SIZE}_${PER}_${PROT}/g" < Sub_Scripts/iridis_sub.sh > sub_${SIZE}_${PER}_${PROT}.sh
            qsub sub_${SIZE}_${PER}_${PROT}.sh
        done
        rm INPUT_${SIZE}_${PER}_X
    done
    rm INPUT_${SIZE}_X_X
done
