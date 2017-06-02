module swap PrgEnv-cray PrgEnv-intel
make clean run

for FIELD in 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3
do
	cp ARC_INPUT.dat ARC_INPUT_"$FIELD"
	sed -i "s/MAGFIELD .*/MAGFIELD ${FIELD}/g" ARC_INPUT_"$FIELD"
done
