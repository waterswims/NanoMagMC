source load_compilers.sh

make run

for SIZE in 1 2 3 4
do
	# sed -i.bu "s/double asd(.*);/double asd($SIZE);/g" Includes/distro.cpp
	# make clean run
	# cp d_run d_run_"$SIZE"
	# sed -i.bu "s/d_run_.*/d_run_${SIZE}/g" dist_sub.sh
	# cp dist_sub.sh dist_sub_"$SIZE"
	# qsub dist_sub_"$SIZE"
	# sleep 2

	cp D_INPUT.dat D_INPUT_"$SIZE"
	sed -i.bu "s/SIZEDEV .*/SIZEDEV ${SIZE}/g" D_INPUT_"$SIZE"
	cp dist_sub.sh dist_sub_"$SIZE"
	sed -i.bu "s/D_INPUT_.*/D_INPUT_${SIZE}/g" dist_sub_"$SIZE"
	echo "Submitting deviation equal to $SIZE"
	qsub dist_sub_"$SIZE"
	sleep 1
done
