source load_compilers.sh

make clean run

for SIZE in 5 15
do
	# sed -i.bu "s/int size(.*);/int size($SIZE);/g" Includes/main.cpp
	# make clean run
	# cp run run_"$SIZE"
	# sed -i.bu "s/run_.*/run_${SIZE}/g" job_sub.sh
	# cp job_sub.sh job_sub_"$SIZE"
	# qsub job_sub_"$SIZE"
	# sleep 2

	cp INPUT.dat INPUT_"$SIZE"
	sed -i.bu "s/SIZE .*/SIZE ${SIZE}/g" INPUT_"$SIZE"
	cp job_sub.sh job_sub_"$SIZE"
	sed -i.bu "s/INPUT_.*/INPUT_${SIZE}/g" job_sub_"$SIZE"
	echo "Submitting size equal to $SIZE"
	qsub job_sub_"$SIZE"
	sleep 1
done
