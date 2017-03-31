source load_compilers.sh

make f_run

for FIELD in 0.18
do
	for SIZE in 20
	do
		cp INPUT.dat F_INPUT_"$SIZE"_"$FIELD"
		sed -i "s/SIZE .*/SIZE ${SIZE}/g" F_INPUT_"$SIZE"_"$FIELD"
		sed -i "s/MAGFIELD .*/MAGFIELD ${FIELD}/g" F_INPUT_"$SIZE"_"$FIELD"

		cp f_sub.sh f_sub_"$SIZE"_"$FIELD"
		sed -i "s/INPUT_.*_.*/INPUT_${SIZE}_${FIELD}/g" f_sub_"$SIZE"_"$FIELD"

		if [ "$SIZE" -eq "25" ]
		then
			sed -i "s/walltime=.*:00:00/walltime=02:00:00/g" f_sub_"$SIZE"_"$FIELD"
		fi

		if [ "$SIZE" -eq "50" ]
		then
			sed -i "s/walltime=.*:00:00/walltime=03:59:00/g" f_sub_"$SIZE"_"$FIELD"
		fi

		if [ "$SIZE" -eq "75" ]
		then
			sed -i "s/walltime=.*:00:00/walltime=08:00:00/g" f_sub_"$SIZE"_"$FIELD"
		fi

		if [ "$SIZE" -eq "100" ]
		then
			sed -i "s/walltime=.*:00:00/walltime=13:00:00/g" f_sub_"$SIZE"_"$FIELD"
		fi

		echo "Submitting size equal to ${SIZE} and field equal to ${FIELD}"
		qsub f_sub_"$SIZE"_"$FIELD"
		sleep 1
	done
done
