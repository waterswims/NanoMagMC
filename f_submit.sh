source load_compilers.sh

make run

for FIELD in -0.01 -0.0095 -0.009 -0.0085 -0.008 -0.0075 -0.007 -0.0065 -0.006 -0.0055 -0.005 -0.0045 -0.004 -0.0035 -0.003 -0.0025 -0.002 -0.0015 -0.001 -0.0005 0.0005 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0055 0.006 0.0065 0.007 0.0075 0.008 0.0085 0.009 0.0095 0.01
	for SIZE in 10
	do
		cp INPUT.dat F_INPUT_"$SIZE"_"$FIELD"
		sed -i "s/SIZE .*/SIZE ${SIZE}/g" F_INPUT_"$SIZE"_"$FIELD"
		sed -i "s/MAGFIELD .*/MAGFIELD ${FIELD}/g" F_INPUT_"$SIZE"_"$FIELD"

		cp f_sub.sh f_sub_"$SIZE"_"$FIELD"
		sed -i "s/INPUT_.*_.*/INPUT_${SIZE}_${FIELD}/g" f_sub_"$SIZE"_"$FIELD"

		if [ "$SIZE" -eq "25" ]
		then
			sed -i "s/walltime=.*:00:00/walltime=03:59:00/g" f_sub_"$SIZE"_"$FIELD"
		fi

		echo "Submitting size equal to ${SIZE} and field equal to ${FIELD}"
		qsub f_sub_"$SIZE"_"$FIELD"
		sleep 1
	done
done
