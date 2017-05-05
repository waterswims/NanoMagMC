source load_compilers.sh

make run

for FIELD in -0.005 -0.00475 -0.0045 -0.00425 -0.004 -0.00375 -0.0035 -0.00325 -0.003 -0.00275 -0.0025 -0.00225 -0.002 -0.00175 -0.0015 -0.00125 -0.001 -0.00075 -0.0005 -0.00025 0.00025 0.0005 0.00075 0.001 0.00125 0.0015 0.00175 0.002 0.00225 0.0025 0.00275 0.003 0.00325 0.0035 0.00375 0.004 0.00425 0.0045 0.00475 0.005
	for SIZE in 25 50
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
