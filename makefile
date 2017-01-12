run:
	cd Includes; make clean all run
	cp Includes/run ./
	# cp Includes/d_run ./

clean:
	-rm run
	# -rm d_run
