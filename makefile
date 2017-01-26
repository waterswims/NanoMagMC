run:
	cd Includes; make run f_run
	cp Includes/run ./
	cp Includes/f_run ./

clean:
	cd Includes; make clean
	-rm run
	-rm f_run
