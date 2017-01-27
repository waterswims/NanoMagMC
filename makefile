run:
	cd Includes; make run f_run
	cp Includes/run ./
	cp Includes/f_run ./

test:
	cd Includes; make test
	cp Includes/test ./

clean:
	cd Includes; make clean
	-rm run
	-rm f_run
	-rm test
