run:
	make -C Includes run f_run
	cp Includes/run ./
	cp Includes/f_run ./

test:
	make -C Includes test
	cp Includes/test ./

clean:
	make -C Includes clean
	-rm run
	-rm f_run
	-rm test
