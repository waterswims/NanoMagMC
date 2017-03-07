RUN_TYPES = run f_run test run_new
.PHONY: $(RUN_TYPES)

$(RUN_TYPES):
	make -C Includes $@
	cp Includes/$@ ./

.PHONY: clean
clean:
	make -C Includes clean
	-rm run
	-rm f_run
	-rm test
	-rm run_new
