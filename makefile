RUN_TYPES = run test
.PHONY: $(RUN_TYPES)

$(RUN_TYPES):
	make -C Includes $@
	cp Includes/$@ ./

.PHONY: clean
clean:
	make -C Includes clean
	-rm run
	-rm test
