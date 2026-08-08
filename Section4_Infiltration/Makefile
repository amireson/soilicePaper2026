RUN_SCRIPTS := $(wildcard run_*.py)
RUN_OUTPUTS := $(patsubst %.py,%.dill,$(RUN_SCRIPTS))

all: $(RUN_OUTPUTS)
	python -B plotResults.py

run_%.dill: run_%.py
	python -B $< 
	mv out.dill $@

.PHONY: all clean

clean:
	rm -f *.png *.dill

