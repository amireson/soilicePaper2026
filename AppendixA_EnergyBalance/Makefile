RUN_SCRIPTS := $(wildcard run_*.py)
RUN_OUTPUTS := $(patsubst %.py,%.dill,$(RUN_SCRIPTS))

all: $(RUN_OUTPUTS)
	python -B reportResult.py

run_%.dill: run_%.py
	python -B $< 
	mv out.dill $@

.PHONY: all clean

clean:
	rm -f *.png *.dill

# all: Kplots.png KFuns.png
# 
# Kplots.png: run.dill $(RUN_OUTPUTS) plots.py
# 	python -B plots.py
# 
# KFuns.png: def_const.txt def_pars.txt Kfuns.py
# 	python -B Kfuns.py
# 
# run.dill: run.py
# 	python -B run.py
