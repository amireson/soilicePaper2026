thawBenchmarks.png: Lunardini_100a.csv plot.py StefanNum1.csv run_Stefan.py run_Lunardini.py analyticalSoln.py
	python -B plot.py
	open thawBenchmarks.png

Lunardini_100a.csv: run_Lunardini.py analyticalSoln.py
	python -B run_Lunardini.py

StefanNum1.csv: run_Stefan.py analyticalSoln.py
	python -B run_Stefan.py

.PHONY: clean

clean:
	rm -rf *.csv *.png

