unfrozenFlow.png: psi_CeliaMiller.csv psi_soiliceMiller.csv theta_soiliceMathias.csv theta_Mathias.csv plot.py
	python -B plot.py
	open unfrozenFlow.png

psi_soiliceMiller.csv: MillerProblem_soilice.py
	python -B MillerProblem_soilice.py

psi_CeliaMiller.csv: MillerProblem_Celia.py
	python -B MillerProblem_Celia.py

theta_soiliceMathias.csv: MathiasProblem.py
	python -B MathiasProblem.py

theta_Mathias.csv: HorizontalInfiltrationModel.m
	matlab -nodesktop -nojvm -nosplash -r "HorizontalInfiltrationModel;exit"

.PHONY: clean

clean:
	rm -f *.csv *.png
