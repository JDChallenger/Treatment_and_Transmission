#Improve all these names
#1. 
scenario1:
	g++ WHM_untreated.cpp RK4.cpp
	./a.out
	python InfectivityPlot_untreated.py

#2.
scenario2:
	date
	g++ WHM_treated.cpp RK4.cpp
	./a.out
	python InfectivityPlot_treated.py
	date

#3. cohort	
scenario3:
	date
	g++ WHM_treated_cohort.cpp RK4.cpp
	./a.out
	date

#4. 
scenario4:
	date
	g++ WHM_treated_ADH.cpp RK4.cpp
	./a.out
	python InfectivityPlot_ADH.py
	date	

#5. Slow!	
scenario5:
	date
	g++ WHM_treated_ADH_cohort.cpp RK4.cpp
	./a.out
	date

# Zip up the folder
tar:
	cd ..; tar -zcvf within_host_model.tar.gz WithinHostModel