# Treatment and Transmission
This code accompanies article 'How delayed treatment and non-adherent treatment contribute to onward transmission of malaria: a modelling study', which will be published in BMJ Global Health. 

The code provided here allows the user to explore five modelling scenarios:
1. An untreated malaria infection, tracking both asexual and sexual parasite densities
   - Running WHM_untreated.cpp simulates the infection, and also calculates the overall infectivity of the malaria episode. This can be visualised with the Python script InfectivityPlot_untreated.py.

2. A single malaria infection treated with Artemether-Lumefantrine (AL)
   - Running WHM_treated.cpp simulates the infection and the impact of drug treatment. The malaria infection and the pharmacokinetic model can be visualised with the Python script InfectivityPlot_treated.py.

3. A cohort of treated patients, calculating the average infectivity of the cohort
   - Running WHM_treated_cohort.cpp simulates a cohort of patients (default size is 500). This code could be modified to examine e.g. the impact of delaying treatment or underdosing.

4. A single treated infection, using dose timings from XYZ.
   - Running WHM_treated_ADH.cpp simulates the infection and the impact of drug treatment. The malaria infection and the pharmacokinetic model can be visualised with the Python script InfectivityPlot_ADH.py.

5. Running the model for all cohort of treated patients, using all adherence profiles in the adherence data
   - Running WHM_treated_ADH_cohort.cpp simulates a treated cohort using the Tanzanian data to determine the doses taken and their timings. Note that the code avoids simulating certain adherence profiles, where overdosing is detected (see details below).

Here is a description of the auxilliary files in the project:

* parameters.h (Header file, storing parameters common to the C++ files)
* RK4.cpp (functions for the pharmacokinetic and pharmacodynamic models)
* Bradley_Infectivity.dat numerical values for the conversation of a gametocyte density to a measure of infectivity (the probability that malaria is transmitted to a feeding mosquito). These values were generated from the study by Bradley et al. https://doi.org/10.7554/eLife.34463.

* makefile (directive which, for Linux or Mac, compiles and runs the C++ model for the five scenarios listed above. In the terminal, navigate to the folder containing the model files. The command **make scenario1** runs scenario 1, the command **make scenario2** runs scenario 2 and so on.

To run the model that utilising the real-world adherence data (described in Bruxvoort et al. https://doi.org/10.1371/journal.pone.0134275) you will need to download the data from [this repository](http://actc.lshtm.ac.uk). Then remove the column headings and save as a text file. We did not use data from all 659 patients: as explained in the article, we removed some individuals who had taken too many pills in one of their doses. This was done because our model does not take toxicity effects into consideration. However, we include all the data here, so that the user can make their own decision for how to deal with this issue.