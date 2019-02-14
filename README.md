# spatialGut
The gut microbiome is a heterogeneous group of microbes spatially distributed along various sections of the intestines and across the mucosa and lumen in each section. 
This program simulates the varying microbiota composition in different sites along the intestinal tract (e.g., proximal, middle and distal small intestine; cecum; proximal, milddle and distal large intestine and in each site the differential microbiota composition in the lumen and on the mucus layer. Given the initial nutrients available, the total microbial biomass on the mucus layer and the oxygen available in each section of the intestines, the spatially differential microbiota composition profile is simulated using a dynamic framework connecting the luminal and mucosal microbiota in each intestinal section as well as connecting consecutive sections embedded with microbiota metabolic models to predict microbial growth and metabolism in each time step.

__Functions__

`spatialGut.m` - the main function to simulate the spatially differential microbiota

`simSpatialGutExample.m` - an example script to setup the parameters for `spatialGut` and simulate using random uptake bounds

`plotSpatialGutResults.m` - to retrieve and plot the data from the simulations. It reproduces Figure 2 in ref. 1 by default



__Data files__

`simulationData.mat` - matlab file containing the metabolic model and nutrients available for simulations

`experimentalData.mat` - matlab file containing the experimental data in ref. 2 for comparison



__References__

1. Chan SHJ, Senftle MN, Friedman ES, Wu GD, Maranas CD. Predicting the longitudinally and radially varying gut microbiota composition using multi-scale microbial metabolic modeling. Submitted.

2. Friedman ES, Bittinger K, Esipova T V., Hou L, Chau L, Jiang J, et al. Microbes vs. chemistry in the origin of the anaerobic gut lumen. Proc Natl Acad Sci. 2018;201718635. 
