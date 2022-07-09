# ImmunoHatch
Immuno-PCR Optimization Algorithm


<!-- TABLE OF CONTENTS -->
## Table of Contents

* [Introduction](#introduction)
  * [Built With](#built-with)
* [The Model](#the-model)
* [Performance](#performance)
* [qPCR Set Up and Data Collection](#qpcr-set-up-and-data-collection)
* [Code Organization](#code-organization)
* [Sample Results](#sample-results)
* [Acknowledgements](#acknowledgements)




<!-- Introduction -->
## Introduction

Immuno-PCRs take a staple proteomics assay, ELISA (Enzyme Linked Immunosorbent Assay), and convert the standard detection method from HRPS + TMB to qPCR. Usually, a properly developed Immuno-PCR will require 4-10x less reagent and yield 2x-1000x better sensitivity with no loss of specificity.

A sandwich immuno-PCR requires a capture antibody immobilized on a PCR plate as well as a detection antibody conjugated to a target DNA strand. Once a workflow is established, the main part of immuno-PCR development is optimizing the concentrations of the capture and detection antibodies to maximize signal while minimizing background. While a myriad of methods exist to determine optimal concentrations, most require a kind of guess-and-check schema. 

Outlined here is a mathematical model and computer software which do away with guess-and-check optimization resulting in faster development and higher quality Immuno-PCRs.

### Built With
R4.0.1
* [stats]
* [nlstools]
* [tidyverse]
* [scatterplot3d]
* [ggpubr]



<!-- The Model -->
## The Model

In general, the signal (Ct value when target antigen present), the background (Ct value when no target present), and the normalized final signal (dCq = background Ct – signal Ct) all follow a rough logistic curve. This model uses a modified Gompertz Equation:

delta Cq= αe^(-e^(βX))+γX+δe^(-e^(θY))+μY+ε  ; ε ~ N(0,σ^2,)

Where α and δ control the upper asymptote for the capture and detection antibodies respectively, β, θ are scaling parameters, γ and μ control the effect of false background signal from the capture and detection antibodies respectively, X and Y are respective concentrations of capture and detection antibodies, e is the natural number, and ε is an error term assumed normal with mean 0 and standard deviation σ^2. However, with resampling techniques (here jackknife) ε need not be assumed normal.


<!-- Performance -->
## Performance

The algorithm has been used to produce over 20 immuno-PCR kits at RayBiotech Life. While lacking formal model testing due to time constraints, it has out preformed 3 existing guess-and-check method immuno-PCRs reducing their LOD by an average of 2.5-fold, and it has matched conditions in 2 existing guess-and-check method immuno-PCRs. Further it thus far has underperformed zero guess-and-check developed immuno-PCRs. The model has only failed or produced undesirable results when low quality qPCR data is fed to it OR when the target outright fails to convert from ELISA format to immuno-PCR format. 

Use of this algorithm has reduced the average number of experiments per target from ~4-5 to ~2-3 and saved subsequent development time and reagent, all the while increasing immuno-PCR kit quality. Further, the algorithm has proven itself useful when new lots require revalidation since new lots oftentimes require extensive redevelopment. 

Overall, this algorithm has proven its value within the workflow it is designed to optimize.


<!-- qPCR Set Up and Data Collection -->
## qPCR Set Up and Data Collection

The development of this methodology uses proprietary antibody pairs and qPCR reagents on a 96 well Quantstudio 5 system. The immuno-PCR plate should be set up to maximize the amount of dCq data points generated, with all signal (Ct when target antigen present) and background (Ct when target antigen NOT present) data points run in at least duplicate. 

Further, one section of the plate should vary detection antibody concentration in a dilution series while holding capture antibody constant, and another section should vary capture antibody holding detection antibody constant. The remainder of the plate must vary both detection and capture in various ways together (holding neither constant). Setting up dilution series is ideal, but the model will work with random concentrations of each antibody as well. It is important to collect as much data across the whole Gompertz surface space as possible. 
Once the experiment is run, the data should be exported to a .csv file. To generate this file, open the “results” tab of the generated Quantstudio excel export. Three columns should be taken from the “results” tab of the generated Quantstudio excel export and saved as a .csv or other R readable file format. The three columns are :

  * The column titled “Sample” (which should specify each well as a signal or background well) 
  * The column titled “Target” (which should specify the two antibody concentrations of the well with a delimiter between the two) Ex: defining the target of a well as “4000x ; 4ug” would indicate a 4000 fold dilution from detection antibody stock and a 4ug/mL capture antibody concentration.
  * The column titled “CT” (which specifies the CT value of the well)

<!-- Code Organization -->
## Code Organization

The first chunk contains four required functions that the algorithm uses to build the model and update predictions. 

  * data_prep: cleans the raw .csv qPCR data to prepare for model fitting. NOTE: this function assumes the dection antibody is specified as a DILUTION factor from a stock of 250ug/ml and assumes the caputure antibody is given as a CONCENTRATION in ug/mL
  * surface_fitter: takes in cleaned data and (currently) fits the Gompertz equation to the data
  * jack_fitter: takes in the parameters from the full data model and undergoes jackknife resampling to estimate bias and create confidence intervals (not yet complete) for the model. NOTE: when 1 or more data points from qPCR data is of low quality or with high standard deviation in dCq the jackknife resampling bias correction is often too extreme and produces poor results. In these cases, removing outliers (not recommended due to data loss) or using the full data model should be done
  * predict_optim: takes in the predictions from the model and scales back the amount of antibodies used to reduce dCq signal down to 95% (default). When background effect is minimal, the model predicts using a maximum concetration of antibodies, despite the fact increasing antibody concentrations has minimal effect. This function ensures real world resource constraints are met 
  * pareto_optim(): takes in the fitted gompertz equation dataframe (table_gompredictior) and picks a random point to run a homebrew pareto optimization front algorithm. The algorithm searches upwards on a normalized dCq axis and leftwards on a normalized total antibody axis to find the pareto front (points where onecannot further optimize one axis without sacrificing the other). This function now runs predict_optim() within it as well.

pareto_optim() then displays a dataframe with 10 rows:
1) Max signal
2) Min overall antibody usage
3) closest to mean of Cab on the optimized front
4) closest to mean of Dab on the optimized front
5) Lowest Cab on pareto front (that isn't lowest overall antibody)
6) Lowest Dab on pareto front (that isn't lowest overall anitbody)
7)  7-9) 3 random other suggestions on the optimized front
10) old optimizer value (highly likely to be near but not on the pareto front)

And returns two objects, one containing plots of the pareto front, one containting the full pareto front data frame (250 points by default) for further analysis if required. 


<!-- SAMPLE RESULTS -->
## Sample Results

Model prediction integrity plots:

![qual_cont](https://user-images.githubusercontent.com/67161057/178116890-9e4a317a-fd1e-4963-877b-3cbf54b89ac2.png)


Visualization of predicted surface (black points) and expiremental data (blue points):

![gomp_surf](https://user-images.githubusercontent.com/67161057/178116916-66842f43-d55a-4c48-bad1-c235bbedaf2f.png)


Pareto boundary visualization:

![pareto_boundary](https://user-images.githubusercontent.com/67161057/178116937-b6332031-40f7-4d9c-a135-607d9fa1431f.png)


Sample suggested results on pareto front:


<img width="253" alt="Screen Shot 2022-07-09 at 1 38 03 PM" src="https://user-images.githubusercontent.com/67161057/178116950-ce0bd177-f6c9-4702-bd47-5ade895438f1.png">


<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
Thanks to RayBiotech Life for providing space and resources to develop and implement this algorithm 

