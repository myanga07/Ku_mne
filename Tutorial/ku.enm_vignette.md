ku.enm vignette
================
Marlon E. Cobos and A. Townsend Peterson
March 15, 2018

-   [Introduction](#introduction)
-   [Getting started](#getting-started)
    -   [Directory structure and necessary data](#directory-structure-and-necessary-data)
    -   [Istalling the package](#istalling-the-package)
    -   [Downloading the example data](#downloading-the-example-data)
-   [Doing the analyses](#doing-the-analyses)
    -   [Making analyses more reproducible](#making-analyses-more-reproducible)
    -   [Calibration models](#calibration-models)
    -   [Evaluation and selection of best models](#evaluation-and-selection-of-best-models)
    -   [Final models creation](#final-models-creation)

<br>

Introduction
------------

**ku.enm** is an R package

This is a brief tutorial for using the functions of the **ku.enm** R package. An example of a disease vector is used in this tutorial to make it more reliable.

Check the functions help to change their arguments according to specific needs.

<br>

Getting started
---------------

### Directory structure and necessary data

Since this package was designed to perform multiple analyses avoiding excessive demands from the computer (especially related to ram memory used for R), it needs certain data and organization in the working directory.Following this structure (Fig.1) will allow working with one or more species in a project and avoid potential problems during the analyses.

Before starting the analyses make sure your working directory has the following components:

-   A folder containing the distinct sets of environmental variables (i.e., M\_variables in Figure 1) that are going to be used (more than one recommended but not mandatory). These variables must represent the area in which your models are going to be calibrated.
-   The maxent.jar application, available at <a href="https://biodiversityinformatics.amnh.org/open_source/maxent/" target="_blank">Maxent repository</a>.
-   A csv file containing all the species occurrences (preferably after cleaning and thinning original data to avoid problems like wrong records or spatial auto-correlation).
-   A csv file containing species occurrences for training models. For obtaining this and the next files the complete set of records must be divided. Occurrences partition can be done in multiple ways, but independence of training and testing data is desired.
-   A csv file containing species occurrences for testing models.

<br>

<img src="Structure.png" alt="Figure 1. Directory structure and necessary data for getting started with the ku.enm R package." style="width:50.0%" />

<br>

### Istalling the package

The **ku.enm** R package is in the CRAN repository and can be installed and/or loaded using the following code.

``` r
if(!require(ku.enm)){
    install.packages("ku.enm")
    library(ku.enm)
}
```

This package can also be installed from a zip file, however all its dependencies should be installed separately.

<br>

### Downloading the example data

Data used as example for testing this package correspond to the turkey tick *Amblyomma americanum*, a vector of various diseases including human monocytotropic ehrlichiosis, canine and human granulocytic ehrlichiosis, tularemia, and southern tick-associated rash illness. This species is distributed in North America and a complete analysis of the risk of its invasion to other areas is being studiend in Raghavan et al. (in prep.).

These data is allready structured as needed for doing the analyses with this package, and can be download (from *url*) and extracted using the code below.

``` r
download.file(url = "url", destfile = "C:/Users/YOUR_USER/Documents/ku_enm_example.zip", 
              quiet = FALSE) #donwload the zipped example folder

unzip(zipfile = "C:/Users/YOUR_USER/Documents/ku_enm_example.zip",
      exdir = "C:/Users/YOUR_USER/Documents") #unzip the example folder

unlink("C:/Users/YOUR_USER/Documents/ku_enm_example.zip") #erase zip file

setwd("C:/Users/YOUR_USER/Documents/ku_enm_example/Species_1") #set the working directory

dir() #check what is in your working directory
```

Your working directory will be structured similar to what is presented in Figure 1.

<br>

Doing the analyses
------------------

### Making analyses more reproducible

Once the working diretory and data are ready, the function *ku.enm.start* (.Rmd) will allow generating an *R Markdown* file as a guide to perform all the analyses that this package includes. By recording all the code chunks that are going to be used during the modeling process this file also helps to make it more reproducible. This file will be written in the working directory.

``` r
help(ku.enm.start)
```

``` r
#Preparing variables to be used in arguments
file_name = "ku_enm_complete_process"
```

``` r
ku.enm.start(file.name = file_name)
```

<br>

### Calibration models

Notice that from this point, the following procedures will be performed in the *R Markdown* file previously created, but only if the *ku.enm.start* function was used.

The function *ku.enm.cal* creates and executes a batch file for generating maxent calibration models that will be written in a sub-directories, named as the parameterizations selected, inside the output directory (Fig. 2 green area). Calibration models will be created with multiple combinations of regularization multipliers, feature classes, and sets of environmental predictors. For each combination this function creates one maxent model with the complete set of occurrences and another with training occurrences only. In some computers the user will be asked if ruining the batch file before starting the modeling process in Maxent.

Maxent will run in command-line interface (**do not close this application**) and its graphic interface will not show up to avoid interfering activities other than the modeling process.

``` r
help(ku.enm.cal)
```

``` r
#Variables with information to be used as arguments
occ_all <- "Sp_all.csv"
occ_tra <- "Sp_cal.csv"
M_var_dir <- "M_vars"
batch_cal <- "ku_enm_calibration_models"
cal_dir <- "Calibration_Models"
reg_mul <- c(seq(0.1,1,0.1),seq(2,6,1),8,10,15,20)
```

``` r
ku.enm.cal(occ.all = occ_all, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
           out.dir = cal_dir, reg.mult = reg_mul, f.clas = "all")
```

<br>

### Evaluation and selection of best models

The function *ku.enm.eval* evaluates models performance based on significance (partial ROC), omission rates (E = 5%), and complexity (AICc), and selects best models based on distinct criteria (see selection.criteria in function help). Partial ROC and omission rates are evaluated on models created with training occurrences, and AICc values are calculated on models created with the full set of occurrences. The outputs will be stored in a folder which will contain: a *CSV* file with the statistics of models meeting various criteria, another with only the selected models based on the chosen criteria, a third one with the performance metrics for all calibration models, a plot *PNG* of the models performance based on the analysed metrics, and an *HTML* file reporting all the results of the model evaluation and selection process designed to guide further interpretations (Fig. 2 purple area).

``` r
help(ku.enm.eval)
```

``` r
#Variables with information to be used as arguments
occ_test <- "Sp_eval.csv"
out_eval <- "ku_enm_evaluation_results"
```

``` r
ku.enm.eval(path = cal_dir, occ.all = occ_all, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
            out.eval = out_eval, omi.val = 5, rand.perc = 50, no.inter = 1000, kept = TRUE, selection = "OR_AICc")
```

<br>

### Final models creation

After selecting parameterizations producing the best models, the next step is creating the final models and, if needed, project them to other areas or scenarios. The *ku.enm.mod* function takes the *CSV* file with the best models result of model selection and writes and executes and batch file for creating final models with the selected parameterizations. Models and projections will be stored in subdirectories inside an output folder and these subdirectories will be named similar to how they are named when calibration models are created. By defining the folder in which scenarios for projections are (i.e. sub-folders with variables of other areas, future or past scenarios, etc.), this function will automatically perform them.

Maxent will run in command-line interface as when creating the calibration models (**again, do not close this application**). However, take into account that the process of creating final models may take considerably more time, especially when projecting, because it will be executing other processes not performed during calibration (e.g., Jackknife analyses).

``` r
help(ku.enm.mod)
```

``` r
#Variables with information to be used as arguments
mod_dir <- "Final_Models"
G_var_dir <- "G_variabless"
##Most of the variables used here as arguments were already created for the previous function
```

``` r
ku.enm.mod(occ.all = occ_all, M.var.dir = M_var_dir, out.eval = out_eval, rep.n = 10, rep.type = "Bootstrap", 
           out.dir = mod_dir, out.format = "logistic", project = TRUE, G.var.dir, 
           ext.type = "all", write.mess = FALSE, write.clamp = FALSE)
```

<br>

At the end of your process you working directory will have the structure and data presented below.

<img src="Structure1.png" alt="Figure 2. Directory structure and data after performing the analyses with the ku.enm R package. Background colors represent necessary data before starting the analyses (blue) and data generated after using the start function (yellow), creating calibration models (green), evaluating calibration models (purple), preparing projection layers (orange), and generating final models and its projections (light grey)." style="width:65.0%" />
