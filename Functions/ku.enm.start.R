ku.enm.start <- function(file.name){
  sink(paste(file.name, ".Rmd", sep = ""))
  cat(
"---
title: \"ku.enm: modeling process\"
---

This R Markdown file is created in the working directory and is intended to make more reproducible the processes of model calibration and final model creation.

Information on using this R Markdown file:

- Try executing code chunks by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.
- Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

A brief tutorial for using functions of the ku.enm R package can be found in https://github.com/manubio13/ku.enm/blob/master/Tutorial/ku.enm_vignette.md. Additionally, function help can be checked to change arguments according to specific needs.

### Calibration models

Calibration models are a large set of candidate models created to respond to the need to test multiple parameter combinations, for example, distinct regularization multiplier values, various feature classes, and different sets of environmental variables. 

Arguments explanation *ku.enm.cal*:

- occ.all (character) is the name of the csv file with all the occurrences, columns must be: species, longitud, latitud
- occ.cal (character) is the name of the csv file with the calibration occurrences, columns equal to occ.all
- M.var.dir (character) is the name of the folder containing other folders with different sets of environmental variables
- batch (character) name of the batch file with the code to create all calibration maxent models
- out.dir (character) name of the folder which will contain all calibration models subfolders
- reg.mult (numeric or numeric vector) regularization multiplier(s) to be evaluated
- Feature clases can be selected from  four different combination sets or manually:
Combination sets are: \"all\", \"basic\", \"no.t.h\", \"no.h\", and \"no.t\". default = \"all\". basic = \"l\", \"lq\", \"lqp\", \"lqpt\", \"lqpth\". 
Manually you can select from the following list:
\"l\", \"q\", \"p\", \"t\", \"h\", \"lq\", \"lp\", \"lt\", \"lh\", \"qp\", \"qt\", \"qh\",
\"pt\", \"ph\", \"th\", \"lqp\", \"lqt\", \"lqh\", \"lpt\", \"lph\", \"qpt\", \"qph\",
\"qth\", \"pth\", \"lqpt\", \"lqph\", \"lqth\", \"lpth\", \"lqpth\"

\```{r}
#Variables with information to be used as arguments
occ_all <- \"Sp_all.csv\"
occ_tra <- \"Sp_cal.csv\"
M_var_dir <- \"M_variables\"
batch_cal <- \"ku_enm_calibration_models\"
cal_dir <- \"Calibration_Models\"
reg_mul <- c(seq(0.1,1,0.1),seq(2,6,1),8,10)
\```

The following is the code for using the function.

\```{r}
ku.enm.cal(occ.all = occ_all, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
out.dir = cal_dir, reg.mult = reg_mul, f.clas = \"all\")
\```


### Evaluation and selection of best models

Evaluation is an important step in model calibration. This step centers on selecting candidate models and their associated parameters to identify the very best models for the purposes of the study. The ku.enm.eval function evaluates candidate models based on three distinct criteria: statistical significance (based on partial ROC), prediction ability (we use omission rates, but other metrics, such as overall correct classification rate, can also be used), and complexity (here evaluated using AICc). 

Arguments explanation for *ku.enm.eval*:

- path is the directory in wich the folders containig calibration models are being or were created 
- occ.all (character) is the name of the csv file with the calibration occurrences, columns must be: species, longitud, latitud
- occ.tra (character) is the name of the csv file with the calibration occurrences, columns must be: species, longitud, latitud
- occ.test (character) is the name of the csv file with the evaluation occurrences, columns must be: species, longitud, latitud
- batch (character) the name of the .bat file created with the KU.ENM.cal function
- out.eval (ccaracter) name of the folder in wich the results of the evaluation will be written
- omi.val (numeric) is the % of omission error allowed (5%)
- rand.perc (numeric) is the percentage of data to be used for the bootstraping process, default (50%) 
- no.inter (numeric) is the number of times that the bootstrap is going to be recalculated, default (100)
- kept (logical) if TRUE all calibration models will be kept after evaluation, default TRUE
- selection (character) model selection criterion, can be \"OR_AICc\", \"AICc\", or \"OR\"; OR = omission rates

\```{r}
#Variables with information to be used as arguments
occ_test <- \"Sp_eval.csv\"
out_eval <- \"ku_enm_evaluation_results\"
##Most of the variables used here as arguments were already created for the previous function
\```

The following code chunk allows evaluating candidate models that were created previously, selecting those with best performance based on the three criteria.

\```{r}
ku.enm.eval(path = cal_dir, occ.all = occ_all, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
out.eval = out_eval, omi.val = 5, rand.perc = 50, no.inter = 1000, kept = TRUE, selection = \"OR_AICc\")
\```

### Final models creation

After selecting parametrizations that produce the best models, the next step is to create the final models, and if needed transfer them to other environmental data sets (e.g., to other periods or geographic regions).

Arguments explanation for *ku.enm.mod*:

- occ.all (character) is the  csv file with all the occurrences, columns must be: species, longitud, latitud
- M.var.dir (character) name of the forlder containing folders in wich calibration (M) variables are
- out.eval (character) name of the folder were evaluation results were written
- rep.n (numeric) number of model replicates
- rep.type (character) is the replicate type, can be: \"Crossvalidate\", \"Bootstrap\", \"Subsample\"
- out.dir (character) name of the output directory to be created and in which all models subdirectories will be created
- out.format (character) is the models output format, can be: \"raw\", \"logistic\", \"cloglog\", \"cumulative\"
- project (logical) if TRUE your models will be projected to the scenarios in G.var.dir, default = FALSE
- G.var.dir (character) if project is TRUE, name of the forlder containing folders in wich variables of your projection scenarios are
- ext.type (character) if project is TRUE, is the extrapolation type of projections, can be: \"all\", \"ext_clam\", \"ext\", and \"no_ext\", default = \"all\"
- write.mess (logical) if TRUE, grids of MESS analysis results will be written, default = FALSE
- write.clamp (logical) if TRUE, a grid of the spatial distribution of clamping will be written, default = FALSE

\```{r, eval=FALSE, include=TRUE}
#Variables with information to be used as arguments
mod_dir <- \"Final_Models\"
G_var_dir <- \"G_variabless\"
##Most of the variables used here as arguments were already created for the previous function
\```

The *ku.enm.mod* function has the following syntax:

\```{r, eval=FALSE, include=TRUE}
ku.enm.mod(occ.all = occ_all, M.var.dir = M_var_dir, out.eval = out_eval, rep.n = 10, rep.type = \"Bootstrap\", 
out.dir = mod_dir, out.format = \"logistic\", project = TRUE, G.var.dir, 
ext.type = \"all\", write.mess = FALSE, write.clamp = FALSE)
\```"
  )
  sink()
  file.edit(paste(file.name, ".Rmd", sep = ""))
}