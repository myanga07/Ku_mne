ku.enm.eval <- function(path, occ.all, occ.tra, occ.test, batch, out.eval, omi.val = 5, 
                        rand.perc = 50, no.inter = 1000, kept = TRUE, selection = "OR_AICc"){
  #####
  #Packages
  if(!require(yaml)){
    install.packages("yaml")
    library(yaml)
  }
  if(!require(ENMeval)){
    install.packages("ENMeval")
    library(ENMeval)
  }
  if(!require(sqldf)){
    install.packages("sqldf")
    library(sqldf)
  }
  if(!require(rmarkdown)){
    install.packages("rmarkdown")
    library(rmarkdown)
  }
  if(!require(knitr)){
    install.packages("knitr")
    library(knitr)
  }
  
  rasterOptions(maxmemory = 1e+15) #giving more power to the raster package
  options(digits = 16) #number of digits that r can recognize
  
  #####
  #Functions
  ##Partial ROC related functions
  PartialROC <- function(PresenceFile=NA, PredictionFile=NA, OmissionVal=NA, RandomPercent=NA, NoOfIteration=NA)
  {
    ##Raw to logistic values
    InRastlog = PredictionFile
    
    ## Currently fixing the number of classes to 100. But later flexibility should be given in the parameter.
    InRast = round((InRastlog/cellStats(InRastlog,max))*1000)
    
    ## This function should be called only once outside the loop. This function generates values for x-axis.
    ## As x-axis is not going to 
    ClassPixels = AreaPredictedPresence(InRast)
    
    Occur = read.csv(PresenceFile)
    Occur = Occur[,-1]
    ExtRast = extract(InRast, Occur)
    
    ## Remove all the occurrences in the class NA. As these points are not used in the calibration.
    OccurTbl = cbind(Occur, ExtRast)
    OccurTbl = OccurTbl[which(is.na(OccurTbl[,3]) == FALSE),]
    
    PointID = seq(1:nrow(OccurTbl))
    OccurTbl = cbind(PointID, OccurTbl)
    names(OccurTbl)= c("PointID", "Longitude", "Latitude", "ClassID")
    # # ## Generate the % points within each class in this table. Write SQL, using sqldf package
    # # OccurINClass = sqldf("Select count(*), ClassID from OccurTbl group by ClassID order by ClassID desc")
    # # OccurINClass = cbind(OccurINClass, cumsum(OccurINClass[,1]), cumsum(OccurINClass[,1]) / nrow(OccurTbl))
    # # names(OccurINClass) = c("OccuCount", "ClassID", "OccuSumBelow", "Percent")
    
    ## Use option cl.cores to choose an appropriate cluster size.
    lapply(X = 1:NoOfIteration,FUN =  function(x){
      ll = sample(nrow(OccurTbl), round(RandomPercent/100 * nrow(OccurTbl)), replace=TRUE)
      OccurTbl1 = OccurTbl[ll,]
      ## Generate the % points within each class in this table. Write SQL, using sqldf package
      OccurINClass = sqldf("Select count(*), ClassID from OccurTbl1 group by ClassID order by ClassID desc")
      OccurINClass = cbind(OccurINClass, cumsum(OccurINClass[,1]), cumsum(OccurINClass[,1]) / nrow(OccurTbl1))
      names(OccurINClass) = c("OccuCount", "ClassID", "OccuSumBelow", "Percent")
      
      #### Raster file will contain all the classes in ClassID column, while occurrences table may not have all the classes.
      #### Somehow we have to make the ClassID same as raster file ClassID. This could be done with SQL command update.
      #### but update is not working, not sure about the reason. So I am running the loop which is very very slow.
      XYTable = GenerateXYTable(ClassPixels,OccurINClass)
      #plot(XYTable[,2], XYTable[,3])
      AreaRow = CalculateAUC(XYTable, OmissionVal, x)
      names(AreaRow) <- c("IterationNo", paste("AUC_at_Value_", OmissionVal, sep = ""), "AUC_at_0.5", "AUC_ratio")
      #AreaRow[1] <- as.integer(AreaRow[1])
      return(AreaRow)
    })
  }
  
  AreaPredictedPresence <- function(InRast)
  {
    ### Now calculate proportionate area predicted under each suitability
    ClassPixels = freq(InRast)
    ### Remove the NA pixels from the table.
    if (is.na(ClassPixels[dim(ClassPixels)[1],1])== TRUE)
    {
      ClassPixels = ClassPixels[-dim(ClassPixels)[1],]
    }
    ClassPixels = ClassPixels[order(nrow(ClassPixels):1),]
    TotPixPerClass = cumsum(ClassPixels[,2])
    PercentPixels = TotPixPerClass / sum(ClassPixels[,2])
    
    ClassPixels = cbind(ClassPixels, TotPixPerClass, PercentPixels)
    ClassPixels = ClassPixels[order(nrow(ClassPixels):1),]
    return(ClassPixels)
  }
  
  ## This function generates the XY coordinate table. Using this table areas is calculated.
  GenerateXYTable<-function(ClassPixels, OccurINClass)
  {
    XYTable = ClassPixels[,c(1,4)]
    XYTable = cbind(XYTable,rep(-1,nrow(XYTable)))
    # names(XYTable) = c("ClassID", "XCoor", "YCoor")
    ## Set the previous value for 1-omission, i.e Y-axis as the value of last
    ## class id in Occurrence table. LAst class id will always smallest
    ## area predicted presence.
    PrevYVal = OccurINClass[1,4]
    for (i in nrow(ClassPixels):1)
    {
      CurClassID = XYTable[i,1]
      YVal = OccurINClass[which(OccurINClass[,2]==CurClassID),4]
      ## print(paste("Length of YVal :",length(YVal), "Current Loop count :", i, "Current value of YVal : ", YVal, sep = " " ))
      
      if (length(YVal) == 0 )
      {
        XYTable[i,3] = PrevYVal
      }
      else
      {
        XYTable[i,3] = YVal
        PrevYVal = YVal
      }
    }
    ## Add A dummy class id in the XYTable with coordinate as 0,0
    XYTable = rbind(XYTable, c(XYTable[nrow(XYTable),1] + 1, 0, 0))
    XYTable = as.data.frame(XYTable)
    names(XYTable) = c("ClassID", "XCoor", "YCoor")
    ### Now calculate the area using trapezoid method.
    return(XYTable)
  }
  
  CalculateAUC <- function(XYTable, OmissionVal, IterationNo)
  {
    ## if OmissionVal is 0, then calculate the complete area under the curve. Otherwise calculate only partial area
    if (OmissionVal > 0)
    {
      PartialXYTable = XYTable[which(XYTable[,3] >= OmissionVal),]
      ### Here calculate the X, Y coordinate for the parallel line to x-axis depending upon the OmissionVal
      ### Get the classid which is bigger than the last row of the XYTable and get the XCor and Ycor for that class
      ### So that slope of the line is calculated and then intersection point between line parallel to x-axis and passing through
      ### ommissionval on Y-axis is calculated.
      PrevXCor = XYTable[which(XYTable[,1]==PartialXYTable[nrow(PartialXYTable),1])+1,2]
      PrevYCor = XYTable[which(XYTable[,1]==PartialXYTable[nrow(PartialXYTable),1])+1,3]
      XCor1 = PartialXYTable[nrow(PartialXYTable),2]
      YCor1 = PartialXYTable[nrow(PartialXYTable),3]
      ## Calculate the point of intersection of line parallel to x-asiz and this line. Use the equation of line
      ## in point-slope form y1 = m(x1-x2)+y2
      Slope = (YCor1 - PrevYCor) / (XCor1 - PrevXCor)
      YCor0 = OmissionVal
      XCor0 = (YCor0 - PrevYCor + (Slope * PrevXCor)) / Slope
      ### Add this coordinate in the PartialXYTable with classid greater than highest class id in the table.
      ### Actually class-id is not that important now, only the place where we add this xcor0 and ycor0 is important.
      ### add this as last row in the table
      PartialXYTable = rbind(PartialXYTable, c(PartialXYTable[nrow(PartialXYTable),1]+1, XCor0, YCor0))
    }
    else
    {
      PartialXYTable = XYTable
    } ### if OmissionVal > 0
    
    ## Now calculate the area under the curve on this table.
    XCor1 = PartialXYTable[nrow(PartialXYTable),2]
    YCor1 = PartialXYTable[nrow(PartialXYTable),3]
    AUCValue = 0
    AUCValueAtRandom = 0
    for (i in (nrow(PartialXYTable)-1):1)
    {
      XCor2 = PartialXYTable[i,2]
      YCor2 = PartialXYTable[i,3]
      
      # This is calculating the AUCArea for 2 point trapezoid.
      TrapArea = (YCor1 * (abs(XCor2 - XCor1))) + (abs(YCor2 - YCor1) * abs(XCor2 - XCor1)) / 2
      AUCValue = AUCValue + TrapArea
      # now caluclate the area below 0.5 line.
      # Find the slope of line which goes to the point
      # Equation of line parallel to Y-axis is X=k and equation of line at 0.5 is y = x
      TrapAreaAtRandom = (XCor1 * (abs(XCor2 - XCor1))) + (abs(XCor2 - XCor1) * abs(XCor2 - XCor1)) / 2
      AUCValueAtRandom = AUCValueAtRandom + TrapAreaAtRandom
      XCor1 = XCor2
      YCor1 = YCor2
    }
    NewRow = c(IterationNo, AUCValue, AUCValueAtRandom, AUCValue/AUCValueAtRandom)
    return(NewRow)
  }
  
  rasToBinary <- function(rasterPre,threshold){
    rDF <- data.frame(rasterToPoints(rasterPre))
    names(rDF) <- c("x","y","binary")
    reclass <- sapply(rDF$binary, function(x){
      if(x >= threshold) return(1)
      else return(0)
    })
    rDF$binary <- reclass
    coordinates(rDF) <- ~ x + y
    # coerce to SpatialPixelsDataFrame
    gridded(rDF) <- TRUE
    # coerce to raster
    binRaster <- raster(rDF)
    return(binRaster)
  }
  
  ##AICc related functions
  ###Function to get number of parameters (Modified from ENMeval)
  n.par <- function(x){
    lambdas <- x[1:(length(x) - 4)]
    countNonZeroParams <- function(x) {
      if (strsplit(x, split = ", ")[[1]][2] != "0.0") 
        1
    }
    no.params <- sum(unlist(sapply(lambdas, countNonZeroParams)))
    return(no.params)
  }
  
  #####
  #Data
  ##Source of initial information for model evaluation order
  bat <- readLines(paste(batch, ".bat", sep = "")) #reading the batch file written to create the calibration models
  
  ###Recognizing the folders names and separating them for different procedures
  fol <- gregexpr("outputd.\\S*", bat)
  fold <- regmatches(bat, fol)
  folde <- unlist(fold)
  extract <- paste("outputdirectory=", path, "\\", sep = "")
  folder <- gsub(extract, "", folde, fixed = T) #names of all the calibration models folders
  
  folder_a <- gregexpr("\\S*all", folder)
  folder_al <- regmatches(folder, folder_a)
  folder_all <- unlist(folder_al) #folders with the models for calculating AICcs
  
  folder_c <- gregexpr("\\S*cal", folder)
  folder_ca <- regmatches(folder, folder_c)
  folder_cal <- unlist(folder_ca) #folder with the models for calculating pROCs and omission rates
  
  ##Models
  ###For AICc
  dir_names <- as.vector(paste(getwd(), "/", path, "/", folder_all, sep = "")) #vector of the subdirectories with the models
  
  ###For pROC and omission rates
  dir_names1 <- as.vector(paste(getwd(), "/", path, "/", folder_cal, sep = "")) #vector of the subdirectories with the models
  
  ###Names of the models to be evaluated
  mod_nam <- as.vector(gsub("_all", "", folder_all, fixed = TRUE)) #names of the models (taken from the folders names)
  
  ##Complete set and calibration and evaluation occurrences
  oc <- read.csv(occ.all) #read all occurrences
  oc <- oc[,-1] #erase species name column
  
  occ <- read.csv(occ.tra) #read calibration occurrences
  occ <- occ[,-1] #erase species name column
  
  occ1 <- read.csv(occ.test) #read test occurrences
  occ1 <- occ1[,-1] #erase species name column
  
  ##Omission value
  OmissionVal <- omi.val / 100
  
  ###Place of the occ.tra with the value to be considered the omi.val
  val <- ceiling(length(occ[,1]) * omi.val / 100) + 1 
  
  #####
  #pROCs, omission rates, and AICcs calculation
  cat("\nPartial ROCs, omission rates, and AICcs calculation, please wait...\n")
  
  AICcs <- list() #empty list of AICc results
  proc_res_med <- list() #empty list of mean AUC values
  proc_res_lt1 <- vector() #empty vector of significance values of the AUCs
  om_rates <- vector() #empty vector of ommision rates 
  
  pb <- winProgressBar(title = "Progress bar", min = 0, max = length(dir_names), width = 300) #progress bar
  
  for(i in 1:length(dir_names)){ #calculating AUC rations in a loop
    Sys.sleep(0.1)
    setWinProgressBar(pb, i, title = paste( round(i / length(dir_names) * 100, 2), "% of the evaluation process has finished"))
    
    #AICc calculation
    suppressWarnings(while (!file.exists(as.vector(list.files(dir_names[i], pattern=".lambdas", full.names=TRUE)))) {
      Sys.sleep(1)
    })
    lbds <- as.vector(list.files(dir_names[i], pattern=".lambdas", full.names=TRUE)) #lambdas file
    lambdas <- readLines(lbds)
    par_num <- n.par(lambdas) #getting the number of parameters for each model
    suppressWarnings(while (!file.exists(dir_names[i], pattern="asc", full.names=TRUE)) {
      Sys.sleep(1)
    })
    mods <- list.files(dir_names[i], pattern="asc", full.names=TRUE) #name of ascii model
    mod <- raster(mods) #reading each ascii model created with the complete set of occurrences
    AICcs[[i]] <- suppressWarnings(calc.aicc(nparam = par_num, occ = oc, predictive.maps = mod)) #calculating AICc for each model
    
    #pROCs calculation
    suppressWarnings(while (!file.exists(dir_names1[i], pattern="asc", full.names=TRUE)) {
      Sys.sleep(1)
    })
    mods1 <- list.files(dir_names1[i], pattern="asc", full.names=TRUE) #list of the ascii models
    mod1 <- raster(mods1) #reading each ascii model created with the calibration occurrences
    proc_res <- PartialROC(PresenceFile = occ.test, PredictionFile = mod1, 
                           OmissionVal = OmissionVal, RandomPercent = rand.perc, 
                           NoOfIteration = no.inter) #Partial ROC analyses for each model
    proc_res <- as.data.frame(do.call(rbind, proc_res)) #converting each list of AUC ratios interations in a table for each model
    proc_res_med[[i]] <- apply(proc_res, 2, mean) #mean of AUC ratios interations for each model
    proc_res_lt1[i] <- sum(proc_res[,4] <= 1)/length(proc_res[,4]) #proportion of AUC ratios <= 1 for each model  
    
    #Omission rates calculation
    suit_val_cal <- na.omit(extract(mod1, occ))
    suit_val_eval <- na.omit(extract(mod1, occ1))
    omi.val.suit <- sort(suit_val_cal)[val]
    om_rates[i] <- as.data.frame(length(suit_val_eval[suit_val_eval < omi.val.suit])/length(suit_val_eval))
    
    #Erasing calibration models after evaluating them if kept = FALSE
    if(kept == FALSE){
      unlink(dir_names[i], recursive = T)
      unlink(dir_names1[i], recursive = T)
    }
  }
  suppressMessages(close(pb))
  n.mod <- i
  
  ##Erasing main folder of calibration models if kept = FALSE
  if(kept == FALSE){
    unlink(path, recursive = T)
    cat("\nAll calibration models were deleted\n")
  }else{
    cat("\nAll calibration models were kept\n")
  }
  
  ##Creating the final tables
  ###From AICc analyses
  AICcs <- do.call(rbind, AICcs) #joining tables
  for (i in 1:length(AICcs[,1])) {
    AICcs[i,2] <- (AICcs[i,1] - min(AICcs[,1], na.rm = TRUE))
    AICcs[i,3] <- (exp(-0.5 * AICcs[i,2])) / (sum(exp(-0.5 * AICcs[,2]), na.rm = TRUE))
  }
  
  ###From pROC analyses
  proc_res_m <- do.call(rbind, proc_res_med)[,-(2:3)] #joining tables of the mean AUC ratios
  proc_res_m[,1] <- mod_nam #declaring model names
  proc_res_m <- data.frame(proc_res_m, proc_res_lt1) #adding a new column with the number of AUC ratios interations < 1
  
  ###Omission rates
  om_rates <- unlist(om_rates)
  
  #####
  #Joining the results
  KU_ENM_eval <- as.data.frame(cbind(proc_res_m, om_rates, AICcs))
  KU_ENM_eval <- data.frame(KU_ENM_eval[,1], as.numeric(levels(KU_ENM_eval[,2]))[KU_ENM_eval[,2]], KU_ENM_eval[,3], 
                            KU_ENM_eval[,4], KU_ENM_eval[,5], KU_ENM_eval[,6], KU_ENM_eval[,7], KU_ENM_eval[,8])
  colnames(KU_ENM_eval) <- c("Model", "Mean_AUC_ratio", "Partial_ROC",#changing column names in the final table
                             paste("Ommission_rate_at_", omi.val, "%", sep = ""), "AICc",
                             "delta_AICc", "W_AICc", "num_parameters")
  
  
  #####
  #Choosing the best models
  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR"){
    if(selection == "OR_AICc"){
      KU_ENM_bes <- KU_ENM_eval[KU_ENM_eval[,3] <= 0.05,]
      KU_ENM_best <- KU_ENM_bes[KU_ENM_bes[,4] <= OmissionVal,]
      if(length(KU_ENM_best[,4]) != 0){
        for (i in 1:length(KU_ENM_best[,1])) {
          KU_ENM_best[i,6] <- (KU_ENM_best[i,5] - min(KU_ENM_best[,5], na.rm = TRUE))
          KU_ENM_best[i,7] <- (exp(-0.5 * KU_ENM_best[i,6])) / (sum(exp(-0.5 * KU_ENM_best[,6]), na.rm = TRUE))
        }
        KU_ENM_best <- na.omit(KU_ENM_best[KU_ENM_best[,6] <= 2,])
        if(length(KU_ENM_best[,6]) > 10){
          KU_ENM_best <- KU_ENM_best[order(KU_ENM_best[,6]),][1:10,] 
        }else{
          KU_ENM_best <- KU_ENM_best[order(KU_ENM_best[,6]),]
        }
      }else{
        mesKU <- "\nNone of your models meets the omission rates criterion,\n
        models with the smallest omission rates will be presented\n"
        KU_ENM_best <- KU_ENM_bes[order(KU_ENM_bes[,4]),][1:100,]
        for (i in 1:length(KU_ENM_best[,1])) {
          KU_ENM_best[i,6] <- (KU_ENM_best[i,5] - min(KU_ENM_best[,5], na.rm = TRUE))
          KU_ENM_best[i,7] <- (exp(-0.5 * KU_ENM_best[i,6])) / (sum(exp(-0.5 * KU_ENM_best[i,6]), na.rm = TRUE))
        }
        KU_ENM_best <- na.omit(KU_ENM_best[KU_ENM_best[,6] <= 2,])
        if(length(KU_ENM_best[,6]) > 10){
          KU_ENM_best <- KU_ENM_best[order(KU_ENM_best[,6]),][1:10,] 
        }else{
          KU_ENM_best <- KU_ENM_best[order(KU_ENM_best[,6]),]
        }
      }
    }
    
    if(selection == "AICc"){
      KU_ENM_best1 <- KU_ENM_eval[KU_ENM_eval[,3] <= 0.05,]
      KU_ENM_best1 <- na.omit(KU_ENM_best1[KU_ENM_best1[,6] <= 2,])
      if(length(KU_ENM_best1[,6]) > 10){
        KU_ENM_best1 <- KU_ENM_best1[order(KU_ENM_best1[,6]),][1:10,] 
      }else{
        KU_ENM_best1 <- KU_ENM_best1[order(KU_ENM_best1[,6]),]
      }
    }
    
    if(selection == "OR"){
      KU_ENM_bes <- KU_ENM_eval[KU_ENM_eval[,3] <= 0.05,]
      KU_ENM_best2 <- KU_ENM_bes[KU_ENM_bes[,4] <= OmissionVal,]
      if(length(KU_ENM_best2[,4]) != 0){
        if(length(KU_ENM_best2[,4]) > 10){
          KU_ENM_best2 <- KU_ENM_best2[order(KU_ENM_best2[,4]),][1:10,] 
        }else{
          KU_ENM_best2 <- KU_ENM_best2[order(KU_ENM_best2[,4]),]
        } 
      }else{
        mesKU <- "\nNone of your models meets the omission rates criterion,\n
        models with the smallest omission rates will be presented\n"
        KU_ENM_best2 <- KU_ENM_bes[order(KU_ENM_bes[,4]),][1:10,]
      }
    }
    }else{
      cat("\nNo valid model selection criterion has been defined,\n
          no file containing the best models will be created.\n
          Select your best models from the complete list.\n")
    }
  
  #####
  #Statistics of the process
  ##Counting
  KU_ENM_best_OR_AICc <- KU_ENM_bes[KU_ENM_bes[,4] <= OmissionVal,]
  if(length(KU_ENM_best_OR_AICc[,4]) != 0){
    for (i in 1:length(KU_ENM_best_OR_AICc[,1])) {
      KU_ENM_best_OR_AICc[i,6] <- (KU_ENM_best_OR_AICc[i,5] - min(KU_ENM_best_OR_AICc[,5], na.rm = TRUE))
      KU_ENM_best_OR_AICc[i,7] <- (exp(-0.5 * KU_ENM_best_OR_AICc[i,6])) / (sum(exp(-0.5 * KU_ENM_best_OR_AICc[,6]), na.rm = TRUE))
    }
    KU_ENM_best_OR_AICc <- na.omit(KU_ENM_best_OR_AICc[KU_ENM_best_OR_AICc[,6] <= 2,])
  }
  
  KU_ENM_best_AICc <- na.omit(KU_ENM_bes[KU_ENM_bes[,6] <= 2,])
  
  KU_ENM_best_OR <- KU_ENM_bes[KU_ENM_bes[,4] <= OmissionVal,]
  
  ##Preparing the table
  r_names <- c("Statistically significant models", "Models meeting OR criteria", 
               "Models meeting AICc critera", "Models meeting OR and AICc criteria")
  statis <- c(length(KU_ENM_bes[,3]),
              length(KU_ENM_best_OR[,4]),
              length(KU_ENM_best_AICc[,6]),
              length(KU_ENM_best_OR_AICc[,2]))
  KU_ENM_stats <- cbind(r_names, statis)
  colnames(KU_ENM_stats) <- c("Criteria", "Number of models")
  
  #####
  #Writing the results
  ##csv files
  cat("\nWriting ku.enm.eval results...\n")
  dnam <- "ku_enm_evaluation_results"
  dir.create(dnam)
  
  name <- paste(dnam, "evaluation_results.csv", sep = "/")
  name0 <- paste(dnam, "evaluation_stats.csv", sep = "/")
  name1 <- paste(dnam, "best_models_OR_AICc.csv", sep = "/")
  name2 <- paste(dnam, "best_models_AICc.csv", sep = "/")
  name3 <- paste(dnam, "best_models_OR.csv", sep = "/")
  
  
  write.csv(KU_ENM_eval, file = name, eol = "\n", na = "NA", row.names = FALSE)
  write.csv(KU_ENM_stats, file = name0, eol = "\n", na = "NA", row.names = FALSE)
  
  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR"){
    if(selection == "OR_AICc"){
      write.csv(KU_ENM_best, file = name1, eol = "\n", na = "NA", row.names = FALSE)
    }
    if(selection == "AICc"){
      write.csv(KU_ENM_best1, file = name2, eol = "\n", na = "NA", row.names = FALSE)
    }
    if(selection == "OR"){
      write.csv(KU_ENM_best2, file = name3, eol = "\n", na = "NA", row.names = FALSE)
    }
  }
  
  ##Plot
  png(paste(dnam, "evaluation_figure.png", sep = "/"), width = 80, height = 80, units = "mm", res = 600)
  par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.5)
  plot(`Ommission_rate_at_5%`~log(AICc), data = na.omit(KU_ENM_eval),
       xlab = "Natural logarithm of AICc", ylab = paste("Omission rates at", 
                                                        paste(omi.val, "%", sep = ""), "threshold value", sep = " "),
       las = 1, col = "gray35")
  points(`Ommission_rate_at_5%`~log(AICc), data = na.omit(KU_ENM_eval[!KU_ENM_eval[,1] %in% KU_ENM_bes[,1],]),
         col = "red1", pch = 19, cex = 1.5)
  
  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR"){
    if(selection == "OR_AICc"){
      points(`Ommission_rate_at_5%`~log(AICc), data = na.omit(KU_ENM_best),
             col = "dodgerblue1", pch = 19, cex = 1.8)
      legend("bottomright", legend = c("Selected models OR & AICc", "Non significant models", "All models"),
             pt.cex = c(1.8, 1.5, 1), pch = c(19, 19, 1), col = c("dodgerblue1", "red1", "gray35"), bty = "n",
             inset = c(0.01, 0))
    }
    if(selection == "AICc"){
      points(`Ommission_rate_at_5%`~log(AICc), data = na.omit(KU_ENM_best1),
             col = "darkorchid1", pch = 19, cex = 1.8)
      legend("bottomright", legend = c("Selected models AICc", "Non significant models", "All models"), 
             pt.cex = c(1.8, 1.5, 1), pch = c(19, 19, 1), col = c("darkorchid1", "red1", "gray35"), bty = "n",
             inset = c(0.01, 0))
    }
    if(selection == "OR"){
      points(`Ommission_rate_at_5%`~log(AICc), data = na.omit(KU_ENM_best2),
             col = "orange2", pch = 19, cex = 1.8)
      legend("bottomright", legend = c("Selected models OR", "Non significant models", "All models"),
             pt.cex = c(1.8, 1.5, 1), pch = c(19, 19, 1), col = c("orange2", "red1", "gray35"), bty = "n",
             inset = c(0.01, 0))
    }
  }
  dev.off()
  
  ##html file
  ###Writing the html file
  sink(paste(dnam, "evaluation_results.Rmd", sep = "/"))
  cat("---
title: \"ku.enm: evaluation results\"
output:
  html_document:
      toc: true
      toc_depth: 4
---

\```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
\```

<br>

### Breif description of the models evaluation and selection process

\```{r, echo=FALSE}
st4 <- read.csv(\"evaluation_results.csv\")
sett <- as.character(st4[,1])
setts <- strsplit(sett, split = \"_\")
rm <- vector()
for (i in 1:length(setts)) {
rm[i] <- setts[[i]][2]
}
f.clas <- vector()
for (i in 1:length(setts)) {
f.clas[i] <- setts[[i]][4]
}
var.di <- vector()
for (i in 1:length(setts)) {
var.di[i] <- paste(setts[[i]][5:length(setts[[i]])], collapse = \"_\")
}
rm1 <- paste(unique(rm), collapse = \", \")
f.clas1 <- paste(unique(f.clas), collapse = \", \")
var.di1 <- paste(unique(var.di), collapse = \", \")
par <- rbind(rm1, f.clas1, var.di1)
\```

This is the final report of the ku.enm.eval function implemented in the ku.enm R package. 

A total of \`r length(st4[,1])\` calibration models with parameters resulted from the combination of \`r length(unique(rm))\` regularization multipliers, \`r length(unique(f.clas))\` feature classes, and \`r length(unique(var.di))\` distinct sets of environmental variables, have been evaluated. Models peformance was evaluated based on statistical significance (Partial_ROC), omission rates, and the Akaike information criterion corrected for small sample sizes (AICc).

\```{r par, echo=FALSE}
colnames(par) <- \"Parameters\"
row.names(par) <- c(\"Regularization multipliers\", \"Feature classes\", \"Sets of predictors\")
knitr::kable(par, digits=c(0,0), row.names = TRUE, caption = \"Table 1. Parameters of the evaluated models.\")
\```

<br>
All the results presented below can be found in the evaluation output folder for further analyses.

<br>
<br>

### Models evaluation statistics

In the following table you are going to find information about how many models meet the four criteria of selection that this function uses.

\```{r, echo=FALSE}
st <- read.csv(\"evaluation_stats.csv\")
colnames(st) <- c(\"Criteria\",	\"Number_of_models\")
knitr::kable(st, digits=c(0,0), caption = \"Table 2. General statistics of models that meet distinct criteria.\")
\```

<br>
<br>

### Best models according to pre-defined criteria

The following table contains the best models selected by your pre-defined criteria.

Notice that if the selection criterion was, models with the best Omission rates and among them those with lower AICc values, the delta_AICc values were recalculated among the new candidate models.

\```{r, echo=FALSE}
best <- list.files(pattern = \"best\")
st1 <- read.csv(best)
colnames(st1) <- c(\"Model\",	\"Mean_AUC_ratio\",	\"Partial_ROC\", \"Ommission_rate_5%\", \"AICc\",	\"delta_AICc\",	\"W_AICc\",	\"num_parameters\")
knitr::kable(st1, digits=c(0,3,3,3,3,3,3,0), caption = \"Table 3. Performance statistics of the best models selected based on the pre-defined critera.\")
\```

<br>
<br>

### Models preformance plot

The figure below shows the position of your selected models in the distribution of all the calibration models according to omission rates and AICc values.

![Figure 1. Distribution of all, non-statistically significant, and selected models according to their omission rates and AICc values.](evaluation_figure.png){width=60%}

<br>
<br>

### Performance evaluation statistics for all models

Following you will find the performance statistics for all the calibration models.

\```{r, echo=FALSE}
st4 <- read.csv(\"evaluation_results.csv\")
colnames(st4) <-  c(\"Model\",	\"Mean_AUC_ratio\",	\"Partial_ROC\", \"Ommission_rate_5%\", \"AICc\",	\"delta_AICc\",	\"W_AICc\",	\"num_parameters\")
knitr::kable(st4, digits=c(0,3,3,3,3,3,3,0), caption = \"Table 4. Performance statistics of all the calibration models.\")
\```
      ")
  sink()
  render(paste(dnam, "evaluation_results.Rmd", sep = "/"), "html_document", quiet = TRUE)
  unlink(paste(dnam, "evaluation_results.Rmd", sep = "/"))
  
  #####
  #Finalizing the function
  cat("\nProcess finished\n")
  cat(paste("A folder containing the results of the evaluation of", n.mod,
            "\ncalibration models has been written\n", sep = " "))
  
  cat(paste("\nThe folder", dnam, "contains:\n", sep = " "))
  cat("   -A html file and its dependencies that sum all the results, check\n")
  cat(paste("    ", "evaluation_results.html\n", sep = ""))
  
  cat("   -Two csv files with the results and stats of the evaluation,\n")
  
  if(selection == "OR_AICc" | selection == "AICc" | selection == "OR"){
    if(selection == "OR_AICc"){
      cat("    and an aditional csv file containing the best models selected by OR and AICc.\n")
    }
    if(selection == "AICc"){
      cat("    and an aditional csv file containing the best models selected by AICc.\n")
    }
    if(selection == "OR"){
      cat("    and an aditional csv file containing the best models selected by OR.\n")
    }
  }
  
  if(exists("mesKU") == TRUE){
    warning(mesKU)
  }
  
  cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
}