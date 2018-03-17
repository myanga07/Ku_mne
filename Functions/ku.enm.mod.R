ku.enm.mod <- function(occ.all, M.var.dir, out.eval, rep.n = 10, rep.type = "Bootstrap", 
                       out.dir, out.format = "logistic", project = FALSE, G.var.dir,  
                       ext.type = "all", write.mess = FALSE, write.clamp = FALSE){
  #####
  #Funtion to get free ram 
  get_free_ram <- function(){
    if(Sys.info()[["sysname"]] == "Windows"){
      x <- system2("wmic", args =  "OS get FreePhysicalMemory /Value", stdout = TRUE)
      x <- x[grepl("FreePhysicalMemory", x)]
      x <- gsub("FreePhysicalMemory=", "", x, fixed = TRUE)
      x <- gsub("\r", "", x, fixed = TRUE)
      as.integer(x)
    } else {
      cat("\nOnly supported on Windows OS\n1.6 Gb of the memory will be used for java runnings\n")
      x <- 4000000
    }
  }
  
  #####
  #Data
  ##Data from best models table
  best <- list.files(path = out.eval, pattern = "best")
  sett <- read.csv(best)
  sett1 <- as.character(sett[,1])
  setts <- strsplit(sett1, split = "_")
  
  ###Regularization multipliers
  rm <- vector()
  for (i in 1:length(setts)) {
    rm[i] <- setts[[i]][2]
  }
  
  ###Feature classes
  f.clas <- vector()
  for (i in 1:length(setts)) {
    f.clas[i] <- setts[[i]][4]
  }
  
  ###Calibration (M) variables
  var.di <- vector()
  for (i in 1:length(setts)) {
    var.di[i] <- setts[[i]][5]
  }
  var.dir <- paste(M.var.dir, var.di, sep = "\\")
  
  #output directory
  dir.create(out.dir)
  
  #Defining maximum ram to be used (50% of free memory)
  ram <- paste("-mx", (round((get_free_ram()/1000)*0.5)), "m", sep = "")
  
  #####
  #Maxent settings
  ##Environmental calibration variables sets
  env <- vector()
  for (i in 1:length(var.dir)) {
    env[i] <- paste("environmentallayers=", var.dir[i], sep = "")
  }
  
  ##Species occurrences
  samp <- paste("samplesfile=", occ.all, sep = "")
  
  ##Feature classes combinations
  fea <- c("linear=true quadratic=false product=false threshold=false hinge=false",
           "linear=false quadratic=true product=false threshold=false hinge=false",
           "linear=false quadratic=false product=true threshold=false hinge=false",
           "linear=false quadratic=false product=false threshold=true hinge=false",
           "linear=false quadratic=false product=false threshold=false hinge=true",
           "linear=true quadratic=true product=false threshold=false hinge=false",
           "linear=true quadratic=false product=true threshold=false hinge=false",
           "linear=true quadratic=false product=false threshold=true hinge=false",
           "linear=true quadratic=false product=false threshold=false hinge=true",
           "linear=false quadratic=true product=true threshold=false hinge=false",
           "linear=false quadratic=true product=false threshold=true hinge=false",
           "linear=false quadratic=true product=false threshold=false hinge=true",
           "linear=false quadratic=false product=true threshold=true hinge=false",
           "linear=false quadratic=false product=true threshold=false hinge=true",
           "linear=false quadratic=false product=false threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=false hinge=false",
           "linear=true quadratic=true product=false threshold=true hinge=false",
           "linear=true quadratic=true product=false threshold=false hinge=true",
           "linear=true quadratic=false product=true threshold=true hinge=false",
           "linear=true quadratic=false product=true threshold=false hinge=true",
           "linear=false quadratic=true product=true threshold=true hinge=false",
           "linear=false quadratic=true product=true threshold=false hinge=true",
           "linear=false quadratic=true product=false threshold=true hinge=true",
           "linear=false quadratic=false product=true threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=true hinge=false",
           "linear=true quadratic=true product=true threshold=false hinge=true",
           "linear=true quadratic=true product=false threshold=true hinge=true",
           "linear=true quadratic=false product=true threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=true hinge=true")
  
  names(fea) <- c("l", "q", "p", "t", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh",
                  "pt", "ph", "th", "lqp", "lqt", "lqh", "lpt", "lph", "qpt", "qph",
                  "qth", "pth", "lqpt", "lqph", "lqth", "lpth", "lqpth")
  
  ###Selected feature classes using data from the best models table
  fea <- fea[f.clas]
  
  ##Projection (G) variables folders and subfolders, extrapolation types, and writting clamp and MESS
  if(project == TRUE){
    G.dir <- paste(G.var.dir, var.di, sep = "/")
    G.dirs <- vector()
    for (i in 1:length(G.dir)) {
      dirs <- dir(G.dir[i])
      dires <- vector()
      for (j in 1:length(dirs)) {
        dires[j] <- paste(G.var.dir, var.di[i], dirs[j], sep = "\\")
      }
      G.dirs[i] <- paste("projectionlayers=", paste(dires, collapse = ","), sep = "")
    }
    
    if(ext.type == "ext_clam"){
      mid.com <- "extrapolate=true doclamp=true responsecurves=true jackknife=true"
      ext.nam <- "_EC"
    }
    if(ext.type == "ext"){
      mid.com <- "extrapolate=true doclamp=false responsecurves=true jackknife=true"
      ext.nam <- "_E"
    }
    if(ext.type == "no_ext"){
      mid.com <- "extrapolate=false doclamp=false responsecurves=true jackknife=true"
      ext.nam <- "_NE"
    }
    if(ext.type == "all"){
      mid.com <- c("extrapolate=true doclamp=true responsecurves=true jackknife=true",
                   "extrapolate=true doclamp=false responsecurves=true jackknife=true",
                   "extrapolate=false doclamp=false responsecurves=true jackknife=true")
      ext.nam <- c("_EC", "_E", "_NE")
    }
    
    if(write.mess == FALSE){
      w.mess <- "writeclampgrid=false"
    }
    if(write.clamp == FALSE){
      w.clamp <- "writemess=false"
    }
  }else{
    mid.com <- "extrapolate=false doclamp=false writeclampgrid=false writemess=false responsecurves=true jackknife=true"
  }

  ##Output format
  out <- paste("outputformat=", out.format, sep = "")
  
  ##Number of replicates
  rep <-  paste("replicates=", rep.n, sep = "")
  
  ##Replicate type
  rept <- paste("replicatetype=", rep.type, sep = "")
  
  #Fixed commands
  ##Intitial command
  in.comm <- paste("java", ram, "-jar maxent.jar", sep = " ")
  
  ##Autofeature
  a.fea <- "autofeature=false"
  
  ##Other maxent settings
  fin.com <- "warnings=false visible=false redoifexists autorun\n"
  
  #####
  #Final code
  if(project == TRUE){
    if(write.clamp == FALSE | write.mess == FALSE){
      if(write.clamp == FALSE & write.mess == TRUE){
        pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
        b.nam <- "ku_enm_final_models.bat"
        sink(b.nam)
        
        for (i in 1:length(sett1)) {
          Sys.sleep(0.1)
          setWinProgressBar(pb, i, title=paste(round(i/length(sett1)*100, 0), "% finished"))
          
          for (j in 1:length(ext.nam)) {
            subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
            dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))
            
            reg.m <- paste("betamultiplier=", rm[i], sep = "")
            cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, out, mid.com[j], w.clamp, fin.com, sep = " "))
          }
        }
      }
      if(write.mess == FALSE & write.clamp == TRUE){
        pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
        b.nam <- "ku_enm_final_models.bat"
        sink(b.nam)
        
        for (i in 1:length(sett1)) {
          Sys.sleep(0.1)
          setWinProgressBar(pb, i, title=paste(round(i/length(sett1)*100, 0), "% finished"))
          
          for (j in 1:length(ext.nam)) {
            subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
            dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))
            
            reg.m <- paste("betamultiplier=", rm[i], sep = "")
            cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, out, mid.com[j], w.mess, fin.com, sep = " "))
          }
        } 
      }
      if(write.clamp == FALSE & write.mess == FALSE){
        pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
        b.nam <- "ku_enm_final_models.bat"
        sink(b.nam)
        
        for (i in 1:length(sett1)) {
          Sys.sleep(0.1)
          setWinProgressBar(pb, i, title=paste(round(i/length(sett1)*100, 0), "% finished"))
          
          for (j in 1:length(ext.nam)) {
            subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
            dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))
            
            reg.m <- paste("betamultiplier=", rm[i], sep = "")
            cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, out, mid.com[j], w.clamp, w.mess, fin.com, sep = " "))
          }
        } 
      }
    }else{
      pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
      b.nam <- "ku_enm_final_models.bat"
      sink(b.nam)
      
      for (i in 1:length(sett1)) {
        Sys.sleep(0.1)
        setWinProgressBar(pb, i, title=paste(round(i/length(sett1)*100, 0), "% finished"))
        
        for (j in 1:length(ext.nam)) {
          subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], ext.nam[j], sep = "")
          dir.create(paste(out.dir, "/", sett1[i], ext.nam[j], sep = ""))
          
          reg.m <- paste("betamultiplier=", rm[i], sep = "")
          cat(paste(in.comm, env[i], samp, G.dirs[i], subfol, reg.m, a.fea, fea[i], rep, rept, out, mid.com[j], fin.com, sep = " "))
        }
      }
    }
    
    sink()
    suppressMessages(close(pb))
    
    cat("\nIf asked, allow runing as administrator.")
    shell.exec(file.path(getwd(), "ku_enm_final_models.bat"))
    
    cat("\nProcess finished\n")
    cat(paste("A maxent batch file for creating", length(sett1) * length(ext.nam), "final models and their projections has been written", sep = " "))
    cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
  }else{
    pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sett1), width = 300) #progress bar
    b.nam <- "ku_enm_final_models.bat"
    sink(b.nam)
    
    for (i in 1:length(sett1)) {
      Sys.sleep(0.1)
      setWinProgressBar(pb, i, title=paste(round(i/length(sett1)*100, 0), "% finished"))
      
      subfol <- paste("outputdirectory=", out.dir, "\\", sett1[i], sep = "")
      dir.create(paste(out.dir, "/", sett1[i], sep = ""))
      
      reg.m <- paste("betamultiplier=", rm[i], sep = "")
      cat(paste(in.comm, env[i], samp, subfol, reg.m, a.fea, fea[i], rep, rept, out, mid.com, fin.com, sep = " "))
    }
    
    sink()
    suppressMessages(close(pb))
    
    cat("\nIf asked, allow runing as administrator.")
    shell.exec(file.path(getwd(), "ku_enm_final_models.bat"))
    
    cat("\nProcess finished\n")
    cat(paste("A maxent batch file for creating", length(sett1), "final models without projections has been written", sep = " "))
    cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
  }
}