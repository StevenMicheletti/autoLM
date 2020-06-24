print("R ALERT: Checking for R dependencies. Will attempt to install automatically if not present.")

#Install packages if they don't exist already
	list.of.packages <- c("data.table", "nlme","MuMIn","lme4")
	new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
	if(length(new.packages)) print("R ALERT: Installing dependencies for first time use....")
	if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

suppressMessages(require(lme4))
suppressMessages(require(nlme))
suppressMessages(require(MuMIn))
suppressMessages(require(data.table))

args <- commandArgs()

fnameo = args[6]
env.input = args[7]
freq.input = args[8]
standardize = args[9]
coords.input = args[10]


if (is.na(fnameo) ==  T){
 print("Not running via shell wrapper, variables must be included in Rscript")
 fnameo = 'autoLM.out'
 env.input = 'env.txt'
 freq.input = 'all_freq.txt'
 standardize = TRUE
 coords.input = 'locations.txt' # set to "nope" if not interested
}



autolm <- function(fnameo,
                  env.input ,
                  freq.input,
                  standardize = TRUE,
                  coords.input = NULL) {

  #standardized env pop variables 
  e.data <-  fread(env.input, sep='\t', header=TRUE)
  e.data <- as.data.frame(e.data)
  if (standardize == TRUE | standardize == "TRUE"){
    e.data <- data.frame(scale(e.data))
  }
  
  # Population, clust, east, north 
  if (file.exists(coords.input)) {
    c.data <- fread(coords.input, header=TRUE)
    if (ncol(c.data) < 4) stop ("Coordinate file has too many columns. 
              Make sure it is tab delimited and includes: Pop Name, genetic cluster, east (m), north (m)")
    c.data <- as.data.frame(c.data)
    colnames(c.data) = c("ID", "clust", "east", "north")
    if (min(c.data$east) < 0 | (min(c.data$north) < 0) ) stop ("Coordinates must be in meters!")   
    east  = c.data$east
    north = c.data$north
    clust = c.data$clust
    v= (as.numeric(length(e.data)))
    Var <- list()
      for (i in 1:v) {
        Var[[i]] <-e.data[, i]
    }
    names(Var) <- colnames(e.data)
    if (length(Var) < 1) stop ("No environmental variables loaded. Check input file")
  }
  
  # population, minor allele fz with row and column headings
  f.data <-  data.frame(fread(freq.input, header=TRUE), row.names=1)
  
  f.data <- as.data.frame(t(f.data))
  f.data[1,] <-as.character(f.data[1,])
  #f.data <- sapply( f.data, as.numeric )

  vpop = as.numeric(nrow(f.data))
  epop = as.numeric(ncol(f.data))
  env = as.numeric(ncol(e.data))
  
  Snp <- list()
  for (j in 1:epop) {
    Snp[[j]] <-f.data[, j]
  }
  
 
  names(Snp) <-colnames(f.data)
  cfile = paste0(getwd(),"/",fnameo)
  
  if (file.exists(fnameo)){
     print(" Warning: Output file exists and is being overwritten!")
  }
  
  out.headerz <-  noquote(paste("snp", "env_var", "pval", "r2", "max_freq_diff", "mean_freq", "score", sep= '\t'))
  write.table(out.headerz, file= fnameo, sep='\t', col.names = F, row.names = F, quote = F)
  

  if (file.exists(coords.input)) {            
    print ("Coordinates and clusters provided as random effects: Running model")
    for (i in 1:env) { 
     test <- as.vector(Var[[i]])
     ci <- names(Var[i])
     for (j in 1:vpop) { 
        tryCatch({
        cn <- row.names(f.data[j,])
        freq <- as.numeric(f.data[j,])
        #model
        mod.gau <- lme(fixed = freq ~ test , random = ~ 1 | clust, method = "ML", 
                       correlation = corGaus (1, form = ~ east + north), control = lmeControl(returnObject = TRUE))
        #Print
        me <-mean(freq)
        mxmi <-  abs((max(freq) - min(freq)))
        mo <- median(max(freq), min(freq))
        mdiff <- mean(abs(diff(freq)))
        spr = mean(as.numeric(diff(quantile(freq))))
        s1 <- summary(mod.gau)$tTable[2,5]
        s3 <- r.squaredLR(mod.gau, null.RE = TRUE)[1]
        scorez <- ((1-s1) + (mxmi) + 1-abs(0.5 - me - 0.5) + (s3*7) + (spr * 4)) / 10.5
        s4 <- cbind (cn,ci,s1,s3, mxmi, me, scorez)
        s5 <-as.data.frame(s4)
        write.table(s5, file=fnameo, sep='\t',
                    col.names = F, row.names = F, quote=FALSE, append=TRUE)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      } 
    }
  }
  

  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
  if (is.null(coords.input) | (!file.exists(coords.input))) {
    print ("No coordinates provided, there will be no random effects: Running model")
    for (i in 1:env) { 
      test <- as.vector(e.data[[i]])
      ci <- colnames(e.data[i])
      for (j in 1:vpop) { 
        tryCatch({
          cn <- row.names(f.data[j,])
          freq <- as.numeric(f.data[j,])
          #model
          mod.gau <- lm(freq ~ test, method = "qr")
          #Print
          me <-mean(freq)
          mxmi <-  abs((max(freq) - min(freq)))
          mo <- median(max(freq), min(freq))
          mdiff <- mean(abs(diff(freq)))
          spr = mean(as.numeric(diff(quantile(freq))))
          s1 <- lmp(mod.gau)
          s3 <- summary(mod.gau)$r.squared
          scorez <- ((1-s1) + (mxmi) + 1-abs(0.5 - me - 0.5) + (s3*7) + (spr * 4)) / 10.5
          s4 <- cbind (cn,ci,s1,s3, mxmi, me, scorez)
          s5 <-as.data.frame(s4)
          write.table(s5, file= fnameo, sep='\t',
                      col.names = F, row.names = F, quote=FALSE, append=TRUE)
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      } 
    }
  }
  }

autolm(fnameo, env.input, freq.input, standardize, coords.input)
