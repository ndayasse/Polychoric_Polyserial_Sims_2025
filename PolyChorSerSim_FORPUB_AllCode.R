# C.D. Ayasse

# Simulations comparing accuracy of the following correlation methods:
#   (A) Polychoric (2-step & ML) for ordinal-ordinal vs spearman's (and/or pearson's)
#   (B) Polyserial (2-step & ML) for ordinal-continuous vs spearman's (and/or pearson's)


# Clear working environment:
rm(list=ls()) # removes all objects from default environment
gc() # 'garbage collection' - free unused R memory & memory usage report


# Load Packages: --------------------------------------------------------------- ####

# For data wrangling:
library(tidyr) #data wrangling
library(dplyr) #data wrangling
#   Note: loading "tidyverse" loads additional not-needed packages

# For simulating data or analyzing data:
library(faux) #for correlated data and distribution transformations (normal, gamma, and others)
library(polycor) #various types of correlations 
library(moments) #to compute (observed) skewness & kurtosis





# Compute Item Intercepts (Functions): ----------------------------------------- ####




# CENTERED-Normal Intercepts: --- #

choose.norm.sds <- 2.5 

#   Function for computing Centered-Normal Intercepts:
ints.center.fntn <- function(val, num.categ, n.digits=2) {
  vec <- seq(from=(-1*val), to=val, length.out=(num.categ+1))
  ints <- round(vec[c(2:(length(vec)-1))], digits=n.digits)
  return(ints)
}

#   Computed based on X SDs (choose.norm.sds):
ints.center <- vector(mode="list", length=length(nresp.cat.vec))
val <- latent.sd*choose.norm.sds
for (cur.num in nresp.cat.vec) {
  ind <- which(nresp.cat.vec==cur.num)
  ints.center[[ind]] <- ints.center.fntn(val=val, num.categ=cur.num, 
                                         n.digits=n.dg)
  names(ints.center)[ind] <- as.character(cur.num)
}




# SKEWED Intercepts - Function to Compute: --- #

#   Note: gamma rate parameter = sqrt(shape)

ints.skewed.fntn <- function(
    norm.val, num.categ, gamma.shape, gamma.rate, skew.direction="negative", 
    min.percentile=0, max.percentile=1, adj.percentile=0.0001, n.digits=2) {
  
  # NOTE: if @min.percentile is 0 and/or @max.percentile is 1 (exactly), each
  #       will be adjusted by @adj.percentile to avoid -Inf and Inf values 
  #       (actual minimum percentile will be @min.percentile + @adj.percentile 
  #        and actual maximum percentile will be 
  #        @max.percentile - @adj.percentile)
  
  
  # Compute values corresponding to percentiles in gamma distribution:
  if (min.percentile==0) {
    min.val <- min.percentile + adj.percentile
  } else { min.val <- min.percentile }
  if (max.percentile==1) {
    max.val <- max.percentile - adj.percentile
  } else { max.val <- max.percentile }
  perc.seq <- seq(from=min.val, to=max.val, length.out=(num.categ+2))
  gamma.vals <- qgamma(p=perc.seq, shape=gamma.shape, rate=gamma.rate)
  
  # Compute the difference between each:
  diff.vals <- diff(gamma.vals)
  
  # Convert the differences to intercepts:
  start.val <- -1*norm.val
  full.vec <- vector(mode="numeric", length=0)
  for (diff.cur in diff.vals) {
    cur.val <- start.val + diff.cur
    start.val <- cur.val # RESET
    full.vec <- c(full.vec, round(cur.val, digits=n.digits))
  }
  
  # Remove 2 elements (keep the middle - remove first & final):
  ints <- full.vec[c(2:(length(full.vec)-1))]
  
  # If skew direction is positive instead of negative:
  #   (DEFAULT results in NEGATIVE skew intercepts)
  if (tolower(skew.direction)!="negative" & tolower(skew.direction)!="neg") {
    ints <- sort(-1*ints) # NOTE: IMPORTANT TO REMEMBER TO SORT!!! (must be in ascending order)
  }
  
  # Return the vector of intercepts:
  return(ints)
  
} # end of ints.skewed.fntn





# Polychotomize Continuous Variable Into Ordinal Item (Function): -------------- ####

# C.D. Ayasse
# 8 December 2023

# Simple function to create ordinal data (polychotomize continuous into ordinal)
#   just using a vector (continuous variable) and intercepts


# Details on Arguments / Inputs:
# @vector : (default is NULL) an inputted numeric vector to use as the threshold;
#                       if supplied, will ignore norm.mean and norm.sd.
#                       Note: if simulated within function, right now only 
#                       supports normal distribution.
# @norm.mean : (default is 0) if @vector is not supplied, will use this numeric
#                       value as the mean of the simulated threshold vector (normal).
# @norm.sd : (default is 1) if @vector is not supplied, will use this nuemric
#                       value as the SD of the simulated threshold vector (normal).
# @sample.size : (default is 100) if @vector is not supplied, will use this
#                       numeric value as the length of the vector to simulate.
# @item.intercepts : (default is NULL) if supplied, will use the vector of
#                       numeric intercepts to cut @vector into ordinal responses.
# @nresp.cat : (default is NULL) will be ignored if @item.intercepts is supplied;
#                       if @item.intercepts is NULL, will use this value to 
#                       choose item intercepts to produce this number of 
#                       response categories.
# @lowest.ord.value : (default is 1) the value at which the ordinal scale will 
#                       start, typically either 1 or 0.
# @center.ord.value : (default is NULL) if supplied, instead of choosing values
#                       for the ordinal scale based on lowest.ord.value, will
#                       instead choose values such that the resulting scale is
#                       centered on (mean is) center.ord.value.
#                       NOTE: If @nresp.cat is EVEN, and @center.ord.value is
#                       supplied, @center.ord.value MUST be an increment of 0.5
#                       (this way all response options are increments of 1).
# @output.df : (default is FALSE) if FALSE, will output just the new ordinal
#                       vector; if TRUE, will output a dataframe with the
#                       inputted vector and the new ordinal vector as columns.
# NOTE: MUST supply ONE OF @item.intercepts or @nresp.cat (both cannot be NULL)


very.simple.ord.item.fntn <- function(
    vector=NULL, norm.mean=0, norm.sd=1, sample.size=100, 
    item.intercepts=NULL, nresp.cat=NULL, 
    lowest.ord.value=1, center.ord.value=NULL, resp.options=NULL, 
    output.df=FALSE) {
  
  library(dplyr)
  
  # Checking whether to prioritize the lowest.ord.value OR the center.ord.value:
  if (is.null(center.ord.value)) {
    prioritize.lowest.ord.val <- TRUE
    center.ord.value <- 0
  } else {
    prioritize.lowest.ord.val <- FALSE
  }
  
  # Setting the number of response categories to default to if neither 
  #   nresp.cat nor item.intercepts are provided:
  default.nresp.cat <- 5
  
  
  # If not supplied, simulate threshold vector:
  if (is.null(vector)) {
    vector <- rnorm(n=sample.size, mean=norm.mean, sd=norm.sd)
  }
  
  # If not supplied, create item intercepts:
  if (is.null(item.intercepts)) {
    if (is.null(nresp.cat)) {
      warning(paste0("Must supply ONE OF @item.intercepts or @nresp.cat. ",
                     "Neither was supplied. Will default to ",
                     default.nresp.cat," response categories."))
      nresp.cat <- default.nresp.cat
    } 
    
    if (nresp.cat < 9) { half.int.range <- 1.7 } else { half.int.range <- 2.0 }
    max.int <- center.ord.value + half.int.range
    min.int <- center.ord.value - half.int.range
    item.intercepts <- sort(seq(
      from=min.int, to=max.int, length.out=(nresp.cat-1)), decreasing=F)
    
  } else { #end of if item.intercepts is NULL
    # if intercepts are not in correct order the polychotomizing won't work well:
    item.intercepts <- sort(item.intercepts, decreasing=F)
  }
  #   If @item.intercepts was supplied and not @nresp.cat, make sure we calc that:
  if (is.null(nresp.cat)) {
    nresp.cat <- length(item.intercepts) + 1
  }
  
  
  # Checking whether resp.options supplied & ensuring it's the right length:
  if (!is.null(resp.options)) {
    if (length(resp.options)!=nresp.cat) {
      warning(paste0(
        "If supplied, @resp.options MUST be a vector with length equal to ",
        "@nresp.cat. Length of supplied @resp.options is ",length(resp.options),
        " and @nresp.cat is ",nresp.cat,". Will revert to default values for ",
        "@resp.options."))
      resp.options <- NULL
    }
  }
  
  
  
  # Create ordinal vector:
  dataframe <- data.frame("Input"=vector)
  
  #   Make sequential response options - IF resp.options has NOT been supplied (or was the wrong length):
  if (is.null(resp.options)) {
    if (prioritize.lowest.ord.val) {
      resp.opt.vec <- seq(from=lowest.ord.value, length.out=nresp.cat, by=1)
    } else {
      if (nresp.cat %% 2 != 0) { #for ODD values of nresp.cat (e.g., 3, 5, 7)
        lowest.ord.value <- center.ord.value - ((nresp.cat-1)/2)
        resp.opt.vec <- seq(from=lowest.ord.value, length.out=nresp.cat, by=1)
        if (mean(resp.opt.vec) != center.ord.value) {
          stop("Check computation of response options!")
        }
      } else { #for EVEN values of nresp.cat (e.g., 4, 6)
        if (center.ord.value %% 0.5 != 0) {stop(paste0(
          "If nresp.cat is an even number and center.ord.value is supplied, ",
          "center.ord.value MUST be an increment of 0.5 (e.g., 0.5)"))}
        lowest.ord.value <- center.ord.value - ((nresp.cat)/2) + 0.5
        resp.opt.vec <- seq(from=lowest.ord.value, length.out=nresp.cat, by=1)
      }
    }
  } else { resp.opt.vec <- resp.options }
  
  dataframe$New.Ordinal <- NA
  for (resp.opt.cur in resp.opt.vec) {
    index.cur <- which(resp.opt.vec==resp.opt.cur)
    
    if (resp.opt.cur == resp.opt.vec[1]) { #first response option:
      
      dataframe[which(dataframe$Input < item.intercepts[[1]]),"New.Ordinal"] <- resp.opt.cur
      temp.sub.vec <- dataframe[which(dataframe$Input < item.intercepts[[1]]),"Input"]
      
    } else if (resp.opt.cur == resp.opt.vec[length(resp.opt.vec)]) { #last response option:
      
      dataframe[which(dataframe$Input >= 
                        item.intercepts[[length(item.intercepts)]]),
                "New.Ordinal"] <- resp.opt.cur
      temp.sub.vec <- dataframe[which(
        dataframe$Input >= item.intercepts[[length(item.intercepts)]]),"Input"]
      
    } else { #all other response options (middle ones):
      
      dataframe[which(dataframe$Input >= item.intercepts[[(index.cur-1)]] &
                        dataframe$Input < item.intercepts[[index.cur]]),
                "New.Ordinal"] <- resp.opt.cur
      temp.sub.vec <- dataframe[which(
        dataframe$Input >= item.intercepts[[(index.cur-1)]] &
          dataframe$Input < item.intercepts[[index.cur]]),"Input"]
      
    }
    
  } #end of loop through response options vector (resp.opt.vec)
  
  
  # Return output, either as dataframe or just vector:
  if (output.df) {
    return(dataframe)
  } else if (!output.df) {
    ordinal.vec <- dataframe$New.Ordinal
    return(ordinal.vec)
  } 
  
} #end of function [very.simple.ord.item.fntn]





# Main Simulation & Computations (Function): ----------------------------------- ####


# Notes on Arguments / Inputs:
#  @sample.size:  the sample size of the simulated sample to produce 
#                   (i.e., the number of people in the simulated dataset)
#  @data.type:    whether to simulate two ordinal items ('polychoric') or one
#                   ordinal and one continuous item ('polyserial')
#  @corr.pop.val: the population correlation strength between the two 
#                   latent variables
#  @latent.mean:  population mean (mu) of the latent variables [default 0]
#  @latent.sd:    population standard deviation (sigma) of the latent variables
#                   [default 1]
#  @gamma.shape:   gamma shape parameter to use for skewed distributions
#  @gamma.shape.norm: gamma shape parameter to use for 'normal' distributions
#  @skew.dir.1:   values of 0, 1, or -1 ; indicates whether item 1 is skewed 
#                   and if so which direction [1 for positive; -1 for negative]
#  @skew.dir.2:   values of 0, 1, or -1 ; indicates whether item 2 is skewed 
#                   and if so which direction [1 for positive; -1 for negative]
#  @int.vec.1:    vector of the item response thresholds [a.k.a. item 
#                   intercepts] to use for the first item
#  @int.vec.2:    vector of the item response thresholds [a.k.a. item 
#                   intercepts] to use for the second item
#  @int.vec:      vector of the item response thresholds [a.k.a. item 
#                   intercepts] to use for the ordinal item for ordinal-
#                   continuous ['polyserial'] data
#  @incl.pears:   indicates whether or not to compute Pearson's R; if TRUE, 
#                   Pearson's R will be computed, if FALSE it will not 
#                   [default FALSE]
#  @run.ml.poly:  indicates whether or not to compute the maximum likelihood 
#                   (ML) version of the polychoric or polyserial correlation, in
#                   addition to the 2-step version which will always be 
#                   calculated; if TRUE, ML version is computed, if FALSE it 
#                   will not [default TRUE]


poly.chor.ser.fntn <- function(
    sample.size, data.type, corr.pop.val, latent.mean=0, latent.sd=1, 
    gamma.shape=0, gamma.shape.norm=999, skew.dir.1=0, skew.dir.2=0, 
    int.vec.1=NULL, int.vec.2=NULL, int.vec=NULL, 
    incl.pears=FALSE, run.ml.poly=TRUE ) { 
  
  
  
  
  # Sim using correlated continuous Latent vars: --- #
  
  poly.dat <- data.frame("USUBJID"=c(1:sample.size))
  
  
  # Start w/ 2 normal latent vars that are correlated - regardless of condition:
  temp.dat <- faux::rnorm_multi(n=sample.size, vars=2, mu=latent.mean, 
                                sd=latent.sd, r=corr.val, empirical=FALSE, 
                                varnames=c("OrigLatent_1","OrigLatent_2"))
  
  # Join latent variables with the USUBJID's:
  poly.dat <- cbind(poly.dat, temp.dat)
  
  
  
  
  # Conduct gamma transformations of latent variables: --- #
  
  #   Compute Gamma Rates:
  gamma.rate.norm <- sqrt(gamma.shape.norm/(latent.sd^2)) # 0 # 
  gamma.rate <- sqrt(gamma.shape/(latent.sd^2))
  
  #   For 'normal' latent vars - Transform using faux package:
  if (skew.val.1==0) {
    
    # Transform:
    poly.dat$Latent_1 <- faux::norm2gamma( 
      x=poly.dat$OrigLatent_1, mu=latent.mean, sd=latent.sd, 
      shape=gamma.shape.norm, rate=gamma.rate.norm)
    
    # Re-center at intended population mean:
    obs1.mean <- mean(poly.dat$Latent_1, na.rm=T)
    obsmexp1.mean <- obs1.mean - latent.mean
    poly.dat$Latent_1 <- poly.dat$Latent_1 - obsmexp1.mean
    
  }
  if (skew.val.2==0) {
    
    # Transform:
    poly.dat$Latent_2 <- faux::norm2gamma( 
      x=poly.dat$OrigLatent_2, mu=latent.mean, sd=latent.sd, 
      shape=gamma.shape.norm, rate=gamma.rate.norm)
    
    # Re-center at intended population mean:
    obs2.mean <- mean(poly.dat$Latent_2, na.rm=T)
    obsmexp2.mean <- obs2.mean - latent.mean
    poly.dat$Latent_2 <- poly.dat$Latent_2 - obsmexp2.mean
    
  }
  
  #   For skewed gamma latent vars - Transform using faux package:
  if (skew.val.1!=0) {
    
    # Transform to POSITIVE skewed (gamma) using faux package:
    poly.dat$Latent_1 <- faux::norm2gamma( 
      x=poly.dat$OrigLatent_1, mu=latent.mean, sd=latent.sd, # pop mean and SD of X (input var)
      shape=gamma.shape, rate=gamma.rate)
    
    # If want NEGATIVE Skewed:
    if (skew.val.1<0) { poly.dat$Latent_1 <- (-1)*poly.dat$Latent_1 }
    
  } 
  if (skew.val.2!=0) {
    
    # Transform to POSITIVE skewed (gamma) using faux package:
    poly.dat$Latent_2 <- faux::norm2gamma( 
      x=poly.dat$OrigLatent_2, mu=latent.mean, sd=latent.sd, # pop mean and SD of X (input var)
      shape=gamma.shape, rate=gamma.rate)
    
    # If want NEGATIVE Skewed:
    if (skew.val.2<0) { poly.dat$Latent_2 <- (-1)*poly.dat$Latent_2 }
    
  } 
  
  
  
  
  # Create the ordinal item(s) from latent vars: --- #
  
  if (distrib.cur=="Polychoric") {
    
    # Item 1 based on Latent_1:
    ord.item1 <- very.simple.ord.item.fntn(
      vector=poly.dat$Latent_1, item.intercepts=int.vec.1, output.df=FALSE)
    
    # Item 2 based on Latent_2:
    ord.item2 <- very.simple.ord.item.fntn(
      vector=poly.dat$Latent_2, item.intercepts=int.vec.2, output.df=FALSE)
    
    # Join both onto the main df:
    poly.dat$Item_1 <- ord.item1
    poly.dat$Item_2 <- ord.item2
    
  } else if (distrib.cur=="Polyserial") {
    
    # Item 1 based on Latent_1:
    ord.item1 <- very.simple.ord.item.fntn(
      vector=poly.dat$Latent_1, item.intercepts=int.vec, output.df=FALSE)
    #   Join onto the main df:
    poly.dat$Item_1 <- ord.item1
    
    # Item 2 IS Latent_2:
    poly.dat$Item_2 <- poly.dat$Latent_2
    
  }
  
  
  
  
  # Compute Correlations (All): --- #
  
  
  
  # Between Latent Variables:
  
  #   Correlations between latent variables underlying ordinal items:
  latent.pears <- cor.test(poly.dat$Latent_1, poly.dat$Latent_2, 
                           method="pearson")
  latent.spear <- cor.test(poly.dat$Latent_1, poly.dat$Latent_2, 
                           method="spearman")
  latent.pears.r <- latent.pears$estimate
  latent.spear.r <- latent.spear$estimate
  
  #   Correlations between *original* continuous-normal latent variables:
  origltnt.pears <- cor.test(poly.dat$OrigLatent_1, poly.dat$OrigLatent_2, 
                             method="pearson")
  origltnt.spear <- cor.test(poly.dat$OrigLatent_1, poly.dat$OrigLatent_2, 
                             method="spearman")
  origltnt.pears.r <- origltnt.pears$estimate
  origltnt.spear.r <- origltnt.spear$estimate
  
  #   Correls b/w *original* latent vars & their transformed versions:
  origvltnt.v1.pears <- cor.test(poly.dat$OrigLatent_1, poly.dat$Latent_1, 
                                 method="pearson")
  origvltnt.v1.spear <- cor.test(poly.dat$OrigLatent_1, poly.dat$Latent_1, 
                                 method="spearman")
  origvltnt.v2.pears <- cor.test(poly.dat$OrigLatent_2, poly.dat$Latent_2, 
                                 method="pearson")
  origvltnt.v2.spear <- cor.test(poly.dat$OrigLatent_2, poly.dat$Latent_2, 
                                 method="spearman")
  origvltnt.v1.pears.r <- origvltnt.v1.pears$estimate
  origvltnt.v1.spear.r <- origvltnt.v1.spear$estimate
  origvltnt.v2.pears.r <- origvltnt.v2.pears$estimate
  origvltnt.v2.spear.r <- origvltnt.v2.spear$estimate
  
  #   Concordance ICC b/w *original* latent vars & their transformed versions:
  
  #     Var 1:
  v1.df <- dplyr::select(
    poly.dat, c("OrigLatent_1","Latent_1")) %>% as.data.frame()
  origvltnt.v1.icc.all <- psych::ICC(x=v1.df)
  origvltnt.v1.icc.c1 <- origvltnt.v1.icc.all$results[which(
    origvltnt.v1.icc.all$results$type=="ICC3"),"ICC"]
  origvltnt.v1.icc.c1.locl <- origvltnt.v1.icc.all$results[which(
    origvltnt.v1.icc.all$results$type=="ICC3"),"lower bound"]
  origvltnt.v1.icc.c1.upcl <- origvltnt.v1.icc.all$results[which(
    origvltnt.v1.icc.all$results$type=="ICC3"),"upper bound"]
  
  #     Var 2:
  v2.df <- dplyr::select(
    poly.dat, c("OrigLatent_2","Latent_2")) %>% as.data.frame()
  origvltnt.v2.icc.all <- psych::ICC(x=v2.df)
  origvltnt.v2.icc.c1 <- origvltnt.v2.icc.all$results[which(
    origvltnt.v2.icc.all$results$type=="ICC3"),"ICC"]
  origvltnt.v2.icc.c1.locl <- origvltnt.v2.icc.all$results[which(
    origvltnt.v2.icc.all$results$type=="ICC3"),"lower bound"]
  origvltnt.v2.icc.c1.upcl <- origvltnt.v2.icc.all$results[which(
    origvltnt.v2.icc.all$results$type=="ICC3"),"upper bound"]
  
  
  
  # Between Items:
  
  
  # Spearman:
  spear.item <- cor.test(poly.dat$Item_1, poly.dat$Item_2, 
                         method="spearman")
  spear.item.r <- spear.item$estimate
  
  
  # Pearson:
  if (incl.pears) {
    pears.item <- cor.test(poly.dat$Item_1, poly.dat$Item_2, 
                           method="pearson")
    pears.item.r <- pears.item$estimate
  } else { pears.item.r <- NA }
  
  
  # Polychoric:
  if (distrib.cur=="Polychoric") {
    
    # Via Polycor Package (fewer issues/errors than Psych package):
    polychor.results <- polycor::polychor(
      poly.dat$Item_1, poly.dat$Item_2, ML=FALSE, std.err=FALSE)
    poly.r.2stp <- polychor.results
    
    # If run.ml.polycor is TRUE, also run via ML:
    if (run.ml.polycor) {
      polychor.results.ml <- polycor::polychor(
        poly.dat$Item_1, poly.dat$Item_2, ML=TRUE, std.err=FALSE)
      poly.r.ml <- polychor.results.ml
    } else { poly.r.ml <- NA }
    
  } #end if polychoric - calc polychoric correlations
  
  
  # Polyserial:
  if (distrib.cur=="Polyserial") {
    
    # Via Polycor Package (fewer issues/errors than Psych package):
    polyser.results <- polycor::polyserial(
      poly.dat$Item_2, poly.dat$Item_1, ML=FALSE, std.err=FALSE)
    poly.r.2stp <- polyser.results
    
    # If run.ml.polycor is TRUE, also run via ML:
    if (run.ml.polycor) {
      polyser.results.ml <- polycor::polyserial(
        poly.dat$Item_2, poly.dat$Item_1, ML=TRUE, std.err=FALSE)
      poly.r.ml <- polyser.results.ml
    } else { poly.r.ml <- NA }
    
    
  } #end if polyserial - calc polyserial correlations
  
  
  
  
  # Calculate Observed Skewness & Kurtosis of Latent Variables & Items: --- #
  
  obs.skew.ltnt1 <- moments::skewness(poly.dat$Latent_1)
  obs.skew.ltnt2 <- moments::skewness(poly.dat$Latent_2)
  obs.skew.item1 <- moments::skewness(poly.dat$Item_1)
  obs.skew.item2 <- moments::skewness(poly.dat$Item_2)
  obs.skew.oltnt1 <- moments::skewness(poly.dat$OrigLatent_1)
  obs.skew.oltnt2 <- moments::skewness(poly.dat$OrigLatent_2)
  
  obs.kurt.ltnt1 <- moments::kurtosis(poly.dat$Latent_1)
  obs.kurt.ltnt2 <- moments::kurtosis(poly.dat$Latent_2)
  obs.kurt.item1 <- moments::kurtosis(poly.dat$Item_1)
  obs.kurt.item2 <- moments::kurtosis(poly.dat$Item_2)
  obs.kurt.oltnt1 <- moments::kurtosis(poly.dat$OrigLatent_1)
  obs.kurt.oltnt2 <- moments::kurtosis(poly.dat$OrigLatent_2)
  
  
  
  
  # Return results: --- #
  
  results <- list(
    
    # Data Type:
    "Data.Type" = data.type, 
    
    # Observed Correlations Results:
    "Latent.Pearson" = latent.pears.r,
    "Latent.Spearman" = latent.spear.r,
    "Item.Pearson" = pears.item.r,
    "Item.Spearman" = spear.item.r,
    "Item.PolyCorr.2step" = poly.r.2stp,
    "Item.PolyCorr.ML" = poly.r.ml,
    #   Additional checking via correlations:
    "LatentOrigNorm.Pearson" = origltnt.pears.r,
    "LatentOrigNorm.Spearman" = origltnt.spear.r,
    "LatentNormVTrsGam.Var1.Pearson" = origvltnt.v1.pears.r, 
    "LatentNormVTrsGam.Var1.Spearman" = origvltnt.v1.spear.r, 
    "LatentNormVTrsGam.Var2.Pearson" = origvltnt.v2.pears.r, 
    "LatentNormVTrsGam.Var2.Spearman" = origvltnt.v2.spear.r, 
    #   Additional checking via ICC:
    "LatentNormVTrsGam.Var1.ICCc1" = origvltnt.v1.icc.c1,
    "LatentNormVTrsGam.Var1.ICCc1.locl" = origvltnt.v1.icc.c1.locl,
    "LatentNormVTrsGam.Var1.ICCc1.upcl" = origvltnt.v1.icc.c1.upcl, 
    "LatentNormVTrsGam.Var2.ICCc1" = origvltnt.v2.icc.c1, 
    "LatentNormVTrsGam.Var2.ICCc1.locl" = origvltnt.v2.icc.c1.locl, 
    "LatentNormVTrsGam.Var2.ICCc1.upcl" = origvltnt.v2.icc.c1.upcl,
    
    # Observed Distribution Skewness & Kurtosis Results:
    #   Skewness:
    "Ltnt1.ObsSkew" = obs.skew.ltnt1,
    "Ltnt2.ObsSkew" = obs.skew.ltnt2,
    "Item1.ObsSkew" = obs.skew.item1,
    "OLtnt1.ObsSkew" = obs.skew.oltnt1,
    "OLtnt2.ObsSkew" = obs.skew.oltnt2,
    #   Kurtosis:
    "Ltnt1.ObsKurt" = obs.kurt.ltnt1,
    "Ltnt2.ObsKurt" = obs.kurt.ltnt2,
    "Item1.ObsKurt" = obs.kurt.item1,
    "Item2.ObsKurt" = obs.kurt.item2,
    "OLtnt1.ObsKurt" = obs.kurt.oltnt1,
    "OLtnt2.ObsKurt" = obs.kurt.oltnt2
  )
  
  return(results)
  
  
} #end function [poly.chor.ser.fntn]




