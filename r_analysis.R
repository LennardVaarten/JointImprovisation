#__________________________________________________________________________


# Loading packages, data, etc. --------------------------------------------

#__________________________________________________________________________


# Installing libraries
install.packages("ggplot2")
install.packages("RTransferEntropy")
install.packages("rEDM")
install.packages("devtools")
install.packages("lme4")
install.packages("lmerTest")
install.packages("MuMIn")
install.packages("glme")
install.packages("future")
install.packages("performance")
install.packages("stringr")

# Importing libraries
library(ggplot2)
library(RTransferEntropy)
library(rEDM)
library(devtools)
library(lme4)
library(lmerTest)
library(MuMIn)
library(glme)
library(future)
library(performance)
library(stringr)
 
# Importing our sliding window-averaged data
rms_timeseries = read.csv("C:/Shortcutsensei/0. Thesis/rms_timeseries.csv")
tonalcentroid_timeseries_raw = read.csv("C:/Shortcutsensei/0. Thesis/tonal_centroid_timeseries.csv")
spectralflatness_timeseries = read.csv("C:/Shortcutsensei/0. Thesis/spectral_flatness_timeseries.csv")

# Importing datasets with info about the trials
musicianDataset = read.csv("C:/Users/lenna/Downloads/study_3_4/study_3_4/endImpro_4R_booths_3_trio_12.csv", sep=";")
trialDataset = read.csv("C:/Shortcutsensei/JointImprovisation/Data/endingDataset_fromThomas_withconditions.csv", sep=";") 
listenerDataset = read.csv("C:/Shortcutsensei/JointImprovisation/Data/Appreciation_Data_N46.csv", sep=";")

View(musicianDataset)
View(trialDataset)
View(listenerDataset)

#_____________________________________________________________________________________

# Further preprocessing of Tonnetz data ----------------------------------------------

#_____________________________________________________________________________________

# We compute a 'tonnetz distance' measure. At each window, it computes the Euclidean distance between
# the current and the previous tonnetz (both of which are 6-dimensional vectors).
tonnetzdistance_timeseries = data.frame("g1_t1_b1" = rep(0, 3000))

for(group in 1:12) {
  for(trial in 1:16) {
    for(booth in 1:3) {
      
      # Create string featuring the name of the current column
      colname = sprintf("g%s_t%s_b%s", group, trial, booth)
      # Initialize empty vector, which will contain distance values between each windows tonnetz TONw and 
      # the previous window's tonnetz, TONw-1.
      distVec = c()
      
      curBooth = tonalcentroid_timeseries_raw[,grepl(colname, names(tonalcentroid_timeseries_raw))]
      
      for(row in 1:3000){
        tonnetz = unname(unlist(curBooth[row,]))
        if(row==1){
          # for the very first window's tonnetz, compute the tonnetz distance as simply the distance from a tonnetz of only 0's
          distance = dist(rbind(tonnetz, rep(0, 6)))[1]
        } else {
          # for all subsequent windows, compute the euclidean distance between the current and previous tonnetz. 
          # This represents a measure of harmonic change from one window to the next.
          prevTonnetz = unname(unlist(curBooth[row-1,]))
          distance = dist(rbind(tonnetz, prevTonnetz))[1]
        }
        # Append to our vector of distances
        distVec = c(distVec, distance)
      }
      # Add distance to dataframe
      tonnetzdistance_timeseries[[colname]] = distVec
    }
  }
  # Logging
  print(sprintf("Preprocess Tonnetz g%s", group))
}

#_____________________________________________________________________________________

# TE for real pairs / trios and random pairs ------------

#_____________________________________________________________________________________

# REAL PAIRS + TRIALS
# Calculate TE values for all 'real pairs' of musicians; i.e., same group, same trial. Also calculate
# TE for trios by computing mean of pairwise TE for each pair in a trio. This will come in handy in RQ2.

# Initialize empty vectors
rms_TEvalues_pertrial = list()
spectralflatness_TEvalues_pertrial = list()
tonnetzdistance_TEvalues_pertrial = list()

rms_TEvalues_perpair = list()
spectralflatness_TEvalues_perpair = list()
tonnetzdistance_TEvalues_perpair = list()

# enable parallel processing for all future transfer_entropy calls
plan(multisession)

for(group in 1:12){
  for(trial in 1:16){
      # Get all possible combinations of booths: booth 1-2, booth 1-3, booth 2-3
      pairs = combn(c(1:3), m=2)
      
      trial_rms_TE = 0
      trial_spectralflatness_TE = 0
      trial_tonnetzdistance_TE = 0
      
      for(pair in 1:3){
        
        colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
        colname2 = sprintf("g%s_t%s_b%s", group, trial, pairs[2,pair])
        
        # Calculate transfer entropy for all possible pairs of RMS amplitude time series within the trial
        rmsTE = transfer_entropy(rms_timeseries[,colname1],
                             rms_timeseries[,colname2], shuffles=100, nboot=1, quiet=TRUE, seed=TRUE)
        rmsTExy = rmsTE$coef[1]
        rmsTEyx = rmsTE$coef[2]
        trial_rms_TE = sum(trial_rms_TE, rmsTExy, rmsTEyx)
        rms_TEvalues_perpair[sprintf("%s-%s", colname1, colname2)] = rmsTExy
        rms_TEvalues_perpair[sprintf("%s-%s", colname2, colname1)] = rmsTEyx

        # Calculate transfer entropy for all possible pairs of tonnetz distance time series within the trial
        tonnetzdistanceTE = transfer_entropy(tonnetzdistance_timeseries[,colname1], 
                                             tonnetzdistance_timeseries[,colname2], shuffles=100, nboot=1, quiet=TRUE, seed=TRUE)
        tonnetzdistanceTExy = tonnetzdistanceTE$coef[1]
        tonnetzdistanceTEyx = tonnetzdistanceTE$coef[2]
        trial_tonnetzdistance_TE = sum(trial_tonnetzdistance_TE, tonnetzdistanceTExy, tonnetzdistanceTEyx)
        tonnetzdistance_TEvalues_perpair[sprintf("%s-%s", colname1, colname2)] = tonnetzdistanceTExy
        tonnetzdistance_TEvalues_perpair[sprintf("%s-%s", colname2, colname1)] = tonnetzdistanceTEyx
        
        # Calculate transfer entropy for all possible pairs of spectral flatness time series within the trial
        spectralflatnessTE = transfer_entropy(spectralflatness_timeseries[,colname1],
                                 spectralflatness_timeseries[,colname2], shuffles=100, nboot=1, quiet=TRUE, seed=TRUE)
        spectralflatnessTExy = spectralflatnessTE$coef[1]
        spectralflatnessTEyx = spectralflatnessTE$coef[2]
        trial_spectralflatness_TE = sum(trial_spectralflatness_TE, spectralflatnessTExy, spectralflatnessTEyx)
        spectralflatness_TEvalues_perpair[sprintf("%s-%s", colname1, colname2)] = spectralflatnessTExy
        spectralflatness_TEvalues_perpair[sprintf("%s-%s", colname2, colname1)] = spectralflatnessTEyx
        
      }
      
      # We have summed 6 TE values per trial (booth 1-2, 1-3, 2-3, 2-1, and so on). Divide by 6 to
      # get a trial-wide average TE value.
      rms_TEvalues_pertrial[sprintf("g%s_t%s", group, trial)] = trial_rms_TE / 6
      tonnetzdistance_TEvalues_pertrial[sprintf("g%s_t%s", group, trial)] = trial_tonnetzdistance_TE / 6
      spectralflatness_TEvalues_pertrial[sprintf("g%s_t%s", group, trial)] = trial_spectralflatness_TE / 6
      
  }
  # Logging
  print(sprintf("TE REAL g%s", group))
}

# RANDOM PAIRS
# Calculate TE values for 'random pairs' of musicians; i.e., same group, different trial.

rms_TEvalues_random = c()
spectralflatness_TEvalues_random = c()
tonnetzdistance_TEvalues_random = c()

for(group in 1:12){
  for(trial in 1:16){
    pairs = combn(c(1:3), m=2)
    for(pair in 1:3){
      
      # Create random pairs by matching g1_t1_b1 to g1_t2_b2, g1_t1_b1 to g1_t2_b3, g1_t1_b2 to g1_t2_b3, and so on.
      trialRandom = if (trial != 16) trial+1 else 1
      colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
      colname2 = sprintf("g%s_t%s_b%s", group, trialRandom, pairs[2,pair])
      
      # Calculate transfer entropy for pairs of RMS amplitude time series
      rmsTE = transfer_entropy(rms_timeseries[,colname1],
                               rms_timeseries[,colname2], shuffles=1, nboot=1, quiet=TRUE, seed=TRUE)
      rmsTExy = rmsTE$coef[1]
      rmsTEyx = rmsTE$coef[2]

      rms_TEvalues_random = c(rms_TEvalues_random, rmsTExy, rmsTEyx)
      
      # Calculate transfer entropy for all possible pairs of tonnetz distance time series within the trial
      tonnetzdistanceTE = transfer_entropy(tonnetzdistance_timeseries[,colname1], 
                                           tonnetzdistance_timeseries[,colname2], shuffles=100, nboot=1, quiet=TRUE, seed=TRUE)
      tonnetzdistanceTExy = tonnetzdistanceTE$coef[1]
      tonnetzdistanceTEyx = tonnetzdistanceTE$coef[2]
      
      tonnetzdistance_TEvalues_random = c(tonnetzdistance_TEvalues_random, tonnetzdistanceTExy, tonnetzdistanceTEyx)

      # Calculate transfer entropy for all possible pairs of spectral flatness time series within the trial
      spectralflatnessTE = transfer_entropy(spectralflatness_timeseries[,colname1],
                                            spectralflatness_timeseries[,colname2], shuffles=1, nboot=1, quiet=TRUE, seed=TRUE)
      spectralflatnessTExy = spectralflatnessTE$coef[1]
      spectralflatnessTEyx = spectralflatnessTE$coef[2]
      
      spectralflatness_TEvalues_random = c(spectralflatness_TEvalues_random, spectralflatnessTExy, spectralflatnessTEyx)
      
    }
  }
  # Logging
  print(sprintf("TE RANDOM g%s", group))
}

# Get just the values from each key-value pair in the list, so that we can perform our statistical tests
rms_TEvalues_raw = unname(unlist(rms_TEvalues_perpair))
tonnetzdistance_TEvalues_raw = unname(unlist(tonnetzdistance_TEvalues_perpair))
spectralflatness_TEvalues_raw = unname(unlist(spectralflatness_TEvalues_perpair))

#_____________________________________________________________________________________

# Rho for real trios and random trios ------------

#_____________________________________________________________________________________

# REAL TRIOS

rms_rhovalues = list()
tonnetzdistance_rhovalues = list()
spectralflatness_rhovalues = list()

# List which contains, for each trial, the window at which at least one of the musicians has
# stopped playing.
endPoints = list()

for(group in 1:12){
  for(trial in 1:16){
    
    colname1 = sprintf("g%s_t%s_b1", group, trial)
    colname2 = sprintf("g%s_t%s_b2", group, trial)
    colname3 = sprintf("g%s_t%s_b3", group, trial)
    
    colnames = c(colname1, colname2, colname3)
    
    endPoint = min(c(which(is.na(rms_timeseries[colname1]))[1], 
                      which(is.na(rms_timeseries[colname2]))[1], 
                      which(is.na(rms_timeseries[colname3]))[1]))
    
    endPoints[sprintf("g%s_t%s", group, trial)] = endPoint
    
    # Calculate rho for all real pairs of RMS amplitude time series
    
    rmsDF = cbind(rms_timeseries[colname1], rms_timeseries[colname2], rms_timeseries[colname3])
    
    trialRMSRho = 0
    
    # For each trial, calculate three values of rho, each with a different booth as the target
    for(colname in colnames){
      targRhos = c()
      # Calculate rho values from 1 lag up to 20 lags
      for(forecast_time in 1:20){
        simplex_output = block_lnlp(rmsDF, lib = c(1, endPoint), pred = c(1, endPoint),
                                     target_column = colname, tp = forecast_time, exclusion_radius=forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        targRhos = c(targRhos, simplex_output$stats$rho$rho)
      }
      trialRMSRho = trialRMSRho + mean(targRhos)
    }
    
    # Calculate rho for all real pairs of spectral flatness time series
    
    spectralFlatnessDF = cbind(spectralflatness_timeseries[colname1], spectralflatness_timeseries[colname2], spectralflatness_timeseries[colname3])
    
    trialSpectralFlatnessRho = 0
    
    for(colname in colnames){
      targRhos = c()
      for(forecast_time in 1:20){
        simplex_output = block_lnlp(spectralFlatnessDF, lib = c(1, endPoint), pred = c(1, endPoint),
                                    target_column = colname, tp = forecast_time, exclusion_radius=forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        targRhos = c(targRhos, simplex_output$stats$rho$rho)
      }
      trialSpectralFlatnessRho = trialSpectralFlatnessRho + mean(targRhos)
    }
    
    # Calculate rho for all real pairs of tonnetz distance time series
    
    tonnetzdistanceDF = cbind(tonnetzdistance_timeseries[colname1], tonnetzdistance_timeseries[colname2], tonnetzdistance_timeseries[colname3])
    
    trialTonnetzDistanceRho = 0
    
    for(colname in colnames){
      targRhos = c()
      for(forecast_time in 1:20){
        simplex_output = block_lnlp(tonnetzdistanceDF, lib = c(1, endPoint), pred = c(1, endPoint),
                                    target_column = colname, tp = forecast_time, exclusion_radius=forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        targRhos = c(targRhos, simplex_output$stats$rho$rho)
      }
      trialTonnetzDistanceRho = trialTonnetzDistanceRho + mean(targRhos)
    }
    
    # The variables trialRMSrho, trialSpectralFlatnessRho and trialTonnetzDistancRho now contain the sum of the rho
    # values for the three different targets. We divide by 3 to get final rho values for the trial.
    rms_rhovalues[sprintf("g%s_t%s", group, trial)] = trialRMSRho / 3
    spectralflatness_rhovalues[sprintf("g%s_t%s", group, trial)] = trialSpectralFlatnessRho / 3
    tonnetzdistance_rhovalues[sprintf("g%s_t%s", group, trial)] = trialTonnetzDistanceRho / 3
  }
  print(sprintf("EDM REAL g%s", group))
}

# RANDOM TRIOS

# Initialize empty vectors
rms_rhovalues_random = c()
spectralflatness_rhovalues_random = c()
tonnetzdistance_rhovalues_random = c()

for(group in 1:12){
  for(trial in 1:16){

    # Same procedure as before for creating random pairs of musicians
    trialRandom1 = if (trial != 16) trial+1 else 1
    trialRandom2 = if (trial != 1) trial-1 else 16
    
    colname1 = sprintf("g%s_t%s_b1", group, trial)
    colname2 = sprintf("g%s_t%s_b2", group, trialRandom1)
    colname3 = sprintf("g%s_t%s_b3", group, trialRandom2)
    
    colnames = c(colname1, colname2, colname3)
    
    # This gives us the point, for this trial, at which at least 1 musician has stopped playing for good.
    endPoint = min(c(which(is.na(rms_timeseries[colname1]))[1], 
                     which(is.na(rms_timeseries[colname2]))[1], 
                     which(is.na(rms_timeseries[colname3]))[1]))
    
    # Calculate rho for random pairs of RMS amplitude time series

    rmsDF = cbind(rms_timeseries[colname1], rms_timeseries[colname2], rms_timeseries[colname3])

    trialRMSRho = 0

    for(colname in colnames){
      targRhos = c()
      for(forecast_time in 1:20){
        simplex_output = block_lnlp(rmsDF, lib = c(1, endPoint), pred = c(1, endPoint),
                                    target_column = colname, tp = forecast_time, exclusion_radius=forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        targRhos = c(targRhos, simplex_output$stats$rho$rho)
      }
    }
    trialRMSRho = trialRMSRho + mean(targRhos)
    
    # Calculate rho for random pairs of tonnetz distance time series
    
    tonnetzdistanceDF = cbind(tonnetzdistance_timeseries[colname1], tonnetzdistance_timeseries[colname2], tonnetzdistance_timeseries[colname3])
    
    trialTonnetzDistanceRho = 0
    
    for(colname in colnames){
      targRhos = c()
      for(forecast_time in 1:20){
        simplex_output = block_lnlp(tonnetzdistanceDF, lib = c(1, endPoint), pred = c(1, endPoint),
                                    target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        targRhos = c(targRhos, simplex_output$stats$rho$rho)
      }
    }
    
    trialTonnetzDistanceRho = trialTonnetzDistanceRho + mean(targRhos)

    # Calculate rho for random pairs of spectral flatness time series

    spectralFlatnessDF = cbind(spectralflatness_timeseries[colname1], spectralflatness_timeseries[colname2], spectralflatness_timeseries[colname3])

    trialSpectralFlatnessRho = 0

    for(colname in colnames){
      targRhos = c()
      for(forecast_time in 1:20){
        simplex_output = block_lnlp(spectralFlatnessDF, lib = c(1, endPoint), pred = c(1, endPoint),
                                     target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        targRhos = c(targRhos, simplex_output$stats$rho$rho)
      }
    }

    trialSpectralFlatnessRho = trialSpectralFlatnessRho + mean(targRhos)
    
    rms_rhovalues_random = c(rms_rhovalues_random, trialRMSRho / 3)
    tonnetzdistance_rhovalues_random = c(tonnetzdistance_rhovalues_random, trialTonnetzDistanceRho / 3)
    spectralflatness_rhovalues_random = c(spectralflatness_rhovalues_random, trialSpectralFlatnessRho / 3)
    
  }
  # logging
  print(sprintf("EDM RANDOM g%s", group))
}

rms_rhovalues_raw = unname(unlist(rms_rhovalues))
tonnetzdistance_rhovalues_raw = unname(unlist(tonnetzdistance_rhovalues))
spectralflatness_rhovalues_raw = unname(unlist(spectralflatness_rhovalues))

#_____________________________________________________________________________________

# RQ1 ------------

#_____________________________________________________________________________________

# RQ1a

# Test for normality
shapiro.test(rms_TEvalues_raw)

mean(rms_TEvalues_raw)
sd(rms_TEvalues_raw)
mean(rms_TEvalues_random)
sd(rms_TEvalues_random)

wilcox.test(rms_TEvalues_raw, rms_TEvalues_random)

RQ1a_df = data.frame(rms_TEvalues_raw, rms_TEvalues_random)
names(RQ1a_df) = c("RMS TE real", "RMS TE random")
RQ1a_df = gather(RQ1a_df)
ggplot(RQ1a_df, aes(x=value, y=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("Transfer Entropy") +
  ylab("real-random") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black")

# RQ1b
shapiro.test(tonnetzdistance_TEvalues_raw)

mean(tonnetzdistance_TEvalues_raw)
sd(tonnetzdistance_TEvalues_raw)
mean(tonnetzdistance_TEvalues_random)
sd(tonnetzdistance_TEvalues_random)

wilcox.test(tonnetzdistance_TEvalues_raw, tonnetzdistance_TEvalues_random)

RQ1b_df = data.frame(tonnetzdistance_TEvalues_raw, tonnetzdistance_TEvalues_random)
names(RQ1b_df) = c("Tonnetz distance TE real", "Tonnetz distance TE random")
RQ1b_df = gather(RQ1b_df)
ggplot(RQ1b_df, aes(x=value, y=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("Transfer Entropy") +
  ylab("real-random") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black")

# RQ1c
shapiro.test(spectralflatness_TEvalues_raw)

mean(spectralflatness_TEvalues_raw)
sd(spectralflatness_TEvalues_raw)
mean(spectralflatness_TEvalues_random)
sd(spectralflatness_TEvalues_random)

wilcox.test(spectralflatness_TEvalues_raw, spectralflatness_TEvalues_random)

RQ1c_df = data.frame(spectralflatness_TEvalues_raw, spectralflatness_TEvalues_random)
names(RQ1c_df) = c("Spectral flatness TE real", "Spectral flatness TE random")
RQ1c_df = gather(RQ1c_df)
ggplot(RQ1c_df, aes(x=value, y=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("Transfer Entropy") +
  ylab("real-random") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black")

# RQ1d
shapiro.test(rms_rhovalues_raw)

mean(rms_rhovalues_raw)
sd(rms_rhovalues_raw)
mean(rms_rhovalues_random)
sd(rms_rhovalues_random)

wilcox.test(rms_rhovalues_raw, rms_rhovalues_random)

RQ1d_df = data.frame(rms_rhovalues_raw, rms_rhovalues_random)
names(RQ1d_df) = c("RMS rho real", "RMS rho random")
RQ1d_df = gather(RQ1d_df)
ggplot(RQ1d_df, aes(x=value, y=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("average rho") +
  ylab("real-random") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black")

# RQ1e
shapiro.test(tonnetzdistance_rhovalues_raw)

mean(tonnetzdistance_rhovalues_raw)
sd(tonnetzdistance_rhovalues_raw)
mean(tonnetzdistance_rhovalues_random)
sd(tonnetzdistance_rhovalues_random)

wilcox.test(tonnetzdistance_rhovalues_raw, tonnetzdistance_rhovalues_random)

RQ1e_df = data.frame(tonnetzdistance_rhovalues_raw, tonnetzdistance_rhovalues_random)
names(RQ1e_df) = c("Tonnetz distance rho real", "Tonnetz distance rho random")
RQ1e_df = gather(RQ1e_df)
ggplot(RQ1e_df, aes(x=value, y=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("average rho") +
  ylab("real-random") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black")

# RQ1f
shapiro.test(spectralflatness_rhovalues_raw)

mean(spectralflatness_rhovalues_raw)
sd(spectralflatness_rhovalues_raw)
mean(spectralflatness_rhovalues_random)
sd(spectralflatness_rhovalues_random)

wilcox.test(spectralflatness_rhovalues_raw, spectralflatness_rhovalues_random)

RQ1f_df = data.frame(spectralflatness_rhovalues_raw, spectralflatness_rhovalues_random)
names(RQ1f_df) = c("Spectral flatness rho real", "Spectral flatness rho random")
RQ1f_df = gather(RQ1f_df)
ggplot(RQ1f_df, aes(x=value, y=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("average rho") +
  ylab("real-random") + 
  stat_summary(fun = "mean",
               geom = "point",
               colour = "black")

#_____________________________________________________________________________________

# RQ2 ------------

#_____________________________________________________________________________________

for (row in 1:nrow(trialDataset)){
  group = trialDataset[row, "trio"]
  trial = trialDataset[row, "take"]
  trialDataset$rho_rms_full[row] = rms_rhovalues[[sprintf("g%s_t%s", group, trial)]]
  trialDataset$TE_rms_full[row] = rms_TEvalues_pertrial[sprintf("g%s_t%s", group, trial)]

  unidirectionality_full = 0
  
  pairs = combn(c(1:3), m=2)
  for (pair in 1:3){
    colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
    colname2 = sprintf("g%s_t%s_b%s", group, trial, pairs[2,pair])
    
    pairTExy = rms_TEvalues_perpair[[sprintf("%s-%s", colname1, colname2)]]
    pairTEyx = rms_TEvalues_perpair[[sprintf("%s-%s", colname2, colname1)]]
    
    unidirectionality_pair = max(c(pairTExy / pairTEyx, pairTEyx / pairTExy))
    unidirectionality_full = unidirectionality_full + unidirectionality_pair
  }
  
  trialDataset$unidirectionality_index_full[row] = unidirectionality_full / 3
  
}

trialDataset$TE_rms_full = unlist(trialDataset$TE_rms_full)
trialDataset$rho_rms_full = unlist(trialDataset$rho_rms_full)
trialDataset$unidirectionality_index_full = unlist(trialDataset$unidirectionality_index_full)

trialDataset$rho_full = unlist(trialDataset$rho_rms_full + trialDataset$rho_tonnetzdistance_full + trialDataset$rho_spectralflatness_full)

# RQ2a

RQ2a = lmer(goodImproAll ~ TE_rms_full + (1|trio),
            trialDataset[trialDataset$take > 4,])
summary(RQ2a)

RQ2a_df = data.frame(trialDataset[trialDataset$take > 4,]$goodImproAll, trialDataset[trialDataset$take > 4,]$TE_rms_full)
RQ2a_df = RQ2a_df[complete.cases(RQ2a_df),]
names(RQ2a_df) = c("Musician_Appreciation", "Transfer_Entropy")

RQ2a_df$predlmer = predict(RQ2a)

ggplot(RQ2a_df, aes(x=Transfer_Entropy, y=Musician_Appreciation)) + 
  geom_point(na.rm=TRUE) +
  geom_smooth(aes(y = predlmer), size = 1) +
  scale_y_continuous(breaks=seq(1,7))

check_model(RQ2a)

# RQ2b
RQ2b = lmer(goodImproAll ~ unidirectionality_index_full + (1|trio),
            trialDataset[trialDataset$take > 4,])
summary(RQ2b)

RQ2b_df = data.frame(trialDataset[trialDataset$take > 4,]$goodImproAll, trialDataset[trialDataset$take > 4,]$unidirectionality_index_full)
RQ2b_df = RQ2b_df[complete.cases(RQ2b_df),]
names(RQ2b_df) = c("Musician_Appreciation", "Unidirectionality_Index")
RQ2b_df$predlmer = predict(RQ2b)

ggplot(RQ2b_df, aes(x=Unidirectionality_Index, y=Musician_Appreciation)) + 
  geom_point(na.rm=TRUE) +
  scale_y_continuous(breaks=seq(1,7)) +
  geom_smooth(aes(y = predlmer), size = 1) 

# RQ2c

RQ2c = lmer(goodImproAll ~ rho_rms_full + (1|trio),
            trialDataset[trialDataset$take > 4,])
summary(RQ2c)

RQ2c_df = data.frame(trialDataset[trialDataset$take > 4,]$goodImproAll, trialDataset[trialDataset$take > 4,]$rho_rms_full)
RQ2c_df = RQ2c_df[complete.cases(RQ2c_df),]
names(RQ2c_df) = c("Musician_Appreciation", "Rho")
RQ2c_df$predlmer = predict(RQ2c)

ggplot(RQ2c_df, aes(x=Rho, y=Musician_Appreciation)) + 
  geom_point(na.rm=TRUE) +
  geom_smooth(aes(y = predlmer), size = 1) +
  scale_y_continuous(breaks=seq(1,7))

check_model(RQ2c)

#_____________________________________________________________________________________

# Compute post-prompt TE, post-prompt rho, post-prompt unidirectionality index

#_____________________________________________________________________________________

windowSizeInSeconds = (1/22050) * 3520

View(trialDataset)

trialDataset$TE_after_prompt = c(rep(NA, nrow(trialDataset)))
trialDataset$unidirectionality_index_after_prompt = c(rep(NA, nrow(trialDataset)))

TE_perpair_after_prompt = list()

for(group in 1:12) {
  for(trial in 1:16) {
    idx = which(trialDataset$trio == group & trialDataset$take == trial)
    promptTime = trialDataset$promptTime[idx]
    promptWindow = ceiling(promptTime / windowSizeInSeconds)
    promptType = trialDataset$type[idx]
    promptNumber = trialDataset$number[idx]
    endPoint = endPoints[[sprintf("g%s_t%s", group, trial)]]
    
    if (!(promptWindow %in% c(0, NA)) & endPoint - promptWindow >= round(12/windowSizeInSeconds)) { # require at least 19 seconds of playing after prompt
      pairs = combn(c(1:3), m=2)
      unidirectionality = 0
      TE_after_prompt = 0
      for(pair in 1:3){
        
        colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
        colname2 = sprintf("g%s_t%s_b%s", group, trial, pairs[2,pair])
        
        rmsTE = transfer_entropy(rms_timeseries[,colname1][promptWindow:endPoint],
                                 rms_timeseries[,colname2][promptWindow:endPoint], shuffles=100, nboot=1, quiet=TRUE, seed=TRUE)
        rmsTExy = rmsTE$coef[1]
        rmsTEyx = rmsTE$coef[2]
        unidirectionality_pair = max(c(rmsTExy / rmsTEyx, rmsTEyx / rmsTExy))
        unidirectionality = unidirectionality + unidirectionality_pair
        TE_after_prompt = TE_after_prompt + rmsTExy + rmsTEyx
        TE_perpair_after_prompt[sprintf("%s-%s", colname1, colname2)] = rmsTExy
        TE_perpair_after_prompt[sprintf("%s-%s", colname2, colname1)] = rmsTEyx
      }
      TE_after_prompt = TE_after_prompt / 6
      unidirectionality = unidirectionality / 3
      trialDataset$TE_after_prompt[idx] = TE_after_prompt
      trialDataset$unidirectionality_index_after_prompt[idx] = unidirectionality
    }
  }
  print(sprintf("TE ENDINGS g%s", group))
}

trialDataset$rho_after_prompt = c(rep(NA, nrow(trialDataset)))
trialDataset$rho_change = c(rep(NA, nrow(trialDataset)))

for (group in 1:12){
  for (trial in 1:16){
    
    idx = which(trialDataset$trio == group & trialDataset$take == trial)
    promptTime = trialDataset$promptTime[idx]
    promptWindow = ceiling(promptTime / windowSizeInSeconds)
    promptType = trialDataset$type[idx]
    promptNumber = trialDataset$number[idx]
    endPoint = endPoints[[sprintf("g%s_t%s", group, trial)]]
    
    colname1 = sprintf("g%s_t%s_b1", group, trial)
    colname2 = sprintf("g%s_t%s_b2", group, trial)
    colname3 = sprintf("g%s_t%s_b3", group, trial)
    colnames = c(colname1, colname2, colname3)
    
    rmsDF = cbind(rms_timeseries[colname1], rms_timeseries[colname2], rms_timeseries[colname3])
    
    trialRhoAfterPrompt = 0
    trialRhoBeforePrompt = 0
    
    if (!(promptWindow %in% c(0, NA)) & endPoint - promptWindow >= round(10/windowSizeInSeconds)) {
      for(colname in colnames){
        targRhosAfterPrompt = c()
        targRhosBeforePrompt = c()
        for(forecast_time in 1:20){
          simplex_before_prompt = block_lnlp(rmsDF, lib = sprintf("1 %s", promptWindow), pred = c(1, promptWindow),
                                            target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
          simplex_after_prompt = block_lnlp(rmsDF, lib = sprintf("1 %s", promptWindow), pred = c(promptWindow+1, endPoint),
                                       target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
          targRhosAfterPrompt = c(targRhos, simplex_after_prompt$stats$rho$rho)
          targRhosBeforePrompt = c(targRhosBeforePrompt, simplex_before_prompt$stats$rho$rho)
        }
        trialRhoAfterPrompt = trialRhoAfterPrompt + mean(targRhosAfterPrompt)
        trialRhoBeforePrompt = trialRhoBeforePrompt + mean(targRhosBeforePrompt)
      }
      trialDataset$rho_after_prompt[idx] = trialRhoAfterPrompt / 3
      trialDataset$rho_change[idx] = (trialRhoAfterPrompt / 3) / (trialRhoBeforePrompt / 3)
    }
  }
  print(sprintf("EDM ENDINGS g%s", group))
}


trialDataset$TE_after_prompt = unlist(trialDataset$TE_after_prompt)
trialDataset$rho_after_prompt = unlist(trialDataset$rho_after_prompt)

#_____________________________________________________________________________________

# RQ3 ------------

#_____________________________________________________________________________________

View(listenerDataset)

trialDataset$type = factor(trialDataset$type)
trialDataset$trio = factor(trialDataset$trio)

# RQ3a
RQ3a.full = lmer(log(TE_after_prompt) ~ number + type + number:type + (1|trio),
              trialDataset[!is.na(trialDataset$TE_after_prompt),])
step(RQ3a.full, direction = "both")
RQ3a.final = lmer(log(TE_after_prompt) ~ number + type + (1|trio), 
                  trialDataset[!is.na(trialDataset$TE_after_prompt),])
summary(RQ3a.final)

# RQ3b
RQ3b.full = lmer(sqrt(rho_after_prompt) ~ number + type + number:type + (1|trio),
                 trialDataset[!is.na(trialDataset$rho_after_prompt),])
step(RQ3b.full, direction = "both")
RQ3b.final = lmer(sqrt(rho_after_prompt) ~ number + (1|trio), 
                  trialDataset[!is.na(trialDataset$rho_after_prompt),])
summary(RQ3b.final)

mean(trialDataset$TE_after_prompt[!is.na(trialDataset$TE_after_prompt) & trialDataset$number == 1])
mean(trialDataset$TE_after_prompt[!is.na(trialDataset$TE_after_prompt) & trialDataset$number == 2])
mean(trialDataset$TE_after_prompt[!is.na(trialDataset$TE_after_prompt) & trialDataset$number == 3])

RQ3b = lmer(rho_after_prompt ~ number + type + number:type + (1|trio),
            trialDataset[!is.na(trialDataset$rho_after_prompt),])

summary(RQ3b)
MuMIn::r.squaredGLMM(RQ3b)
confint(RQ3b)

mean(trialDataset$rho_after_prompt[!is.na(trialDataset$rho_after_prompt) & trialDataset$number == 1])
mean(trialDataset$rho_after_prompt[!is.na(trialDataset$rho_after_prompt) & trialDataset$number == 2])
mean(trialDataset$rho_after_prompt[!is.na(trialDataset$rho_after_prompt) & trialDataset$number == 3])

# Individual Predictability Data ------------------------------------------

#_____________________________________________________________________________________

# RQ4 ------------

#_____________________________________________________________________________________

# RQ4a


# RQ4b

# RQ4c

from = c()
to = c()
from_prompt = c()
to_prompt = c()
TE_pair_full = c()
TE_pair_afterprompt = c()
directionality_ratio_afterprompt = c()

musicianDataset$booth = factor(musicianDataset$booth)
levels(musicianDataset$booth) = c(1,2,3)

for(group in 1:12) {
  for(trial in 1:16) {
    pairs = combn(c(1:3), m=2)
    for(pair in 1:3){
      
      prompt_heard_1 = musicianDataset[musicianDataset$trio == group & musicianDataset$take == trial & musicianDataset$booth == pairs[1,pair],]$prompt_heard
      prompt_heard_2 = musicianDataset[musicianDataset$trio == group & musicianDataset$take == trial & musicianDataset$booth == pairs[2,pair],]$prompt_heard
      
      if (length(prompt_heard_1) != 0 & length(prompt_heard_2) != 0){
        colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
        colname2 = sprintf("g%s_t%s_b%s", group, trial, pairs[2,pair])
        from = c(from, colname1, colname2)
        to = c(to, colname2, colname1)
        from_prompt = c(from_prompt, prompt_heard_1, prompt_heard_2)
        to_prompt = c(to_prompt, prompt_heard_2, prompt_heard_1)
        TE_pair_full = c(TE_pair_full,
                         rms_TEvalues_perpair[sprintf("%s-%s", colname1, colname2)], 
                         rms_TEvalues_perpair[sprintf("%s-%s", colname2, colname1)])
        TE_pair_afterprompt = c(TE_pair_afterprompt, 
                                 TE_perpair_after_prompt[sprintf("%s-%s", colname1, colname2)],
                                 TE_perpair_after_prompt[sprintf("%s-%s", colname2, colname1)])
      }
    }
  }
}

RQ4_pair_df = data.frame(from = from)

RQ4_pair_df$to = to
RQ4_pair_df$TE_pair_full = TE_pair_full
RQ4_pair_df$TE_pair_afterprompt = TE_pair_afterprompt # a[sapply(a, is.null)] <- NA
RQ4_pair_df$from_prompt = from_prompt
RQ4_pair_df$to_prompt = to_prompt

RQ4_pair_df$directionality_ratio_after_prompt = rep(NA, nrow(RQ4_pair_df))

for(row in 1:nrow(RQ4_pair_df)){
  if(row%%2==1){
    if(length(unname(unlist(RQ4_pair_df$TE_pair_afterprompt[row])))!=0){
      RQ4_pair_df$directionality_ratio_after_prompt[row] = unname(unlist(RQ4_pair_df$TE_pair_afterprompt[row])) / unname(unlist(RQ4_pair_df$TE_pair_afterprompt[row+1]))
  }
  }
}

for(col in names(RQ4_pair_df)){
  RQ4_pair_df[col] = unlist(RQ4_pair_df[col])
}

#RQ4_pair_df$individual_predictability_change = rep(NA, nrow(RQ4_pair_df))
#RQ4_pair_df$individual_predictability_full = rep(NA, nrow(RQ4_pair_df))

for (group in 2:3){
  for (trial in 1:16){
    idx = which(trialDataset$trio == group & trialDataset$take == trial)
    promptTime = trialDataset$promptTime[idx]
    promptWindow = ceiling(promptTime / windowSizeInSeconds)
    endPoint = endPoints[[sprintf("g%s_t%s", group, trial)]]
    
    colname1 = sprintf("g%s_t%s_b1", group, trial)
    colname2 = sprintf("g%s_t%s_b2", group, trial)
    colname3 = sprintf("g%s_t%s_b3", group, trial)
    colnames = c(colname1, colname2, colname3)
    
    rmsDF = cbind(rms_timeseries[colname1], rms_timeseries[colname2], rms_timeseries[colname3])
    
    for(colname in colnames){
      total_rho = c()
      for(forecast_time in 1:20){
        simplex_full = block_lnlp(unlist(rms_timeseries[colname]), tp=forecast_time)
        total_rho = c(total_rho, simplex_full$const_pred_rho$rho)
      }
      RQ4_pair_df[match(colname, RQ4_pair_df$from),]$individual_predictability_full = mean(total_rho)
    }
    
    if (!(promptWindow %in% c(0, NA)) & endPoint - promptWindow >= round(10/windowSizeInSeconds)) {
      for(colname in colnames){
        total_rho_after_prompt = c()
        total_rho_before_prompt = c()
        for(forecast_time in 1:20){
          simplex_after_prompt = block_lnlp(unlist(rms_timeseries[colname]), lib = sprintf("1 %s", promptWindow), pred = c(promptWindow+1, endPoint), tp=forecast_time)
          simplex_before_prompt = block_lnlp(unlist(rms_timeseries[colname]), lib = sprintf("1 %s", promptWindow), pred = c(1, promptWindow), tp=forecast_time, exclusion_radius = forecast_time)
          total_rho_after_prompt = c(total_rho_after_prompt, simplex_after_prompt$const_pred_rho$rho)
          total_rho_before_prompt = c(total_rho_before_prompt, simplex_before_prompt$const_pred_rho$rho)
        }
        individual_predictability_change = mean(total_rho_after_prompt) / mean(total_rho_before_prompt)
        RQ4_pair_df[match(colname, RQ4_pair_df$from),]$individual_predictability_change = individual_predictability_change
      }
    }
  }
  print(sprintf("Rho calculated for %s", group))
}

RQ4_pair_df$from_prompt = factor(RQ4_pair_df$from_prompt)
RQ4_pair_df$to_prompt = factor(RQ4_pair_df$to_prompt)
RQ4_pair_df$trio = rep(1:12, each=nrow(RQ4_pair_df)/12)

RQ4_pair_df

# RQ4a

# RQ4b

# RQ4c
RQ4c = lm(directionality_ratio_after_prompt ~ relevel(from_prompt, ref="No-Goal") + relevel(to_prompt, ref="No-Goal"), data=RQ4_pair_df)
summary(RQ4c)

# IMPORTANT
# RQ4_pair_df[match(unique(RQ4_pair_df$from), RQ4_pair_df$from),]




RQ4_pair_df[RQ4_pair_df$from=="g1_t3_b1",]$individual_predictability_full

