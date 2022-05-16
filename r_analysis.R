#__________________________________________________________________________


# Loading packages, data, etc. --------------------------------------------

#__________________________________________________________________________

# Turn off scientific notation
options(scipen=999)

dev.off()

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
install.packages("tidyr")
install.packages("sjPlot")
install.packages("effects")

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
library(tidyr)
library(sjPlot)
library(effects)

# Importing libraries
# Importing our sliding window-averaged data
rms_timeseries = read.csv("C:/Shortcutsensei/0. Thesis/rms_timeseries.csv")
tonalcentroid_timeseries_raw = read.csv("C:/Shortcutsensei/0. Thesis/tonal_centroid_timeseries.csv")
spectralflatness_timeseries = read.csv("C:/Shortcutsensei/0. Thesis/spectral_flatness_timeseries.csv")

# Importing datasets with info about the trials
musicianDataset = read.csv("C:/Users/lenna/Downloads/study_3_4/study_3_4/endImpro_4R_booths_3_trio_12.csv", sep=";")
trialDataset = read.csv("C:/Shortcutsensei/JointImprovisation/Data/endingDataset_fromThomas_withconditions.csv", sep=";") 
listenerDataset = read.csv("C:/Shortcutsensei/JointImprovisation/Data/Appreciation_Data_N46.csv", sep=";")

windowSizeInSeconds = (1/22050) * 3520

#_____________________________________________________________________________________

# Further preprocessing of Tonnetz data ----------------------------------------------

#_____________________________________________________________________________________

# We compute a 'tonnetz distance' measure. At each window, it computes the Euclidean distance between
# the current and the previous tonnetz (both of which are 6-dimensional vectors).
tonnetzdistance_timeseries = data.frame("g1_t1_b1" = rep(0, 3000))

for(group in 1:12) {
  for(trial in 1:16) {
    for(booth in 1:3) {
      
      # Create string featuring the identifier of the current time series
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
          # for all subsequent windows, compute the Euclidean distance between the current and previous tonnetz. 
          # This represents a measure of harmonic change from one window to the next.
          prevTonnetz = unname(unlist(curBooth[row-1,]))
          distance = dist(rbind(tonnetz, prevTonnetz))[1]
        }
        # Append to our vector of tonnetz distances
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

# Plotting example of acoustic feature time series ------------

#_____________________________________________________________________________________

# For each acoustic feature, make a line plot of the 3 musicians in trio 1, trial 1.
# See figure 1 in text.

# RMS amplitude
rms_timeseries_example_b1 = rms_timeseries$g1_t1_b1[!is.na(rms_timeseries$g1_t1_b1)]
rms_timeseries_example_b2 = rms_timeseries$g1_t1_b2[!is.na(rms_timeseries$g1_t1_b1)]
rms_timeseries_example_b3 = rms_timeseries$g1_t1_b3[!is.na(rms_timeseries$g1_t1_b1)]
rms_timeseries_example = data.frame(seconds=seq(length(rms_timeseries_example_b1)) * windowSizeInSeconds,
                                                b1=rms_timeseries_example_b1,
                                                b2=rms_timeseries_example_b2,
                                                b3=rms_timeseries_example_b3)[0:60/windowSizeInSeconds,]

ggplot(data=rms_timeseries_example, aes(x=seconds)) +
  geom_line(aes(y=b1), color="steelblue") +
  geom_line(aes(y=b2), color="red") +
  geom_line(aes(y=b3), color="black") +
  ylab("RMS amplitude") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 22))

# Tonnetz distance
tonnetzdistance_timeseries_example_b1 = tonnetzdistance_timeseries$g1_t1_b1[!is.na(tonnetzdistance_timeseries$g1_t1_b1)]
tonnetzdistance_timeseries_example_b2 = tonnetzdistance_timeseries$g1_t1_b2[!is.na(tonnetzdistance_timeseries$g1_t1_b1)]
tonnetzdistance_timeseries_example_b3 = tonnetzdistance_timeseries$g1_t1_b3[!is.na(tonnetzdistance_timeseries$g1_t1_b1)]
tonnetzdistance_timeseries_example = data.frame(seconds=seq(length(tonnetzdistance_timeseries_example_b1)) * windowSizeInSeconds,
                                    b1=tonnetzdistance_timeseries_example_b1,
                                    b2=tonnetzdistance_timeseries_example_b2,
                                    b3=tonnetzdistance_timeseries_example_b3)[0:60/windowSizeInSeconds,]

ggplot(data=tonnetzdistance_timeseries_example, aes(x=seconds)) +
  geom_line(aes(y=b1), color="steelblue") +
  geom_line(aes(y=b2), color="red") +
  geom_line(aes(y=b3), color="black") +
  ylab("Tonnetz distance") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 22))

# Spectral flatness
spectralflatness_timeseries_example_b1 = spectralflatness_timeseries$g1_t1_b1[!is.na(spectralflatness_timeseries$g1_t1_b1)]
spectralflatness_timeseries_example_b2 = spectralflatness_timeseries$g1_t1_b2[!is.na(spectralflatness_timeseries$g1_t1_b1)]
spectralflatness_timeseries_example_b3 = spectralflatness_timeseries$g1_t1_b3[!is.na(spectralflatness_timeseries$g1_t1_b1)]
spectralflatness_timeseries_example = data.frame(seconds=seq(length(spectralflatness_timeseries_example_b1)) * windowSizeInSeconds,
                                                b1=spectralflatness_timeseries_example_b1,
                                                b2=spectralflatness_timeseries_example_b2,
                                                b3=spectralflatness_timeseries_example_b3)[0:60/windowSizeInSeconds,]

ggplot(data=spectralflatness_timeseries_example, aes(x=seconds)) +
  geom_line(aes(y=b1), color="steelblue") +
  geom_line(aes(y=b2), color="red") +
  geom_line(aes(y=b3), color="black") +
  ylab("Spectral flatness") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 22))

#_____________________________________________________________________________________

# Find the characteristic time step (Markov order) for ETE calculations ------------

#_____________________________________________________________________________________

# Enable parallel processing for all calc_ete() calls
plan(multisession)

# char_time = rep(0, 20)

char_time_rho = rep(0, 20)

for(group in 1:12){
  for(pair in 1:3){
    # colname1 = sprintf("g%s_t%s_b%s", group, 6, pairs[1,pair])
    # colname2 = sprintf("g%s_t%s_b%s", group, 6, pairs[2,pair])

    # # Calculate transfer entropy for all pairs of RMS amplitude time series within the trial

    # Run calc_ete() with Markov orders from 1 up to 20, to mitigate the issue of
    # having a fixed time delay. calc_ete() calculates effective transfer entropy.
    
    # for(markov_order in 1:20){
    #   TExy = calc_ete(rms_timeseries[,colname1],
    #                   rms_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
    #   TEyx = calc_ete(rms_timeseries[,colname2],
    #                   rms_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
    #   char_time[markov_order] = char_time[markov_order] + TExy + TEyx
    #   
    #   TExy = calc_ete(tonnetzdistance_timeseries[,colname1],
    #                   tonnetzdistance_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
    #   TEyx = calc_ete(tonnetzdistance_timeseries[,colname2],
    #                   tonnetzdistance_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
    #   char_time[markov_order] = char_time[markov_order] + TExy + TEyx
    #   
    #   TExy = calc_ete(spectralflatness_timeseries[,colname1],
    #                   spectralflatness_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
    #   TEyx = calc_ete(spectralflatness_timeseries[,colname2],
    #                   spectralflatness_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
    #   char_time[markov_order] = char_time[markov_order] + TExy + TEyx
    # }
    
    print("pass")
  }
  
  colnames = c(sprintf("g%s_t6_b1", group), 
               sprintf("g%s_t6_b2", group), 
               sprintf("g%s_t6_b3", group))
  
  # Get the ending of the performance (i.e. window at which the first musician stops playing)
  endPoint = min(c(which(is.na(rms_timeseries[colnames[1]]))[1], 
                   which(is.na(rms_timeseries[colnames[2]]))[1], 
                   which(is.na(rms_timeseries[colnames[3]]))[1]))
  
  rmsDF = cbind(rms_timeseries[colnames[1]], rms_timeseries[colnames[2]], rms_timeseries[colnames[3]])
  tonnetzdistanceDF = cbind(tonnetzdistance_timeseries[colnames[1]], tonnetzdistance_timeseries[colnames[2]], tonnetzdistance_timeseries[colnames[3]])
  spectralFlatnessDF = cbind(spectralflatness_timeseries[colnames[1]], spectralflatness_timeseries[colnames[2]], spectralflatness_timeseries[colnames[3]])
  
  # For each trial, calculate three values of rho, each with a different booth as the target
  for(colname in colnames){
    # Calculate rho values from 1 lag up to 20 lags
    for(forecast_time in 1:20){
      # First, we use the first half of the performance to construct the attractor manifold
      # and use the second half to predict
      simplex_output1 = block_lnlp(rmsDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                   target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
      # Then, we turn it around; we use the second half to construct the manifold and the first half
      # for prediction purposes
      simplex_output2 = block_lnlp(rmsDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                   target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
      char_time_rho[forecast_time] = char_time_rho[forecast_time] + mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho)
      
      simplex_output1 = block_lnlp(tonnetzdistanceDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                   target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
      simplex_output2 = block_lnlp(tonnetzdistanceDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                   target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
      char_time_rho[forecast_time] = char_time_rho[forecast_time] + mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho)
      
      simplex_output1 = block_lnlp(spectralFlatnessDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                   target_column = colname, tp = forecast_time,  stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
      simplex_output2 = block_lnlp(spectralFlatnessDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                   target_column = colname, tp = forecast_time,  stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
      char_time_rho[forecast_time] = char_time_rho[forecast_time] + mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho)
    }
  }
  print("EDM")
  
  print(sprintf("GROUP: %s", group))
} 

# char_time = char_time / (3 * 12 * 6)
char_time_rho = char_time_rho / (3 * 12)

char_time_rho

plot(char_time_rho, type='b', pch=19, ylab=
       "Rho", xlab="Forecast time", main =
       "Average Rho for forecast times 1-20")

#_____________________________________________________________________________________

# TE for real pairs / trials and random pairs ------------

#_____________________________________________________________________________________

# REAL PAIRS + TRIALS
# Calculate TE values for all 'real pairs' of musicians; i.e., same group, same trial. Also calculate
# TE for trials by computing mean of all pairwise TEs a trial. This will come in handy in RQ2.

# Initialize empty vectors

#rms_TEvalues_perpair = list()
#spectralflatness_TEvalues_perpair = list()
#tonnetzdistance_TEvalues_perpair = list()

#rms_TEvalues_pertrial = list()
#spectralflatness_TEvalues_pertrial = list()
#tonnetzdistance_TEvalues_pertrial = list()

format(Sys.time(), "%X")

for(group in c(6)){
  for(trial in 13:16){
      # Get all combinations of booths: booth 1-2, booth 1-3, booth 2-3
      pairs = combn(c(1:3), m=2)
      
      trial_rms_TE = c()
      trial_spectralflatness_TE = c()
      trial_tonnetzdistance_TE = c()
      
      for(pair in 1:3){
        
        colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
        colname2 = sprintf("g%s_t%s_b%s", group, trial, pairs[2,pair])
        
        # Calculate transfer entropy for all pairs of RMS amplitude time series within the trial
        
        rmsTExy = 0
        rmsTEyx = 0
        
        # Run calc_ete() with Markov orders from 1 up to 20, to mitigate the issue of
        # having a fixed time delay. calc_ete() calculates effective transfer entropy.
        for(markov_order in 1:20){
          TExy = calc_ete(rms_timeseries[,colname1],
                          rms_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
          TEyx = calc_ete(rms_timeseries[,colname2],
                          rms_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
          rmsTExy = rmsTExy + TExy
          rmsTEyx = rmsTEyx + TEyx
        }
        
        # We have calculated 20 TE values both ways for this pair. Divide by 20 to average over
        # the different markov order settings.
        rmsTExy = rmsTExy / 20
        rmsTEyx = rmsTEyx / 20
        
        # Store results in list
        rms_TEvalues_perpair[sprintf("%s-%s", colname1, colname2)] = rmsTExy
        rms_TEvalues_perpair[sprintf("%s-%s", colname2, colname1)] = rmsTEyx
        
        # We will later take the mean of this vector to get one final TE value for the whole trial
        trial_rms_TE = c(trial_rms_TE, rmsTExy, rmsTEyx)

        # Calculate transfer entropy for all possible pairs of tonnetz distance time series within the trial
        
        tonnetzdistanceTExy = 0
        tonnetzdistanceTEyx = 0
        
        for(markov_order in 1:20){
          TExy = calc_ete(tonnetzdistance_timeseries[,colname1],
                                      tonnetzdistance_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
          TEyx = calc_ete(tonnetzdistance_timeseries[,colname2],
                                        tonnetzdistance_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
          tonnetzdistanceTExy = tonnetzdistanceTExy + TExy
          tonnetzdistanceTEyx = tonnetzdistanceTEyx + TEyx
        }
        
        tonnetzdistanceTExy = tonnetzdistanceTExy / 20
        tonnetzdistanceTEyx = tonnetzdistanceTEyx / 20
        
        tonnetzdistance_TEvalues_perpair[sprintf("%s-%s", colname1, colname2)] = tonnetzdistanceTExy
        tonnetzdistance_TEvalues_perpair[sprintf("%s-%s", colname2, colname1)] = tonnetzdistanceTEyx
        
        trial_tonnetzdistance_TE = c(trial_tonnetzdistance_TE, tonnetzdistanceTExy, tonnetzdistanceTEyx)
        
        # Calculate transfer entropy for all possible pairs of spectral flatness time series within the trial
        
        spectralflatnessTExy = 0
        spectralflatnessTEyx = 0
        
        for(markov_order in 1:20){
          TExy = calc_ete(spectralflatness_timeseries[,colname1],
                         spectralflatness_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
          TEyx = calc_ete(spectralflatness_timeseries[,colname2],
                         spectralflatness_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
          spectralflatnessTExy = spectralflatnessTExy + TExy
          spectralflatnessTEyx = spectralflatnessTEyx + TEyx
        }
        
        spectralflatnessTExy = spectralflatnessTExy / 20
        spectralflatnessTEyx = spectralflatnessTEyx / 20
        
        spectralflatness_TEvalues_perpair[sprintf("%s-%s", colname1, colname2)] = spectralflatnessTExy
        spectralflatness_TEvalues_perpair[sprintf("%s-%s", colname2, colname1)] = spectralflatnessTEyx
        
        trial_spectralflatness_TE = c(trial_spectralflatness_TE, spectralflatnessTExy, spectralflatnessTEyx)
      }
      
      print(format(Sys.time(), "%X"))
      
      # Get trial-wide TE values by averaging over the 6 pairwise TE values for each trial.
      # Store in list.
      rms_TEvalues_pertrial[sprintf("g%s_t%s", group, trial)] = mean(trial_rms_TE)
      tonnetzdistance_TEvalues_pertrial[sprintf("g%s_t%s", group, trial)] = mean(trial_tonnetzdistance_TE)
      spectralflatness_TEvalues_pertrial[sprintf("g%s_t%s", group, trial)] = mean(trial_spectralflatness_TE)
      
      # Logging
      print(sprintf("TE REAL g%s t%s", group, trial)) 
  }
  #beep(sound=9)
}

lapply(rms_TEvalues_perpair, write, "C:/Shortcutsensei/JointImprovisation/scripts/rms_TEvalues_perpair.txt", append=TRUE)
lapply(rms_TEvalues_pertrial, write, "C:/Shortcutsensei/JointImprovisation/scripts/rms_TEvalues_pertrial.txt", append=TRUE)
lapply(tonnetzdistance_TEvalues_perpair, write, "C:/Shortcutsensei/JointImprovisation/scripts/tonnetzdistance_TEvalues_perpair.txt", append=TRUE)
lapply(tonnetzdistance_TEvalues_pertrial, write, "C:/Shortcutsensei/JointImprovisation/scripts/tonnetzdistance_TEvalues_pertrial.txt", append=TRUE)
lapply(spectralflatness_TEvalues_perpair, write, "C:/Shortcutsensei/JointImprovisation/scripts/spectralflatness_TEvalues_perpair.txt", append=TRUE)
lapply(spectralflatness_TEvalues_pertrial, write, "C:/Shortcutsensei/JointImprovisation/scripts/spectralflatness_TEvalues_pertrial.txt", append=TRUE)

# RANDOM PAIRS
# Calculate TE values for 'random pairs' of musicians; i.e., same group, different trial.

# rms_TEvalues_random = c()
# spectralflatness_TEvalues_random = c()
# tonnetzdistance_TEvalues_random = c()
# 
# rms_TEvalues_random_pertrial = c()
# spectralflatness_TEvalues_random_pertrial = c()
# tonnetzdistance_TEvalues_random_pertrial = c()

for(group in 11:12){
  for(trial in 1:16){
    pairs = combn(c(1:3), m=2)
    trial_rms_TE = c()
    trial_tonnetzdistance_TE = c()
    trial_spectralflatness_TE = c()
    for(pair in 1:3){
      
      # Create random pairs by matching g1_t1_b1 to g1_t2_b2, g1_t1_b1 to g1_t2_b3, g1_t1_b2 to g1_t2_b3, and so on.
      trialRandom = if (trial != 16) trial+1 else 1
      colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
      colname2 = sprintf("g%s_t%s_b%s", group, trialRandom, pairs[2,pair])
      
      # Calculate effective transfer entropy for random pairs of RMS amplitude time series
      rmsTExy = 0
      rmsTEyx = 0
      
      for(markov_order in 1:20){
        TExy = calc_ete(rms_timeseries[,colname1],
                       rms_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
        TEyx = calc_ete(rms_timeseries[,colname2],
                       rms_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
        rmsTExy = rmsTExy + TExy
        rmsTEyx = rmsTEyx + TEyx
      }
      
      rmsTExy = rmsTExy / 20
      rmsTEyx = rmsTEyx / 20
      
      # add all random-pair TE values to a vector
      rms_TEvalues_random = c(rms_TEvalues_random, rmsTExy, rmsTEyx)
      trial_rms_TE = c(trial_rms_TE, rmsTExy, rmsTEyx)
      
      # Calculate transfer entropy for random pairs of tonnetz distance time series
      tonnetzdistanceTExy = 0
      tonnetzdistanceTEyx = 0
      
      for(markov_order in 1:20){
        TExy = calc_ete(tonnetzdistance_timeseries[,colname1],
                       tonnetzdistance_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
        TEyx = calc_ete(tonnetzdistance_timeseries[,colname2],
                       tonnetzdistance_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
        tonnetzdistanceTExy = tonnetzdistanceTExy + TExy
        tonnetzdistanceTEyx = tonnetzdistanceTEyx + TEyx
      }
      
      tonnetzdistanceTExy = tonnetzdistanceTExy / 20
      tonnetzdistanceTEyx = tonnetzdistanceTEyx / 20
      
      tonnetzdistance_TEvalues_random = c(tonnetzdistance_TEvalues_random, tonnetzdistanceTExy, tonnetzdistanceTEyx)
      trial_tonnetzdistance_TE = c(trial_tonnetzdistance_TE, tonnetzdistanceTExy, tonnetzdistanceTEyx)
      
      # Calculate transfer entropy for random pairs of spectral flatness time series

      spectralflatnessTExy = 0
      spectralflatnessTEyx = 0
      
      for(markov_order in 1:20){
        TExy = calc_ete(spectralflatness_timeseries[,colname1],
                       spectralflatness_timeseries[,colname2], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
        TEyx = calc_ete(spectralflatness_timeseries[,colname2],
                       spectralflatness_timeseries[,colname1], lx=markov_order, ly=markov_order, shuffles=20, seed=TRUE)
        spectralflatnessTExy = spectralflatnessTExy + TExy
        spectralflatnessTEyx = spectralflatnessTEyx + TEyx
      }
      
      spectralflatnessTExy = spectralflatnessTExy / 20
      spectralflatnessTEyx = spectralflatnessTEyx / 20
      
      spectralflatness_TEvalues_random = c(spectralflatness_TEvalues_random, spectralflatnessTExy, spectralflatnessTEyx)
      trial_spectralflatness_TE = c(trial_spectralflatness_TE, spectralflatnessTExy, spectralflatnessTEyx)
      
    }
    # Logging
    print(sprintf("TE RANDOM g%s t%s", group, trial))
    
    print(format(Sys.time(), "%X"))
    
    rms_TEvalues_random_pertrial = c(rms_TEvalues_random_pertrial, mean(trial_rms_TE))
    tonnetzdistance_TEvalues_random_pertrial = c(tonnetzdistance_TEvalues_random_pertrial, mean(trial_tonnetzdistance_TE))
    spectralflatness_TEvalues_random_pertrial = c(spectralflatness_TEvalues_random_pertrial, mean(trial_spectralflatness_TE))
  }
}

rms_TEvalues_random_pertrial = unique(spectralflatness_TEvalues_random_pertrial)
tonnetzdistance_TEvalues_random_pertrial = unique(tonnetzdistance_TEvalues_random_pertrial)
spectralflatness_TEvalues_random_pertrial = unique(spectralflatness_TEvalues_random_pertrial)

lapply(rms_TEvalues_random, write, "C:/Shortcutsensei/JointImprovisation/scripts/rms_TEvalues_random.txt", append=TRUE)
lapply(rms_TEvalues_random_pertrial, write, "C:/Shortcutsensei/JointImprovisation/scripts/rms_TEvalues_random_pertrial.txt", append=TRUE)
lapply(tonnetzdistance_TEvalues_random, write, "C:/Shortcutsensei/JointImprovisation/scripts/tonnetzdistance_TEvalues_random.txt", append=TRUE)
lapply(tonnetzdistance_TEvalues_random_pertrial, write, "C:/Shortcutsensei/JointImprovisation/scripts/tonnetzdistance_TEvalues_random_pertrial.txt", append=TRUE)
lapply(spectralflatness_TEvalues_random, write, "C:/Shortcutsensei/JointImprovisation/scripts/spectralflatness_TEvalues_random.txt", append=TRUE)
lapply(spectralflatness_TEvalues_random_pertrial, write, "C:/Shortcutsensei/JointImprovisation/scripts/spectralflatness_TEvalues_random_pertrial.txt", append=TRUE)

rms_TEvalues_random_pertrial

# Get just the values from each key-value pair in the list, so that we can perform our statistical tests
rms_TEvalues_raw = unname(unlist(rms_TEvalues_perpair))
tonnetzdistance_TEvalues_raw = unname(unlist(tonnetzdistance_TEvalues_perpair))
spectralflatness_TEvalues_raw = unname(unlist(spectralflatness_TEvalues_perpair))

length(rms_TEvalues_pertrial)

#_____________________________________________________________________________________

# Rho for real trios and random trios ------------

#_____________________________________________________________________________________

# REAL TRIOS

# Initialize empty lists
#rms_rhovalues = list()
#tonnetzdistance_rhovalues = list()
#spectralflatness_rhovalues = list()

# This list will contain, for each trial, the window at which at least one of the musicians has
# stopped playing.
endPoints = list()

for(group in 5:12){
  for(trial in 1:16){
    
    # Get the three recordings from a trial
    colnames = c(sprintf("g%s_t%s_b1", group, trial), 
                 sprintf("g%s_t%s_b2", group, trial), 
                 sprintf("g%s_t%s_b3", group, trial))
    
    # Get the ending of the performance (i.e. window at which the first musician stops playing)
    endPoint = min(c(which(is.na(rms_timeseries[colnames[1]]))[1], 
                      which(is.na(rms_timeseries[colnames[2]]))[1], 
                      which(is.na(rms_timeseries[colnames[3]]))[1]))
    
    endPoints[sprintf("g%s_t%s", group, trial)] = endPoint
    
    # Calculate rho for all real pairs of RMS amplitude time series
    
    rmsDF = cbind(rms_timeseries[colnames[1]], rms_timeseries[colnames[2]], rms_timeseries[colnames[3]])
    
    trialRMSRho = c()
    
    # For each trial, calculate three values of rho, each with a different booth as the target
    for(colname in colnames){
      rhos_forecast_times = c()
      # Calculate rho values from 1 lag up to 20 lags
      for(forecast_time in 1:20){
        # First, we use the first half of the performance to construct the attractor manifold
        # and use the second half to predict
        simplex_output1 = block_lnlp(rmsDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                     target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        # Then, we turn it around; we use the second half to construct the manifold and the first half
        # for prediction purposes
        simplex_output2 = block_lnlp(rmsDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                     target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        # Average the two rho values we got out of this and append to rhos_forecast_times vector
        rhos_forecast_times = c(rhos_forecast_times, mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho))
      }
      # Average rhos_forecast_times (a vector of length 20); append to trialRMSRho,
      # which we will later average to get the group-level predictability
      trialRMSRho = c(trialRMSRho, mean(rhos_forecast_times))
    }
    
    # Calculate rho for all real pairs of tonnetz distance time series
    
    tonnetzdistanceDF = cbind(tonnetzdistance_timeseries[colnames[1]], tonnetzdistance_timeseries[colnames[2]], tonnetzdistance_timeseries[colnames[3]])

    trialTonnetzDistanceRho = c()

    for(colname in colnames){
       rhos_forecast_times = c()
       for(forecast_time in 1:20){
       simplex_output1 = block_lnlp(tonnetzdistanceDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                    target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
       simplex_output2 = block_lnlp(tonnetzdistanceDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                    target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
       rhos_forecast_times = c(rhos_forecast_times, mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho))
       }
     trialTonnetzDistanceRho = c(trialTonnetzDistanceRho, mean(rhos_forecast_times))
    }
    
    # Calculate rho for all real pairs of spectral flatness time series
    
    spectralFlatnessDF = cbind(spectralflatness_timeseries[colnames[1]], spectralflatness_timeseries[colnames[2]], spectralflatness_timeseries[colnames[3]])
    
    trialSpectralFlatnessRho = c()
    
    for(colname in colnames){
      rhos_forecast_times = c()
      for(forecast_time in 1:20){
        simplex_output1 = block_lnlp(spectralFlatnessDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                   target_column = colname, tp = forecast_time,  stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        simplex_output2 = block_lnlp(spectralFlatnessDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                    target_column = colname, tp = forecast_time,  stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        rhos_forecast_times = c(rhos_forecast_times, mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho))
       }
       trialSpectralFlatnessRho = c(trialSpectralFlatnessRho, mean(rhos_forecast_times))
    }
    
    
    # The vectors trialRMSrho, trialSpectralFlatnessRho and trialTonnetzDistanceRho each contain 3 rho values;
    # one for musician (i.e. target). Average them to get group-level predictability.
    rms_rhovalues[sprintf("g%s_t%s", group, trial)] = mean(trialRMSRho)
    spectralflatness_rhovalues[sprintf("g%s_t%s", group, trial)] = mean(trialSpectralFlatnessRho)
    tonnetzdistance_rhovalues[sprintf("g%s_t%s", group, trial)] = mean(trialTonnetzDistanceRho)
    print(trial)
    print(format(Sys.time(), "%X"))
  }
  # Logging
  print(sprintf("EDM REAL g%s", group))
}

# Get raw values to make data amenable to statistical tests
rms_rhovalues_raw = unname(unlist(rms_rhovalues))
tonnetzdistance_rhovalues_raw = unname(unlist(tonnetzdistance_rhovalues))
spectralflatness_rhovalues_raw = unname(unlist(spectralflatness_rhovalues))

# RANDOM TRIOS

# Initialize empty vectors
# rms_rhovalues_random = c()
# spectralflatness_rhovalues_random = c()
# tonnetzdistance_rhovalues_random = c()

for(group in 1:12){
  for(trial in 1:16){

    # Same procedure as before for creating random pairs of musicians
    trialRandom1 = if (trial != 16) trial+1 else 1
    trialRandom2 = if (trial != 1) trial-1 else 16
    
    colnames = c(sprintf("g%s_t%s_b1", group, trial),
                 sprintf("g%s_t%s_b2", group, trialRandom1),
                 sprintf("g%s_t%s_b3", group, trialRandom2))
    
    # This gives us the point, for this random trio, at which at least 1 musician has stopped playing for good.
    endPoint = min(c(which(is.na(rms_timeseries[colnames[1]]))[1], 
                     which(is.na(rms_timeseries[colnames[2]]))[1], 
                     which(is.na(rms_timeseries[colnames[3]]))[1]))
    
    # Calculate rho for random pairs of RMS amplitude time series

    rmsDF = cbind(rms_timeseries[colnames[1]], rms_timeseries[colnames[2]], rms_timeseries[colnames[3]])
    
    trialRMSRho = c()

    for(colname in colnames){
      rhos_forecast_times = c()
      for(forecast_time in 1:20){
        simplex_output1 = block_lnlp(rmsDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                     target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        simplex_output2 = block_lnlp(rmsDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                     target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
        rhos_forecast_times = c(rhos_forecast_times, mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho))
      }
      trialRMSRho = c(trialRMSRho, mean(rhos_forecast_times))
    }
    
    # Calculate rho for random pairs of tonnetz distance time series
    
    tonnetzdistanceDF = cbind(tonnetzdistance_timeseries[colnames[1]], tonnetzdistance_timeseries[colnames[2]], tonnetzdistance_timeseries[colnames[3]])
    
    trialTonnetzDistanceRho = c()
     
    for(colname in colnames){
       rhos_forecast_times = c()
       for(forecast_time in 1:20){
       simplex_output1 = block_lnlp(tonnetzdistanceDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                    target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
       simplex_output2 = block_lnlp(tonnetzdistanceDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                    target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
       rhos_forecast_times = c(rhos_forecast_times, mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho))
       }
     trialTonnetzDistanceRho = c(trialTonnetzDistanceRho, mean(rhos_forecast_times))
    }
    

    # Calculate rho for random pairs of spectral flatness time series

    spectralFlatnessDF = cbind(spectralflatness_timeseries[colnames[1]], spectralflatness_timeseries[colnames[2]], spectralflatness_timeseries[colnames[3]])
     
    trialSpectralFlatnessRho = c()
     
    for(colname in colnames){
       rhos_forecast_times = c()
       for(forecast_time in 1:20){
         simplex_output1 = block_lnlp(spectralFlatnessDF, lib = c(1, floor(endPoint / 2)), pred = c(floor(endPoint/2)+1, endPoint),
                                    target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
         simplex_output2 = block_lnlp(spectralFlatnessDF, lib = c(floor(endPoint/2)+1, endPoint), pred = c(1, floor(endPoint / 2)),
                                      target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
         rhos_forecast_times = c(rhos_forecast_times, mean(simplex_output1$stats$rho$rho, simplex_output2$stats$rho$rho))
       }
      trialSpectralFlatnessRho = c(trialSpectralFlatnessRho, mean(rhos_forecast_times))
    }
    
    rms_rhovalues_random = c(rms_rhovalues_random, mean(trialRMSRho))
    tonnetzdistance_rhovalues_random = c(tonnetzdistance_rhovalues_random, mean(trialTonnetzDistanceRho))
    spectralflatness_rhovalues_random = c(spectralflatness_rhovalues_random, mean(trialSpectralFlatnessRho))
    
  }
  # Logging
  print(sprintf("EDM RANDOM g%s", group))
}

#_____________________________________________________________________________________

# RQ1 ------------

#_____________________________________________________________________________________

# RQ1a

# Test for normality
qqnorm(rms_TEvalues_raw)

# Descriptive statistics
mean(rms_TEvalues_raw)
sd(rms_TEvalues_raw)
mean(rms_TEvalues_random_pertrial)
sd(rms_TEvalues_random_pertrial)

# Test for significant difference between TE values for RMS amplitude time series
# of real vs random pairs
wilcox.test(rms_TEvalues_raw, rms_TEvalues_random_pertrial)

# Plot the relationship
# 773 * 657

RQ1a.1_df = data.frame(rms_TEvalues_raw, rms_TEvalues_random_pertrial)
names(RQ1a.1_df) = c("RMS ETE real", "RMS ETE random")
RQ1a.1_df = gather(RQ1a.1_df)
ggplot(RQ1a.1_df, aes(x=value, y=key, fill=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("ETE") +
  ylab("Random-real") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  theme(legend.position="none", text = element_text(size = 15))

# Tonnetz distance TE
qqnorm(tonnetzdistance_TEvalues_raw)

mean(tonnetzdistance_TEvalues_raw)
sd(tonnetzdistance_TEvalues_raw)
mean(tonnetzdistance_TEvalues_random_pertrial)
sd(tonnetzdistance_TEvalues_random)

wilcox.test(tonnetzdistance_TEvalues_raw, tonnetzdistance_TEvalues_random_pertrial)

RQ1a.2_df = data.frame(tonnetzdistance_TEvalues_raw, tonnetzdistance_TEvalues_random_pertrial)
names(RQ1a.2_df) = c("Tonnetz dist. ETE real", "Tonnetz dist. ETE random")
RQ1a.2_df = gather(RQ1a.2_df)
ggplot(RQ1a.2_df, aes(x=value, y=key, fill=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("ETE") +
  ylab("Random-real") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  theme(legend.position="none", text = element_text(size = 17))

# Spectral flatness TE
qqnorm(spectralflatness_TEvalues_raw)

mean(spectralflatness_TEvalues_raw)
sd(spectralflatness_TEvalues_raw)
mean(spectralflatness_TEvalues_random_pertrial)
sd(spectralflatness_TEvalues_random_pertrial)

wilcox.test(spectralflatness_TEvalues_raw, spectralflatness_TEvalues_random_pertrial)

RQ1a.3_df = data.frame(spectralflatness_TEvalues_raw, spectralflatness_TEvalues_random_pertrial)
names(RQ1a.3_df) = c("Spectral fltn. ETE real", "Spectral fltn. ETE random")
RQ1a.3_df = gather(RQ1a.3_df)
ggplot(RQ1a.3_df, aes(x=value, y=key, fill=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("ETE") +
  ylab("Random-real") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  theme(legend.position="none", text = element_text(size = 17))

# RQ1b - same as before, but now with rho values

# RMS amplitude rho
qqnorm(rms_rhovalues_raw)

mean(rms_rhovalues_raw)
sd(rms_rhovalues_raw)
mean(rms_rhovalues_random)
sd(rms_rhovalues_random)

wilcox.test(rms_rhovalues_raw, rms_rhovalues_random)

RQ1b.1_df = data.frame(rms_rhovalues_raw, rms_rhovalues_random)
names(RQ1b.1_df) = c("RMS rho real", "RMS rho random")
RQ1b.1_df = gather(RQ1b.1_df)
ggplot(RQ1b.1_df, aes(x=value, y=key, fill=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("Rho") +
  ylab("Random-real") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  theme(legend.position="none", text = element_text(size = 17))

# Tonnetz distance rho
qqnorm(tonnetzdistance_rhovalues_raw)

mean(tonnetzdistance_rhovalues_raw)
sd(tonnetzdistance_rhovalues_raw)
mean(tonnetzdistance_rhovalues_random)
sd(tonnetzdistance_rhovalues_random)

wilcox.test(tonnetzdistance_rhovalues_raw, tonnetzdistance_rhovalues_random)

RQ1b.2_df = data.frame(tonnetzdistance_rhovalues_raw, tonnetzdistance_rhovalues_random)
names(RQ1b.2_df) = c("Tonnetz dist. rho real", "Tonnetz dist. rho random")
RQ1b.2_df = gather(RQ1b.2_df)
ggplot(RQ1b.2_df, aes(x=value, y=key, fill=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("Rho") +
  ylab("Random-real") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black") +
  theme(legend.position="none", text = element_text(size = 17))

# Spectral flatness rho
qqnorm(spectralflatness_rhovalues_raw)

mean(spectralflatness_rhovalues_raw)
sd(spectralflatness_rhovalues_raw)
mean(spectralflatness_rhovalues_random)
sd(spectralflatness_rhovalues_random)

wilcox.test(spectralflatness_rhovalues_raw, spectralflatness_rhovalues_random)

RQ1b.3_df = data.frame(spectralflatness_rhovalues_raw, spectralflatness_rhovalues_random)
names(RQ1b.3_df) = c("Spectral fltn. rho real", "Spectral fltn. rho random")
RQ1b.3_df = gather(RQ1b.3_df)
ggplot(RQ1b.3_df, aes(x=value, y=key, fill=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("Rho") +
  ylab("Random-real") + 
  stat_summary(fun = "mean",
               geom = "point",
               colour = "black") +
  theme(legend.position="none", text = element_text(size = 17))

#_____________________________________________________________________________________

# RQ2 ------------

#_____________________________________________________________________________________

rms_TEvalues_perpair[[sprintf("%s-%s", "g2_t1_b3", "g2_t1_b1")]]

smooth = min(unlist(rms_TEvalues_perpair)[unlist(rms_TEvalues_perpair) > 0])

options(scipen=999)

for (row in 1:nrow(trialDataset)){
  # Add trial-wide rho and TE values to dataframe
  group = trialDataset[row, "trio"]
  trial = trialDataset[row, "take"]
  trialDataset$rho_full[row] = rms_rhovalues[[sprintf("g%s_t%s", group, trial)]] + tonnetzdistance_rhovalues[[sprintf("g%s_t%s", group, trial)]] + spectralflatness_rhovalues[[sprintf("g%s_t%s", group, trial)]]
  trialDataset$TE_full[row] = rms_TEvalues_pertrial[[sprintf("g%s_t%s", group, trial)]] + tonnetzdistance_TEvalues_pertrial[[sprintf("g%s_t%s", group, trial)]] + spectralflatness_TEvalues_pertrial[[sprintf("g%s_t%s", group, trial)]]
  trialDataset$TE_rms_full[row] = rms_TEvalues_pertrial[[sprintf("g%s_t%s", group, trial)]]
  trialDataset$rho_rms_full[row] = rms_rhovalues[[sprintf("g%s_t%s", group, trial)]]
  
  unidirectionality_indices = c()
  
  # Get all pairs within trial
  pairs = combn(c(1:3), m=2)
  for (pair in 1:3){
    colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
    colname2 = sprintf("g%s_t%s_b%s", group, trial, pairs[2,pair])
    
    # Retrieve pairwise TE values from list
    pairTExy = rms_TEvalues_perpair[[sprintf("%s-%s", colname1, colname2)]] + smooth
    pairTEyx = rms_TEvalues_perpair[[sprintf("%s-%s", colname2, colname1)]] + smooth 
    
    # Calculate the unidirectionality index of a pair. Unidirectionality index is 
    # always a value of at least 1, with 1 signifying perfect bidirectionality.
    unidirectionality_index_pair = max(c(pairTExy / pairTEyx, pairTEyx / pairTExy))
    if(unidirectionality_index_pair > 10){
      unidirectionality_index_pair = 10
    }
    # Append pair unidirectionality index to vector of unidirectionality indices in the trial.
    unidirectionality_indices = c(unidirectionality_indices, unidirectionality_index_pair)
  }
  
  # Take average of unidirectionality indices within the trial to get the unidirectionality
  # index for the whole trial. 'full' here indicates that we take TE over the whole performance,
  # not just the part after the prompt.
  trialDataset$unidirectionality_index_full[row] = mean(unidirectionality_indices)
}

# Get data in suitable format for statistical testing
trialDataset$TE_full = unlist(trialDataset$TE_full)
trialDataset$rho_full = unlist(trialDataset$rho_full)
trialDataset$TE_rms_full = unlist(trialDataset$TE_rms_full)
trialDataset$rho_rms_full = unlist(trialDataset$rho_rms_full)
trialDataset$unidirectionality_index_full = unlist(trialDataset$unidirectionality_index_full)

trialDataset$rho_rms_full_scaled = scale(trialDataset$rho_rms_full)
trialDataset$TE_rms_full_scaled = scale(trialDataset$TE_rms_full)

# RQ2a
# Does the amount of information flow predict subjective quality of improvisations?
RQ2.1 = lmer(goodImproAll ~ TE_rms_full + rho_rms_full + unidirectionality_index_full + (1|trio),
            trialDataset[trialDataset$take > 4,])

summary(RQ2.1)
check_model(RQ2.1)

RQ2.2 = lmer(goodImproAll ~ rho_rms_full + TE_rms_full + (1|trio),
            trialDataset[trialDataset$take > 4,])
summary(RQ2.2)

plot(compare_performance(RQ2.1, RQ2.2))

RQ2.3 = lmer(goodImproAll ~ rho_rms_full + (1|trio),
             trialDataset[trialDataset$take > 4,])
summary(RQ2.3)

compare_performance(RQ2.1, RQ2.2, RQ2.3)

anova(RQ2.3, RQ2.4)
plot(compare_performance(RQ2.3, RQ2.4))

compare_performance(RQ2.3, RQ2.4)

View(trialData)
# Plot the relationship
RQ2.1_df = data.frame(trialDataset[trialDataset$take > 4,]$goodImproAll, trialDataset[trialDataset$take > 4,]$unidirectionality_index_full_bin)
RQ2.1_df = RQ2a_df[complete.cases(RQ2.1_df),]
names(RQ2.1_df) = c("Subjective_Quality", "Transfer_Entropy")
RQ2.1_df$predlmer = predict(RQ2.1)
ggplot(RQ2.1_df, aes(x=Transfer_Entropy, y=Subjective_Quality)) + 
  geom_point(na.rm=TRUE) +
  geom_smooth(aes(y = predlmer), size = 1) +
  scale_y_continuous(breaks=seq(1,7))


RQ2TE = lmer(goodImproAll ~ TE_rms_full + (1|trio),
             trialDataset[trialDataset$take > 4,])
RQ2TE_df = data.frame(trialDataset[trialDataset$take > 4,]$goodImproAll, trialDataset[trialDataset$take > 4,]$TE_rms_full)
RQ2TE_df = RQ2TE_df[complete.cases(RQ2TE_df),]
names(RQ2TE_df) = c("Subjective_Quality", "ETE")
RQ2TE_df$predlmer = predict(RQ2TE)

ggplot(RQ2TE_df, aes(x=ETE, y=Subjective_Quality)) + 
  geom_point(na.rm=TRUE) +
  geom_smooth(aes(y = predlmer), size = 1) +
  scale_y_continuous(breaks=seq(1,7)) + 
  theme(text = element_text(size = 22))

RQ2Uni = lmer(goodImproAll ~ unidirectionality_index_full + (1|trio),
             trialDataset[trialDataset$take > 4,])
RQ2Uni_df = data.frame(trialDataset[trialDataset$take > 4,]$goodImproAll, trialDataset[trialDataset$take > 4,]$unidirectionality_index_full)
RQ2Uni_df = RQ2Uni_df[complete.cases(RQ2Uni_df),]
names(RQ2Uni_df) = c("Subjective_Quality", "Unidirectionality_Index")
RQ2Uni_df$predlmer = predict(RQ2Uni)

ggplot(RQ2Uni_df, aes(x=Unidirectionality_Index, y=Subjective_Quality)) + 
  geom_point(na.rm=TRUE) +
  geom_smooth(aes(y = predlmer), size = 1) +
  scale_y_continuous(breaks=seq(1,7)) + 
  theme(text = element_text(size = 22))

RQ2Rho_df = data.frame(trialDataset[trialDataset$take > 4,]$goodImproAll, trialDataset[trialDataset$take > 4,]$rho_rms_full)
RQ2Rho_df = RQ2Rho_df[complete.cases(RQ2Rho_df),]
names(RQ2Rho_df) = c("Subjective_Quality", "Rho")
RQ2Rho_df$predlmer = predict(RQ2.3)

ggplot(RQ2Rho_df, aes(x=Rho, y=Subjective_Quality)) + 
  geom_point(na.rm=TRUE) +
  geom_smooth(aes(y = predlmer), size = 1) +
  scale_y_continuous(breaks=seq(1,7)) + 
  theme(text = element_text(size = 22))

#_____________________________________________________________________________________

# Compute post-prompt TE, post-prompt rho --------

#_____________________________________________________________________________________

# Add empty columns to dataframe
trialDataset$TE_after_prompt = c(rep(NA, nrow(trialDataset)))
trialDataset$TE_before_prompt = c(rep(NA, nrow(trialDataset)))
trialDataset$unidirectionality_index_after_prompt = c(rep(NA, nrow(trialDataset)))

# Initialize empty list of pairwise post-prompt TE values
TE_perpair_after_prompt = list()
TE_perpair_before_prompt = list()

for(group in 1:12) {
  for(trial in 5:16) {
    # Get info about the current trial
    idx = which(trialDataset$trio == group & trialDataset$take == trial)
    promptTime = trialDataset$promptTime[idx]
    promptWindow = ceiling(promptTime / windowSizeInSeconds)
    promptType = trialDataset$type[idx]
    promptNumber = trialDataset$number[idx]
    endPoint = endPoints[[sprintf("g%s_t%s", group, trial)]]
    
    # Get TE of the part of the improvisation that occurred after the prompt.
    # Only include trials where all 3 musicians continued playing for at least 10 seconds after the prompt
    if (!(promptWindow %in% c(0, NA)) & endPoint - promptWindow >= round(10/windowSizeInSeconds)) {
      pairs = combn(c(1:3), m=2)
      TE_after_prompt = c()
      TE_before_prompt = c()
      for(pair in 1:3){
        
        colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
        colname2 = sprintf("g%s_t%s_b%s", group, trial, pairs[2,pair])
        
        TExyPair = 0
        TEyxPair = 0
        TExyPairBeforePrompt = 0
        TEyxPairBeforePrompt = 0
        
        for(markov_order in 1:20){
          TExy = calc_ete(rms_timeseries[,colname1][promptWindow:endPoint],
                         rms_timeseries[,colname2][promptWindow:endPoint], lx=markov_order, ly=markov_order, seed=TRUE)
          TEyx = calc_ete(rms_timeseries[,colname2][promptWindow:endPoint],
                         rms_timeseries[,colname1][promptWindow:endPoint], lx=markov_order, ly=markov_order, seed=TRUE)
          
          TExyBeforePrompt = calc_ete(rms_timeseries[,colname1][0:promptWindow],
                          rms_timeseries[,colname2][0:promptWindow], lx=markov_order, ly=markov_order, seed=TRUE)
          TEyxBeforePrompt = calc_ete(rms_timeseries[,colname2][0:promptWindow],
                          rms_timeseries[,colname1][0:promptWindow], lx=markov_order, ly=markov_order, seed=TRUE)
          TExyPair = TExyPair + TExy
          TEyxPair = TEyxPair + TEyx
          TExyPairBeforePrompt = TExyPairBeforePrompt + TExyBeforePrompt
          TEyxPairBeforePrompt = TEyxPairBeforePrompt + TEyxBeforePrompt
        }
        
        TExyPair = TExy / 20
        TEyxPair = TEyx / 20
      
        # Save pairwise post-prompt TE values
        TE_perpair_after_prompt[sprintf("%s-%s", colname1, colname2)] = TExyPair
        TE_perpair_after_prompt[sprintf("%s-%s", colname2, colname1)] = TEyxPair
        TE_perpair_before_prompt[sprintf("%s-%s", colname1, colname2)] = TExyPairBeforePrompt
        TE_perpair_before_prompt[sprintf("%s-%s", colname2, colname1)] = TEyxPairBeforePrompt
        
        TE_after_prompt = c(TE_after_prompt, TExyPair, TEyxPair)
        TE_before_prompt = c(TE_after_prompt, TExyPairBeforePrompt, TEyxPairBeforePrompt)
      }
      # Save trial-wide post-prompt TE value
      trialDataset$TE_after_prompt[idx] = mean(TE_after_prompt)
      trialDataset$TE_before_prompt[idx] = mean(TE_before_prompt)
    }
  }
  # Logging
  print(sprintf("TE ENDINGS g%s", group))
}

trialDataset$rho_before_prompt = c(rep(NA, nrow(trialDataset)))
trialDataset$rho_after_prompt = c(rep(NA, nrow(trialDataset)))
trialDataset$rho_change = c(rep(NA, nrow(trialDataset)))

for (group in 1:12){
  for (trial in 5:16){
    # Get info about the current trial
    idx = which(trialDataset$trio == group & trialDataset$take == trial)
    promptTime = trialDataset$promptTime[idx]
    promptWindow = ceiling(promptTime / windowSizeInSeconds)
    promptType = trialDataset$type[idx]
    promptNumber = trialDataset$number[idx]
    endPoint = endPoints[[sprintf("g%s_t%s", group, trial)]]
    
    colnames = c(sprintf("g%s_t%s_b1", group, trial), 
                 sprintf("g%s_t%s_b2", group, trial), 
                 sprintf("g%s_t%s_b3", group, trial))
    
    rmsDF = cbind(rms_timeseries[colnames[1]], rms_timeseries[colnames[2]], rms_timeseries[colnames[3]])
    
    # again, select only those trials with at least 10 s of group performance after the prompt
    if (!(promptWindow %in% c(0, NA)) & endPoint - promptWindow >= round(10/windowSizeInSeconds)) {
      for(colname in colnames){
        rhosAfterPrompt = c()
        rhosBeforePrompt = c()
        for(forecast_time in 1:20){
          
          # calculate the pre-prompt Rho using Simplex algorithm twice: first predicting 2nd half of performance based
          # on 1st half, then predicting 1st half based on second half
          simplex_before_prompt1 = block_lnlp(rmsDF, lib = c(1, floor(promptWindow/2)), pred = c(floor(promptWindow/2)+1, promptWindow),
                                            target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
          simplex_before_prompt2 = block_lnlp(rmsDF, lib = c(floor(promptWindow/2)+1, promptWindow), pred = c(1, floor(promptWindow/2)),
                                             target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
          # calculate post-prompt Rho using everything before the prompt as the manifold
          simplex_after_prompt = block_lnlp(rmsDF, lib = c(1, promptWindow), pred = c(promptWindow+1, endPoint),
                                       target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
          rhosAfterPrompt = c(rhosAfterPrompt, simplex_after_prompt$stats$rho$rho)
          rhosBeforePrompt = c(rhosBeforePrompt, mean(simplex_before_prompt1$stats$rho$rho, simplex_before_prompt2$stats$rho$rho))
        }
        trialRhoAfterPrompt = mean(rhosAfterPrompt)
        trialRhoBeforePrompt = mean(rhosBeforePrompt)
      }
      # add to DF
      trialDataset$rho_after_prompt[idx] = trialRhoAfterPrompt
      trialDataset$rho_before_prompt[idx] = trialRhoBeforePrompt
    }
    print(sprintf("g%s_t%s", group, trial))
  }
  # Logging
  print(sprintf("EDM ENDINGS g%s", group))
}

# Get raw values
trialDataset$TE_after_prompt = unlist(trialDataset$TE_after_prompt)
trialDataset$rho_after_prompt = unlist(trialDataset$rho_after_prompt)

#_____________________________________________________________________________________

# RQ3 ------------

#_____________________________________________________________________________________

# Turn into factor variables
trialDataset$type = factor(trialDataset$type)
trialDataset$trio = factor(trialDataset$trio)

# RQ3a
RQ3a.1 = lmer(sqrt(TE_after_prompt) ~ number + type + number:type + (1|trio),
              trialDataset[!is.na(trialDataset$TE_after_prompt),])
summary(RQ3a.1)

RQ3a.2 = lmer(sqrt(TE_after_prompt) ~ type + number:type + (1|trio),
              trialDataset[!is.na(trialDataset$TE_after_prompt),])
summary(RQ3a.2)


RQ3a.3 = lmer(sqrt(TE_after_prompt) ~ type + (1|trio),
              trialDataset[!is.na(trialDataset$TE_after_prompt),])
summary(RQ3a.3)

compare_performance(RQ3a.1, RQ3a.2, RQ3a.3)


trialDataset[trialDataset$type==1, "TE_after_prompt"]

RQ3a_df = data.frame(rms_TEvalues_raw, rms_TEvalues_random)
names(RQ1a_df) = c("RMS TE real", "RMS TE random")
RQ1a_df = gather(RQ1a_df)
ggplot(RQ1a_df, aes(x=value, y=key)) + 
  geom_violin() +
  coord_flip() +
  xlab("Transfer Entropy") +
  ylab("Random-real") +
  stat_summary(fun = "mean",
               geom = "point", 
               colour = "black")

effects_RQ3a.3 <- as.data.frame(effect(term="type", mod=RQ3a.3))
urchin_plot <- ggplot() + 
  geom_point(data=trialDataset[!is.na(trialDataset$TE_after_prompt),], aes(type, sqrt(TE_after_prompt))) + 
  geom_point(data=effects_RQ3a.3, aes(x=type, y=fit), color="blue") +
  geom_line(data=effects_RQ3a.3, aes(x=type, y=fit), color="blue") +
  geom_ribbon(data= effects_RQ3a.3, aes(x=type, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  labs(x="Urchins (centered & scaled)", y="Coral Cover")
urchin_plot

step(RQ3a.full, direction = "both")
# Save best fitting model
RQ3a.final = lmer(TE_after_prompt ~ type + (1|trio), 
                  trialDataset[!is.na(trialDataset$TE_after_prompt),])
summary(RQ3a.final)

compare_performance(RQ3a.full, RQ3a.final)

View(musicianDataset)

# RQ3b
RQ3b.1 = lmer(rho_after_prompt ~ number + type + number:type + (1|trio),
                 trialDataset[!is.na(trialDataset$rho_after_prompt),])
summary(RQ3b.1)

qqnorm(residuals(RQ3b.1))

RQ3b.2 = lmer(rho_after_prompt ~ number + type + (1|trio), 
              trialDataset[!is.na(trialDataset$rho_after_prompt),])
summary(RQ3b.2)

RQ3b.3 = lmer(rho_after_prompt ~ number + (1|trio), 
              trialDataset[!is.na(trialDataset$rho_after_prompt),])
summary(RQ3b.3)

compare_performance(RQ3b.1, RQ3b.2, RQ3b.3)

# Compare means

MuMIn::r.squaredGLMM(RQ3b.final)
confint(RQ3b.final)

mean(trialDataset$rho_after_prompt[!is.na(trialDataset$rho_after_prompt) & trialDataset$number == 1])
mean(trialDataset$rho_after_prompt[!is.na(trialDataset$rho_after_prompt) & trialDataset$number == 2])
mean(trialDataset$rho_after_prompt[!is.na(trialDataset$rho_after_prompt) & trialDataset$number == 3])

#_____________________________________________________________________________________

# RQ4 ------------

#_____________________________________________________________________________________

# Initialize empty vectors
from = c()
to = c()
from_prompt = c()
to_prompt = c()
TE_pair_full = c()
TE_pair_afterprompt = c()
directionality_afterprompt = c()

# Cast to factor variable
musicianDataset$booth = factor(musicianDataset$booth)
levels(musicianDataset$booth) = c(1,2,3)

for(group in 1:12) {
  for(trial in 1:16) {
    pairs = combn(c(1:3), m=2)
    for(pair in 1:3){
      
      # What prompt did each of the musicians in the pair hear?
      prompt_heard_1 = musicianDataset[musicianDataset$trio == group & musicianDataset$take == trial & musicianDataset$booth == pairs[1,pair],]$prompt_heard
      prompt_heard_2 = musicianDataset[musicianDataset$trio == group & musicianDataset$take == trial & musicianDataset$booth == pairs[2,pair],]$prompt_heard
      
      # CHECK LATER
      if (length(prompt_heard_1) != 0 & length(prompt_heard_2) != 0){
        colname1 = sprintf("g%s_t%s_b%s", group, trial, pairs[1,pair])
        colname2 = sprintf("g%s_t%s_b%s", group, trial, pairs[2,pair])
        from = c(from, colname1, colname2)
        to = c(to, colname2, colname1)
        from_prompt = c(from_prompt, prompt_heard_1, prompt_heard_2)
        to_prompt = c(to_prompt, prompt_heard_2, prompt_heard_1)
        # Retrieve pairwise TE for the full performance, for RQ4c
        TE_pair_full = c(TE_pair_full,
                         rms_TEvalues_perpair[sprintf("%s-%s", colname1, colname2)], 
                         rms_TEvalues_perpair[sprintf("%s-%s", colname2, colname1)])
        # Retrieve post-prompt pairwise TE, for use in RQ4a
        TE_pair_afterprompt = c(TE_pair_afterprompt, 
                                 TE_perpair_after_prompt[sprintf("%s-%s", colname1, colname2)],
                                 TE_perpair_after_prompt[sprintf("%s-%s", colname2, colname1)])
      }
    }
  }
}

# Create dataframe with the new variables
RQ4_pair_df = data.frame(from = from)
RQ4_pair_df$to = to
RQ4_pair_df$TE_pair_full = TE_pair_full
RQ4_pair_df$TE_pair_afterprompt = TE_pair_afterprompt 
RQ4_pair_df$from_prompt = from_prompt
RQ4_pair_df$to_prompt = to_prompt
RQ4_pair_df$directionality_ratio_after_prompt = rep(NA, nrow(RQ4_pair_df))

for(row in 1:nrow(RQ4_pair_df)){
  # Only add post-prompt directionality ratio to DF once per pair. Adding this twice for one pair
  # twice for a pair (so g1_t1_b1-->g1_t1_b2 AND g1_t1_b2-->g1_t1_b1) would not add any information.
  if(row%%2==1){
    if(length(unname(unlist(RQ4_pair_df$TE_pair_afterprompt[row])))!=0){
      RQ4_pair_df$directionality_ratio_after_prompt[row] = unname(unlist(RQ4_pair_df$TE_pair_afterprompt[row])) / unname(unlist(RQ4_pair_df$TE_pair_afterprompt[row+1]))
    }
  }
}

# get raw numeric values
for(col in names(RQ4_pair_df)){
  RQ4_pair_df[col] = unlist(RQ4_pair_df[col])
}

# Add columns for change in individual (i.e. per musician) predictability following the prompt,
# and for individual predictability during the full trial
RQ4_pair_df$individual_predictability_change = rep(NA, nrow(RQ4_pair_df))
RQ4_pair_df$individual_predictability_full = rep(NA, nrow(RQ4_pair_df))

for(group in c(1:12)){
  for (trial in c(1:16)){
    # Retrieve info from trial
    idx = which(trialDataset$trio == group & trialDataset$take == trial)
    promptTime = trialDataset$promptTime[idx]
    promptWindow = ceiling(promptTime / windowSizeInSeconds)
    endPoint = endPoints[[sprintf("g%s_t%s", group, trial)]]
    
    colnames = c(sprintf("g%s_t%s_b1", group, trial), 
                 sprintf("g%s_t%s_b2", group, trial), 
                 sprintf("g%s_t%s_b3", group, trial))
    
    rmsDF = cbind(rms_timeseries[colnames[1]], rms_timeseries[colnames[2]], rms_timeseries[colnames[3]])
    
    if(nrow(RQ4_pair_df[RQ4_pair_df$from == colname1,]) > 0){
      # Calculate individual predictability for each musician, using only their own playing as the
      # state space history
      for(colname in colnames){
        total_rho = c()
        rmsTS = unlist(rms_timeseries[colname])
        
        for(forecast_time in 1:20){
          simplex_output1 = block_lnlp(rmsTS, lib = c(1, floor(length(ts) / 2)), pred = c(floor(length(ts)/2)+1, length(ts)),
                                       target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
          simplex_output2 = block_lnlp(rmsDF, lib = c(floor(length(ts)/2)+1, length(ts)), pred = c(1, floor(length(ts) / 2)),
                                       target_column = colname, tp = forecast_time, stats_only = FALSE, first_column_time = FALSE, silent = TRUE)
          
          total_rho = c(total_rho, mean(simplex_output1$const_pred_rho$rho, simplex_output2$const_pred_rho$rho))
        }
        RQ4_pair_df[match(colname, RQ4_pair_df$from),]$individual_predictability_full = mean(total_rho)
      }
      
      # For musicians in trials that continued for at least 10 seconds after prompt, calculate
      # the factor with which their individual predictability changed following the prompt. 
      # For this, we need to calculate their pre-prompt Rho and post-prompt Rho, and then calculate the 
      # difference between these two as a percentage.
      if (!(promptWindow %in% c(0, NA)) & endPoint - promptWindow >= round(10/windowSizeInSeconds)) {
        for(colname in colnames){
          total_rho_after_prompt = c()
          total_rho_before_prompt = c()
          for(forecast_time in 1:20){
            rhoAfterPrompt_forecast_time = 0
            rhoBeforePrompt_forecast_time = 0
            
            simplex_after_prompt = block_lnlp(unlist(rms_timeseries[colname]), lib = c(1, promptWindow), pred = c(promptWindow+1, endPoint), tp=forecast_time)
            simplex_before_prompt1 = block_lnlp(unlist(rms_timeseries[colname]), lib = c(1, floor(promptWindow/2)), pred = c(floor(promptWindow/2)+1, promptWindow), tp=forecast_time)
            simplex_before_prompt2 = block_lnlp(unlist(rms_timeseries[colname]), lib = c(floor(promptWindow/2)+1, promptWindow), pred = c(1, floor(promptWindow/2)), tp=forecast_time)
            rhoAfterPrompt_forecast_time = rhoAfterPrompt_forecast_time + simplex_after_prompt$const_pred_rho$rho
            rhoBeforePrompt_forecast_time = rhoBeforePrompt_forecast_time + mean(simplex_before_prompt1$const_pred_rho$rho, simplex_before_prompt2$const_pred_rho$rho)
            
            total_rho_after_prompt = c(total_rho_after_prompt, rhoAfterPrompt_forecast_time)
            total_rho_before_prompt = c(total_rho_before_prompt, rhoBeforePrompt_forecast_time)
          }
          individual_predictability_change = (mean(total_rho_after_prompt) - mean(total_rho_before_prompt)) / mean(total_rho_before_prompt)
          RQ4_pair_df[match(colname, RQ4_pair_df$from),]$individual_predictability_change = individual_predictability_change
        }
      }
    }
  }
  # Logging
  print(sprintf("Rho calculated for %s", group))
}

# Turn into factors
RQ4_pair_df$from_prompt = factor(RQ4_pair_df$from_prompt)
RQ4_pair_df$to_prompt = factor(RQ4_pair_df$to_prompt)
# Add trio to DF (CHECK LATER)
RQ4_pair_df$trio = rep(1:12, each=nrow(RQ4_pair_df)/12)

# RQ4a
# Use square root transformation, because residuals highly non-normal (CHECK LATER)
RQ4a = lm(sqrt(directionality_ratio_after_prompt) ~ relevel(from_prompt, ref="No-Goal") + relevel(to_prompt, ref="No-Goal"), RQ4_pair_df[!is.na(RQ4_pair_df$directionality_ratio_after_prompt),])
summary(step(RQ4a, direction = "both"))
summary(RQ4a)

hist(sqrt(RQ4_pair_df$directionality_ratio_after_prompt))

RQ4_pair_df[RQ4_pair_df$from=="g1_t3_b1",]$individual_predictability_full

# RQ4b
RQ4b = lm(individual_predictability_change ~ relevel(from_prompt, ref="No-Goal"), RQ4_pair_df[!is.na(RQ4_pair_df$individual_predictability_change),])
summary(step(RQ4b, direction = "both"))
summary(RQ4b)

# RQ4c
RQ4c = lm(TE_pair_full ~ individual_predictability_full, RQ4_pair_df)
summary(step(RQ4c, direction = "both"))
summary(RQ4c)




















