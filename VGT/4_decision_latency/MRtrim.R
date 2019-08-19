# Eve Limbrick-Oldfield Feb 2018

#Funtion to implement van Selst modified recursive RT trimming, with option to NOT remove outliers below
# the min cut off.

# R code adapted from https://figshare.com/articles/RT_Trimming_ToolBox_zip/717189.

# Input is a data frame with four columns: Trial, Participant, condition, RT

# min_crop 0 or 1: If 0, values below the min threshold are not removed
# min_sampleSize: choose not to remove trials of rare events (due to slot machine data having no exp. control)

# Output is object containding data_trimmed (same as input), and number of trials rmeoved per participant/condition, and removed trials.

MRtrim<-function(RT_data,min_crop,minsampleSize){
  
  criterion <- data.frame(sampleSize=c(1:100),modifiedRecursive=c(8.000, 6.200, 5.300, 4.800, 4.475, 
                                                                4.250, 4.110, 4.000, 3.920, 3.850, 
                                                                3.800, 3.750, 3.736, 3.723, 3.709, 
                                                                3.700, 3.681, 3.668, 3.654, 3.640, 
                                                                3.631, 3.622, 3.613, 3.604, 3.595, 
                                                                3.586, 3.577, 3.568, 3.559, 3.550, 
                                                                3.548, 3.546, 3.544, 3.542, 3.540, 
                                                                3.538, 3.536, 3.534, 3.532, 3.530, 
                                                                3.528, 3.526, 3.524, 3.522, 3.520, 
                                                                3.520, 3.516, 3.514, 3.512, 3.510, 
                                                                3.510, 3.510, 3.509, 3.509, 3.509, 
                                                                3.509, 3.509, 3.509, 3.508, 3.508, 
                                                                3.508, 3.508, 3.507, 3.507, 3.507, 
                                                                3.507, 3.507, 3.506, 3.506, 3.506, 
                                                                3.506, 3.506, 3.505, 3.505, 3.505,
                                                                3.505, 3.505, 3.504, 3.504, 3.504, 
                                                                3.504, 3.504, 3.503, 3.503, 3.503, 
                                                                3.503, 3.503, 3.502, 3.502, 3.502,
                                                                3.502, 3.502, 3.501, 3.501, 3.501, 
                                                                3.501, 3.501, 3.500, 3.500, 3.500))
 
  conditions <- subset(RT_data$condition, !duplicated(RT_data$condition))
  numConditions <- length(conditions)
  conditionNames <- vector(length = numConditions)
    
  for (a in 1:numConditions){
    conditionNames[a] <- toString(conditions[a])
  }

  Participants <- subset(RT_data$Participant, !duplicated(RT_data$Participant))
  if (min_crop==1){
    Count_frame <- data.frame(Participant=integer(),condition=integer(),mincount=integer(),maxcount=integer(),count=integer(),ntrials=integer())
  }
  if (min_crop==0){
    Count_frame <- data.frame(Participant=integer(),condition=integer(),count=integer(),ntrials=integer())
  }
  q <- numeric(0)
  rounding <- 0
  mod_RT <- numeric(0)
  p_list <- numeric(0)
  cond_list <- numeric(0)
  trial <- numeric(0)
  count <- 0

  for (i in Participants){ 
    for (cond in conditionNames) {#for each condition
      q = q+1
      mincount<-0
      maxcount<-0
      tempDatasubset <- subset(RT_data, Participant == i & condition == cond)
      restrictedData <- tempDatasubset
      sampleSize <- nrow(restrictedData)
      #if sample is greater than 100, use SDs for 100.
      if(sampleSize > 100){sampleSize <- 100}
      stDev <- criterion$modifiedRecursive[sampleSize]    
  
      if (nrow(restrictedData)>=minsampleSize){
      repeat{  
        x <- max(restrictedData$RT, na.rm = TRUE)#find the largest value in the data structure
        tempData <- restrictedData$RT[restrictedData$RT !=x] #temporarily remove largest value
        sdVal <- sd(tempData,  na.rm = TRUE) #find SD of tecondmporary data    sdMax <- stDev * sdVal #find desired SDs of temporary data
        sdMax <- stDev * sdVal #find desired SDs of temporary data
        maxCutoff <- sdMax + mean(tempData, na.rm = TRUE) #find maximum cutoff value of main data
        minCutoff <- mean(tempData, na.rm = TRUE) - sdMax #find minimumc cutoff value of main data
        
        x <- max(restrictedData$RT,na.rm = TRUE)#find the largest value in the main data structure
        y <- min(restrictedData$RT,na.rm = TRUE)#find the smallest value in the main data structure
        
        removedTrials <-0   
        
        if (sampleSize>=minsampleSize){
          #if there is a data point above the cutoff, remove it
          if(x>maxCutoff ){
            restrictedData$RT[ restrictedData$RT==x ] <- NA
            removedTrials <- 1
            maxcount = maxcount + 1
          } 
          if (min_crop==1){#if there is a data poin below the cutoff, remove it
            if(y<minCutoff){
            restrictedData$RT[restrictedData$RT==y] <- NA
            removedTrials <- 1
            mincount = mincount + 1
            }
          }
          #when there are no trials removed on the current iteration, break the loop.
        }
        if(removedTrials == 0){break}  
        
      }
        
      }
      ntrials<-length(restrictedData$RT)
      trial<-c(trial,restrictedData$Trial)
      mod_RT <- c(mod_RT,restrictedData$RT)
      p_temp<-c(rep(i, length(restrictedData$RT)))
      p_list<-c(p_list,p_temp)
  
      if (min_crop == 1){
        count <-maxcount+mincount
        if (nrow(restrictedData)<minsampleSize){count=NA}
        temp_frame<-data.frame(i,cond,mincount,maxcount,count,ntrials)
      }
      if (min_crop == 0){
        count <-maxcount
        if (nrow(restrictedData)<minsampleSize){count=NA}
        temp_frame<-data.frame(i,cond,count,ntrials)
      }
      Count_frame<-rbind(Count_frame,temp_frame)
    
    }
  }
   
  data_output<-data.frame(trial,p_list,mod_RT)
  data_output <- data_output[order(data_output$p_list,data_output$trial ),] 
  removed_trials<- data_output[is.na(data_output$mod_RT),]
  toReturn<- list("trimmed_data_withNAs"=data_output,"removedCount"=Count_frame,"removedTrials"=removed_trials) 
  return(toReturn) 

}
