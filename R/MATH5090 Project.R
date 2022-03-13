#Time the whole process
ptm <- proc.time()

#Expected Sample Size, Type I and Type II error code
Type_I_II_error_and_Expected_Sample_Size = function(lambda, gamma, n1, n2, theta, a0, b0){
  
  #All possible outcomes in stage 1, stage 2 and combined  
  y_1s = 0:n1
  y_2squiggles = 0:(n2-n1)
  y_2s = 0:n2
  
  #Progression threshold for stage 1
  C1 = 1 - lambda * (n1 / n2)^gamma
  
  #Decision for each possible case in stage 1
  does_it_stop_1 = pbeta(0.5, a0 + y_1s, b0 + n1 - y_1s) > C1
  
  #Probabilities of each outcome in stage 1
  probabilities_stage_1 = dbinom(y_1s, n1, theta)
  
  #Probability of stopping in stage 1
  s1 = sum(probabilities_stage_1*does_it_stop_1)
  
  #Probabilities of each outcome in stage 2 ONLY (not including stage 1 events)
  probabilities_stage_2_only = dbinom(y_2squiggles, n2-n1, theta)
  
  #Matrix to store probabilities of every combination of events (i.e. P(y1 = i)*P(y2squiggle = j))
  prob = matrix(nrow = n1 + 1, ncol = n2 - n1 + 1)
  
  #Calculating the probabilities of every combination of events (i.e. P(y1 = i)*P(y2squiggle = j)) and removing the events when stage 2 is not reached
  for (i in 1:(n1+1)){
    for(j in 1:(n2-n1+1)){
      prob[i,j] = probabilities_stage_1[i]*probabilities_stage_2_only[j]*(1-does_it_stop_1[i])
    }
  }
  
  #Creating vector to store probabilities P(y1 + y2squiggle = i)
  probabilities_stage_2 = rep(0, n2+1)
  
  #Calculating the probabilities P(y1 + y2squiggle = i)
  for(i in 1:(n2+1)){
    probabilities_stage_2[i] = sum(prob[c(row(prob) + col(prob) == i + 1)]) 
  }
  
  #Calculating P(y1 + y2squiggle = i)/P(stage 2 reached) i.e. conditional probability of i successes given stage 2 is reached
  probabilities_stage_2 = probabilities_stage_2/sum(probabilities_stage_1*(1-does_it_stop_1))
  
  #Progression threshold for stage 2
  C2 = 1 - lambda * (n2 / n2)^gamma
  
  #Decision for each possible case in stage 2
  does_it_stop_2 = pbeta(0.5, a0 + y_2s, b0 + n2 - y_2s) > C2
  
  #Probability of stopping in stage 2
  s2 = sum(probabilities_stage_2*does_it_stop_2)
  
  #Returns Type I and II error along with expected sample size. (if else statement used in event when no cases reach stage 2)
  return(c("Probability of declaring successful result (Type I error under the null hypothesis)" = ifelse(C1 == 0, 0, 1-(s1 + (1-s1)*s2)),"Probability of failing to declare a successful result (Type II error under the alternative hypothesis)" = ifelse(C1 == 0, 1, s1 + (1-s1)*s2), "Expected Sample Size" = n1*s1+n2*(1-s1)))
}

#Limits for gamma and lambda in grid search
lam_min = 0
lam_max = 1
gam_min = 0
gam_max = 10

#Matrices to store Type I and II errors along with Expected sample size for each gamma and lambda
TypeI = matrix(nrow = 101, ncol = 101)
TypeII = matrix(nrow = 101, ncol = 101)
SampleSize = matrix(nrow = 101, ncol = 101)

#Chosen values for sample size in each stage
n_1 = 25
n_2 = 80

#Creates sequence of lambda and gamma and names columms and rows using these sequences
lamb = seq(from = lam_min, to = lam_max, by = (lam_max  - lam_min)/100)
gam = seq(from = gam_min, to = gam_max, by = (gam_max - gam_min)/100)
colnames(SampleSize) = colnames(TypeII) = colnames(TypeI) = gam
rownames(SampleSize) = rownames(TypeII) = rownames(TypeI) = lamb

#Double loop to simulate Type I and Type II errors, if passes criteria expected sample size also 
for (i in 1:101){
  for (j in 1:101){
    TypeI[i,j] = Type_I_II_error_and_Expected_Sample_Size(lambda = lamb[i], gamma = gam[j], n1 = n_1, n2 = n_2, theta = 0.5, a0 = 0.5, b0 = 0.5)[1]
    TypeII[i,j] = Type_I_II_error_and_Expected_Sample_Size(lambda = lamb[i], gamma = gam[j], n1 = n_1, n2 = n_2, theta = 0.7, a0 = 0.5, b0 = 0.5)[2]
    if (TypeII[i,j] <= 0.2 & TypeI[i,j] <= 0.05){
      SampleSize[i,j] = Type_I_II_error_and_Expected_Sample_Size(lambda = lamb[i], gamma = gam[j], n1 = n_1, n2 = n_2, theta = 0.5, a0 = 0.5, b0 = 0.5)[3]
    }
  }
}

#Outputs the gamma and lambda that minimise the expected sample size, along with the expected sample size itself
c("Gamma" = gam[ceiling(which.min(SampleSize)/101)], "Lambda" = lamb[which.min(SampleSize)%%101], "Min" = min(SampleSize, na.rm = TRUE))

#Time the whole process
proc.time() - ptm