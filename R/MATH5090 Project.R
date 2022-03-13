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
  probabilities_stage_2 = dbinom(y_2squiggles, n2-n1, theta)
  
  #Matrix to calculate probabilities of each event happening
  prob = matrix(nrow = n1 + 1, ncol = n2 - n1 + 1)
  
  #Calculating the probabilities of every combination of events (i.e. P(y1 = i)*P(y2squiggle = j)) and removing the events when stage 2 is not reached
  for (i in 1:(n1+1)){
    for(j in 1:(n2-n1+1)){
      prob[i,j] = probabilities_stage_1[i]*probabilities_stage_2[j]*(1-does_it_stop_1[i])
    }
  }
  
  #Creating vector to store probabilities
  Probabilities = rep(0, n2+1)
  
  #Calculating the probabilities P(y1 + y2squiggle = i)
  for(i in 1:(n2+1)){
    Probabilities[i] = sum(prob[c(row(prob) + col(prob) == i + 1)]) 
  }
  
  #Calculating P(y1 + y2squiggle = i)/P(stage 2 reached) i.e. conditional probability of i successes given stage 2 is reached
  Probabilities = Probabilities/sum(probabilities_stage_1*(1-does_it_stop_1))
  
  #Progression threshold for stage 2
  C2 = 1 - lambda * (n2 / n2)^gamma
  
  #Decision for each possible case in stage 2
  does_it_stop_2 = pbeta(0.5, a0 + y_2s, b0 + n2 - y_2s) > C2
  
  #Probability of stopping in stage 2
  s2 = sum(Probabilities*does_it_stop_2)
  
  #Returns Type I and II error along with expected sample size. (if else statement used in event when no cases reach stage 2)
  return(c("Probability of declaring successful result (Type I error under the null hypothesis)" = ifelse(C1 == 0, 0, 1-(s1 + (1-s1)*s2)),"Probability of failing to declare a successful result (Type II error under the alternative hypothesis)" = ifelse(C1 == 0, 1, s1 + (1-s1)*s2), "Expected Sample Size" = n1*s1+n2*(1-s1)))
}