
library(tidyverse)

sample_size = function(r){
  m0_m1 = 99 ## it means m0/m1=99, m0/(m0+m1)=0.99
  alpha = .8 * .05 / (1 + (1 - .05) * .99)
  
  Za = qnorm(1 - alpha / 2)
  Zb = qnorm(.2)
  
  C = 0.5 * log((1+r)/(1-r))
  N = ((Za+Zb)/C)**2 + 3
  return(N)
}

data.frame("Correlation" = seq(.15,1, by=.01),
           "Samples" = unlist(lapply(seq(.15,1, by=.01), function(x){sample_size(x)}))) %>% 
  ggplot() + geom_line(aes(x=Correlation, y=Samples), size=1.5) + theme_bw() + 
  labs(title="Correlation sample size calculations", x= "Correlation Coefficient")