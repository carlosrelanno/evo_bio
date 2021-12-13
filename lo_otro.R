library(ggplot2)

rm(list=ls())
replicates = 100
N_values = c(100, 300, 1000, 3000)
bottleneck = c(1, 0.5, 0.25)
days = 10000
n_loc = 1

cosos = matrix(0,0,4)
colnames(cosos) = c('N', 'bn', 't_fix', 'p_fix')
for (bn in bottleneck) {
  for (N in N_values) {
    tfix = matrix(0, replicates, 1)
    # Contenedor para grafica
    prop_mut_matrix = matrix(0, replicates, days/200)
    
    for (r in 1:replicates) {
      
      pop = matrix(0, n_loc, N)
      
      mutated_ind = round(runif(60, 1, N))
      pop[mutated_ind] = pop[mutated_ind] + 1
      
      for (day in 1:days) {
        # Reproduction
        offspring = sample(1:N, N * bn, replace = TRUE) # select offspring
        pop = matrix(rep(pop[offspring], length.out=N), 1, N)
        # print(sum(pop))
        
        if (sum(pop) == 0) {
          # record day of extinction
          break
        }
        
        if (sum(pop) == N) {
          tfix[r, 1] = day
          # record day of fixation
          break
        }
        if (day %in% seq(0, days, by = 200)) {
          prop_mut_matrix[r,day/200] = sum(pop)/N
        }
      } 
    }
    #print(tfix)
    mean_tfix = mean(tfix[tfix[,1] != 0])
    p_fix = length(tfix[tfix[,1] != 0])/replicates
    cosos = rbind(cosos, c(N, bn, mean_tfix, p_fix))
    
    ggplot(data=prop_mut_matrix, aes(x = days/200 )
    )
  }
}
}

print(cosos)

data = as.data.frame(cosos)
data$N = factor(data$N)
data$bn = factor(data$bn)

ggplot(data, aes(x = data$N,y = data$t_fix),fill = data$bn)
