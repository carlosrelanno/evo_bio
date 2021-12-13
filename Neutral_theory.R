# Neutral theory exp simulation

N = 100 # Number of individuals
n_loc = 5 # Number of alleles
u = 0.0001 # Mutation rate
days = 5000

pop = matrix(0, n_loc, N) # Population matrix
S1 = matrix(0, n_loc, days/200) # Container to record single substitution
S2 = matrix(0, n_loc, days/200)
S3 = matrix(0, n_loc, days/200)
S4 = matrix(0, n_loc, days/200)

obs = 0 # Initialize number of records

# Algorithm for growth and mutation
for (day in 1:days) {
  # Reproduction
  offspring = sample(1:N, N/2, replace = TRUE) # Choose N/2 individuals among the population
  pop = cbind(pop[,offspring], pop[,offspring], deparse.level = 0)
  
  # Mutation
  mutations = rpois(1, N*u)
  ind = round(runif(mutations, 1, N)) # Allocate mutations across individuals
  gen = round(runif(mutations, 1, n_loc)) # Allocate mutations across genes
  coords = rbind(cbind(c(gen), c(ind))) # Create an indexing matrix to introduce mutations
  pop[coords] = pop[coords] + 1 # Add 1 change to the selected genes and individualst

  # Data retrieving
  if (day %in% seq(1, 50000, by = 200)) {
    obs =  obs + 1
    for (i in (1:n_loc)) {  # For each gene
      S1[i,obs] = sum(pop[i,]==1)
      S2[i,obs] = sum(pop[i,]==2)
      S3[i,obs] = sum(pop[i,]==3)
      S4[i,obs] = sum(pop[i,]==4)
    }
    print(day)
  }
}
# Visualization
print(rbind(S1[,obs],S2[,obs],S3[,obs],S4[,obs]))
plot(1:obs, S1[1:obs], type = 'l', ylim = c(1, N))
for (i in 1:n_loc) {
  lines(1:obs, S1[i,1:obs], col = i, lwd=2.5)
  lines(1:obs, S2[i,1:obs], col = i, lwd=3.5)
  lines(1:obs, S3[i,1:obs], col = i, lwd=4.5)
  lines(1:obs, S4[i,1:obs], col = i, lwd=5.5)
}
legend(1, 95,legend=c(1:n_loc),
       col=c(1:n_loc), lty=1, cex=0.8,
       title="Locus", text.font=4)
