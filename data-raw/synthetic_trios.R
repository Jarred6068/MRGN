
## this script generates data using the simData.from.graph function to simulate one of the 5 model topologies
## with confounder, intermediate and common child variables

#Simulate All Data and Parameters:
library("MRGN")

#set.seed(234)
#pre-allocate simulation conditions
#model.types=c("model0","model1")
model.types=c("model0","model1","model2","model3","model4")
number.of.datasets=100
#preallocate the sim parameters for later saving::
params=as.data.frame(matrix(0, nrow = number.of.datasets, ncol = 6))
colnames(params)=c("model", "SD", "minor.freq", "b.snp","b.med",
                   "num.confounder", 'num.intermediate','num.common.child')

#preallocate list for all simulated datasets
sim.datasets=list()

print("Simulating Data...")
#simulate all datasets
resamp.vec = NULL
pc.name.list = list()
for(j in 1:number.of.datasets){
  init=0
  k=1
  resamples=NULL
  #simulate trio under parameters
  #use while statement to resample until all genotypes are represented
  #simulate parameters
  params$model[j]=sample(model.types, 1)
  #simulate minor allele frequency between 1 and 50%
  params$minor.freq[j]=sample(seq(0.01, 0.5, 0.01),1)
  #simulate snp effects and mediation effects
  params$b.snp[j]=sample(seq(0.5, 1.5, 0.1), 1)
  params$b.med[j]=sample(seq(0.5, 1, 0.1), 1)
  #simulate Std error of residuals (noise) to be between 1/3 and 1 + 1/2 times mediation signal
  params$SD[j]=(sample(seq(0.3, 1.5, 0.1), 1)*params$b.med[j])
  #generate the number of each type of confounding variable so that
  params$num.confounder[j] = sample(1:5, 1)
  params$num.intermediate[j] = 1
  params$num.common.child[j] = 1

  #checkpoint
  print(paste0("printing params for sim trio ", j))
  print(params[j,])

  while(length(init)<=2){




    #simulate data under parameters
    X = simData.from.graph(model = params$model[j],
                           theta = params$minor.freq[j],
                           b0.1 = 0,
                           b.snp = params$b.snp[j],
                           b.med = params$b.med[j],
                           sd.1 = params$SD[j],
                           conf.num.vec = c(K = 0, U = params$num.confounder[j],
                                            W = params$num.intermediate[j],
                                            Z = params$num.common.child[j]),
                           simulate.confs = TRUE,
                           sample.size = 500,
                           plot.graph = FALSE,
                           conf.coef.ranges = list(K = c(0, 0),
                                                   U = c(0.15, 0.5),
                                                   W = c(0.15, 0.5),
                                                   Z = c(1, 1.5)))

    init=unique(X$data$V1)
    resamples[k]=k
    k = k+1
  }

  print(paste0("resampled ", k, " times before 3 genotypes were represented"))
  colnames(X$data) = paste0(colnames(X$data),".", j)
  sim.datasets[[j]]=X
}

####post-processing####
sim.datasets2 = lapply(sim.datasets, function(x) x$data)
#extract only confounding vars
confs.only.list = lapply(sim.datasets2, function(x) x[,-c(1:3)])
#bind into confounder matrix
conf.mat = do.call("cbind", confs.only.list)
print("Done!...Saving...")

#save files in MRGN_extra
save(sim.datasets, file = "C:/Users/Bruin/Documents/GitHub/MRGN_extra/synthetic_trios_data.RData")
save(conf.mat, file = "C:/Users/Bruin/Documents/GitHub/MRGN_extra/synthetic_trios_confs.RData")
save(params, file = "C:/Users/Bruin/Documents/GitHub/MRGN_extra/synthetic_trios_params.RData")

print("Done!")
