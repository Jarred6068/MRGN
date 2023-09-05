load(file = "C:/Users/Bruin/Documents/GitHub/MRGN_extra/exampleWBtrios.RData")
load(file = "C:/Users/Bruin/Documents/GitHub/MRGN_extra/WBpcs.RData")
load(file = "C:/Users/Bruin/Documents/GitHub/MRGN_extra/synthetic_trios_data.RData")
load(file = "C:/Users/Bruin/Documents/GitHub/MRGN_extra/synthetic_trios_confs.RData")
load(file = "C:/Users/Bruin/Documents/GitHub/MRGN_extra/synthetic_trios_params.RData")
#Bandita's CNA trios
CNAtrios=readRDS("C:/Users/Bruin/Documents/GitHub/MRGN_extra/100triosdata.Rdata")
#M1 example trio
M1trio=MRPC::simu_data_M1
#WBtrios
WBtrios=wbtriolist[1:100]
WBtrios=lapply(WBtrios, function(x){row.names(x)=c(1:nrow(x));return(x)})
WBscores=PCs.matrix
triomat = do.call('cbind', WBtrios)
snp.idx = seq(1, dim(triomat)[2]-2, 3)
WBsnps = triomat[,snp.idx]
WBgenes = triomat[,-snp.idx]
row.names(WBscores)=c(1:nrow(WBscores))
#synthetic trios -- rename
synData = sim.datasets
synTrios = lapply(sim.datasets, function(x) x$data[,1:3])
synConfs = conf.mat
synParams = params
usethis::use_data(WBtrios, overwrite = T)
usethis::use_data(WBscores, overwrite = T)
usethis::use_data(M1trio, overwrite = T)
usethis::use_data(CNAtrios, overwrite = T)
usethis::use_data(WBgenes, overwrite = T)
usethis::use_data(WBsnps, overwrite = T)
usethis::use_data(synData, overwrite = T)
usethis::use_data(synTrios, overwrite = T)
usethis::use_data(synConfs, overwrite = T)
usethis::use_data(synParams, overwrite = T)
