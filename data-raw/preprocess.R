load("C:/Users/Bruin/Documents/GitHub/MRGN_extra/exampleWBtrios.RData")
load("C:/Users/Bruin/Documents/GitHub/MRGN_extra/WBpcs.RData")
CNAtrios=readRDS("C:/Users/Bruin/Documents/GitHub/MRGN_extra/100triosdata.Rdata")
M1trio=MRPC::simu_data_M1
WBtrios=wbtriolist[1:100]
WBtrios=lapply(WBtrios, function(x){row.names(x)=c(1:nrow(x));return(x)})
WBscores=PCs.matrix
row.names(WBscores)=c(1:nrow(WBscores))
usethis::use_data(WBtrios, overwrite = T)
usethis::use_data(WBscores, overwrite = T)
usethis::use_data(M1trio, overwrite = T)
usethis::use_data(CNAtrios, overwrite = T)
