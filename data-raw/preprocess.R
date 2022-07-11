load("C:/Users/Bruin/Documents/GitHub/MRGN_extra/exampleWBtrios.RData")
load("C:/Users/Bruin/Documents/GitHub/MRGN_extra/WBpcs.RData")
M1trio=MRPC::simu_data_M1
WBtrios=wbtriolist
WBtrios=lapply(WBtrios, function(x){row.names(x)=c(1:nrow(x));return(x)})
WBscores=PCs.matrix
row.names(WBscores)=c(1:nrow(WBscores))
usethis::use_data(WBtrios[1:200], overwrite = T)
usethis::use_data(WBscores, overwrite = T)
usethis::use_data(M1trio, overwrite = T)
