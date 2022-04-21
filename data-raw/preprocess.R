load("C:/Users/Bruin/Documents/GitHub/MRGN_extra/exampleWBtrios.RData")
load("C:/Users/Bruin/Documents/GitHub/MRGN_extra/WBpcs.RData")
M1trio=MRPC::simu_data_M1
WBtrios=wbtriolist
WBscores=PCs.matrix
usethis::use_data(WBtrios, overwrite = T)
usethis::use_data(WBscores, overwrite = T)
usethis::use_data(M1trio, overwrite = T)
