# infer.triplet
# Internal function
#
infer.triplet <- function (data,
                           type,
                           alpha = 0.01,
                           half1 = FALSE,
                           verbose = FALSE) {
  #preform the regressions in step 1
  lmodel <- stats::lm(data[,3]~., data = data[,-3])
  if(verbose){print(summary(lmodel))}
  coefs <- as.data.frame(summary(lmodel)$coefficients)
  pval <- coefs$`Pr(>|t|)`[3]
  if (half1)
    return(pval)
  rej <- pval < alpha
  model <- if (!rej) {"M1.1"}
  else if (type == 2) {"M2.1"}
  else {"Other"}
  data.frame(rej = rej, pval = pval, Inferred.Model = model)
}

# Half infer.triplet
infer.triplet.half1 <- function (data,
                                 verbose = FALSE) {
  #preform the regressions in step 1
  lmodel <- stats::lm(data[,3]~., data = data[,-3])
  if(verbose){print(summary(lmodel))}
  coefs <- as.data.frame(summary(lmodel)$coefficients)
  pval <- coefs$`Pr(>|t|)`[3]
  return(pval)
}
