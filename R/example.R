

#' infer a single trio
#' result=infer.trio(M1trio)
#' get the predicted topology
#' which.model=class.vec(result)
#' print(which.model)

#' infer a set of trios:
#' result2 = sapply(exampleTrios[sample(1:10000, 10)], infer.trio)
#' get the predicted topology
#' models = apply(result2, 2, class.vec)
#' print(models)
