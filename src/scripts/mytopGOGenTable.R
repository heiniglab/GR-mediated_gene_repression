# This function is identical with the GenTable function within topGO
# EXCEPT for the formatting of p.values where eps was set to "eps = 1e-30" leading to some downstream trouble
# we set this eps cutoff for very small pvalues to FALSE

    mytopGOGenTable <- function (object, ..., orderBy = 1, ranksOf = 2, 
                        topNodes = 10, numChar = 40, format.FUN = format.pval, 
                        decreasing = FALSE, useLevels = FALSE) 
    {
      resList <- list(...)
      if (!all(sapply(resList, is, "topGOresult"))) 
        stop("Use: topGOdata, topGOresult_1, topGOresult_2, ..., \"parameters\".")
      if (is.null(names(resList))) 
        names(resList) <- paste("result", 1:length(resList), 
                                sep = "")
      resList <- lapply(resList, score)
      if (length(resList) == 1) {
        orderBy <- ranksOf <- 1
        l <- data.frame(resList)
        names(l) <- ifelse(is.null(names(resList)), "", names(resList))
      }
      else {
        l <- topGO:::.sigAllMethods(resList)
      }
      index <- order(l[, orderBy], decreasing = decreasing)
      l <- l[index, , drop = FALSE]
      if (decreasing) 
        rr <- rank(-l[, ranksOf], ties.method = "first")
      else rr <- rank(l[, ranksOf], ties.method = "first")
      whichTerms <- rownames(l)[1:topNodes]
      l <- l[whichTerms, , drop = FALSE]
      rr <- as.integer(rr[1:topNodes])
      shortNames <- topGO:::.getTermsDefinition(whichTerms, ontology(object), 
                                        numChar = numChar)
      infoMat <- data.frame(`GO ID` = whichTerms, Term = shortNames, 
                            stringsAsFactors = FALSE)
      if (useLevels) {
        nodeLevel <- buildLevels(graph(object), leafs2root = TRUE)
        nodeLevel <- unlist(mget(whichTerms, envir = nodeLevel$nodes2level))
        infoMat <- data.frame(infoMat, Level = as.integer(nodeLevel))
      }
      annoStat <- termStat(object, whichTerms)
      if (ranksOf != orderBy) {
        dim(rr) <- c(length(rr), 1)
        colnames(rr) <- paste("Rank in ", ifelse(is.character(ranksOf), 
                                                 ranksOf, colnames(l)[ranksOf]), sep = "")
        infoMat <- data.frame(infoMat, annoStat, rr, apply(l, 
                                                           2, format.FUN, dig = 2, eps=FALSE), check.names = FALSE, 
                              stringsAsFactors = FALSE)
      }
      else {
        infoMat <- data.frame(infoMat, annoStat, apply(l, 
                                                       2, format.FUN, dig = 2, eps=FALSE), check.names = FALSE, 
                              stringsAsFactors = FALSE)
      }
      rownames(infoMat) <- 1:length(whichTerms)
      return(infoMat)
    }
