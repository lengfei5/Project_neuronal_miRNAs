#library(BSgenome.Dmelanogaster.UCSC.dm3)
library(parallel)
#options(scipen=999)


as.num <- function(..)
  {
    as.numeric(as.character(..))
  }

as.df <- function(d)
  {
    if (is.null(dim(d)))
      {
        d <- t(data.frame(d))
      }
    return(d)
  }

isoCount <- function(x, accuracy=0, aboveIsoCnt = 0)
  {
    n <- as.character(x[1, "gname"])
    type <- paste0(sort(unique(as.character(x[, "type"]))), collapse = "|")
    
    isoCnt <- sapply(split(as.num(x[, "count"]), paste(as.character(x[,"chr"]), as.character(x[,"p5"]), as.character(x[,"strand"]), sep = "#")), sum)
    isoSeed <- sapply(split(rep("NNNNNNN", nrow(x)), paste(as.character(x[,"chr"]), as.character(x[,"p5"]), as.character(x[,"strand"]), sep = "#")), unique)
#    isoSeed <- sapply(split(substr(as.character(x[, "sequence"]), 1, 7), paste(as.character(x[,"chr"]), as.character(x[,"p5"]), as.character(x[,"strand"]), sep = "#")), unique)
    m <- max(isoCnt)/sum(isoCnt)
    f <- isoCnt/sum(isoCnt)
    if ( f == 1 && length(unlist(isoSeed)) > 1) #mismatches in seed
      {
        isoSeed <- unlist(isoSeed)[1]
      }
    f <- cbind(isoSeed=isoSeed, isoCnt=isoCnt, isoFrac=f)
    if (length(f[,1]) != 1) #else data frame converted to vector if only one row
      {
        f <- f[order(-as.numeric(as.character(f[, "isoFrac"]))), ]
      }
    f <- cbind(f, isoform="secondary")
    if (aboveIsoCnt == 0)
      {
        f[f[,"isoFrac"] == m, "isoform"] <- "primary"
      } else {
        f[as.num(f[,"isoCnt"]) > aboveIsoCnt, "isoform"] <- "primary"
      }

    n <- cbind(do.call(rbind, strsplit(rownames(f), '#', fixed = T)), n)
    n <- cbind(n, type)
    colnames(n) <- c("chr", "p5", "strand", "gname", "type")
    f <- cbind(n, f)
    
    if (accuracy != 0)
      {
        isoform.primary <- f[,"isoform"]
        names(isoform.primary) <- rownames(f)
        isoform.primary <- names(isoform.primary[isoform.primary=="primary"])
        p5.acc <- sapply(isoform.primary, function(i)
                         {
                           i <- unlist(strsplit(i, "#", fixed = T))
                           f.close <- as.df(f[(f[,"chr"] == i[1] & f[,"strand"] == i[3]), c("p5", "isoCnt")])
                           p5 <- as.num(i[2])
                           f.close <- as.df(f.close[as.num(f.close[,"p5"]) > p5-accuracy & as.num(f.close[,"p5"]) < p5+accuracy, ])
                           acc <- sum(abs(as.num(f.close[,"p5"]) - p5) * as.num(f.close[,"isoCnt"]))/sum(as.num(f.close[,"isoCnt"]))
                           return(acc)
                         })
      } else {
        p5.acc <- NA
      }

    f <- cbind(as.df(f[isoform.primary,]), p5.acc=unlist(p5.acc))
    rownames(f) <- isoform.primary
    
    return(f)
  }

#mergeEntriesIntoOne (NULL = no merge; negative value = merge all; positive value = merge top N hits
isoCountFeature <- function(x, feat="p3", egalitarian=FALSE, onlyPrimary=TRUE, accuracy=0, composition=FALSE, featuresAsCol = NULL, mergeEntriesIntoOne = NULL, sum = FALSE)
  {
    n <- as.character(x[1, "gname"])

    if (egalitarian == TRUE)
      {
        x[,"count"] <- 1
        x <- unique(x[,c("chr", "p5", "strand", "count", feat)])
      }
    
    featCnt <- split(x, paste(as.character(x[,"chr"]), as.character(x[,"p5"]), as.character(x[,"strand"]), sep = "#"))
    featCnt <- lapply(featCnt, function(y)
           {
             if (composition == TRUE)
               {
                 y <- y[y[,feat] != "",]
                 f.col <- as.character(y[, feat])
                 A <- sum((nchar(f.col) - nchar(gsub("A", "", f.col)))*as.num(y[, "count"]))
                 T <- sum((nchar(f.col) - nchar(gsub("T", "", f.col)))*as.num(y[, "count"]))
                 C <- sum((nchar(f.col) - nchar(gsub("C", "", f.col)))*as.num(y[, "count"]))
                 G <- sum((nchar(f.col) - nchar(gsub("G", "", f.col)))*as.num(y[, "count"]))
                 isoCnt <- c(A, T, C, G)
                 names(isoCnt) <- c("A", "T", "C", "G")
               } else {
                 isoCnt <- sapply(split(as.numeric(as.character(y[, "count"])), list(as.character(y[,feat]))), sum)
               }

             if (is.null(featuresAsCol))
               {
             
                 m <- max(isoCnt)/sum(isoCnt)
                 f <- isoCnt/sum(isoCnt)
                 f <- cbind(isoCnt=isoCnt, isoFrac=f)
                 if (length(f[,1]) != 1) #else data frame converted to vector if only one row
                   {
                     f <- f[order(-as.numeric(as.character(f[, "isoFrac"]))), ]
                   }
                 f <- cbind(f, isoStatus="secondary")
                 f[f[,"isoFrac"] == m, "isoStatus"] <- "primary"
                 f <- cbind(isoValue=rownames(f), f)

                 if (accuracy != 0)
                   {
                     p3.acc <- NULL
                     for (p3 in as.num(f[,"isoValue"]))
                       {
                         f.close <- as.df(f[,c("isoValue", "isoCnt")])
                         f.close <- as.df(f.close[as.num(f.close[,"isoValue"]) > p3-accuracy & as.num(f.close[,"isoValue"]) < p3+accuracy, ])

                         #∑( abs(distance) x reads )/∑(reads)                
                         acc <- sum(abs(as.num(f.close[,"isoValue"]) - p3) * as.num(f.close[,"isoCnt"]))/sum(as.num(f.close[,"isoCnt"]))
                         p3.acc <- c(p3.acc, acc)
                       }
                     f <- cbind(f, isoAccuracy=p3.acc)
                   }
           
             
                 if (onlyPrimary == TRUE)
                   {
                     f <- as.df(f[f[,"isoStatus"] == "primary", ])
                   }
                 colnames(f) <- gsub("iso", feat, colnames(f))

                 if (! is.null(mergeEntriesIntoOne))
                   {
                     if (length(f[,1]) > 1)
                       {
                         if (mergeEntriesIntoOne < 0 || length(f[,1]) < mergeEntriesIntoOne)
                           {
                             mergeEntriesIntoOne <- length(f[,1])
                           }
                         
                         for (i in 1:(mergeEntriesIntoOne-1))
                           {
                             if (i == 1)
                               {
                                 tmp <- paste(f[1,], f[i+1,], sep = "|")
                               } else {
                                 tmp <- paste(tmp, f[i+1,], sep = "|")
                               }
                           }
                         cname <- colnames(f)
                         f <- as.df(tmp)
                         colnames(f) <- cname
                       }
                   }

                 
               } else { #use feature values as separate colnames
                 isoCnt.col <- rep(0, length(featuresAsCol))
                 names(isoCnt.col) <- featuresAsCol
                 cname <- intersect(names(isoCnt.col), names(isoCnt))
                 isoCnt.col[cname] <- isoCnt[cname]
                 names(isoCnt.col) <- paste(feat, names(isoCnt.col), sep = ".")
                 f <- as.df(isoCnt.col)

                 if (typeof(featuresAsCol) == "integer")
                   {
                     isoCnt.n0 <- isoCnt[names(isoCnt) != 0]
                     m <- sum(as.num(names(isoCnt.n0)) * isoCnt.n0) / sum(isoCnt.n0)
                     f <- cbind(f, isoMeanEx0=m)
                     colnames(f) <- gsub("iso", feat, colnames(f))
                   }
               }
             
             return(f)
           })

    f <- NULL
    for (isoform in names(featCnt))
      {
        rows <- length(featCnt[[isoform]][,1])
        iso <- cbind(do.call(rbind, strsplit(isoform, '#', fixed = T)), n)
        iso <- matrix(rep(iso, each=rows), nrow=rows)
        colnames(iso) <- c("chr", "p5", "strand", "gname")
        f <- rbind(f, cbind(iso, featCnt[[isoform]]))
      }

    return(f)
  }

add5p3p <- function(seqCnt)
  {
    #p5 & p3 correspond to the start/end in BED. Hence p5 is 0-based for + and 1-based for - strand. To reconstruct BED positions just add/subtract length => 0-based start & 1-based end.
    tmp <- t(apply(seqCnt, 1, function(l)
                   {
                     if (as.character(l["strand"]) == "-")
                       {
                         p5 <- as.num(l["end"])
                         l["end"] <- as.num(l["start"])
                         l["start"] <- p5
                       }
                     return(c(p5=as.num(l["start"]), p3=as.num(l["end"])))
                   }))
    seqCnt <- cbind(seqCnt, tmp)
    seqCnt$p3 <- format(seqCnt$p3, scientific=FALSE)
    #seqCnt$p5 <- format(seqCnt$p5, scientific=FALSE)
    seqCnt$p5 <- 0
    return(seqCnt)
  }

isoCntPrimary <- function(seqCnt, aboveIsoCnt)
  {
    print("isoCntPrimary ...")
    tmp <- split(seqCnt, as.character(seqCnt$gname))
    
    #most abundant 5'
    tmp.p5 <- lapply(tmp, function(x)
                     {
                       isoCnt <- isoCount(x, accuracy=5, aboveIsoCnt=aboveIsoCnt)
                     })

    tmp.primary <- do.call(rbind, tmp.p5)
    tmp.primary <- tmp.primary[,c("chr", "p5", "strand", "gname", "isoSeed", "type", "isoFrac", "p5.acc"), drop = FALSE]
    return(tmp.primary)
  }

reduceToPrimary <- function(tmp.primary, seqCnt)
  {
    print("reduceToPrimary ...")
    cn <- colnames(seqCnt)
    tmp <- tmp.primary[,c("chr", "p5", "strand", "gname"), drop = FALSE]
    tmp <- merge(seqCnt, tmp, by = c("chr", "p5", "strand", "gname"))
    tmp <- tmp[,cn]
    tmp <- split(tmp, as.character(tmp$gname))

    return(tmp)
  }

isoCntAbundance <- function(tmp.primary, tmp)
  {
    print("isoCntAbundance ...")
    tmp.tail<- lapply(tmp, function(x)
                      {
                        x <- x[, grep("align", colnames(x), invert = TRUE, value = TRUE)] 
                        tailed <- NULL
                        tailed[as.num(x$tailLen) > 0] <- "PM"
                        tailed[as.num(x$tailLen) == 0] <- "GM"
                        x <- cbind(x, align=tailed)
                        isoCnt.tailed <- isoCountFeature(x, feat="align", featuresAsCol=c("GM", "PM"))
                      })
    tmp.tail <- do.call(rbind, tmp.tail)

    if (is.null(tmp.tail))
    {
        return(NULL)
    } else {
        tmp.primary <- merge(tmp.primary, tmp.tail, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
    }
    return(tmp.primary)
  }

isoCntP3Accuracy <- function(tmp.primary, tmp)
  {
    print("isoCntP3Accuracy ...")
    tmp.p3 <- lapply(tmp, function(x)
                     {
                       isoCnt.p3 <- isoCountFeature(x, accuracy = 5, mergeEntriesIntoOne=-1)
                     })

    tmp.p3 <- do.call(rbind, tmp.p3)
    tmp.p3 <- tmp.p3[,c("chr", "p5", "strand", "gname", "p3Accuracy"), drop=FALSE]
    tmp.p3 <- unique(tmp.p3)
    
    tmp.primary <- merge(tmp.primary, tmp.p3, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
    return(tmp.primary)
  }

isoCntTailValue <- function(tmp.primary, tmp)
  {
    print("isoCntTailValue ...")
    tmp.tailSeq <- lapply(tmp, function(x)
                          {
                            x <- x[! is.na(x[,"tail"]), ]
                            x <- x[as.character(x[,"tail"]) != "", ]
                            if (length(x[,1]) == 0)
                              {
                                isoCnt.tailSeq <- NULL
                              } else {
                                isoCnt.tailSeq <- isoCountFeature(x, feat="tail", mergeEntriesIntoOne=3, onlyPrimary = FALSE)
                              }
                          })

    tmp.tailSeq <- do.call(rbind, tmp.tailSeq)


    if (is.null(tmp.tailSeq))
      {
        tmp.primary <- cbind(tmp.primary, tailValue=NA)
      } else {
          tmp.tailSeq <- tmp.tailSeq[,c("chr", "p5", "strand", "gname", "tailValue"), drop=FALSE]
          tmp.primary <- merge(tmp.primary, tmp.tailSeq, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
      }
    return(tmp.primary)
  }

isoCntUMIValue <- function(tmp.primary, tmp)
  {
    print("isoCntUMIValue ...")
    tmp.umiSeq <- lapply(tmp, function(x)
                          {
                            x <- x[! is.na(x[,"UMInum"]), ]
                            x <- x[as.character(x[,"UMInum"]) != "", ]
                            if (length(x[,1]) == 0)
                              {
                                isoCnt.umiSeq <- NULL
                              } else {
                                  isoCnt.umiSeq <- isoCountFeature(x, feat="UMInum", mergeEntriesIntoOne=3, onlyPrimary = FALSE)
                              }
                          })

    tmp.umiSeq <- do.call(rbind, tmp.umiSeq)


    if (is.null(tmp.umiSeq))
      {
        tmp.primary <- cbind(tmp.primary, UMInumValue=NA)
      } else {
          tmp.umiSeq <- tmp.umiSeq[,c("chr", "p5", "strand", "gname", "UMInumValue"), drop=FALSE]
          tmp.primary <- merge(tmp.primary, tmp.umiSeq, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
      }
    return(tmp.primary)
  }


isoCntLast8mer <- function(tmp.primary, tmp)
  {
    print("isoCntLast8mer ...")
    tmp.last8mer <- lapply(tmp, function(x)
                          {
                            last8mer <- apply(x, 1, function(y)
                                              {
                                                y <- gsub(paste(as.character(y["tail"]), "$", sep = ""), "", as.character(y["sequence"]), perl = TRUE)
                                                y <- substr(y, nchar(y)-8, nchar(y))
                                                return(y)
                                              })
                            x <- cbind(x, last8mer=last8mer)
                            isoCnt.last8mer <- isoCountFeature(x, feat="last8mer", mergeEntriesIntoOne=3, onlyPrimary = FALSE)
                          })

    tmp.last8mer <- do.call(rbind, tmp.last8mer)
#    tmp.last8mer <- tmp.last8mer[,c("chr", "p5", "strand", "gname", "last8mer")]

#    if (is.null(tmp.last8mer))
#      {
#        tmp.primary <- cbind(tmp.primary, last8mer=NA)
#      } else {
        tmp.primary <- merge(tmp.primary, tmp.last8mer, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
#      }
    return(tmp.primary)
  }


isoCntComposition <- function(tmp.primary, tmp)
  {
    print("isoCntComposition ...")
    tmp.comp <- lapply(tmp, function(x)
                       {
                         isoCnt.tailPos <- isoCountFeature(x, feat="tail", onlyPrimary = F, composition=TRUE, featuresAsCol=c("A", "T", "C", "G"))
                       })
    tmp.comp <- do.call(rbind, tmp.comp)

    tmp.primary <- merge(tmp.primary, tmp.comp, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
    return(tmp.primary)
  }


isoCntMonotail <- function(tmp.primary, tmp, monoTail = c("A", "T"))
  {
    print("isoCntMonotail ...")

    tmp.monoTail <- lapply(tmp, function(x)
                           {
                             x <- x[as.character(x[,"tail"]) != "", ]                         
                             
                             if (length(x[,1]) == 0)
                               {
                                 isoCnt.monoTail <- NULL
                               } else {
                                 tailed <- rep("None", nrow(x))
                                 for (m in monoTail)
                                   {
                                     tailed[unlist(lapply(strsplit(as.character(x[,"tail"]), "|"), function(nuc) { all(nuc == m) }))] <- m
                                   }

                                 x <- cbind(x, PMmonotail=tailed)
                                 isoCnt.monoTail <- isoCountFeature(x, feat="PMmonotail", featuresAsCol=monoTail)
                               
                               }
                         })
    tmp.monoTail <- do.call(rbind, tmp.monoTail)

    if (is.null(tmp.monoTail))
      {
          tmp.primary <- cbind(tmp.primary, monoTail=NA)
      } else {
          tmp.primary <- merge(tmp.primary, tmp.monoTail, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
      }
    
    return(tmp.primary)
  }

#length + taillength + monotail
isoCntLTM <- function(tmp.primary, tmp, monoTail = c("A", "T"), lenDis)
  {
    print("isoCntLTM ...")

    tmp.LTM <- lapply(tmp, function(x)
                           {
                             tailed <- rep("None", nrow(x))
                             for (m in monoTail)
                               {
                                 tailed[unlist(lapply(strsplit(as.character(x[,"tail"]), "|"), function(nuc) { all(nuc == m) }))] <- m
                               }
                             x <- cbind(x, PMmonotail=tailed)
                             x <- x[as.character(x[,"PMmonotail"]) != "None", ]
                             
                             if (length(x[,1]) == 0)
                               {
                                 isoCnt.LTM <- NULL
                               } else {

                                 x <- cbind(x, lenDisPMeT=(nchar(as.character(x[, "sequence"]))-as.num(x[,"tailLen"])))
                                 x <- cbind(x, LTM=apply(x, 1, function(y) {paste(y[c("PMmonotail", "lenDisPMeT", "tailLen")], collapse="_")}))

                                 LTM <- unlist(lapply(unlist(lapply(monoTail, paste, lenDis, sep="_")), paste, 1:4, sep="_"))
                                 isoCnt.LTM <- isoCountFeature(x, feat="LTM", featuresAsCol=LTM)
                               
                               }
                         })
    tmp.LTM <- do.call(rbind, tmp.LTM)

    tmp.primary <- merge(tmp.primary, tmp.LTM, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
    return(tmp.primary)
  }



isoCntLastNT <- function(tmp.primary, tmp, monoTail = NULL)
  {
    print("isoCntLastNT ...")
    tmp.lastNT <- lapply(tmp, function(x)
                         {
                           x <- x[as.character(x[,"tail"]) != "", ]

                           if (! is.null(monoTail))
                               {
                                 x <- x[unlist(lapply(strsplit(as.character(x[,"tail"]), "|"), function(nuc) { all(nuc == monoTail) })), ]
                               }
                           
                           if (length(x[,1]) == 0)
                             {
                               isoCnt.lastNT <- NULL
                             } else {
                               lastNT <- apply(x, 1, function(y)
                                               {
                                                 y <- gsub(paste(as.character(y["tail"]), "$", sep = ""), "", as.character(y["sequence"]), perl = TRUE)
                                                 y <- substr(y, nchar(y), nchar(y))
                                                 return(y)
                                               })
                               x <- cbind(x, lastNT=lastNT)
                               if (! is.null(monoTail))
                                 {
                                   colnames(x) <- gsub("lastNT", paste0("lastNT", monoTail,"tail"), colnames(x))
                                   isoCnt.lastNT <- isoCountFeature(x, feat=paste0("lastNT", monoTail,"tail"), onlyPrimary = F, composition=TRUE, featuresAsCol=c("A", "T", "C", "G"))
                                 } else {
                                   isoCnt.lastNT <- isoCountFeature(x, feat="lastNT", onlyPrimary = F, composition=TRUE, featuresAsCol=c("A", "T", "C", "G"))
                                 }
                             }
                         })
    tmp.lastNT <- do.call(rbind, tmp.lastNT)

    if (is.null(tmp.lastNT))
      {
        if (! is.null(monoTail))
          {
            cn <- paste0(paste0("lastNT", monoTail,"tail."), c("A", "T", "C", "G"))
          } else {
            cn <- paste0("lastNT.", c("A", "T", "C", "G"))
          }
        tmp.lastNT <- matrix(0, nrow=nrow(tmp.primary), ncol=4)
        colnames(tmp.lastNT) <- cn
        tmp.primary <- cbind(tmp.primary, tmp.lastNT)
      } else {
        tmp.primary <- merge(tmp.primary, tmp.lastNT, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
      }
    return(tmp.primary)
  }

isoCntLastDiNT <- function(tmp.primary, tmp, monoTail = NULL)
  {
    print("isoCntLastDiNT ...")
    dinuc <- NULL
    for (n1 in c("A", "T", "C", "G"))
      {
        for (n2 in c("A", "T", "C", "G"))
          {
            dinuc <- c(dinuc, paste0(n1, n2))
          }
      }
    
    tmp.lastDiNT <- lapply(tmp, function(x)
                         {
                           x <- x[as.character(x[,"tail"]) != "", ]

                           if (! is.null(monoTail))
                               {
                                 x <- x[unlist(lapply(strsplit(as.character(x[,"tail"]), "|"), function(nuc) { all(nuc == monoTail) })), ]
                               }
                           
                           if (length(x[,1]) == 0)
                             {
                               isoCnt.lastDiNT <- NULL
                             } else {
                               lastDiNT <- apply(x, 1, function(y)
                                               {
                                                 y <- gsub(paste(as.character(y["tail"]), "$", sep = ""), "", as.character(y["sequence"]), perl = TRUE)
                                                 y <- substr(y, nchar(y)-1, nchar(y))
                                                 return(y)
                                               })
                               x <- cbind(x, lastDiNT=lastDiNT)
                               if (! is.null(monoTail))
                                 {
                                   colnames(x) <- gsub("lastDiNT", paste0("lastDiNT", monoTail,"tail"), colnames(x))
                                   isoCnt.lastDiNT <- isoCountFeature(x, feat=paste0("lastDiNT", monoTail,"tail"), onlyPrimary = F, composition=FALSE, featuresAsCol=dinuc)
                                 } else {
                                   isoCnt.lastDiNT <- isoCountFeature(x, feat="lastDiNT", onlyPrimary = F, composition=FALSE, featuresAsCol=dinuc)
                                 }
                             }
                         })
    tmp.lastDiNT <- do.call(rbind, tmp.lastDiNT)

    if (is.null(tmp.lastDiNT))
      {
        if (! is.null(monoTail))
          {
            cn <- paste0(paste0("lastDiNT", monoTail,"tail."), dinuc)
          } else {
            cn <- paste0("lastDiNT.", dinuc)
          }
        tmp.lastDiNT <- matrix(0, nrow=nrow(tmp.primary), ncol=16)
        colnames(tmp.lastDiNT) <- cn
        tmp.primary <- cbind(tmp.primary, tmp.lastDiNT)
      } else {
        tmp.primary <- merge(tmp.primary, tmp.lastDiNT, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
      }
    return(tmp.primary)
  }


isoCntLastNTGM <- function(tmp.primary, tmp)
  {
    print("isoCntLastNTGM ...")
    tmp.lastNTGM <- lapply(tmp, function(x)
                           {
                             x <- x[as.character(x[,"tail"]) == "", ]
                             if (length(x[,1]) == 0)
                               {
                                 isoCnt.lastNTGM <- NULL
                               } else {
                                 lastNTGM <- apply(x, 1, function(y)
                                                   {
                                                     y <- as.character(y["sequence"])
                                                     y <- substr(y, nchar(y), nchar(y))
                                                     return(y)
                                                   })
                                 x <- cbind(x, lastNTGM=lastNTGM)
                                 isoCnt.lastNTGM <- isoCountFeature(x, feat="lastNTGM", onlyPrimary = F, composition=TRUE, featuresAsCol=c("A", "T", "C", "G"))
                               }
                           })
    tmp.lastNTGM <- do.call(rbind, tmp.lastNTGM)

    tmp.primary <- merge(tmp.primary, tmp.lastNTGM, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
    return(tmp.primary)
  }

isoCntTailLen <- function(tmp.primary, tmp)
  {
    print("isoCntTailLen ...")
    tmp.tailLen <- lapply(tmp, function(x)
                          {
                            isoCnt.tailLen <- isoCountFeature(x, feat="tailLen", featuresAsCol=1:6)
                          })
    tmp.tailLen <- do.call(rbind, tmp.tailLen)

    tmp.primary <- merge(tmp.primary, tmp.tailLen, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
    return(tmp.primary)
  }


isoCntTailLenMonotail <- function(tmp.primary, tmp, monoTail=c("A", "T"))
{
        print("isoCntTailLenMonotail ...")
        tmp.TM <- lapply(tmp, function(x)
                           {
                             tailed <- rep("None", nrow(x))
                             for (m in monoTail)
                               {
                                 tailed[unlist(lapply(strsplit(as.character(x[,"tail"]), "|"), function(nuc) { all(nuc == m) }))] <- m
                               }
                             x <- cbind(x, PMmonotail=tailed)
                             x <- x[as.character(x[,"PMmonotail"]) != "None", ]
                             
                             if (length(x[,1]) == 0)
                               {
                                 isoCnt.LTM <- NULL
                               } else {

                             x <- cbind(x, TM=apply(x, 1, function(y) {paste(y[c("PMmonotail", "tailLen")], collapse="_")}))
                             TM <- unlist(lapply(monoTail, paste, 1:10, sep="_"))
                             isoCnt.TM <- isoCountFeature(x, feat="TM", featuresAsCol=TM)
                               
                               }
                         })
    tmp.TM <- do.call(rbind, tmp.TM)

    tmp.primary <- merge(tmp.primary, tmp.TM, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
    return(tmp.primary)

  }


isoCntLenDisGM <- function(tmp.primary, tmp, lenDis)
  {
    print("isoCntLenDisGM ...")
    tmp.lenDis <- lapply(tmp, function(x)
                         {
                           x <- x[as.num(x[,"tailLen"]) == 0, ]
                           if (length(x[,1]) == 0)
                             {
                               isoCnt.lenDis <- NULL
                             } else {
                               x <- cbind(x, lenDisGM=nchar(as.character(x[, "sequence"])))
                               isoCnt.lenDis <- isoCountFeature(x, feat="lenDisGM", featuresAsCol=lenDis)
                             }
                         })
    tmp.lenDis <- do.call(rbind, tmp.lenDis)

    tmp.primary <- merge(tmp.primary, tmp.lenDis, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
    return(tmp.primary)
  }

isoCntLenDisPM <- function(tmp.primary, tmp, lenDis)
  {
    print("isoCntLenDisPM ...")
    tmp.lenDisPM <- lapply(tmp, function(x)
                           {
                             x <- x[as.num(x[,"tailLen"]) > 0, ]
                             if (length(x[,1]) == 0)
                               {
                                 isoCnt.lenDis <- NULL
                               } else {
                                 x <- cbind(x, lenDisPM=nchar(as.character(x[, "sequence"])))
                                 isoCnt.lenDis <- isoCountFeature(x, feat="lenDisPM", featuresAsCol=lenDis)
                               }
                           })
    tmp.lenDisPM <- do.call(rbind, tmp.lenDisPM)
    if (is.null(tmp.lenDisPM))
      {
        tmp.lenDisPM <- matrix(0, nrow=nrow(tmp.primary), ncol=length(lenDis))
        colnames(tmp.lenDisPM) <- paste0("lenDisPM.", lenDis)
        tmp.primary <- cbind(tmp.primary, tmp.lenDisPM)
      } else {
        tmp.primary <- merge(tmp.primary, tmp.lenDisPM, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
      }
    return(tmp.primary)
  }

isoCntLenDisPMeT <- function(tmp.primary, tmp, lenDis)
  {
    print("isoCntLenDisPMeT ...")
    tmp.lenDisPMeT <- lapply(tmp, function(x)
                             {
                               x <- x[as.num(x[,"tailLen"]) > 0, ]
                               if (length(x[,1]) == 0)
                                 {
                                   isoCnt.lenDis <- NULL
                                 } else {
                                   x <- cbind(x, lenDisPMeT=(nchar(as.character(x[, "sequence"]))-as.num(x[,"tailLen"])))
                                   isoCnt.lenDis <- isoCountFeature(x, feat="lenDisPMeT", featuresAsCol=lenDis)
                                 }
                             })
    tmp.lenDisPMeT <- do.call(rbind, tmp.lenDisPMeT)

    if (is.null(tmp.lenDisPMeT))
      {
        tmp.lenDisPMeT <- matrix(0, nrow=nrow(tmp.primary), ncol=length(lenDis))
        colnames(tmp.lenDisPMeT) <- paste0("lenDisPMeT.", lenDis)
        tmp.primary <- cbind(tmp.primary, tmp.lenDisPMeT)
      } else {
        tmp.primary <- merge(tmp.primary, tmp.lenDisPMeT, all.x = TRUE, by = c("chr", "p5", "strand", "gname"))
      }
    return(tmp.primary)
  }

miRmiRStar <- function(tmp.primary, seqCnt)
  {
    tmp.p5 <- seqCnt[grep("-5p", as.character(seqCnt$gname)), c("gname", "count")]
    tmp.p5 <- sapply(split(as.num(tmp.p5$count), as.character(tmp.p5$gname)), sum)
    tmp.p5 <- cbind(gname = gsub("-5p", "", names(tmp.p5)), rpm.p5 = tmp.p5)
    tmp.p3 <- seqCnt[grep("-3p", as.character(seqCnt$gname)), c("gname", "count")]
    tmp.p3 <- sapply(split(as.num(tmp.p3$count), as.character(tmp.p3$gname)), sum)
    tmp.p3 <- cbind(gname = gsub("-3p", "", names(tmp.p3)), rpm.p3 = tmp.p3)

    if (nrow(tmp.p5) == 0 && nrow(tmp.p3) == 0)
      {
        return(cbind(tmp.primary, miR=NA))
      }
    
    tmp.mir <- as.matrix(merge(tmp.p5, tmp.p3, by = "gname", all = TRUE))
    tmp.mir[is.na(tmp.mir)] <- 0
    tmp.mir <- cbind(tmp.mir, p5 = (as.num(tmp.mir[,"rpm.p5"]) >= as.num(tmp.mir[,"rpm.p3"])))
    tmp.mir <- cbind(tmp.mir, p3 = (as.num(tmp.mir[,"rpm.p5"]) < as.num(tmp.mir[,"rpm.p3"])))

    tmp.p5 <- data.frame(tmp.mir[,c("gname", "p5")])
    tmp.p5$gname <- paste0(tmp.p5$gname, "-5p")
    colnames(tmp.p5) <- c("gname", "miR")
    tmp.p3 <- data.frame(tmp.mir[,c("gname", "p3")])
    tmp.p3$gname <- paste0(tmp.p3$gname, "-3p")
    colnames(tmp.p3) <- c("gname", "miR")
    tmp.mir <- rbind(tmp.p5, tmp.p3)

    tmp.primary <- merge(tmp.primary, tmp.mir, by = "gname", all.x = TRUE)
    return(tmp.primary)
  }
    
isoform5pAllInOne <- function(seqCnt, lenDis, aboveIsoCnt, last8mer = FALSE)
  {
    #add 5p & 3p
    seqCnt <- add5p3p(seqCnt)
    
    #get primary
    tmp.primary <- isoCntPrimary(seqCnt, aboveIsoCnt)

    #miR & miR-star
    tmp.primary <- miRmiRStar(tmp.primary, seqCnt)
    
    #Reduce data set to primary
    tmp <- reduceToPrimary(tmp.primary, seqCnt)

    tmp.primary <- isoCntAbundance(tmp.primary,tmp)                 #is aligned GM or PM
    if (is.null(tmp.primary)) { return(NULL)}

      tmp.primary$align.total <- apply(tmp.primary[,c("align.GM", "align.PM")], 1, function(x) {sum(as.numeric(as.character(x)))})
      
    tmp.primary <- tmp.primary[,c("chr", "p5", "strand", "gname", "isoSeed", "type", "isoFrac", "miR", "align.GM", "align.PM", "align.total", "p5.acc")]     #5' heterogeneity
    tmp.primary <- isoCntP3Accuracy(tmp.primary, tmp)               #3' heterogeneity
      tmp.primary <- isoCntTailValue(tmp.primary, tmp)                #top3 sequences
      tmp.primary <- isoCntUMIValue(tmp.primary, tmp)                #values of top3 most abundant seqs
    if ( last8mer == TRUE)
      {
        tmp.primary <- isoCntLast8mer(tmp.primary, tmp)            #top3 aligned ends
      }
    tmp.primary <- isoCntMonotail(tmp.primary, tmp)                 #PMmonotail
    tmp.primary <- isoCntComposition(tmp.primary, tmp)              #composition
    tmp.primary <- isoCntLastNT(tmp.primary, tmp)                   #last NT before tail
    tmp.primary <- isoCntLastNT(tmp.primary, tmp, monoTail = "A")   #last NT before tail
    tmp.primary <- isoCntLastNT(tmp.primary, tmp, monoTail = "C")   #last NT before tail
    tmp.primary <- isoCntLastNT(tmp.primary, tmp, monoTail = "G")   #last NT before tail
    tmp.primary <- isoCntLastNT(tmp.primary, tmp, monoTail = "T")   #last NT before tail
    tmp.primary <- isoCntLastNTGM(tmp.primary, tmp)                 #last NT of GM
    tmp.primary <- isoCntTailLen(tmp.primary, tmp)                  #tailLen
    tmp.primary <- isoCntLenDisGM(tmp.primary, tmp, lenDis)         #lenDist GM
    tmp.primary <- isoCntLenDisPM(tmp.primary, tmp, lenDis)         #lenDist PM
    tmp.primary <- isoCntLenDisPMeT(tmp.primary, tmp, lenDis)       #lenDist PM excl Tail

    return(tmp.primary)
  }


isoform5pAllInOne.modular <- function(seqCnt, lenDis, aboveIsoCnt, last8mer = FALSE, module=module)
  {
    #add 5p & 3p
    seqCnt <- add5p3p(seqCnt)
    
    #get primary
    tmp.primary <- isoCntPrimary(seqCnt, aboveIsoCnt)

    #miR & miR-star
    tmp.primary <- miRmiRStar(tmp.primary, seqCnt)
    
    #Reduce data set to primary
    tmp <- reduceToPrimary(tmp.primary, seqCnt)

    tmp.primary <- isoCntAbundance(tmp.primary,tmp)                 #is aligned GM or PM
      tmp.primary$align.total <- apply(tmp.primary[,c("align.GM", "align.PM")], 1, function(x) {sum(as.numeric(as.character(x)))})
      tmp.primary <- tmp.primary[,c("chr", "p5", "strand", "gname", "isoSeed", "type", "isoFrac", "miR", "align.GM", "align.PM", "align.total")]
     
    if (any(module == "isoCntLastDiNT"))
      {
        tmp.primary <- isoCntLastDiNT(tmp.primary, tmp)                 #last DiNT before tail
      }

    if (any(module == "isoCntMonotail"))
      {
        tmp.primary <- isoCntMonotail(tmp.primary, tmp)                 #PMmonotail
      }
      
    if (any(module == "isoCntLTM"))
      {
        tmp.primary <- isoCntLTM(tmp.primary, tmp, lenDis=lenDis)                 #length + taillength + monotail
      }

    if (any(module == "isoCntTailLenMonotail"))
    {
        tmp.primary <- isoCntTailLenMonotail(tmp.primary, tmp)
    }
      
    return(tmp.primary)
  }


getSeqCnt <- function(seqCnt.file, file.suffixPrimary, file.suffixSecondary, normalizeCnt, tailFilter, includePrimary, includeType, excludeId, excludeType, spliceSite = FALSE, annotWithAS = TRUE, minLength=minLength, countColumn = "count")
  {

#    print(paste("seqCnt.file", seqCnt.file))
#    print(paste("file.suffixPrimary", file.suffixPrimary))
#    print(paste("file.suffixSecondary", file.suffixSecondary))
#    print(paste("normalizeCnt", normalizeCnt))
#    print(paste("tailFilter", tailFilter))
#    print(paste("includePrimary", includePrimary))
#    print(paste("includeType", includeType))
#    print(paste("excludeId", excludeId))
#    print(paste("excludeType", excludeType))
#    print(paste("spliceSite", spliceSite))
#    print(paste("annotWithAS", annotWithAS))
    
    if (annotWithAS)
      {
        includePrimary <- paste0(includePrimary, ";")
        includeType <- paste0(includeType, ";")
      }

    if (is.null(file.suffixSecondary) || file.suffixPrimary != file.suffixSecondary)
      {
          seqCnt <- read.delim(seqCnt.file)
        seqCnt <- seqCnt[grep(includePrimary, seqCnt$type), ]
        seqCnt[, "count"] <- seqCnt[,countColumn]
      } else {
          seqCnt <- read.delim(pipe(paste0("zcat ", seqCnt.file, " | head -n 1")))
          seqCnt[, "count"] <- seqCnt[,countColumn]
      }
    
    if (!(is.null(file.suffixSecondary)))
      {
        file.sec <- sub(file.suffixPrimary, file.suffixSecondary, seqCnt.file)
        #seqCnt.sec <- read.delim(file.sec)
        #prefilter already during loading to avoid massive memory
        if (is.null(excludeType))
          {
              seqCnt.sec <- read.delim(pipe(paste0("zcat ", file.sec, " | egrep '", paste(includeType, collapse = "|"), "'")), header = FALSE)
          } else {
            if (annotWithAS)
              {
                excludeType <- paste0(excludeType, ";")
              }
            seqCnt.sec <- read.delim(pipe(paste0("zcat ", file.sec, " | egrep '", paste(includeType, collapse = "|"), "'" , " | egrep -v '", paste(excludeType, collapse = "|"), "'")), header = FALSE)
          }
        colnames(seqCnt.sec)<- as.character(unlist(read.delim(pipe(paste0("zcat ", file.sec, " | head -n 1")), header = FALSE)))
        seqCnt.sec[, "count"] <- seqCnt.sec[,countColumn]
        colnames(seqCnt.sec) <- sub("repeat", "repeat.", colnames(seqCnt.sec), fixed = TRUE)
        seqCnt.sec <- seqCnt.sec[grep(includePrimary, seqCnt.sec$type, invert=TRUE), ]

        for (t in includeType)
          {
            seqCnt.tmp <- seqCnt.sec[grep(t, seqCnt.sec$type), ]
            if ( t == "intron" )
              {
                seqCnt.tmp$gname <- sapply(strsplit(as.character(seqCnt.tmp$gname), "|", fixed = TRUE), function(gname)
                                           {
                                             gname <- gname[grep("-in$", gname)]
                                             gname <- paste(gname, collapse = "|")
                                             return(gname)
                                           })
                if (spliceSite == TRUE)
                  {
                    last2NT <- apply(seqCnt.tmp, 1, function(y)
                                     {
                                       y <- gsub(paste(as.character(y["tail"]), "$", sep = ""), "", as.character(y["sequence"]), perl = TRUE)
                                       y <- substr(y, nchar(y)-1, nchar(y))
                                       return(y)
                                     })
                    seqCnt.tmp <- subset(seqCnt.tmp, last2NT == "AG")                                     
                  }
              } else if (t == "esiRNA") {
                seqCnt.tmp$gname <- sapply(strsplit(as.character(seqCnt.tmp$gname), "|", fixed = TRUE), function(gname)
                                           {
                                             gname <- gname[grep("esiRNA", gname)]
                                             gname <- paste(gname, collapse = "|")
                                             return(gname)
                                       })
                seqCnt.tmp$gname[seqCnt.tmp$strand == "+"] <- gsub("esiRNA", "esiRNAp", as.character(seqCnt.tmp$gname[seqCnt.tmp$strand == "+"]))
                seqCnt.tmp$gname[seqCnt.tmp$strand == "-"] <- gsub("esiRNA", "esiRNAm", as.character(seqCnt.tmp$gname[seqCnt.tmp$strand == "-"]))
              }
            
            seqCnt <- rbind(seqCnt, seqCnt.tmp)
          }
        seqCnt <- unique(seqCnt)
      }

    ###############
    #hierarchy type
    hierarchy <- c("mapped_reads", "mito", "rRNA", "rRNA_AS", "tRNA", "tRNA_AS", "good_reads", "pre_miRNA_5p3p", "pre_miRNA_5p3p_AS", "pre_miRNA", "pre_miRNA_AS", "transposable_element", "transposable_element_AS", "__no_feature", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS", "pseudogene", "pseudogene_AS", "ncRNA", "ncRNA_AS", "snoRNA", "snoRNA_AS", "snRNA", "snRNA_AS", "mRNA", "mRNA_AS", "intron", "intron_AS")

    levels(seqCnt$type) <- c(levels(seqCnt$type), setdiff(hierarchy, levels(seqCnt$type)))
    levels(seqCnt$type) <- c(levels(seqCnt$type), setdiff(c("mito;"), levels(seqCnt$type)))
    if (any(seqCnt$chr == "dmel_mitochondrion_genome"))
      {
        seqCnt[seqCnt$chr == "dmel_mitochondrion_genome", "type"] <- "mito;"
      } else if (any(seqCnt$chr == "MtDNA"))
        {
          seqCnt[seqCnt$chr == "MtDNA", "type"] <- "mito;"
        }
    
    for (t in hierarchy)
      {
        index <- grep(paste0(t, ";"), as.character(seqCnt$type))
        if (length(index) != 0)
          {
            seqCnt[index, "type"] <- t
          }
      }
    seqCnt$type <- droplevels(seqCnt$type)
    #hierarchy type - END
    #####################
    
    if (minLength != -1)
      {
       seqCnt <- seqCnt[nchar(as.character(seqCnt$sequence)) > minLength, ]
      } 
    
    if (tailFilter < 1)
      {
        seqCnt <- seqCnt[nchar(as.character(seqCnt$tail))/nchar(as.character(seqCnt$sequence)) <= tailFilter, ]
      }

    print(head(seqCnt))
    
    if (normalizeCnt == TRUE)
      {
        print(paste("Total (File =", seqCnt.file, ", TF =", tailFilter, "): ", sum(as.num(seqCnt[,"count"]))))
        system(paste("echo '", "Total (File =", seqCnt.file, ", TF =", tailFilter, "): ", sum(as.num(seqCnt[,"count"])), "' >> cntStat.miRNA_used_for_normalization.txt"))
        total.pm <- sum(as.num(seqCnt[,"count"]))/1000000

        print(paste("Within feature normalization:", total.pm))
        
        seqCnt[,"count"] <- as.num(seqCnt[,"count"])/total.pm
      } else if (normalizeCnt == "all")
        {
            file.allCnt <- sub(file.suffixPrimary, "bam.seqCnt.txt.gz", seqCnt.file)
            if (countColumn == "count")
            {
                totalAligned <- sum(read.delim(pipe(paste0("zcat ", file.allCnt, " | cut -f 12")))[,"count"])
            } else if (countColumn == "UMInum")
            {
                totalAligned <- sum(read.delim(pipe(paste0("zcat ", file.allCnt, " | cut -f 13")))[,"UMInum"])
            } else if (countColumn == "UMIfr")
            {
                totalAligned <- sum(read.delim(pipe(paste0("zcat ", file.allCnt, " | cut -f 14")))[,"UMIfr"])
            }                
                
          print(paste("Total (File =", file.allCnt, "): ", totalAligned))
          totalAligned <- totalAligned/1000000

          print(paste("All reads normalization:", totalAligned))
          
          seqCnt[,"count"] <- as.num(seqCnt[,"count"])/totalAligned
        } else if (normalizeCnt != FALSE)
        {
            if (countColumn == "count")
            {
                hierarchy <- read.delim("cnt.typeHierarchy.txt")
            } else {
                hierarchy <- read.delim(paste0("cnt.typeHierarchy.", countColumn, ".txt"))
            }
            hierarchy <- hierarchy[as.character(hierarchy$type) %in% normalizeCnt, ]
            rownames(hierarchy) <- as.character(hierarchy$type)
            hierarchy <- hierarchy[,grep("count", colnames(hierarchy)),drop=FALSE]
            total <- colSums(hierarchy, na.rm = TRUE)

            sampleID <- seqCnt.file
            sampleID <- gsub(".*#(\\d+)_.*", "\\1", sampleID,  perl = TRUE)
            sampleID <- gsub("(\\d+)_.*", "\\1", sampleID,  perl = TRUE)
            sampleID <- gsub(".trimmed.*", "", sampleID,  perl = TRUE)

            total <- total[paste0("count.", sampleID)]
            total <- total/1000000

            print(paste("Standard normalization:", total))
            
            seqCnt[,"count"] <- as.num(seqCnt[,"count"])/total
          }
    
    return(seqCnt)
  }


analyseSeqCnt <- function(seqCnt, seqCnt.file, lenDis, aboveIsoCnt, last8mer)
  {
      d <- isoform5pAllInOne(seqCnt=seqCnt, lenDis = lenDis, aboveIsoCnt=aboveIsoCnt, last8mer=last8mer)

      cname <- colnames(d)
      sampleID <- seqCnt.file
      sampleID <- gsub(".*#(\\d+)_.*", "\\1", sampleID,  perl = TRUE)
      sampleID <- gsub("(\\d+)_.*", "\\1", sampleID,  perl = TRUE)
      sampleID <- gsub(".trimmed.*", "", sampleID,  perl = TRUE)
      cname[6:length(d)] <- gsub("$", paste(".", sampleID, sep = ""), cname[6:length(d)], perl = TRUE)
      colnames(d) <- cname

      return(d)
  }

########################################################
########################################################
# Section: here is the function to use for piRNA normalization
########################################################
########################################################

analysisAll <- function(file.suffixPrimary, output.middleName, file.suffixSecondary = NULL, 
                        normalizeCnt = c("pre_miRNA_5p3p", "pre_miRNA", "transposable_element", "transposable_element_AS", 
                                         "__no_feature", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS",
                                         "mRNA", "mRNA_AS", "intron", "intron_AS"), 
                        tailFilter = 0.12, includePrimary = "pre_miRNA_5p3p", includeType = NULL, excludeId = NULL, excludeType = NULL, lenDis = 18:30, aboveIsoCnt = 0, processes = 1, spliceSite = FALSE, last8mer = FALSE, annotWithAS = TRUE, minLength=-1, countColumn=countColumn)
{
  seqCnt.all <- list.files(".", paste0(file.suffixPrimary, "$"))
  
  if (processes == 1)
    {
      miRNA.list <- lapply(seqCnt.all, function(seqCnt.file)
                           {
                             show(seqCnt.file)
                             seqCnt <- getSeqCnt(seqCnt.file=seqCnt.file, file.suffixPrimary=file.suffixPrimary, file.suffixSecondary=file.suffixSecondary, normalizeCnt=normalizeCnt, tailFilter=tailFilter, includePrimary=includePrimary, includeType=includeType, excludeId=excludeId, excludeType=excludeType, spliceSite=spliceSite, annotWithAS=annotWithAS, minLength=minLength, countColumn=countColumn)
                             return(analyseSeqCnt(seqCnt, seqCnt.file, lenDis, aboveIsoCnt, last8mer))
                           })
    } else {
      cl <- makeCluster(rep("localhost", (processes-1)),
                        type = "SOCK"
                        )

      #export all functions to workers
      ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
      clusterExport(cl, ex)

      miRNA.list <- parLapply(cl, seqCnt.all, function(seqCnt.file)
                              {
                                seqCnt <- getSeqCnt(seqCnt.file=seqCnt.file, file.suffixPrimary=file.suffixPrimary, file.suffixSecondary=file.suffixSecondary, normalizeCnt=normalizeCnt, tailFilter=tailFilter, includePrimary=includePrimary, includeType=includeType, excludeId=excludeId, excludeType=excludeType, spliceSite=spliceSite, annotWithAS=annotWithAS, minLength=minLength, countColumn=countColumn)
                                return(analyseSeqCnt(seqCnt, seqCnt.file, lenDis, aboveIsoCnt, last8mer))
                              })
      stopCluster(cl)
    }

  miRNA.list <- Filter(Negate(is.null), miRNA.list)
  miRNA <- Reduce(function(x,y) merge(x, y, all = TRUE, by = c("chr", "p5", "strand", "gname", "isoSeed")), miRNA.list)
  for ( i in 1:ncol(miRNA))
    {
      if (is.list(miRNA[,i]))
        {
          miRNA[,i] <- vapply(miRNA[,i], paste, collapse = ", ", character(1L)) #if one column is a list ...
        }
    }
  
  if (tailFilter == 1)
    {
      if (aboveIsoCnt > 0)
        {
          write.table(miRNA, file = paste0("allProperties.", output.middleName,".tf1.ic",aboveIsoCnt,".",countColumn,".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        } else {
          write.table(miRNA, file = paste0("allProperties.", output.middleName,".tf1.",countColumn, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        }
    } else {
      if (aboveIsoCnt > 0)
        {
          write.table(miRNA, file = paste0("allProperties.", output.middleName,".tf", tailFilter, ".ic",aboveIsoCnt, ".",countColumn, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        } else {
          write.table(miRNA, file = paste0("allProperties.", output.middleName,".tf", tailFilter, ".",countColumn, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        }
    }
}


##
analyseSeqCnt.modular <- function(seqCnt, seqCnt.file, lenDis, aboveIsoCnt, last8mer, module)
  {
    if (is.null(module))
      {
        d <- isoform5pAllInOne(seqCnt=seqCnt, lenDis = lenDis, aboveIsoCnt=aboveIsoCnt, last8mer=last8mer)
      } else {
        d <- isoform5pAllInOne.modular(seqCnt=seqCnt, lenDis = lenDis, aboveIsoCnt=aboveIsoCnt, last8mer=last8mer, module = module)
      }

      cname <- colnames(d)
      sampleID <- seqCnt.file
      sampleID <- gsub(".*#(\\d+)_.*", "\\1", sampleID,  perl = TRUE)
      sampleID <- gsub("(\\d+)_.*", "\\1", sampleID,  perl = TRUE)
      sampleID <- gsub(".trimmed.*", "", sampleID,  perl = TRUE)
      cname[6:length(d)] <- gsub("$", paste(".", sampleID, sep = ""), cname[6:length(d)], perl = TRUE)
      colnames(d) <- cname

      return(d)
  }

analysisAll.modular <- function(file.suffixPrimary, output.middleName, file.suffixSecondary = NULL, normalizeCnt = c("pre_miRNA_5p3p", "pre_miRNA", "transposable_element", "transposable_element_AS", "__no_feature", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS", "mRNA", "mRNA_AS", "intron", "intron_AS"), tailFilter = 0.12, includePrimary = "pre_miRNA_5p3p", includeType = NULL, excludeId = NULL, excludeType = NULL, lenDis = 18:30, aboveIsoCnt = 0, processes = 1, spliceSite = FALSE, last8mer = FALSE, annotWithAS=TRUE, module = NULL, minLength=-1, countColumn=countColumn)
{
  seqCnt.all <- list.files(".", paste0(file.suffixPrimary, "$"))
  
  if (processes == 1)
    {
      miRNA.list <- lapply(seqCnt.all, function(seqCnt.file)
                           {
                             show(seqCnt.file)
                             seqCnt <- getSeqCnt(seqCnt.file=seqCnt.file, file.suffixPrimary=file.suffixPrimary, file.suffixSecondary=file.suffixSecondary, normalizeCnt=normalizeCnt, tailFilter=tailFilter, includePrimary=includePrimary, includeType=includeType, excludeId=excludeId, excludeType=excludeType, spliceSite=spliceSite, annotWithAS=annotWithAS, minLength=minLength, countColumn=countColumn)
                             return(analyseSeqCnt.modular(seqCnt, seqCnt.file, lenDis, aboveIsoCnt, last8mer, module=module))
                           })
    } else {
      cl <- makeCluster(rep("localhost", (processes-1)),
                        type = "SOCK"
                        )

      #export all functions to workers
      ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
      clusterExport(cl, ex)

      miRNA.list <- parLapply(cl, seqCnt.all, function(seqCnt.file)
                              {
                                seqCnt <- getSeqCnt(seqCnt.file=seqCnt.file, file.suffixPrimary=file.suffixPrimary, file.suffixSecondary=file.suffixSecondary, normalizeCnt=normalizeCnt, tailFilter=tailFilter, includePrimary=includePrimary, includeType=includeType, excludeId=excludeId, excludeType=excludeType, spliceSite=spliceSite, annotWithAS=annotWithAS, minLength=minLength, countColumn=countColumn)
                                return(analyseSeqCnt.modular(seqCnt, seqCnt.file, lenDis, aboveIsoCnt, last8mer, module=module))
                              })
      stopCluster(cl)
    }

  miRNA <- Reduce(function(x,y) merge(x, y, all = TRUE, by = c("chr", "p5", "strand", "gname", "isoSeed")), miRNA.list)
  for ( i in 1:ncol(miRNA))
    {
      if (is.list(miRNA[,i]))
        {
          miRNA[,i] <- vapply(miRNA[,i], paste, collapse = ", ", character(1L)) #if one column is a list ...
        }
    }
  
  if (tailFilter == 1)
    {
      if (aboveIsoCnt > 0)
        {
          write.table(miRNA, file = paste0("allProperties.", paste(module, collapse="."), ".", output.middleName,".tf1.ic",aboveIsoCnt,".",countColumn,".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        } else {
          write.table(miRNA, file = paste0("allProperties.", paste(module, collapse="."), ".", output.middleName,".tf1.",countColumn,".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        }
    } else {
      if (aboveIsoCnt > 0)
        {
          write.table(miRNA, file = paste0("allProperties.", paste(module, collapse="."), ".", output.middleName,".tf", tailFilter, ".ic",aboveIsoCnt,".",countColumn, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        } else {
          write.table(miRNA, file = paste0("allProperties.", paste(module, collapse="."), ".", output.middleName,".tf", tailFilter, ".",countColumn, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        }
    }
}

###############################
# the function used in the original pipeline
# the last step to quantify the read counts and normalized read counts (cpm)
###############################
analysisAllmiRNA <- function(file.suffix = "bam.seqCntMirna.txt.gz", tailFilter = 0.12, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffix, output.middleName = "miRNA", tailFilter = tailFilter, processes = processes, annotWithAS=annotWithAS)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("~/shared/Ameres.GRP/Ivica.Sowemimo/miRNA/head/")
#R: setwd("~/shared/Ameres.GRP/Ivica.Sowemimo/miRNA/muscle/")
#R: setwd("~/shared/Ameres.GRP/Ivica.Sowemimo/miRNA/ovary/")
#R: analysisAllWithinMiRNA(processes = 6)

analysisAllWithinMiRNA <- function(file.suffix = "bam.seqCntMirna.txt.gz", tailFilter = 0.12, processes = 1, annotWithAS=TRUE, normalizeCnt = TRUE)
{
  analysisAll(file.suffixPrimary = file.suffix, output.middleName = "WithinMiRNA", tailFilter = tailFilter, processes = processes, annotWithAS=annotWithAS, normalizeCnt=normalizeCnt)
}


#in mirBase the pre_miRNA_5p3p is currently pre_miRNA... & normalizeCnt = TRUE (within miRNA)
#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("~/shared/Cochella.GRP/Chiara.Alberti/miRNA_R3289/")
#R: setwd("~/clustertmp/Chiara.Alberti/")
#R: analysisAllmiRNABase(processes = 6, aboveIsoCnt = 1)
#R: setwd("/clustertmp/bioinfo/thomas/Ivica.Sowemimo")
#R: analysisAllmiRNABase(file.suffix = "tailor.bam.seqCntMirna.txt.gz", processes = 6, aboveIsoCnt = 1)
#R: system("mv allProperties.miRNA.tf0.12.ic1.txt allProperties.miRNA.tf0.12.ic1.dm.txt")
#R: analysisAllmiRNABase(file.suffix = "mm10.bam.seqCntMirna.txt.gz", processes = 6, aboveIsoCnt = 1)
#R: system("mv allProperties.miRNA.tf0.12.ic1.txt allProperties.miRNA.tf0.12.ic1.mm.txt")
#R: analysisAllmiRNABase(file.suffix = "tailor.bam.seqCntMirna.txt.gz", processes = 6, aboveIsoCnt = 0)
#R: system("mv allProperties.miRNA.tf0.12.txt allProperties.miRNA.tf0.12.dm.txt")
#R: analysisAllmiRNABase(file.suffix = "mm10.bam.seqCntMirna.txt.gz", processes = 6, aboveIsoCnt = 0)
#R: system("mv allProperties.miRNA.tf0.12.txt allProperties.miRNA.tf0.12.mm.txt")
#R: analysisAllmiRNABase(file.suffix = "PlusDm3.seqCntMirna.txt.gz", processes = 6, aboveIsoCnt = 1)
#R: analysisAllmiRNABase(file.suffix = "PlusDm3.seqCntMirna.txt.gz", processes = 6, aboveIsoCnt = 0)
analysisAllmiRNABase <- function(file.suffix = "bam.seqCntMirna.txt.gz", tailFilter = 0.12, processes = 1, 
                                 annotWithAS=TRUE, normalizeCnt = TRUE, aboveIsoCnt = 0)
{
  analysisAll(file.suffixPrimary = file.suffix, output.middleName = "miRNA", includePrimary = "pre_miRNA", 
              tailFilter = tailFilter, processes = processes, annotWithAS=annotWithAS, normalizeCnt = normalizeCnt, aboveIsoCnt = aboveIsoCnt)
}


#source("~/development/smallRNA/functions.miRNASummarize.noFixation.R");
#setwd("/clustertmp/bioinfo/thomas/Chiara.Alberti_R3859")
#setwd("/clustertmp/bioinfo/thomas/miRNA_R4042")
#analysisAllmiRNABase.modular(processes = 6, normalizeCnt=FALSE)
#system("mv allProperties.total.miRNA.tf0.12.count.txt allProperties.total.miRNA.tf0.12.countRaw.txt")
#analysisAllmiRNABase.modular(processes = 6)
#analysisAllmiRNABase.modular(processes = 6, normalizeCnt=FALSE, countColumn="UMIfr")
#system("mv allProperties.total.miRNA.tf0.12.UMIfr.txt allProperties.total.miRNA.tf0.12.UMIfrRaw.txt")
#analysisAllmiRNABase.modular(processes = 6, countColumn="UMIfr")
#analysisAllmiRNABase.modular(processes = 6, normalizeCnt=FALSE, countColumn="UMInum")
#system("mv allProperties.total.miRNA.tf0.12.UMInum.txt allProperties.total.miRNA.tf0.12.UMInumRaw.txt")
#analysisAllmiRNABase.modular(processes = 6, countColumn="UMInum")

analysisAllmiRNABase.modular <- function(file.suffix = "bam.seqCntMirna.txt.gz", tailFilter = 0.12, processes = 1, annotWithAS=TRUE, normalizeCnt = TRUE, aboveIsoCnt = 0, countColumn = "count")
{
  analysisAll.modular(file.suffixPrimary = file.suffix, output.middleName = "miRNA", includePrimary = "pre_miRNA", tailFilter = tailFilter, processes = processes, annotWithAS=annotWithAS, normalizeCnt = normalizeCnt, aboveIsoCnt = aboveIsoCnt, countColumn=countColumn, module="total")
}





#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("/clustertmp/bioinfo/thomas/Exp5_preMiRNA/")
analysisAllPreMiRNA <- function(file.suffix="bam.seqCntPreMirna.txt.gz", tailFilter = 0.05, processes = 1, annotWithAS=TRUE, normalizeCnt = TRUE)
{
  analysisAll(file.suffixPrimary = file.suffix, output.middleName = "premiRNA", includePrimary = c("pre_miRNA"), tailFilter = tailFilter, lenDis = 45:70, processes = processes, annotWithAS=annotWithAS, normalizeCnt = normalizeCnt)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/miRNA/vivo_5.57")
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/miRNA/S2_19636_19660_new_r5.57/")
#R: analysisAllmiRNA.modular(processes=16, tailFilter=0.12, module = c("isoCntLastDiNT"))
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/miRNA/exp5/miRNA_allRep")
#R: analysisAllmiRNA.modular(processes=16, tailFilter=0.12, module = c("isoCntLTM"), annotWithAS = TRUE)
#R: analysisAllmiRNA(processes=16, tailFilter=1) #needed for lenDisGM
#R: analysisAllmiRNA.modular(processes=16, tailFilter=1, module = c("isoCntLTM"), annotWithAS = TRUE)
#R" setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/miRNA/exp5/miRNA_allRep/monoAT")
#R: analysisAllmiRNA.modular(processes=16, tailFilter=1, module = c("isoCntMonotail", "isoCntTailLenMonotail"), annotWithAS = TRUE)
analysisAllmiRNA.modular <- function(file.suffix = "bam.seqCntMirna.txt.gz", tailFilter = 0.12, processes = 1, annotWithAS=TRUE, module = NULL)
{
  analysisAll.modular(file.suffixPrimary = file.suffix, output.middleName = "miRNA", tailFilter = tailFilter, processes = processes, annotWithAS=annotWithAS, module=module)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/miRNA/ath_R4120/monoAT")
#R: analysisAllmiRNAWithin.modular(processes=16, tailFilter=1, module = c("isoCntMonotail", "isoCntTailLenMonotail"), annotWithAS = TRUE, includePrimary = c("pre_miRNA"))
analysisAllmiRNAWithin.modular <- function(file.suffix = "bam.seqCntMirna.txt.gz", tailFilter = 0.12, processes = 1, annotWithAS=TRUE, module = NULL, normalizeCnt = TRUE, includePrimary = c("pre_miRNA_5p3p"))
{
  analysisAll.modular(file.suffixPrimary = file.suffix, output.middleName = "miRNA", tailFilter = tailFilter, processes = processes, annotWithAS=annotWithAS, module=module, normalizeCnt=normalizeCnt, includePrimary = includePrimary)
}


#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/pre_miRNA/S2_20128ff/")
#R: analysisAllPreMiRNA.modular(processes=16, tailFilter=0.05, module = c("isoCntLastDiNT"), annotWithAS = FALSE)
#R: setwd("/clustertmp/bioinfo/thomas/Exp5_preMiRNA/")
#R: analysisAllPreMiRNA.modular(processes=16, tailFilter=0.05, module = c("isoCntLTM"), annotWithAS = TRUE)
#R: analysisAllPreMiRNA(processes=16, tailFilter=1) #needed for lenDisGM
#R: analysisAllPreMiRNA.modular(processes=16, tailFilter=1, module = c("isoCntLTM"), annotWithAS = TRUE)
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/pre_miRNA/exp5/monoAT/")
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/miRNA/exp5/pre_miRNA_mm/monoAT")
#R: analysisAllPreMiRNA.modular(processes=16, tailFilter=1, module = c("isoCntMonotail", "isoCntTailLenMonotail"), annotWithAS = TRUE)


analysisAllPreMiRNA.modular <- function(file.suffix="bam.seqCntPreMirna.txt.gz", tailFilter = 0.05, processes = 1, annotWithAS=TRUE, module=NULL, normalizeCnt = TRUE)
{
  analysisAll.modular(file.suffixPrimary = file.suffix, output.middleName = "premiRNA", includePrimary = c("pre_miRNA"), tailFilter = tailFilter, lenDis = 45:70, processes = processes, annotWithAS=annotWithAS, module=module, normalizeCnt = normalizeCnt)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllmiRNAPlusEsi(processes=16)
#R: analysisAllmiRNAPlusEsi( aboveIsoCnt = 10, processes=16)
analysisAllmiRNAPlusEsi <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("pre_miRNA_5p3p", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "miRNAPlusEsi", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}

analysisAllGenomePlus <- function(file.suffix = "bam.seqCnt.txt.gz", tailFilter = 0.12, processes = 1, annotWithAS=TRUE, normalizeCnt = TRUE, aboveIsoCnt = -1)
{
  analysisAll(file.suffixPrimary = file.suffix, output.middleName = "miRNA", includePrimary = "pre_miRNA", tailFilter = tailFilter, processes = processes, annotWithAS=annotWithAS, normalizeCnt = normalizeCnt, aboveIsoCnt = aboveIsoCnt)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllsmallRNA(aboveIsoCnt = 1, processes=6) #ATTENTION: fast out of memory
#R: analysisAllsmallRNA( aboveIsoCnt = 10, processes=6)
analysisAllsmallRNA <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", 
                                includeType = c("pre_miRNA_5p3p", "transposable_element", "transposable_element_AS", 
                                                "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "smallRNA", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllsmallRNANoMiRNA(aboveIsoCnt = 1, processes=6) #ATTENTION: fast out of memory
#R: analysisAllsmallRNANoMiRNA( aboveIsoCnt = 10, processes=6)
analysisAllsmallRNANoMiRNA <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("transposable_element", "transposable_element_AS", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "smallRNA.NoMiRNA", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}

#Raphael.Manzenreither 3' GFP
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither")
#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllGFP(processes=6)
#R: analysisAllGFP(processes=6, includeType=c("virus"))
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither/GFP/")
#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllGFP(processes=6, includeType=c("pre_miRNA"), construct="GFP")
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither/TR4/")
#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllGFP(processes=6, includeType=c("pre_miRNA"), construct="TR4")
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither/TR57/")
#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllGFP(processes=6, includeType=c("pre_miRNA"), construct="TR57")
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither/GFP")
#R: analysisAllGFP(processes=6, includeType=c("pre_miRNA"), construct="GFP")
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither/ESI21")
#R: analysisAllGFP(processes=6, includeType=c("pre_miRNA"), construct="ESI21")

analysisAllGFP <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("gfp"), tailFilter = 0.07, aboveIsoCnt = -1, processes = 1, annotWithAS=FALSE, normalizeCnt = TRUE, construct=NULL)
{
    analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "3prime", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS, normalizeCnt=normalizeCnt, lenDis = 29:46)

    if (is.null(construct))
    {
        construct <- includeType
    }
    
    ap <- read.delim(paste0("allProperties.3prime.tf", tailFilter, ".txt"))
    ap$strand <- chartr("+-", "-+", ap$strand)
    cn <- colnames(ap)
    cn <- gsub("p3", "pTHREE", cn)
    cn <- gsub("p5", "p3", cn)
    cn <- gsub("pTHREE", "p5", cn)
    colnames(ap) <- cn
    library("Biostrings")
    gfp <- readDNAStringSet(paste0(construct, ".fa"))
    siRNA <- apply(ap[,c("chr", "p3", "strand")], 1, function(l)
    {
        if (l["strand"] == "+")
        {
            start=as.numeric(l["p3"]) - 10
            end=as.numeric(l["p3"])+10
            if (start < 1)
                start <- 1
            if (end > width(gfp)[names(gfp) == l["chr"]])
                end <- width(gfp)[names(gfp) == l["chr"]]
            return(as.character(reverseComplement(subseq(gfp[l["chr"]], start=start, end=end))))
        } else {
            start=as.numeric(l["p3"]) - 9
            end=as.numeric(l["p3"])+11
            if (start < 1)
                start <- 1
            if (end > width(gfp)[names(gfp) == l["chr"]])
                end <- width(gfp)[names(gfp) == l["chr"]]
            return(as.character(subseq(gfp[l["chr"]], start=start, end=end)))
        }
    })
    ap <- cbind(siRNA, ap)
    write.table(ap, file=paste0("allProperties.3prime.tf", tailFilter, ".txt"), quote=FALSE, row.names=FALSE, sep="\t")
}


#Raphael.Manzenreither 3' GFP
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither")
#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllGFPmonoT(processes=6)
#R: analysisAllGFPmonoT(processes=6, includeType=c("virus"))
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Raphael.Manzenreither/constructs_R4592/GFP")
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither/GFP/")
#R: analysisAllGFPmonoT(processes=6, includeType=c("pre_miRNA"))
analysisAllGFPmonoT <- function(file.suffixPrimary = "bam.monoT.seqCnt.txt.gz", file.suffixSecondary="bam.monoT.seqCnt.txt.gz", includePrimary="NONE", includeType = c("gfp"), tailFilter = 1, aboveIsoCnt = -1, processes = 1, annotWithAS=FALSE, normalizeCnt = TRUE)
{
    analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "3primeMonoT", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS, normalizeCnt=normalizeCnt, lenDis = 29:46)

    ap <- read.delim(paste0("allProperties.3primeMonoT.tf", tailFilter, ".txt"))
    ap$strand <- chartr("+-", "-+", ap$strand)
    cn <- colnames(ap)
    cn <- gsub("p3", "pTHREE", cn)
    cn <- gsub("p5", "p3", cn)
    cn <- gsub("pTHREE", "p5", cn)
    colnames(ap) <- cn
    library("Biostrings")
    gfp <- readDNAStringSet(paste0(includeType, ".fa"))
    siRNA <- apply(ap[,c("chr", "p3", "strand")], 1, function(l)
    {
        if (l["strand"] == "+")
        {
            start=as.numeric(l["p3"]) - 10
            end=as.numeric(l["p3"])+10
            if (start < 1)
                start <- 1
            if (end > width(gfp)[names(gfp) == l["chr"]])
                end <- width(gfp)[names(gfp) == l["chr"]]
            return(as.character(reverseComplement(subseq(gfp[l["chr"]], start=start, end=end))))
        } else {
            start=as.numeric(l["p3"]) - 9
            end=as.numeric(l["p3"])+11
            if (start < 1)
                start <- 1
            if (end > width(gfp)[names(gfp) == l["chr"]])
                end <- width(gfp)[names(gfp) == l["chr"]]
            return(as.character(subseq(gfp[l["chr"]], start=start, end=end)))
        }
    })
    ap <- cbind(siRNA, ap)
    write.table(ap, file=paste0("allProperties.3primeMonoT.tf", tailFilter, ".txt"), quote=FALSE, row.names=FALSE, sep="\t")
}


#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllTransposable(aboveIsoCnt=10, processes=6)
#R: analysisAllTransposable(aboveIsoCnt=1, processes=6) #ATTENTION: fast out of memory
analysisAllTransposable <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("transposable_element", "transposable_element_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "transposable", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllcisNAT(processes=16)
#R: analysisAllcisNAT(aboveIsoCnt=10, processes=16)
analysisAllcisNAT <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("cisNATs", "cisNATs_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "cisNAT", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllStructured(processes=16)
#R: analysisAllStructured(aboveIsoCnt=10, processes=16)
analysisAllStructured <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("structured_loci", "structured_loci_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "structured", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllVirus(processes=2)
analysisAllVirus <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("virus", "virus_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE, normalizeCnt=FALSE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "viruses", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS, normalizeCnt=normalizeCnt)
}


#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("/clustertmp/bioinfo/thomas/CG16940_small")
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/CG16940/CG16940_small")
#R: analysisAllFeat( aboveIsoCnt = 10, processes=2)
#R: setwd("/clustertmp/bioinfo/thomas/CG16940_large")
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/CG16940/CG16940_30plus")
#R: analysisAllFeat( aboveIsoCnt = 10, tailFilter=0.05, lenDis=26:50, processes=1)
#R: setwd("/clustertmp/bioinfo/thomas/CG16940_large2")
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/CG16940/CG16940_large40")
#R: analysisAllFeat( aboveIsoCnt = 10, tailFilter=0.05, lenDis=36:50, processes=1)
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Madalena.Pinto/CG16940/16940_IP_25Aug15/")
#R: analysisAllFeat( aboveIsoCnt = 10, tailFilter=0.05, lenDis=26:50, processes=1)

analysisAllFeat <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("mapped_reads", "mito", "rRNA", "rRNA_AS", "tRNA", "tRNA_AS", "good_reads", "pre_miRNA_5p3p", "pre_miRNA_5p3p_AS", "pre_miRNA", "pre_miRNA_AS", "transposable_element", "transposable_element_AS", "__no_feature", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS", "pseudogene", "pseudogene_AS", "ncRNA", "ncRNA_AS", "snoRNA", "snoRNA_AS", "snRNA", "snRNA_AS", "mRNA", "mRNA_AS", "intron", "intron_AS"), tailFilter = 0.12, aboveIsoCnt = 10, processes = 1, annotWithAS=TRUE, normalizeCnt="all", lenDis = 18:30)
{
  
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "allFeat", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS, normalizeCnt=normalizeCnt, lenDis=lenDis)
}


#Pseudo seqCnt file produced with ~/development/scripts/gfp-tailing.pl. Strand "-" not true but it will take 3p as 5p
#R: setwd("/groups/bioinfo/shared/Ameres.GRP/Stefan.Ameres/GFP-tailing")
#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllsiRNAGFP(processes=14)
analysisAllsiRNAGFP <- function(file.suffix = "seqCnt.txt.gz", includePrimary="siRNA", tailFilter = 0.12, processes = 1, normalizeCnt=TRUE, lenDis = 15:45, annotWithAS=TRUE, minLength = 10)
{
    analysisAll(file.suffixPrimary = file.suffix, output.middleName = "siRNA.minLength10", tailFilter = tailFilter, processes = processes, normalizeCnt=normalizeCnt, annotWithAS=annotWithAS, lenDis=lenDis, includePrimary=includePrimary, minLength=minLength)
}


##Raphael.Manzenreither - START

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("/clustertmp/bioinfo/thomas/Raphael.Manzenreither")
#R: analysisAllmiRNARaphael(processes=18)
#R: analysisAllesiRNAdsRNARaphael(processes=18)
#R: analysisAllesiRNARaphael(processes=18)
#R: analysisAlldsRNARaphael(processes=18)

analysisAllmiRNARaphael <- function(file.suffix = "bam.seqCntMirna.txt.gz", tailFilter = 0.12, processes = 1, annotWithAS=TRUE, normalizeCnt = c("pre_miRNA_5p3p", "pre_miRNA", "transposable_element", "transposable_element_AS", "__no_feature", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS", "mRNA", "mRNA_AS", "intron", "intron_AS", "dsRNA", "dsRNA_AS"))
{
  analysisAll(file.suffixPrimary = file.suffix, output.middleName = "miRNA", tailFilter = tailFilter, processes = processes, annotWithAS = annotWithAS, normalizeCnt = normalizeCnt)
}

analysisAllesiRNAdsRNARaphael <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("structured_loci", "structured_loci_AS", "dsRNA", "dsRNA_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE, normalizeCnt = c("pre_miRNA_5p3p", "pre_miRNA", "transposable_element", "transposable_element_AS", "__no_feature", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS", "mRNA", "mRNA_AS", "intron", "intron_AS", "dsRNA", "dsRNA_AS"))
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "esiRNA.dsRNA", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS, normalizeCnt = normalizeCnt)
}

analysisAllesiRNARaphael <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("structured_loci", "structured_loci_AS"), tailFilter = 0.12, aboveIsoCnt = 10, processes = 1, annotWithAS=TRUE, normalizeCnt = c("pre_miRNA_5p3p", "pre_miRNA", "transposable_element", "transposable_element_AS", "__no_feature", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS", "mRNA", "mRNA_AS", "intron", "intron_AS", "dsRNA", "dsRNA_AS"))
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "esiRNA", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS, normalizeCnt = normalizeCnt)
}

analysisAlldsRNARaphael <- function(file.suffixPrimary = "bam.seqCnt.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", includePrimary="NONE", includeType = c("dsRNA", "dsRNA_AS"), tailFilter = 0.12, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE, normalizeCnt = c("pre_miRNA_5p3p", "pre_miRNA", "transposable_element", "transposable_element_AS", "__no_feature", "cisNATs", "cisNATs_AS", "structured_loci", "structured_loci_AS", "mRNA", "mRNA_AS", "intron", "intron_AS", "dsRNA", "dsRNA_AS"))
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "dsRNA", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS, normalizeCnt = normalizeCnt)
}

##Raphael.Manzenreither - END


#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllpremiRNAPlusIntron(processes=16)
#R: analysisAllpremiRNAPlusIntron( aboveIsoCnt = 10, processes=16)
analysisAllpremiRNAPlusIntron <- function(file.suffixPrimary="bam.seqCntPreMirna.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", tailFilter = 0.05, includePrimary="pre_miRNA_5p3p", includeType = c("intron"), excludeId = NULL, aboveIsoCnt = 0, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "premiRNAPlusIntron", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, excludeId = excludeId, lenDis = 45:70, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}              

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllIntronSpliceSite(processes=16) #aboveIsoCnt = -1 to include all
#R: analysisAllIntronSpliceSite( aboveIsoCnt = 10, processes=16)
analysisAllIntronSpliceSite <- function(file.suffixPrimary="bam.seqCntPreMirna.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", tailFilter = 0.05, includePrimary="NONE", includeType = c("intron"), excludeId = NULL, aboveIsoCnt = -1, processes = 1, normalizeCnt = "all", spliceSite = TRUE, last8mer = TRUE, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "intronSpliceSite", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, excludeId = excludeId, lenDis = 45:70, aboveIsoCnt = aboveIsoCnt, processes = processes, normalizeCnt = normalizeCnt, spliceSite = spliceSite, last8mer = last8mer, annotWithAS=annotWithAS)
}              


#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllConstructEsi2(processes=16) #aboveIsoCnt = -1 to include all
#R: analysisAllConstructEsi2( aboveIsoCnt = 1, processes=16)
#R: analysisAllConstructEsi2( aboveIsoCnt = 10, processes=16)
analysisAllConstructEsi <- function(file.suffixPrimary="bam.seqCntMirna.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", tailFilter = 0.12, includePrimary="construct", includeType = c("ago2_sRNA", "cis_NAT", "endo_siRNA"), excludeId = NULL, aboveIsoCnt = -1, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "constructEsi", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, excludeId = excludeId, lenDis = 18:30, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}              


#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: analysisAllpremiRNAPlusIntronExcl( aboveIsoCnt = 10, processes=16)
analysisAllpremiRNAPlusIntronExcl <- function(file.suffixPrimary="bam.seqCntPreMirna.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", tailFilter = 0.05, includePrimary="pre_miRNA_5p3p", includeType = c("intron"), excludeId = NULL, excludeType = c("snoRNA", "tRNA"), aboveIsoCnt = 0, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "premiRNAPlusIntronExcl", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, excludeId = excludeId, excludeType = excludeType, lenDis = 45:70, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}              



#source("~/development/smallRNA/functions.miRNASummarize.R")
#setwd("~/shared/Ameres.GRP/Madalena.Pinto/miRNA/S2_19636_19660_new/")
#analysisAllmiRNATest(processes=5)
analysisAllmiRNATest <- function(file.suffix = "bam.seqCntMirna.txt.gz", tailFilter = 0.12, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffix, output.middleName = "miRNATest", tailFilter = tailFilter, processes = processes, annotWithAS=annotWithAS)
}

#source("~/development/smallRNA/functions.miRNASummarize.R")
#setwd("~/shared/Ameres.GRP/Madalena.Pinto/pre_miRNA/S2_20128ff/")
#analysisAllPreMiRNATest(processes=5)
analysisAllPreMiRNATest <- function(file.suffix="bam.seqCntPreMirna.txt.gz", tailFilter = 0.05, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffix, output.middleName = "premiRNATest", tailFilter = tailFilter, lenDis = 45:70, processes = processes, annotWithAS=annotWithAS)
}


#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("~/shared/Ameres.GRP/Madalena.Pinto/miRNA/S2_19636_19660_new/")#
#R: analysisAllmiRNAPlusEsiTest(processes=16)
#R: analysisAllmiRNAPlusEsiTest( aboveIsoCnt = 10, processes=16)
analysisAllmiRNAPlusEsiTest <- function(file.suffixPrimary="bam.seqCntMirna.txt.gz", file.suffixSecondary="bam.seqCnt.txt.gz", tailFilter = 0.12, includePrimary="pre_miRNA_5p3p", includeType = c("ago2_sRNA", "cis_NAT", "endo_siRNA"), excludeId = NULL, aboveIsoCnt = 0, processes = 1, annotWithAS=TRUE)
{
  analysisAll(file.suffixPrimary = file.suffixPrimary, output.middleName = "miRNAPlusEsiTest", file.suffixSecondary = file.suffixSecondary, tailFilter = tailFilter, includePrimary = includePrimary, includeType = includeType, excludeId = excludeId, lenDis = 18:30, aboveIsoCnt = aboveIsoCnt, processes = processes, annotWithAS=annotWithAS)
}              


#setwd("/home/imp/burkard/shared/Ameres.GRP/Madalena.Pinto/pre_miRNA/S2_20128ff/")
#allPropMean(allReplicates=list(WT=c("20128", "20141", "20142"),
#                          CG1091=c("20129"),
#                          CG16940=c("20130"),
#                          CG1091CG16940=c("20131"),
#                          Exp5=c("20143"),
#                          Exp5CG1091=c("20144"),
#                          Exp5CG16940=c("20145")),
##                        na.rm = TRUE,
##                        na.rm = FALSE,
#                        processes=processes)

#setwd("/clustertmp/bioinfo/thomas/Exp5_preMiRNA")
#allPropMean(allReplicates=list(UT=c("30259", "30269", "30279"),
#              dsGFP=c("30260", "30270", "30280"),
#              dsTrf4=c("30261", "30271", "30281"),
#              dsExp5=c("30262", "30272", "30282"),
#              dsExp5dsTrf4=c("30263", "30273", "30283"),
#              dsExp5dsGFP=c("30264", "30274", "30284"),
#              dsMtr4=c("30265", "30275", "30285"),
#              dsRrp6=c("30266", "30276", "30286"),
#              dsExp5dsMtr4=c("30267", "30277", "30287"),
#              dsExp5dsRrp6=c("30268", "30278", "30288")),
##                        na.rm = TRUE,
##                        na.rm = FALSE,
#                        processes=processes)
            

#setwd("/groups/bioinfo/shared/Ameres.GRP/Raphael.Manzenreither/miRNA/M2196/")
#allPropMean(allReplicates=list(S2SLuc=c("26141"),
#                S2SCG1091=c("26142"),
#                S2STrf4.1=c("26143"),
#                S2STrf4.2=c("26144"),
#                S2SCG11418=c("26145"),
#                S2SMkg=c("26146"),
#                Hen1KOLuc=c("26147"),
#                Hen1KOCG1091=c("26148"),                
#                Hen1KOTrf4.1=c("26149"),
#                Hen1KOTrf4.2=c("26150"),                
#                Hen1KOCG11418=c("26151"),
#                Hen1KO_Mkg=c("26152")),
#                na.rm = TRUE,
#                processes=1)

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("~/shared/Ameres.GRP/Madalena.Pinto/miRNA/S2_19636_19660_new/")
#R: allPropMean()
#R: allPropMean(na.rm = FALSE)
allPropMean <- function(allProp.files = NULL, allReplicates=list(UT=c("19636", "19645", "19654"),
                          GFP=c("19637", "19638", "19646", "19647"),
                          CG1091=c("19639", "19648", "19655"),
                          CG16940=c("19640", "19649", "19656"),
                          CG1091CG16940=c("19641", "19650", "19657"),
                          Exp5=c("19642", "19651", "19658"),
                          Exp5CG1091=c("19643", "19652", "19659"),
                          Exp5CG16940=c("19644", "19653", "19660")),
                        na.rm = TRUE,
                        processes = 1
                        )
  {


    if (is.null(allProp.files))
      {
        allProp.files <- list.files(".", "allProperties\\..*txt$")
      }
    
    for (f in allProp.files)
      {
        allProp <- read.delim(f)
        if (na.rm)
          {
            allProp <- allProp[apply(allProp[,grep("align.GM|align.PM", colnames(allProp))], 1, function(l) { all(!is.na(l))}), ]
#            allProp <- allProp[apply(allProp, 1, function(l) { all(!is.na(l))}), ]
          }
        
        out <- allProp[,c("chr", "p5", "strand", "gname", "isoSeed")]

        for (n in names(allReplicates))
          {
            out <- cbind(out, meanTable(allProp, allReplicates[[n]], n))
          }
        if (na.rm)
          {
            write.table(out, file = sub("allProperties", "allPropertiesMeanNoNA", f), quote = F, row.names = F, sep = "\t")
          } else {
            write.table(out, file = sub("allProperties", "allPropertiesMean", f), quote = F, row.names = F, sep = "\t")
          }
      }
  }

meanTable <- function(allProp, replicate, repName)
  {
    cn <- gsub(paste0(replicate[1], "$"), "", grep(paste0(".", replicate[1], "$"), colnames(allProp), value = TRUE))
    tmp <- lapply(cn, function(n)
           {
             if ( n == "tailValue." )
               {
                 r <- as.data.frame(apply(as.data.frame(allProp[,paste0(n, replicate)]), 1, paste, collapse = ";"))
                 colnames(r) <- "tailValue"
                 return(r)
               } else if ( n == "miR." )
                 {
                   r <- as.data.frame(apply(as.data.frame(allProp[,paste0(n, replicate)]), 1, paste, collapse = ";"))
                   colnames(r) <- "miR"
                   return(r)
                 } else if ( n == "last8merValue." )
                   {
                     r <- as.data.frame(apply(as.data.frame(allProp[,paste0(n, replicate)]), 1, paste, collapse = ";"))
                     colnames(r) <- "last8merValue"
                     return(r)
                   } else {
                     r <- cbind(apply(as.data.frame(allProp[,paste0(n, replicate)]), 1, function(x) { mean(as.num(x), na.rm = TRUE)}),
                                apply(as.data.frame(allProp[,paste0(n, replicate)]), 1, function(x) { sd(as.num(x), na.rm = TRUE)}))
                     colnames(r) <- paste0(n, c("mean", "sd"))
                     return(r)
                   }
           })
    tmp <- do.call(cbind, tmp)
    colnames(tmp) <- paste(colnames(tmp), repName, sep = ".")

    return(tmp)
    
  }

pairs.rnaseq <- function(data, probs = 0.999, main = NULL)
  {
    m <- as.matrix(data)
    t <- apply(m, 2, quantile, probs=probs)
    for (i in 1:length(m[1,]))
    {
      m[ m[,i] > t[i] , i] <- NA
    }
    pairs(m, lower.panel=panel.smooth.rnaseq, upper.panel=panel.cor.rnaseq, main = main)
  }

panel.smooth.rnaseq <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, rm.max = FALSE, ...) 
{
  
  if (rm.max)
    {
      threshold.x <- quantile(x, probs = probs)
      threshold.y <- quantile(y, probs = probs)
      ok <- (x > 0 | y > 0) & x < threshold.x & y < threshold.y & is.finite(x) & is.finite(y)
    } else {
      ok <- (x > 0 | y > 0) & is.finite(x) & is.finite(y)
    }
    
#  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  points(x[ok], y[ok], pch = ".", col = col, bg = bg, cex = 2, ...)

  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

panel.cor.rnaseq <- function(x, y, digits=2, prefix="", cex.cor, rm.max = FALSE, ...)
  {

    if (rm.max)
      {
        threshold.x <- quantile(x, probs = probs)
        threshold.y <- quantile(x, probs = probs)
        ok <- (x > 0 | y > 0) & x < threshold.x & y < threshold.y & is.finite(x) & is.finite(y)
      } else {
        ok <- (x > 0 | y > 0) & is.finite(x) & is.finite(y)
      }
    
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r = (cor(x[ok], y[ok], use = "complete"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor))
      cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex * abs(r))
  }


plotCor <- function(allProp.files = NULL, allReplicates=allReplicates,
                        na.rm = TRUE,
                        processes = 1
                        )
  {

    if (is.null(allProp.files))
      {
        allProp.files <- list.files(".", "allProperties\\..*txt$")
      }
    
    for (f in allProp.files)
      {
        allProp <- read.delim(f)
        if (na.rm)
          {
            allProp <- allProp[apply(allProp[,grep("align.GM|align.PM", colnames(allProp))], 1, function(l) { all(!is.na(l))}), ]
          }

        gm <- allProp[,grep("align.GM", colnames(allProp))]
        pm <- allProp[,grep("align.PM", colnames(allProp))]


        pdf(paste0(f, ".cor.pdf"))
        for (n in names(allReplicates))
          {
            pairs.rnaseq(gm[,paste0("align.GM.", allReplicates[[n]])], probs = 1, main = paste(n, "(GM, q = 1, noNA)"))
            pairs.rnaseq(gm[,paste0("align.GM.", allReplicates[[n]])], probs = 0.95, main = paste(n, "(GM, q = 0.95, noNA)"))

            pairs.rnaseq(pm[,paste0("align.PM.", allReplicates[[n]])], probs = 1, main = paste(n, "(PM, q = 1, noNA)"))
            pairs.rnaseq(pm[,paste0("align.PM.", allReplicates[[n]])], probs = 0.95, main = paste(n, "(PM, q = 0.95, noNA)"))
          }
        dev.off()
      }
  }


tailStatisticsByLengthFrac <- function(file.suffix=".*seqCnt.txt", features = c("tRNA", "pre_miRNA_5p3p"))
  {
    for (type in features)
      {
        tailStatistics <- NULL
        for (f in list.files(".", paste0(file.suffix, "$")))
          {
            feat <- read.delim(f)


            feat.type <- feat[grep(type, as.character(feat$type)), ]
            cntTotal <- sum(feat.type$count)
            
            for (maxFracTail in c(1, 0.3, 0.25, 0.2, 0.17, 0.15, 0.12, 0.1, 0.05))
              {
                feat.f <- feat.type[nchar(as.character(feat.type$tail))/nchar(as.character(feat.type$sequence)) <= maxFracTail, ]
                cntTotal.f <- sum(feat.f$count)
                
                feat.tailed <- feat.f[feat.f$tail != "", ]
                
                cnt <- sum(feat.tailed$count)
                fractLoss <- (cntTotal-cntTotal.f)/cntTotal
                firstA <- sum(feat.tailed[grep("^A", as.character(feat.tailed$tail), perl = TRUE), "count"])/cnt
                firstC <- sum(feat.tailed[grep("^C", as.character(feat.tailed$tail), perl = TRUE), "count"])/cnt
                firstG <- sum(feat.tailed[grep("^G", as.character(feat.tailed$tail), perl = TRUE), "count"])/cnt
                firstT <- sum(feat.tailed[grep("^T", as.character(feat.tailed$tail), perl = TRUE), "count"])/cnt
                
                tailStatistics <- rbind(tailStatistics, c(file=f, type=type, maxTailLength=maxFracTail, fractLoss=fractLoss,cntAll=cntTotal.f, cntTail=cnt, firstA=firstA, firstC=firstC, firstG=firstG, firstT=firstT))
              }
          }
        write.table(tailStatistics, file=paste("tailStatisticsByLengthFrac", type, "txt", sep = "."), row.names = FALSE, sep = "\t", quote = FALSE)
      }

  }


plotLD <- function(allProp)
  {
    isoSeed <- as.character(allProp["isoSeed"])
    gname <- as.character(allProp["gname"])
    samples <- unique(gsub(".*\\.", "", grep("lenDisGM\\.", names(allProp), value = TRUE)))

    if (ceiling(length(samples)/3) > 5)
      {
        pdf(paste(gname, isoSeed, "pdf", sep = "."), width = 21, height = 30)
      } else {
        pdf(paste(gname, isoSeed, "pdf", sep = "."), width = 21, height = ceiling(length(samples)/3)*6)
      }

    #allProp["chr"] <- gsub("chr", "", allProp["chr"])
    allProp["chr"] <- gsub("chr", "", allProp["chr"], ignore.case=T)    

    if (as.character(allProp["strand"] == "+"))
          {
            #p5 is 0-based for + strand (=start of BED region)
            sequence <-getSeq(Dmelanogaster,
                              name=paste("chr", as.character(allProp["chr"]), sep = ""),
                              start=(as.num(allProp["p5"])+1),
                              strand="+",
                              width=100)
          } else {
            #p5 is 1-based for - strand (=end of BED region)
            sequence <-getSeq(Dmelanogaster,
                              name=paste("chr", as.character(allProp["chr"]), sep = ""),
                              end=(as.num(allProp["p5"])),
                              strand="-",
                              width=100)
          }
        sequence <- unlist(strsplit(as.character(sequence), split=""))
        
        for (lenDis.names in list(c("lenDisGM", "lenDisPM"), c("lenDisPM", "lenDisPMeT")))
          {
#            samples <- unique(gsub(".*\\.", "", grep(paste(lenDis.names[1], "\\.", sep = ""), colnames(allProp), value = TRUE)))

            yAxis <- NULL
            for (l in lenDis.names)
              {
                if (is.null(yAxis))
                  {
                    yAxis <- as.numeric(allProp[c(grep(paste(l, "\\..*", sep = ""), names(allProp), value = TRUE))])
                  } else {
                    yAxis <-  rbind(yAxis, as.numeric(allProp[c(grep(paste(l, "\\..*", sep = ""), names(allProp), value = TRUE))]))
                  }
              }
            yAxis[is.na(yAxis)] <- 0
            yAxis <- max(colSums(yAxis), na.rm = TRUE)
            yAxis <- ceiling(yAxis/10)*10
        
            par(mfrow=c(ceiling(length(samples)/3),3))
            legend <- TRUE
        
            for (s in samples)
              {

                lenDis <- NULL
                for (l in lenDis.names)
                  {
                    allProp.ld <- as.numeric(allProp[c(grep(paste(l, "\\..*", s, sep = ""), names(allProp), value = TRUE))])
                    names(allProp.ld) <- gsub(paste(l, "\\.(.*)\\.", s, sep = ""), "\\1", c(grep(paste(l, "\\..*", s, sep = ""), names(allProp), value = TRUE)))
                    allProp.ld[is.na(allProp.ld)] <- 0
                    if (is.null(lenDis))
                      {
                        lenDis <- allProp.ld
                      } else {
                        lenDis <- cbind(lenDis, allProp.ld)
                      }
                  }
                colnames(lenDis) <- lenDis.names

                barplot(t(lenDis), main = s, ylim=c(0,yAxis), ylab = "RPM", names.arg=sequence[as.numeric(rownames(lenDis))], xlab = paste(c("miRNA: ", sequence[1:max(as.numeric(rownames(lenDis)))]), collapse = ""), col = c("gray20", "gray90"))
                if (legend == TRUE)
                  {
                    legend("topright", legend = colnames(lenDis), fill = c("gray20", "gray90"))
                    legend <- FALSE
                  }
              }
          }
        dev.off()
  }

plotLenDis <- function(allProp.files, processes = 1)
  {
    allProp <- NULL
    for (f in allProp.files)
      {
        if (is.null(allProp))
          {
            allProp <- read.delim(f)
          } else {
            tmp <- read.delim(f)
            tmp <- merge(allProp, tmp, all = TRUE, by = c("chr", "p5", "strand", "gname", "isoSeed"))
          }
      }        

    system("mkdir lenDis")
    wd <- getwd()
    setwd(paste(wd, "lenDis", sep = "/"))


    if (processes == 1)
      {
        apply(allProp, 1, plotLD)
      } else {
        cl <- makeCluster(rep("localhost", (processes-1)),
                          type = "SOCK"
                          )

        #export all functions to workers
        ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
        clusterExport(cl, ex)
        clusterEvalQ(cl, {library("BSgenome.Dmelanogaster.UCSC.dm3")})
        parApply(cl, allProp, 1, plotLD)
        stopCluster(cl)
      }
    
    setwd(wd)
  }


seqCntLenDis <- function(seqCnt, file.suffix = "bam.seqCnt.txt.gz", seqSeed = "")
  {
#    setwd("~/shared/Ameres.GRP/Raphael.Manzenreither/miRNA/gfp-bamH1-luc/")
#    setwd("~/shared/Ameres.GRP/Madalena.Pinto/miRNA/S2_19636_19660_new/")
#    setwd("~/shared/Ameres.GRP/Madalena.Pinto/pre_miRNA/S2_20128ff/")

    seqCnt.all <- list.files(".", paste0(file.suffix, "$"))

    if (seqSeed == "")
      {
        pdf("seqCnt.lenDis.pdf", width = 20, height = 10)
      } else {
        pdf(paste0("seqCnt.lenDis.", seqSeed, ".pdf"), width = 20, height = 10)
      }
    
    for (seqCnt.file in seqCnt.all)
      {
    
        seqCnt <- read.delim(seqCnt.file)
        seqCnt <- cbind(seqCnt, len=nchar(as.character(seqCnt$sequence)))
        seqCnt <- seqCnt[grep(paste0("^", seqSeed), as.character(seqCnt$sequence)), ]

        seqCnt$count <- (seqCnt$count * 1000000) / sum(seqCnt$count)

        seqCnt.file <- gsub(".*#", "", seqCnt.file)
        seqCnt.file <- gsub("\\..*", "", seqCnt.file)
        seqCnt.file <- gsub("_.*", "", seqCnt.file)
        
        s.gm <- subset(seqCnt, seqCnt$tailLen == 0)
        s.pm <- subset(seqCnt, seqCnt$tailLen > 0)
        
        par(mfcol=c(1,3))
        if (nrow(seqCnt) > 0)
          barplot(sapply(split(seqCnt$count, seqCnt$len), sum), main = paste("ALL", seqCnt.file))
        if (nrow(s.gm) > 0)
          barplot(sapply(split(s.gm$count, s.gm$len), sum), main = paste("GM", seqCnt.file))
        if (nrow(s.pm) > 0)
          barplot(sapply(split(s.pm$count, s.pm$len), sum), main = paste("PM", seqCnt.file))
      }
    dev.off()
  }

#setwd("~/shared/Ameres.GRP/Raphael.Manzenreither/miRNA/gfp-bamH1-luc/")
#allProp <- read.delim("allProperties.miRNAPlusEsi.tf0.12.ic10.txt")
#allProp.file <- "allProperties.miRNAPlusEsi.tf0.12.ic10.txt"

guidePassengerPair <- function(allProp.file, l.siRNA = 21, l.3pOverhang = 2) #21nt - 2nt 3' overhang (total sequence covered 23nt)
  {

    seq.GBL <- "GGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAggatccTCATCCCTGATCTGATCGGAATGGGTAAGTCCGGCAAGAGCGGGAATGGCTCATATCGCCTCCTGGATCACTACAAGTACCTCACCGCTTGGTTCGAGCTGCTGAACCTTCCAAAGAAAATCATCTTTGTGGGCCACGACTGGGGGGCTTGTCTGGCCTTTCACTACTCCTACGAGCACCAAGACAAGATCAAGGCCATCGTCCATGCTGAGAGTGTCGTGGACGTGATCGAGTCCTGGGACGAGTGGCC"
    
    distance <- l.siRNA - l.3pOverhang
    
    allProp <- read.delim(allProp.file)
    
    p5.1base <- as.data.frame(list(p5=1:nchar(seq.GBL)))
    allProp.1base <- split(allProp, as.character(allProp$strand))
    allProp.1base[["+"]]$p5 <- allProp.1base[["+"]]$p5 + 1
    allProp.1base[["+"]] <- merge(p5.1base, allProp.1base[["+"]], all.x = TRUE)
    allProp.1base[["-"]] <- merge(p5.1base, allProp.1base[["-"]], all.x = TRUE)

    allProp.1base[["+"]] <- allProp.1base[["+"]][order(as.num(allProp.1base[["+"]]$p5)), ]
    allProp.1base[["-"]] <- allProp.1base[["-"]][order(as.num(allProp.1base[["-"]]$p5)), ]
    
    pdf("barplot.GM.1based.pdf", width = 20, height = 10)
    barplot(allProp.1base[["+"]][, "align.GM.21461"], ylim = c(-27000, 27000), main = "GM.21461")
    barplot(-allProp.1base[["-"]][, "align.GM.21461"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.1base[["+"]][, "align.GM.21462"], ylim = c(-27000, 27000), main = "GM.21462")
    barplot(-allProp.1base[["-"]][, "align.GM.21462"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.1base[["+"]][, "align.GM.21463"], ylim = c(-27000, 27000), main = "GM.21463")
    barplot(-allProp.1base[["-"]][, "align.GM.21463"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.1base[["+"]][, "align.GM.21464"], ylim = c(-27000, 27000), main = "GM.21464")
    barplot(-allProp.1base[["-"]][, "align.GM.21464"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.1base[["+"]][, "align.GM.21465"], ylim = c(-27000, 27000), main = "GM.21465")
    barplot(-allProp.1base[["-"]][, "align.GM.21465"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.1base[["+"]][, "align.GM.21466"], ylim = c(-27000, 27000), main = "GM.21466")
    barplot(-allProp.1base[["-"]][, "align.GM.21466"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.1base[["+"]][, "align.GM.21467"], ylim = c(-27000, 27000), main = "GM.21467")
    barplot(-allProp.1base[["-"]][, "align.GM.21467"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.1base[["+"]][, "align.GM.21468"], ylim = c(-27000, 27000), main = "GM.21468")
    barplot(-allProp.1base[["-"]][, "align.GM.21468"], ylim = c(-27000, 27000), add = TRUE)

    dev.off()

    write.table(allProp.1base[["+"]], file = sub(".txt", ".1basedPlus.txt", allProp.file), quote = FALSE, row.names = FALSE, sep = "\t")
    write.table(allProp.1base[["-"]], file = sub(".txt", ".1basedMinus.txt", allProp.file), quote = FALSE, row.names = FALSE, sep = "\t")
    
    
    n <- colnames(allProp)
    cn <- c(n[1:5], n[grep("align.GM", n)])
    allProp <- allProp[,cn]

    allProp.strand <- split(allProp, as.character(allProp$strand))
    

    if (all(allProp.strand[["+"]]$gname == "GBLp"))
      {
        allProp.strand[["+"]] <- cbind(sequence=substring(seq.GBL,
                                         (allProp.strand[["+"]]$p5 - l.3pOverhang + 1),
                                         (allProp.strand[["+"]]$p5 + l.siRNA)
                                      ), allProp.strand[["+"]])
      }

    
    
    allProp.strand[["+"]]$gname <- paste0(sub("p$", "", as.character(allProp.strand[["+"]]$gname)), as.character(allProp.strand[["+"]]$p5))
    allProp.strand[["-"]]$gname <- paste0(sub("m$", "", as.character(allProp.strand[["-"]]$gname)), (as.numeric(as.character(allProp.strand[["-"]]$p5))-distance))
    allProp.gp <- merge(allProp.strand[["+"]], allProp.strand[["-"]], by = "gname", all = TRUE)

    allProp.gp <- allProp.gp[order(as.numeric(as.character(allProp.gp$p5.x))), ]
    write.table(allProp.gp, file = sub("txt$", "guidePass.txt", allProp.file), quote = FALSE, row.names = FALSE, sep = "\t")


    allProp.heat <- allProp.gp[, grep("align.GM", colnames(allProp.gp))]
    allProp.heat <- data.matrix(allProp.heat)
    allProp.heat[is.na(allProp.heat)] <- 0

    allProp.heat <- allProp.heat[,grep("21459", colnames(allProp.heat), invert = TRUE)]
    allProp.heat <- allProp.heat[,grep("21460", colnames(allProp.heat), invert = TRUE)]
    
    library(gplots)

#    heatmapcolors <-  colorRampPalette(c("#00007F", "#007FFF", "#00FFFF","#7FFF7F", "#FFFF00","#FF7F00", "#FF0000","#7F0000"), space="rgb")(30)
    library("RColorBrewer")
    heatmapcolors<-colorRampPalette(brewer.pal(n=9,name="Blues"), space="rgb")(15)
    
    pdf("heatmap.GM.pdf")
#    heatmap.2(allProp.heat, Rowv = NA, Colv = NA)
#    heatmap.2(allProp.heat, Colv = NA)
    heatmap.2(log2(1+allProp.heat), col = heatmapcolors)
    heatmap.2(log2(1+allProp.heat), col = heatmapcolors, Rowv = NA, Colv = NA)

    heatmap(allProp.heat, Rowv = NA, Colv = NA)
    heatmap(allProp.heat, Colv = NA)

    dev.off()

    pdf("barplot.GM.pdf", width = 20, height = 10)
    barplot(allProp.gp[, "align.GM.21461.x"], ylim = c(-27000, 27000))
    barplot(-allProp.gp[, "align.GM.21461.y"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.gp[, "align.GM.21462.x"], ylim = c(-27000, 27000))
    barplot(-allProp.gp[, "align.GM.21462.y"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.gp[, "align.GM.21463.x"], ylim = c(-27000, 27000))
    barplot(-allProp.gp[, "align.GM.21463.y"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.gp[, "align.GM.21464.x"], ylim = c(-27000, 27000))
    barplot(-allProp.gp[, "align.GM.21464.y"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.gp[, "align.GM.21465.x"], ylim = c(-27000, 27000))
    barplot(-allProp.gp[, "align.GM.21465.y"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.gp[, "align.GM.21466.x"], ylim = c(-27000, 27000))
    barplot(-allProp.gp[, "align.GM.21466.y"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.gp[, "align.GM.21467.x"], ylim = c(-27000, 27000))
    barplot(-allProp.gp[, "align.GM.21467.y"], ylim = c(-27000, 27000), add = TRUE)

    barplot(allProp.gp[, "align.GM.21468.x"], ylim = c(-27000, 27000))
    barplot(-allProp.gp[, "align.GM.21468.y"], ylim = c(-27000, 27000), add = TRUE)
    dev.off()

    
    plot(allProp.gp$p5.x, allProp.gp$align.GM.21461.x, type = "l", ylim = c(-27000, 27000))
    lines(allProp.gp$p5.x, -allProp.gp$align.GM.21461.y)

    plot(allProp.gp$p5.x, allProp.gp$align.GM.21462.x, type = "l", ylim = c(-27000, 27000))
    lines(allProp.gp$p5.x, -allProp.gp$align.GM.21462.y)
    
    plot(allProp.gp$p5.x, allProp.gp$align.GM.21463.x, type = "l", ylim = c(-27000, 27000))
    lines(allProp.gp$p5.x, -allProp.gp$align.GM.21463.y)

    plot(allProp.gp$p5.x, allProp.gp$align.GM.21464.x, type = "l", ylim = c(-27000, 27000))
    lines(allProp.gp$p5.x, -allProp.gp$align.GM.21464.y)

    



    
#    plot(allProp.heat[,"align.GM.21461.x"], allProp.heat[,"align.GM.21461.y"])
  }

#R: source("~/development/smallRNA/functions.miRNASummarize.R")
#R: setwd("~/shared/Ameres.GRP/Madalena.Pinto/pre_miRNA/exp5/")
#R: plotBalloon()
plotBalloon <- function()
  {
    library(ggplot2)

#setwd("/clustertmp/bioinfo/thomas/Exp5_preMiRNA/")
files <- list.files(pattern="isoCntLTM.*txt$")
for (f in files)
  {
    allProp.modular <- read.delim(f)
    allProp <- read.delim(sub("isoCntLTM.", "", f))
    allProp <- merge(allProp, allProp.modular, by = c("chr", "p5", "strand", "gname", "isoSeed"))

    for (inclGM in c(TRUE, FALSE))
      {
        if (inclGM)
          {
            pdf(sub(".txt$", ".balloon.pdf", f), height = 15, width = 10)
          } else {
            pdf(sub(".txt$", ".balloonPM.pdf", f), height = 15, width = 10)
          }
        
        apply(allProp, 1, function(ltm)
              {
#ltm <- as.character(t(allProp)[,1])
#names(ltm) <- colnames(allProp)
                gname <- ltm["gname"]

                if (as.character(ltm["strand"]) == "+")
                  {
                    #p5 is 0-based for + strand (=start of BED region)
                    sequence <-getSeq(Dmelanogaster,
                                      name=paste("chr", as.character(ltm["chr"]), sep = ""),
                                      start=(as.num(ltm["p5"])+1),
                                      strand="+",
                                      width=100)
                  } else {
                    #p5 is 1-based for - strand (=end of BED region)
                    sequence <-getSeq(Dmelanogaster,
                                      name=paste("chr", as.character(ltm["chr"]), sep = ""),
                                      end=(as.num(ltm["p5"])),
                                      strand="-",
                                      width=100)
                  }
                sequence <- unlist(strsplit(as.character(sequence), split=""))

                
                if (inclGM)
                  {
                    names(ltm) <- sub("lenDisGM\\.(\\d+)", "LTM.GM_\\1_0", names(ltm), perl = TRUE)
                  }
                names(ltm) <- sub("\\.mean\\.", ".", names(ltm), perl = TRUE)
            
                ltm <- ltm[grep("LTM", names(ltm))]
                ltm <- ltm[grep("\\.sd\\.", names(ltm), invert=TRUE, perl = TRUE)]
                ltm[is.na(ltm)] <- 0

                ltm <- cbind(do.call(rbind, strsplit(sub("LTM.(.*?)\\.(.*)", "\\1_\\2", names(ltm), perl = TRUE), "_")), ltm)

                colnames(ltm) <- c("monotail", "alignLength", "taillength", "condition", "RPM")
#        ltm <- ltm[as.numeric(as.character(ltm[,"ltm"])) > 0, ]
                if (any(as.numeric(as.character(ltm[,"RPM"])) > 0)) {
                  ltm <- as.data.frame(ltm)
                  ltm[,"RPM"] <- as.numeric(as.character(ltm[,"RPM"]))
              
                  p <- ggplot(ltm, aes(x=alignLength, y=taillength, weight=RPM, colour=monotail, size=RPM)) +
                    geom_point( alpha=0.8, guide="none") +
                      facet_grid(condition ~ .) +
                        scale_colour_brewer(palette="Set1", type="qual", name="Monotail") +
                          guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                            ggtitle(gname) + 
                              scale_size_area(max_size = 6) +
                                theme_bw() +
                                  scale_x_discrete( labels=sequence[min(as.num(ltm[,"alignLength"])):max(as.num(ltm[,"alignLength"]))])
                  
                  show(p)
                  return(0)
                } else {
                  return(1)
                }
              })
        dev.off()
      }
  }
}

plotTypePie <- function(pdf=TRUE)
  {
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
#setwd("/clustertmp/bioinfo/thomas/Exp5_preMiRNA/")

    type <- read.delim("cnt.typeHierarchy.txt")
    type[is.na(type)] <- 0
    type <- type[grep("_read", as.character(type[,1]), invert=TRUE), ]
    type.good <- type[grep("mito", as.character(type[,1]), invert=TRUE), ]
    type.good <- type.good[grep("rRNA", as.character(type.good[,1]), invert=TRUE), ]
    type.good <- type.good[grep("tRNA", as.character(type.good[,1]), invert=TRUE), ]
    type.m <- melt(type, id.vars="type")
    type.good.m <- melt(type.good, id.vars="type")

    if (pdf)
        {
            pdf("cnt.typeHierarchy.pie.pdf", width=10, height=15)
        }
    
    colourCount = length(unique(type.m$type))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    
    vis1 = ggplot(data=type.m, aes(x=factor(1), y=value, fill=type)) +
      geom_bar(stat="identity", position="fill") +
        coord_polar(theta="y") +
          facet_wrap(~variable, ncol=3) +
            scale_fill_manual(values = getPalette(colourCount))


    colourCount = length(unique(type.good.m$type))
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))

    vis2 = ggplot(data=type.good.m, aes(x=factor(1), y=value, fill=type)) +
      geom_bar(stat="identity", position="fill") +
        coord_polar(theta="y") +
          facet_wrap(~variable, ncol=3) +
            scale_fill_manual(values = getPalette(colourCount))


    if (pdf)
    {
        show(vis1)
        show(vis2)
        dev.off()
    } else {
        show(vis1)
        show(vis2)
    }
  }
