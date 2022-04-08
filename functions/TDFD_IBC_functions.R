# ====================================================================== #
#                                                                        #
# Functions to measure and test the alpha and beta diversity             #
# Alpha functional diversity is defined as:                              # 
# - the Rao's quadratic entropy (Rao 1982)                               #
# Beta taxonomic diversity is decomposed according to:                   #
# - the framework of Baselga [Global Ecol. Biogeogr. 19, 134-143 (2010)] #
#                                                                        #
# ====================================================================== #


# ======================================== TABLE OF CONTENTS ======================================== #  
#                                                                                                     #
# 01. Rao.CH: Function to calculate the FD from the Rao's quadratic entropy                           #
# 02. betadiv: Function to decompose taxonomic beta-diversity and test significance from null models  #
#                                                                                                     #
# =================================================================================================== #


#--------------------------------------------------------------------#
## 01. Function to calculate the functional diversity from the Rao's quandratic entropy
##     Formule of Rao's index from Champely & Chessel [Environ. Ecol. Stat. 9, 167-177 (2002)]
## WARNING ##
## This function gives similar results that function divc in ade4
## WARNING ##
##
# The function 'Rao.CH' runs the analysis
# It uses three elements:
#' @param taxa. A community matrix (C x S) with C samples (or communities) as rows, and S species as columns
#' @param trait. A matrix or data frame (S x T) with S species as rows, and T traits as columns, must have species names as rownames
#' @param method. whether "euclidean" method is used, dissimilarity matrix is computes from euclidean distance
#
# It return one element:
#' @return Rao.FD. A data frame (C x 1) listing for each C sample the value of functional diversity
# - low value of Rao.FD: species are functionally similar
# - high value of Rao.FD: species are functionally distinct
##
Rao.CH <- function(taxa, trait, method = "euclidean", scale = TRUE) {
  # calculating relative abundances
  Tabun <- taxa / rowSums(taxa)
  Tabun <- drop(as.matrix(Tabun))
  # matrix of functional distance among species
  Fdist <- vegdist(trait, method = method)
  Fdist <- as.matrix(Fdist)
  # observed rao's diversity (looping on the n sites)
  Rao.FD <- as.data.frame(rep(0, nrow(Tabun)))
  names(Rao.FD) <- "Rao_CH"
  for (i in 1:nrow(Tabun)) {
    Rao.FD[i, ] <- (t(Tabun[i, ]) %*% (Fdist^2) %*% Tabun[i, ]) / 2
  }
  if (scale == TRUE) {
    raomax <- divcmax(as.dist(Fdist))$value
    Rao.FD <- Rao.FD/raomax
  }
  return(Rao.FD)
}


#--------------------------------------------------------------------#
## 02. Function to decompose taxonomic beta-diversity and test significance from null models 
##
# The function 'betadiv' runs the analysis
# It uses four elements:
#' @param taxa. A community matrix (C x S) with C samples (or communities) as rows, and S species as columns
#' @param site. A vector corresponding to the names of sites
#' @param method. Whether "BR" method is used, Bray Curtis dissimilarity index is computes from abundance data; and whether "SOR" is used, sorensen dissimilarity index is computes from incidence data
#' @param nrepet. The number of repetition of the null model
#
# It return one element:
#' @return res.beta. A data frame with in row the group site and in columns observed and expected (null model) value of the three component of beta diversity (i.e., total beta diversity, turnover, and nestedness)
##
betadiv <- function(taxa, site, method = "BR", nrepet = 999) {
  require(betapart)
  require(ade4)
  require(vegan)
  f1 <- function(x) {
    X <- x[x$X1 <= 4 & x$X2 >= 18, ]
    X <- apply(X, 2, function(i) mean(i))[3]
    return(X)
  }
  f2 <- function(x, keep1) {
    X <- x[rownames(x) %in% keep1, ]
    X <- apply(X, 2, function(i) mean(i))
    return(X)
  }
  beta.null <- vector(mode = "list", length = nrepet)
  beta.null_BC <- as.data.frame(matrix(nrow = (((1 + length(site[[1]])) * length(site[[1]])) / 2) - length(site[[1]]), ncol = nrepet))
  beta.null_BC.BAL <- as.data.frame(matrix(nrow = (((1 + length(site[[1]])) * length(site[[1]])) / 2) - length(site[[1]]), ncol = nrepet))
  beta.null_BC.GRA <- as.data.frame(matrix(nrow = (((1 + length(site[[1]])) * length(site[[1]])) / 2) - length(site[[1]]), ncol = nrepet))
  if (method == "BR") {
    beta.obs <- betapart::beta.pair.abund(taxa[taxa$cd_site %in% site[[1]], -(1)], index.family = "bray")
    beta.obs_BC <- reshape::melt(as.matrix(beta.obs$beta.bray))[reshape::melt(upper.tri(as.matrix(beta.obs$beta.bray), diag = FALSE))$value, ]
    beta.obs_BC.BAL <- reshape::melt(as.matrix(beta.obs$beta.bray.bal))[reshape::melt(upper.tri(as.matrix(beta.obs$beta.bray.bal), diag = FALSE))$value, ]
    beta.obs_BC.GRA <- reshape::melt(as.matrix(beta.obs$beta.bray.gra))[reshape::melt(upper.tri(as.matrix(beta.obs$beta.bray.gra), diag = FALSE))$value, ]
    beta.samp <- permatfull(taxa[taxa$cd_site %in% site[[1]], -(1)], fixedmar = "both", mtype = "count", times = nrepet)
    for (j in 1:length(beta.samp$perm)) {
      beta.null[[j]] <- betapart::beta.pair.abund(beta.samp$perm[[j]], index.family = "bray")
      beta.null_BC[, j] <- reshape::melt(as.matrix(beta.null[[j]]$beta.bray))[reshape::melt(upper.tri(as.matrix(beta.null[[j]]$beta.bray), 
                                                                                                      diag = FALSE))$value, ]$value
      beta.null_BC.BAL[, j] <- reshape::melt(as.matrix(beta.null[[j]]$beta.bray.bal))[reshape::melt(upper.tri(as.matrix(beta.null[[j]]$beta.bray.bal), 
                                                                                                              diag = FALSE))$value, ]$value
      beta.null_BC.GRA[, j] <- reshape::melt(as.matrix(beta.null[[j]]$beta.bray.gra))[reshape::melt(upper.tri(as.matrix(beta.null[[j]]$beta.bray.gra), 
                                                                                                              diag = FALSE))$value, ]$value
      rownames(beta.null_BC) <- rownames(beta.null_BC.BAL) <- rownames(beta.null_BC.GRA) <- rownames(reshape::melt(as.matrix(beta.null[[j]]$beta.bray))[reshape::melt(upper.tri(as.matrix(beta.null[[j]]$beta.bray),
                                                                                                                                                                                diag = FALSE))$value, ])
    }
    keep.pairwise <- rownames(beta.obs_BC[beta.obs_BC$X1 <= 4 & beta.obs_BC$X2 >= 18, ])
    beta.obs_BC <- f1(x = beta.obs_BC)
    beta.obs_BC.BAL <- f1(x = beta.obs_BC.BAL)
    beta.obs_BC.GRA <- f1(x = beta.obs_BC.GRA)
    beta.null_BC <- f2(x = beta.null_BC, keep1 = keep.pairwise)
    beta.null_BC.BAL <- f2(x = beta.null_BC.BAL, keep1 = keep.pairwise)
    beta.null_BC.GRA <- f2(x = beta.null_BC.GRA, keep1 = keep.pairwise)
    test.SOR <- as.randtest(sim = beta.null_BC, obs = beta.obs_BC, alter = "greater")
    test.SIM <- as.randtest(sim = beta.null_BC.BAL, obs = beta.obs_BC.BAL, alter = "greater")
    test.SNE <- as.randtest(sim = beta.null_BC.GRA, obs = beta.obs_BC.GRA, alter = "greater")
  }
  if (method == "SOR") {
    taxa1 <- taxa[taxa$cd_site %in% site[[1]], -(1)]
    taxa1 <- ifelse(taxa1 > 0, 1, 0)
    beta.obs <- betapart::beta.pair(taxa1, index.family = "sorensen")
    beta.obs_SOR <- reshape::melt(as.matrix(beta.obs$beta.sor))[reshape::melt(upper.tri(as.matrix(beta.obs$beta.sor), diag = FALSE))$value, ]
    beta.obs_SIM <- reshape::melt(as.matrix(beta.obs$beta.sim))[reshape::melt(upper.tri(as.matrix(beta.obs$beta.sim), diag = FALSE))$value, ]
    beta.obs_SNE <- reshape::melt(as.matrix(beta.obs$beta.sne))[reshape::melt(upper.tri(as.matrix(beta.obs$beta.sne), diag = FALSE))$value, ]
    beta.samp <- permatfull(taxa1, fixedmar = "both", mtype = "prab", times = nrepet)
    for (j in 1:length(beta.samp$perm)) {
      beta.null[[j]] <- betapart::beta.pair(beta.samp$perm[[j]], index.family = "sorensen")
      beta.null_BC[, j] <- reshape::melt(as.matrix(beta.null[[j]]$beta.sor))[reshape::melt(upper.tri(as.matrix(beta.null[[j]]$beta.sor), 
                                                                                                     diag = FALSE))$value, ]$value
      beta.null_BC.BAL[, j] <- reshape::melt(as.matrix(beta.null[[j]]$beta.sim))[reshape::melt(upper.tri(as.matrix(beta.null[[j]]$beta.sim), 
                                                                                                         diag = FALSE))$value, ]$value
      beta.null_BC.GRA[, j] <- reshape::melt(as.matrix(beta.null[[j]]$beta.sne))[reshape::melt(upper.tri(as.matrix(beta.null[[j]]$beta.sne), 
                                                                                                         diag = FALSE))$value, ]$value
      rownames(beta.null_BC) <- rownames(beta.null_BC.BAL) <- rownames(beta.null_BC.GRA) <- rownames(reshape::melt(as.matrix(beta.null[[j]]$beta.sor))[reshape::melt(upper.tri(as.matrix(beta.null[[j]]$beta.sor),
                                                                                                                                                                               diag = FALSE))$value, ])
    }
    keep.pairwise <- rownames(beta.obs_SOR[beta.obs_SOR$X1 <= 4 & beta.obs_SOR$X2 >= 18, ])
    beta.obs_SOR <- f1(x = beta.obs_SOR)
    beta.obs_SIM <- f1(x = beta.obs_SIM)
    beta.obs_SNE <- f1(x = beta.obs_SNE)
    beta.null_BC <- f2(x = beta.null_BC, keep1 = keep.pairwise)
    beta.null_BC.BAL <- f2(x = beta.null_BC.BAL, keep1 = keep.pairwise)
    beta.null_BC.GRA <- f2(x = beta.null_BC.GRA, keep1 = keep.pairwise)
    test.SOR <- as.randtest(sim = beta.null_BC, obs = beta.obs_SOR, alter = "greater")
    test.SIM <- as.randtest(sim = beta.null_BC.BAL, obs = beta.obs_SIM, alter = "greater")
    test.SNE <- as.randtest(sim = beta.null_BC.GRA, obs = beta.obs_SNE, alter = "greater")
  }
  res.beta <- data.frame(obs.SOR = test.SOR$obs, exp.SOR = test.SOR$expvar[[2]], 
                         std.obs.SOR = test.SOR$expvar[[1]], var.SOR = test.SOR$expvar[[3]],
                         p.value.SOR = test.SOR$pvalue,
                         obs.SIM = test.SIM$obs, exp.SIM = test.SIM$expvar[[2]], 
                         std.obs.SIM = test.SIM$expvar[[1]], var.SIM = test.SIM$expvar[[3]],
                         p.value.SIM = test.SIM$pvalue,
                         obs.SNE = test.SNE$obs, exp.SNE = test.SNE$expvar[[2]], 
                         std.obs.SNE = test.SNE$expvar[[1]], var.SNE = test.SNE$expvar[[3]],
                         p.value.SNE = test.SNE$pvalue)
  res.beta$group <- names(site)[1]
  if (method == "BR") {
    colnames(res.beta) <- c("obs.BC", "exp.BC", "std.obs.BC", "var.BC", "p.value.BC",
                            "obs.BC.BAL", "exp.BC.BAL", "std.obs.BC.BAL", "var.BC.BAL", "p.value.BC.BAL",
                            "obs.BC.GRA", "exp.BC.GRA", "std.obs.BC.GRA", "var.BC.GRA", "p.value.BC.GRA", "group")
    rownames(res.beta) <- names(site)[1]
  }
  if (method == "SOR") {
    colnames(res.beta) <- c("obs.SOR", "exp.SOR", "std.obs.SOR", "var.SOR", "p.value.SOR",
                            "obs.SIM", "exp.SIM", "std.obs.SIM", "var.SIM", "p.value.SIM",
                            "obs.SNE", "exp.SNE", "std.obs.SNE", "var.SNE", "p.value.SNE", "group")
    rownames(res.beta) <- names(site)[1]
  }
  return(res.beta)
}


# END OF SCRIPT #
