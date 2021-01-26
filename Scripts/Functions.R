# ------------------------------------------------------------------------------
# Project: Evolutionary shifts in social mating systems drive habitat shifts in 
#          birds
# Title: All functions produced to make the coding simpler in script the script 
#        "Complete_process.R"
# Author: Thilina De Silva et al.
# Date: 25/01/2021
# ------------------------------------------------------------------------------

# R packages needed ------------------------------------------------------------
library(phytools)
library(diversitree)
library(scales)
library(hisse)
# ------------------------------------------------------------------------------


# Function for ML reconstructions ----------------------------------------------
# trees: names list of trees to be used
# char_table: data.frame with the characters to be reconstructed. 0, 1, and 0.5 = unknown
# mdir: main directory for results, will be created inside function
# colors: colors for two options of reconstructions for each character in the table
# legends: legend for two options of reconstructions for each character in the table

ml_rec <- function(trees, char_table, mdir, colors, legends) {
  dir.create(mdir)
  
  ## subdirectories
  treenames <- names(trees)
  idir <- paste0(mdir, "/", treenames)
  
  ## names of characters
  char_names <- colnames(char_table)
  
  # running all analyses for all trees in a loop
  results_all <- lapply(1:length(trees), function(w) {
    dir.create(idir[w])
    dir.create(paste0(idir[w], "/Figures"))
    wt <- trees[[w]]
    xcor <- -0.9 
    ycor <- 0.8 
    lo <- ifelse(w == 1, 0.006, 0.02)
    to <- ifelse(w == 1, 0.003, 0.01)
    
    # for each tree: running all analyses in a loop (sapply)
    results <- lapply(1:ncol(char_table), function(i) {
      unkn <- rownames(char_table)[char_table[, i] == 0.5] # find unknowns (0.5 in this case)
      
      character <- setNames(char_table[, i], rownames(char_table)) # setting names to character
      character <- ifelse(character == 0.5, 0, character) # changing unknown by 0 
      character <- as.factor(character) # defining levels of analyses, as.factor
      
      # defining matrix with unknown data
      mat_unk <- to.matrix(character[wt$tip.label], levels(character)) 
      nunkn <- rownames(mat_unk[rownames(mat_unk) %in% unkn, ]) 
      for (j in 1:length(nunkn)) {mat_unk[nunkn[j], ] <- c(0.5, 0.5)} 
      
      write.csv(mat_unk, paste0(idir[w], "/", char_names[i], "_initial_matrix.csv"))
      
      #ancestral reconstruction models ER and ARD
      fitER <- rerootingMethod(wt, mat_unk, model = "ER", tips = TRUE)
      fitARD <- rerootingMethod(wt, mat_unk, model = "ARD", tips = TRUE)
      
      # gathering results
      ## ER
      loglikER <- fitER$loglik # log likelihood
      aicER <- -2 * loglikER + 1 * 2 # AIC calculated as in AIC.ace
      write.csv(fitER$Q, paste0(idir[w], "/", char_names[i], "_Q_ML_ER.csv")) # Q matrix
      write.csv(fitER$marginal.anc, paste0(idir[w], "/", char_names[i], 
                                           "_reconstructions_ML_ER.csv")) # marg anc recons
      
      ## ARD
      loglikARD <- fitARD$loglik # log likelihood
      aicARD <- -2 * loglikARD + 1 * 2 # AIC calculated as in AIC.ace
      write.csv(fitARD$Q, paste0(idir[w], "/", char_names[i], "_Q_ML_ARD.csv")) # Q matrix
      write.csv(fitARD$marginal.anc, paste0(idir[w], "/", char_names[i], 
                                            "_reconstructions_ML_ARD.csv")) # marg anc recons
      
      ## AIC
      loglik <- setNames(c(loglikER, loglikARD), c("loglik_ER", "loglik_ARD"))
      aic <- setNames(c(aicER, aicARD), c("AIC_ER", "AIC_ARD"))
      waic <- aic.w(aic)
      names(waic) <- paste0("W_", names(waic))
      
      
      # plotting tree
      cols <- setNames(colors[[i]], levels(character)) # defining colors 
      
      ## saving the tree in pdf ER
      pdf(paste0(idir[w], "/Figures/", char_names[i], "_recons_ML_ER.pdf"), 
          width = 7, height = 20) # open pdf
      plot(wt, no.margin = T, edge.width = 1, label.offset = lo) # plot tree
      tiplabels(pie = mat_unk, piecol = cols, cex = 0.4, 
                offset = to) # add tip labels
      nodelabels(pie = fitER$marginal[as.character(1:wt$Nnode +length(wt$tip)),], 
                 piecol = cols, cex = 0.4) # add node labels
      xleg <- xcor * par()$usr[1] # x legend coordinates
      yleg <- ycor * par()$usr[4] # y legend coordinates
      add.simmap.legend(colors = cols, prompt = FALSE, leg = legends[[i]], 
                        x = xleg, y = yleg, fsize = 1) # add legend
      dev.off() # close pdf
      
      ## saving the tree in pdf ARD
      pdf(paste0(idir[w], "/Figures/", char_names[i], "_recons_ML_ARD.pdf"), 
          width = 7, height = 20) # open pdf
      plot(wt, no.margin = T, edge.width = 1, label.offset = lo) # plot tree
      tiplabels(pie = mat_unk, piecol = cols, cex = 0.4, 
                offset = to) # add tip labels
      nodelabels(pie = fitARD$marginal[as.character(1:wt$Nnode +length(wt$tip)),], 
                 piecol = cols, cex = 0.4) # add node labels
      xleg <- xcor * par()$usr[1] # x legend coordinates
      yleg <- ycor * par()$usr[4] # y legend coordinates
      add.simmap.legend(colors = cols, prompt = FALSE, leg = legends[[i]], 
                        x = xleg, y = yleg, fsize = 1) # add legend
      dev.off() # close pdf
      
      cat("\tCharacter", i, "of", ncol(X), "finished\n") # indicator of advance
      
      return(list(ER_rec = fitER, ARD_rec = fitARD, loglik = loglik,
                  AIC = aic, WAIC = waic)) 
    })
    
    # log likelihood results
    idirname <- paste0(idir[w], "/Model_comparison_all_characters_ML.csv")
    LogLik <- do.call(rbind, lapply(results, function(z) {z$loglik}))
    AICs <- do.call(rbind, lapply(results, function(z) {z$AIC}))
    WAICs <- do.call(rbind, lapply(results, function(z) {z$WAIC}))
    
    # table for all characters 
    table_results <- data.frame(Character = char_names, LogLik, AICs, WAICs) 
    write.csv(table_results, idirname, row.names = FALSE) 
    
    names(results) <- char_names
    
    cat("Tree", w, "of", length(trees), "finished\n") # indicator of advance
    
    return(results)
  })
  
  names(results_all) <- treenames
  
  return(results_all)
}
# ------------------------------------------------------------------------------


# Function for SCM reconstructions ---------------------------------------------
# trees: names list of trees to be used
# char_table: data.frame with the characters to be reconstructed. 0, 1, and 0.5 = unknown
# mdir: main directory for results, will be created inside function
# colors: colors for two options of reconstructions for each character in the table
# legends: legend for two options of reconstructions for each character in the table
# transition_legend: legends to show transition colors

scm_rec <- function(trees, char_table, mdir, colors, legends, transition_legend) {
  dir.create(mdir)
  
  ## subdirectories
  treenames <- names(trees)
  idir <- paste0(mdir, "/", treenames)
  
  ## names of characters
  char_names <- colnames(char_table)
  
  # running all analyses for all trees in a loop
  results_all <- lapply(1:length(trees), function(w) {
    dir.create(idir[w])
    dir.create(paste0(idir[w], "/Figures"))
    wt <- weaver_tree[[w]]
    xcor <- -0.9 
    ycor <- 0.8 
    lo <- ifelse(w == 1, 0.006, 0.02)
    to <- ifelse(w == 1, 0.003, 0.01)
    xl <- ifelse(w == 1, 0.13, 0.21)
    
    # for each tree: running all analyses in a loop (sapply)
    results <- lapply(1:ncol(char_table), function(i) {
      unkn <- rownames(char_table)[char_table[, i] == 0.5] # find unknowns (0.5 in this case)
      
      character <- setNames(char_table[, i], rownames(char_table)) # setting names to character
      character <- ifelse(character == 0.5, 0, character) # changing unknown by 0 only for defining levels
      character <- as.factor(character) # defining levels of analyses, as.factor
      character1 <- as.factor(setNames(ifelse(character == "1", legends[[i]][2], 
                                              legends[[i]][1]), rownames(char_table))) # levels with names
      
      # defining matrix with unknown data
      mat_unk <- to.matrix(character[wt$tip.label], levels(character)) #  matrix containing the prior probability 
      nunkn <- rownames(mat_unk[rownames(mat_unk) %in% unkn, ]) # finding unknowns in matrix
      for (j in 1:length(nunkn)) {mat_unk[nunkn[j], ] <- c(0.5, 0.5)} # replacing values for unknowns by 0.5 - 0.5
      
      write.csv(mat_unk, paste0(idir[w], "/", char_names[i], "_initial_matrix.csv")) # matrix for reconstructions
      
      #ancestral reconstruction models ER and ARD
      fitER <- make.simmap(wt, mat_unk, nsim = 1000, model = "ER")
      post_probER <- describe.simmap(fitER, plot = FALSE) # posterior probabilities
      objER <- densityMap(fitER, states = levels(character1), plot = FALSE) # density of distinct predictions
      n <- length(objER$cols)
      objER$cols[1:n] <- colorRampPalette(colors[[i]], space = "Lab")(n)
      objER$states <- objER$states[2:1]
      
      fitARD <- make.simmap(wt, mat_unk, nsim = 1000, model = "ARD")
      post_probARD <- describe.simmap(fitARD, plot = FALSE) # posterior probabilities
      objARD <- densityMap(fitARD, states = levels(character1), plot = FALSE) # density of distinct predictions
      objARD$cols[1:n] <- colorRampPalette(colors[[i]], space = "Lab")(n)
      objARD$states <- objARD$states[2:1]
      
      # gathering results
      ## ER
      sink(file = paste0(idir[w], "/", char_names[i], "_summary_SM_ER.txt"))
      print(post_probER); sink()
      write.csv(rbind(post_probER$tips, post_probER$ace), 
                paste0(idir[w], "/", char_names[i], "_reconstructions_SM_ER.csv")) # marg anc recons
      
      ## ARD
      sink(file = paste0(idir[w], "/", char_names[i], "_summary_SM_ARD.txt"))
      print(post_probARD); sink()
      write.csv(rbind(post_probARD$tips, post_probARD$ace), 
                paste0(idir[w], "/", char_names[i], "_reconstructions_SM_ARD.csv")) # marg anc recons
      
      
      # plotting tree
      cols <- setNames(colors[[i]], levels(character)) # defining colors using previous list
      cols1 <- setNames(colors[[i]], levels(character1)) # defining colors using previous list
      
      ## saving the tree in pdf
      ## ER
      pdf(paste0(idir[w], "/Figures/", char_names[i], "_recons_SM_ER.pdf"), 
          width = 7, height = 20) # open pdf
      plot(wt, no.margin = T, edge.width = 1, label.offset = lo) # plot tree
      nodelabels(node = 1:post_probER$tree[[1]]$Nnode + Ntip(post_probER$tree[[1]]), # add node labels
                 pie = post_probER$ace, piecol = cols, cex = 0.4)
      tiplabels(pie = mat_unk, piecol = cols, cex = 0.4, offset = to)# add tip labels
      xleg <- xcor * par()$usr[1] # x legend coordinates
      yleg <- ycor * par()$usr[4] # y legend coordinates
      add.simmap.legend(colors = cols, prompt = FALSE, leg = legends[[i]], # add legend
                        x = xleg, y = yleg, fsize = 1)
      dev.off() # close pdf
      
      pdf(paste0(idir[w], "/Figures/", char_names[i], "_recons_density_SM_ER.pdf"),
          width = 7, height = 20) # open pdf
      plot(objER, no.margin = TRUE, edge.width = 1, label.offset = lo) # plot tree density map
      dev.off() # close pdf
      
      pdf(paste0(idir[w], "/Figures/", char_names[i], "_recons_changes_SM_ER.pdf"), 
          width = 7, height = 20) # open pdf
      dotTree(wt, character1, colors = "transparent", legend = FALSE)
      tiplabels(pie = mat_unk, piecol = cols, cex = 0.4, 
                offset = ifelse(w == 1, lo + lo/5, lo/1.7))# add tip labels
      add.simmap.legend(x = 0, y = -4, colors = cols1, prompt = FALSE, 
                        vertical = FALSE)
      nulo <- sapply(fitER, markChanges, sapply(cols, make.transparent, 0.1))
      add.simmap.legend(colors = sapply(setNames(cols1[2:1], transition_legend[[i]]), 
                                        make.transparent, 0.3), 
                        prompt = FALSE, x = xl, y = -4, vertical = FALSE)
      dev.off() # close pdf
      
      ## ARD
      pdf(paste0(idir[w], "/Figures/", char_names[i], "_recons_SM_ARD.pdf"), 
          width = 7, height = 20) # open pdf
      plot(wt, no.margin = T, edge.width = 1, label.offset = lo) # plot tree
      nodelabels(node = 1:post_probARD$tree[[1]]$Nnode + Ntip(post_probARD$tree[[1]]), # add node labels
                 pie = post_probARD$ace, piecol = cols, cex = 0.4)
      tiplabels(pie = mat_unk, piecol = cols, cex = 0.4, offset = to)# add tip labels
      xleg <- xcor * par()$usr[1] # x legend coordinates
      yleg <- ycor * par()$usr[4] # y legend coordinates
      add.simmap.legend(colors = cols, prompt = FALSE, leg = legends[[i]], # add legend
                        x = xleg, y = yleg, fsize = 1)
      dev.off() # close pdf
      
      pdf(paste0(idir[w], "/Figures/", char_names[i], "_recons_density_SM_ARD.pdf"),
          width = 7, height = 20) # open pdf
      plot(objARD, no.margin = TRUE, edge.width = 1, label.offset = lo) # plot tree density map
      dev.off() # close pdf
      
      pdf(paste0(idir[w], "/Figures/", char_names[i], "_recons_changes_SM_ARD.pdf"), 
          width = 7, height = 20) # open pdf
      dotTree(wt, character1, colors = "transparent", legend = FALSE)
      tiplabels(pie = mat_unk, piecol = cols, cex = 0.4, 
                offset = ifelse(w == 1, lo + lo/5, lo/1.7))# add tip labels
      add.simmap.legend(x = 0, y = -4, colors = cols1, prompt = FALSE, 
                        vertical = FALSE)
      nulo <- sapply(fitARD, markChanges, sapply(cols, make.transparent, 0.1))
      add.simmap.legend(colors = sapply(setNames(cols1[2:1], transition_legend[[i]]), 
                                        make.transparent, 0.3), 
                        prompt = FALSE, x = xl, y = -4, vertical = FALSE)
      dev.off() # close pdf
      
      
      cat("\n\tCharacter", i, "of", ncol(X), "finished\n\n") # indicator of advance
      
      return(list(ER_rec = fitER, ARD_rec = fitARD, pProER = post_probER, 
                  pProARD = post_probARD, DensMapER = objER, DensMapARD = objARD)) # returning log lik to be saved later for all characters
    })
    
    names(results) <- char_names
    
    cat("Tree", w, "of", length(weaver_tree), "finished\n") # indicator of advance
    
    return(results)
  })
  
  names(results_all) <- treenames
  
  return(results_all)
}
# ------------------------------------------------------------------------------


# Function for performing BiSSE analyses ---------------------------------------
# tree: an ultrametric phylogenetic tree
# char_table: data.frame with the characters to be reconstructed. 0, 1, and 0.5 = unknown
# sampling_f: proportion of extant species in state 0 and 1 that are included 
#             in the phylogeny for each character. If null, sampling fraction 
#             method wont be used
# mdir: main directory for results, will be created inside function
# colors: colors for two options of reconstructions for each character in the table
  
bisse_eval <- function(tree, char_table, sampling_f = NULL, mdir, colors) {
  dir.create(mdir)

  ## names of characters
  char_names <- colnames(char_table)
  
  # running the analyses for all characters of interest
  dresults <- lapply(1:length(char_names), function(x) {
  # getting tip states
  char <- char_table[, x]
  char <- ifelse(char == 0.5, NA, char)
  names(char) <- row.names(char_table)
  char <- char[order(match(names(char), tree$tip.label))]
  
  #####
  # making initial function for BISSE analyses
  char_bisse <- make.bisse(tree, char)
  
  # initial parameters
  p <- starting.point.bisse(tree)
  
  # ML search
  mle_p <- find.mle(char_bisse, p)
  
  # test the hypothesis that the speciation rates are different
  ## speciation rates to be equal across character states
  char_bisse_1 <- constrain(char_bisse, lambda1 ~ lambda0) ###############
  
  ## start the ML search again (dropping the lambda1)
  mle_p_1 <- find.mle(char_bisse_1, p[argnames(char_bisse_1)])
  
  # comparing results
  ## anova
  anovaml <- anova(mle_p, equal.l = mle_p_1)
  
  ## MCMC analyses for comparison
  prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))
  
  set.seed(1)
  tmp <- mcmc(char_bisse, mle_p$par, nsteps = 100, prior = prior, lower = 0,
              w = rep(1, 6), print.every = 0)
  w <- diff(sapply(tmp[2:7], range))
  
  ### chain for 10,000 steps
  samples <- mcmc(char_bisse, mle_p$par, nsteps = 10000,  w = w, prior = prior, 
                  lower = 0, print.every = 0)
  
  # plotting initial comparison
  ## folder and colors
  dfdir <- paste0(ddir, "/Figures_", char_names[x])
  dir.create(dfdir)
  col <- colores[[x]]
  
  ## diversification rates
  divs <- cbind(d0 = apply(samples[c("mu0", "lambda0")], 1, diff), 
                d1 = apply(samples[c("mu1", "lambda1")], 1, diff))
  
  ## pdf plot
  pdf(paste0(dfdir, "/", char_names[x], "_MCMC.pdf"), 
      width = 7, height = 1.5) # open pdf
  par(mar = c(4.1, 4.1, 0.3, 0.3), mfrow = c(1, 4))
  par(cex = 0.55)
  
  profiles.plot(samples[c("lambda0", "lambda1")], col.line = col, las = 1, 
                xlab = "Speciation rate estimate", ylab = "Probability density")
  legend("topright", c("lambda0", "lambda1"), border = col, cex = 0.7, 
         fill = alpha(col, 0.5), bty = "n")
  
  profiles.plot(samples[c("mu0", "mu1")], col.line = col, las = 1, 
                xlab = "Extintion rate estimate", ylab = "")
  legend("topright", c("mu0", "mu1"), border = col, cex = 0.7, 
         fill = alpha(col, 0.5), bty = "n")
  
  profiles.plot(samples[c("q01", "q10")], col.line = col, las = 1, 
                xlab = "Transition rate estimate", ylab = "")
  legend("topright", c("q01", "q10"), border = col, cex = 0.7, 
         fill = alpha(col, 0.5), bty = "n")
  
  profiles.plot(divs, col.line = col, las = 1, xlab = "Diversification rate estimate", 
                ylab = "")
  legend("topright", c("d0", "d1"), border = col, fill = alpha(col, 0.5), 
         cex = 0.7, bty = "n")
  
  dev.off()
  
  if (is.null(sampling_f)) {
    # saving results
    save(p, mle_p, mle_p_1, anovaml, samples, divs, 
         file = paste0(ddir, "/", char_names[x], "_diver_results.RData"))
  } else {
    #####
    # incomplete taxonomic sampling
    # calculate what the sampling fraction is for this tree
    char_bisse_p <- make.bisse(tree, char, sampling.f = sampling_f[[x]])
    
    #This can then be optimized, as before:
    p_2 <- starting.point.bisse(tree)
    mle_p_2 <- find.mle(char_bisse_p, p_2)
    
    # test the hypothesis that the speciation rates are different
    ## speciation rates to be equal across character states
    char_bisse_21 <- constrain(char_bisse_p, lambda1 ~ lambda0)
    
    ## start the ML search again (dropping the lambda1)
    mle_p_21 <- find.mle(char_bisse_21, p[argnames(char_bisse_21)])
    
    # comparing results
    ## anova
    anovaml_2 <- anova(mle_p_2, equal.l = mle_p_21)
    
    ## MCMC analyses for comparison
    prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))
    
    set.seed(1)
    tmp <- mcmc(char_bisse_p, mle_p_2$par, nsteps = 100, prior = prior, lower = 0,
                w = rep(1, 6), print.every = 0)
    w <- diff(sapply(tmp[2:7], range))
    
    ### chain for 10,000 steps
    samples_2 <- mcmc(char_bisse_p, mle_p_2$par, nsteps = 10000,  w = w, 
                      prior = prior, lower = 0, print.every = 0)
    
    # plotting initial comparison
    ## diversification rates
    divs_2 <- cbind(d0 = apply(samples_2[c("mu0", "lambda0")], 1, diff), 
                    d1 = apply(samples_2[c("mu1", "lambda1")], 1, diff))
    
    ## pdf plot
    pdf(paste0(dfdir, "/", char_names[x], "_MCMC_excluded_taxa.pdf"), 
        width = 7, height = 1.5) # open pdf
    par(mar = c(4.1, 4.1, 0.3, 0.3), mfrow = c(1, 4))
    par(cex = 0.55)
    
    profiles.plot(samples_2[c("lambda0", "lambda1")], col.line = col, las = 1, 
                  xlab = "Speciation rate estimate", ylab = "Probability density")
    legend("topright", c("lambda0", "lambda1"), border = col, cex = 0.7, 
           fill = alpha(col, 0.5), bty = "n")
    
    profiles.plot(samples_2[c("mu0", "mu1")], col.line = col, las = 1, 
                  xlab = "Extintion rate estimate", ylab = "")
    legend("topright", c("mu0", "mu1"), border = col, cex = 0.7, 
           fill = alpha(col, 0.5), bty = "n")
    
    profiles.plot(samples_2[c("q01", "q10")], col.line = col, las = 1, 
                  xlab = "Transition rate estimate", ylab = "")
    legend("topright", c("q01", "q10"), border = col, cex = 0.7, 
           fill = alpha(col, 0.5), bty = "n")
    
    profiles.plot(divs_2, col.line = col, las = 1, 
                  xlab = "Diversification rate estimate", ylab = "")
    legend("topright", c("d0", "d1"), border = col, fill = alpha(col, 0.5), 
           cex = 0.7, bty = "n")
    
    dev.off()
    
    
    # saving results
    save(p, mle_p, mle_p_1, anovaml, samples, divs, 
         p_2, mle_p_2, mle_p_21, anovaml_2, samples_2, divs_2, 
         file = paste0(ddir, "/", char_names[x], "_diver_results.RData")) 
  }
  
  message(x, " of ", length(char_names), " characters")
})
}


# Function for performing HiSSE analyses ---------------------------------------
# tree: an ultrametric phylogenetic tree
# char_table: data.frame with the characters to be reconstructed. 0, 1, and 0.5 = unknown
# sampling_f: proportion of extant species in state 0 and 1 that are included 
#             in the phylogeny for each character.
# broots: a list of vectors (for each character) indicating fixed root state 
#         probabilities for the two states. If NULL, uncertain.
# hroots: a list of vectors (for each character) indicating fixed root state 
#         probabilities for four states (hiden model) See details in help for the 
#         function hisse. If NULL, uncertain.
# mdir: main directory for results, will be created inside function

hisse_eval <- function(tree, char_table, sampling_f, broots, hroots, mdir) {
  dir.create(mdir)
  
  ## names of characters
  char_names <- colnames(char_table)
  
  dresults <- lapply(1:length(char_names), function(x) {
    # getting tip states
    char <- char_table[, char_names[x]]
    char <- ifelse(char == 0.5, 2, char)
    names(char) <- row.names(char_table)
    char <- char[order(match(names(char), tree$tip.label))]
    dat <- data.frame(Species = tree$tip.label, Character = char)
    rownames(dat) <- NULL
    
    # bisse
    ## transition matrices
    trans_bisse <- TransMatMaker(hidden.states = FALSE)
    trans_bisseet <- trans_bisse
    trans_bisseet[!is.na(trans_bisseet)] <- 1
    
    ## turnovers
    turnbn <- c(1, 1, 0, 0)
    extbn <- c(1, 1, 0, 0)
    
    turnbf <- c(1, 2, 0, 0)
    extbf <- c(1, 2, 0, 0)
    
    ## bisse null
    bn <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = FALSE, 
                turnover.anc = turnbn, eps.anc = extbn, trans.rate = trans_bisse, 
                root.p = broots[[x]])
    
    ## bisse null ET
    bnet <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = FALSE, 
                  turnover.anc = turnbn, eps.anc = extbn, trans.rate = trans_bisseet, 
                  root.p = broots[[x]])
    
    ## bisse full
    bf <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = FALSE, 
                turnover.anc = turnbf, eps.anc = extbf, trans.rate = trans_bisse, 
                root.p = broots[[x]])
    
    ## bisse full ET
    bfet <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = FALSE, 
                  turnover.anc = turnbf, eps.anc = extbf, trans.rate = trans_bisseet, 
                  root.p = broots[[x]])
    
    # hisse 
    ## hisse null
    ### transition rates
    trans_hisse <- TransMatMaker(hidden.states = TRUE)
    trans_hisse <- ParDrop(trans_hisse, c(3, 5, 8, 10))
    trans_hisseet <- trans_hisse
    
    ### turnover
    turnhn <- c(1, 1, 2, 2)
    exthn <- c(1, 1, 2, 2)
    
    ### null ET
    trans_hisseet[!is.na(trans_hisseet) & !trans_hisseet == 0] <- 1
    
    hnet <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = TRUE, 
                  turnover.anc = turnhn, eps.anc = exthn, trans.rate = trans_hisseet, 
                  root.p = hroots[[x]])
    
    ### null no ET
    #### transitions from 0->1 governed by a single rate
    to_change <- cbind(c(1, 3), c(2, 4))
    trans_hisse[to_change] <- 1
    
    #### transitions from 1->0 governed by a single rate:
    to_change <- cbind(c(2, 4), c(1, 3))
    trans_hisse[to_change] <- 2
    
    #### transitions between the hidden state to be a single rate 
    to_change <- cbind(c(1, 3, 2, 4), c(3, 1, 4, 2))
    trans_hisse[to_change] <- 3
    
    hn <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = TRUE, 
                turnover.anc = turnhn, eps.anc = exthn, trans.rate = trans_hisse, 
                root.p = hroots[[x]])
    
    
    ## hisse hidden
    ### transition rates
    trans_hisse <- TransMatMaker(hidden.states = TRUE)
    trans_hisse <- ParDrop(trans_hisse, c(2, 3, 5, 7, 8, 9, 10, 12))
    trans_hisseet <- trans_hisse
    trans_hisseet[!is.na(trans_hisseet) & !trans_hisseet == 0] <- 1
    
    ### turnover
    turnhh2 <- c(1, 2, 0, 3)
    exthh2 <- c(1, 2, 0, 3)
    
    ### hidden ET
    hh2et <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = TRUE, 
                   turnover.anc = turnhh2, eps.anc = exthh2, trans.rate = trans_hisseet, 
                   root.p = hroots[[x]])
    
    ### hidden no ET
    hh2 <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = TRUE, 
                 turnover.anc = turnhh2, eps.anc = exthh2, trans.rate = trans_hisse, 
                 root.p = hroots[[x]])
    
    
    ## hisse full
    ### transition rates
    trans_hisse <- TransMatMaker(hidden.states = TRUE)
    trans_hisse <- ParDrop(trans_hisse, c(3, 5, 8, 10))
    trans_hisseet <- trans_hisse
    trans_hisseet[!is.na(trans_hisseet) & !trans_hisseet == 0] <- 1
    
    ### turnover
    turnhf <- c(1, 2, 3, 4)
    exthf <- c(1, 2, 3, 4)
    
    ### hisse full ET
    hf <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = TRUE, 
                turnover.anc = turnhf, eps.anc = exthf, trans.rate = trans_hisseet, 
                root.p = hroots[[x]])
    
    ### hisse full no ET
    hfet <- hisse(tree, dat, f = sampling_p[[x]], hidden.states = TRUE, 
                  turnover.anc = turnhf, eps.anc = exthf, trans.rate = trans_hisse, 
                  root.p = hroots[[x]])
    
    
    
    # comparison matrix
    turn <- c(paste(turnbn, collapse = " "), paste(turnbn, collapse = " "), 
              paste(turnbf, collapse = " "), paste(turnbf, collapse = " "),
              paste(turnhn, collapse = " "), paste(turnhn, collapse = " "), 
              paste(turnhh2, collapse = " "), paste(turnhh2, collapse = " "), 
              paste(turnhf, collapse = " "), paste(turnhf, collapse = " "))
    
    cmat <- data.frame(Turnover = turn, 
                       rbind(unlist(bn[1:3]), unlist(bnet[1:3]),
                             unlist(bf[1:3]), unlist(bfet[1:3]), 
                             unlist(hn[1:3]), unlist(hnet[1:3]),
                             unlist(hh2[1:3]), unlist(hh2et[1:3]), 
                             unlist(hf[1:3]), unlist(hfet[1:3])))
    
    rownames(cmat) <- c("BiSSE_null", "BiSSE_null_ET", 
                        "BiSSE_full", "BiSSE_full_ET", 
                        "HiSSE_null", "HiSSE_null_ET",
                        "HiSSE_hidden", "HiSSE_hidden_ET", 
                        "HiSSE_full", "HiSSE_full_ET")
    
    cmat <- cmat[order(cmat[, 4]), ]
    
    write.csv(cmat, paste0("HiSSE_diversification/", char_names[x], "_model_comparison.csv"))
    
    return(list(BiSSE_null = bn, BiSSE_null_ET = bnet, 
                BiSSE_full = bf, BiSSE_full_ET = bfet,
                HiSSE_null = hn, HiSSE_null_ET = hnet,
                HiSSE_hidden = hh2, HiSSE_hidden_ET = hh2et,
                HiSSE_full = hf, HiSSE_full_ET = hfet, 
                model_comparison = cmat))
  })
  
  names(dresults) <- char_names
  
  return(dresults)
}
# ------------------------------------------------------------------------------


# Function to perform Mantel test ----------------------------------------------
# tree: an ultrametric phylogenetic tree
# traits: matrix, first column = name of tips, second column = character value 
#         of trait
# runs: number of permutations to be performed
# unknown: character that was used in traits to represent unknown state
# sqrtPhylo: whether to use square root transformation of phylogenetic distance,
#            default = FALSE

mantel_runs <- function(tree, traits, runs = 999, unknown, sqrtPhylo = FALSE) {
  # excluding unknown 
  if (unknown %in% traits[,2]) {
    missingData <- as.character(traits[traits[, 2] == unknown, 1])
    tree <- drop.tip(tree, missingData)
    traits <- traits[traits[, 2] != unknown, ]
  }
  
  # phylogenetic distance
  phylo.dist <- cophenetic(tree)
  if (sqrtPhylo) {
    phylo.dist <- sqrt(phylo.dist)
  }

  # preparing distance matrix for mantel test
  t <- strsplit(paste(traits[, 2], collapse = ""), "")[[1]]
  s <- sort(unique(t))
  m <- sapply(s, function(x){grepl(x, traits[, 2])})
  rownames(m) <- traits[, 1]
  bc.dist <- vegdist(m)
  
  # mantel test
  res.mantel <- mantel(phylo.dist, bc.dist, permutations = runs)

  # preparing results
  pval1.NULL <- (sum(ifelse(res.mantel$perm >= res.mantel$statistic, 1, 
                            0)) + 1) / (runs + 1)
  pval2.NULL <- (sum(ifelse(res.mantel$perm <= res.mantel$statistic, 1, 
                            0)) + 1) / (runs + 1)
  pval3.NULL <- (sum(ifelse(abs(res.mantel$perm) >= abs(res.mantel$statistic), 1, 
                            0)) + 1) / (runs + 1)
  r.Mantel <- res.mantel$statistic
  
  RES <- list(perm.NULL = res.mantel$perm, r.Mantel = r.Mantel,  
              pval1.NULL = pval1.NULL, pval2.NULL = pval2.NULL, 
              pval3.NULL = pval3.NULL)
  
  return(RES)
}
# ------------------------------------------------------------------------------
