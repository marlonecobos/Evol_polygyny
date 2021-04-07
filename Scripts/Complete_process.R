# ------------------------------------------------------------------------------
# Project: Evolution of social mating systems and nesting sociality shapes 
#          historical habitat shifts in birds
# Title: All R code used to perform analyses.
# Author: Thilina De Silva et al.
# Date: 25/01/2021
# ------------------------------------------------------------------------------


# R packages needed ------------------------------------------------------------
# Install them first, if needed, using function install.packages()
library(phytools)
library(diversitree)
library(scales)
library(hisse)
library(geiger)
library(vegan)
# ------------------------------------------------------------------------------


# Functions --------------------------------------------------------------------
source("https://raw.githubusercontent.com/marlonecobos/Evol_polygyny/main/Scripts/Functions.R")
# ------------------------------------------------------------------------------


# Data -------------------------------------------------------------------------
# read character table
data <- read.csv("Characters/all_characters.csv", row.names = 1) 

# names of characters based on columns
char_names <- colnames(data) 

# read trees
ml_tree <- read.nexus("Trees/raxml108.tre") 
by_tree <- read.nexus("Trees/BEASTMCC.tre")
# ------------------------------------------------------------------------------


# Preparing data for reconstructions -------------------------------------------
# list of the two threes
weaver_tree <- list(ML = ml_tree, BY = by_tree)

# colors (list to be used in function below lapply)
colores <- list(polygyny = c("gray75", "brown2"), 
                coloniality = c("dimgray", "darkorange1"),
                diet = c("#998ec3", "#f1a340"), 
                habitat = c("#E8BE4E", "green4"))

# names for the legends (list to be used below)
legends <- list(c("Monogamy", "Polygyny"), c("Non-colonial", "Colonial"), 
                c("Herbivory", "Faunivory"), c("Open habitat", "Dense forest"))
# ------------------------------------------------------------------------------


# Maximum likelihood reconstructions -------------------------------------------
# directory for results
mdir <- "Ancestral_reconstructions_ML" 

# reconstructing all characters with all trees
results_all <- ml_rec(trees = weaver_tree, char_table = data, mdir = mdir, 
                      colors = colores, legends = legends)

# saving results as RData
save(results_all, file = "Ancestral_reconstructions_ML/All_results_ML.RData")
# ------------------------------------------------------------------------------


# Stochastic character mapping reconstructions ---------------------------------
# directory for results
mdir1 <- "Ancestral_reconstructions_SCM"

# names for the legends (list to be used below)
legs1 <- list(c("M -> P", "P -> M"), c("NC -> C", "C -> NC"),
              c("H -> F", "F -> H"), c("OH -> DF", "DF -> OH"))

# reconstructing all characters with all trees
results_allSCM <- scm_rec(trees = weaver_tree, char_table = data, mdir = mdir1, 
                          colors = colores, legends = legends, 
                          transition_legend = legs1)
# saving results as RData
save(results_allSCM, file = "Ancestral_reconstructions_SCM/All_results_SM.RData")
# ------------------------------------------------------------------------------


# Plotting SCM results for polygyny --------------------------------------------
poly <- results_allSCM$BY$Polygyny$DensMapARD

pdf("Ancestral_reconstructions_SCM/BY/Figures/poly_fan.pdf", width = 7, 
    height = 7)
par(cex = 0.5)
plot(poly, type = "fan", no.margin = TRUE, edge.width = 1, label.offset = 0.005)
dev.off()
# ------------------------------------------------------------------------------


# Diversification rates (BiSSE) including sampling fraction --------------------
## fix tolerance for ultrametry
is.ultrametric(by_tree)
by_tree <- force.ultrametric(by_tree)
is.ultrametric(by_tree)

# new directory
ddir <- "Diversification_rates"

# values for sampling.f 
sampling_p <- list(polygyny = c(0.88, 0.96),
                   coloniality = c(0.94, 0.90),
                   herbivory = c(0.90, 0.93),
                   habitat = c(0.93, 0.91))

# running the analyses for all characters of interest
dresults <- bisse_eval(tree = by_tree, char_table = data, sampling_f = sampling_p, 
                       mdir = ddir, colors = colores)
# ------------------------------------------------------------------------------


# Diversification rates (BiSSE) unresolved tip method --------------------------
# adding missing tips to tree
## tree for reference only
by_tree117 <- read.nexus("Trees/117speciesultra.tre") 

## names to be added
names <- by_tree117$tip.label[c(2, 6:8, 11, 56, 71, 100, 106)]

## references to position
tnam <- c("Amblyospiza_albifrons", "Plocepasser_mahali", "Plocepasser_mahali",
          "Plocepasser_mahali", "Pseudonigrita_arnaudi", "Malimbus_racheliae",
          "Ploceus_tricolor", "Ploceus_baglafecht", "multiple")

## additions
### 1
names[1]
which(by_tree$tip.label == "Amblyospiza_albifrons")
cbind(by_tree$edge, by_tree$edge.length)
len <- by_tree$edge.length[213] * 0.999

tree1 <- bind.tip(tree = by_tree, tip.label = names[1], edge.length = len, 
                  where = 1)

cbind(tree1$edge, tree1$edge.length)
tree1$edge.length[1] <- by_tree$edge.length[213] - len
tree1$edge.length[c(2, 3)] <- c(len, len)

### 2
names[2]
which(tree1$tip.label == "Plocepasser_mahali")
cbind(tree1$edge, tree1$edge.length)
len <- tree1$edge.length[10] / 2

tree2 <- bind.tip(tree = tree1, tip.label = names[2], edge.length = len, 
                  where = 5)

cbind(tree2$edge, tree2$edge.length)
tree2$edge.length[10] <- len
tree2$edge.length[c(11, 12)] <- c(len, len)

### 3
names[3]
which(tree2$tip.label == names[2])
cbind(tree2$edge, tree2$edge.length)
len <- tree2$edge.length[12] * 0.999

tree3 <- bind.tip(tree = tree2, tip.label = names[3], edge.length = len, 
                  where = 6)

cbind(tree3$edge, tree3$edge.length)
tree3$edge.length[12] <- tree2$edge.length[12] - len
tree3$edge.length[c(13, 14)] <- c(len, len)

### 4
names[4]
which(tree3$tip.label == "Plocepasser_mahali")
cbind(tree3$edge, tree3$edge.length)
len <- tree3$edge.length[11] * 0.999

tree4 <- bind.tip(tree = tree3, tip.label = names[4], edge.length = len, 
                  where = 5)

cbind(tree4$edge, tree4$edge.length)
tree4$edge.length[11] <- tree3$edge.length[11] - len
tree4$edge.length[c(12, 13)] <- c(len, len)

### 5
names[5]
which(tree4$tip.label == 'Pseudonigrita_arnaudi')
cbind(tree4$edge, tree4$edge.length)
len <- tree4$edge.length[19] / 2

tree5 <- bind.tip(tree = tree4, tip.label = names[5], edge.length = len, 
                  where = 10)

which(tree5$tip.label == 'Pseudonigrita_arnaudi')
cbind(tree5$edge, tree5$edge.length)
tree5$edge.length[19] <- len
tree5$edge.length[c(20, 21)] <- c(len, len)

### 6
names[6]
which(tree5$tip.label == "Malimbus_racheliae")
cbind(tree5$edge, tree5$edge.length)
len <- tree5$edge.length[111] / 2

tree6 <- bind.tip(tree = tree5, tip.label = names[6], edge.length = len, 
                  where = 55)

which(tree6$tip.label == "Malimbus_racheliae")
cbind(tree6$edge, tree6$edge.length)
tree6$edge.length[111] <- len
tree6$edge.length[c(112, 113)] <- c(len, len)

### 7
names[7]
which(tree6$tip.label == "Ploceus_tricolor")
cbind(tree6$edge, tree6$edge.length)
len <- tree6$edge.length[140] / 2

tree7 <- bind.tip(tree = tree6, tip.label = names[7], edge.length = len, 
                  where = 70)

which(tree7$tip.label == "Ploceus_tricolor")
cbind(tree7$edge, tree7$edge.length)
tree7$edge.length[140] <- len
tree7$edge.length[c(141, 142)] <- c(len, len)

### 8
names[8]
which(tree7$tip.label == "Ploceus_baglafecht")
cbind(tree7$edge, tree7$edge.length)
len <- tree7$edge.length[198] / 2

tree8 <- bind.tip(tree = tree7, tip.label = names[8], edge.length = len, 
                  where = 99)

which(tree8$tip.label == "Ploceus_baglafecht")
cbind(tree8$edge, tree8$edge.length)
tree8$edge.length[198] <- len
tree8$edge.length[c(199, 200)] <- c(len, len)

### 9
names[9]
which(tree8$tip.label == "Ploceus_grandis")
cbind(tree8$edge, tree8$edge.length)
len <- tree8$edge.length[209] * 0.99

tree9 <- bind.tip(tree = tree8, tip.label = names[9], edge.length = len, 
                  where = 105)

cbind(tree9$edge, tree9$edge.length)
tree9$edge.length[209] <- tree8$edge.length[209] - len
tree9$edge.length[c(210, 211)] <- c(len, len)

## saving final tree
write.tree(tree9,  file = "Trees/117_beast_Rmodified.tre")

# diversification rate analysis
# read data with additional species
data117 <- read.csv("Characters/all_characters_117.csv", row.names = 1)

# new directory
ddir1 <- "Diversification_rates_117"

# running the analyses for all characters of interest
dresults117 <- bisse_eval(tree = tree9, char_table = data117, mdir = ddir1, 
                          colors = colores)
# ------------------------------------------------------------------------------


# Diversification rates (HiSSE) ------------------------------------------------
# directory
dirr2 <- "HiSSE_diversification"

# values for sampling.f 
sampling_p <- list(polygyny = c(0.88, 0.96),
                   coloniality = c(0.94, 0.90),
                   herbivory = c(0.90, 0.93),
                   habitat = c(0.93, 0.91))

# roots
## BISSE roots
broots <- list(polygyny = c(0, 1),
               coloniality = NULL,
               herbivory = c(1, 0),
               habitat = c(1, 0))

## HISSE roots
hroots <- list(polygyny = c(0, 0.5, 0, 0.5),
               coloniality = NULL,
               herbivory = c(0.5, 0, 0.5, 0),
               habitat = c(0.5, 0, 0.5, 0))

# running the analyses for all characters of interest
dresults2 <- hisse_eval(tree = by_tree, char_table = data, sampling_f = sampling_p, 
                        broots = broots, hroots = hroots, mdir = dirr2)

save(dresults2, file = "HiSSE_diversification/All_results.RData")
# ------------------------------------------------------------------------------


# Mantel test ------------------------------------------------------------------
# read data
data1 <- read.csv("Characters/all_characters.csv") # read csv file
char_names <- colnames(data1)[-1] # names of characters based on columns

# change character states 
data1 <- as.matrix(data1)
data1[data1 == "0.0"] <- "A"
data1[data1 == "1.0"] <- "B"
data1[data1 == "0.5"] <- "AB"

# rates for calculations
rates <- list(
  polygyny = c(as.matrix(read.csv("Ancestral_reconstructions_ML/BY/Polygyny_Q_ML_ARD.csv", 
                                  row.names = 1)))[2:3],
  coloniality = c(as.matrix(read.csv("Ancestral_reconstructions_ML/BY/Coloniality_Q_ML_ARD.csv", 
                                     row.names = 1)))[2:3],
  diet = c(as.matrix(read.csv("Ancestral_reconstructions_ML/BY/Diet_Q_ML_ARD.csv", 
                              row.names = 1)))[2:3],
  habitat = c(as.matrix(read.csv("Ancestral_reconstructions_ML/BY/Habitat_Q_ML_ARD.csv", 
                                 row.names = 1)))[2:3]
)

# running in loop for various characters
mante_results <- lapply(1:length(char_names), function(x) {
  # get character of interest
  traits <- data1[, c(1, x + 1)]
  
  # run the Mantel test
  mtel <- mantel_runs(tree = by_tree, traits = traits, unknown = "AB", 
                      sqrtPhylo = FALSE)
  mtel1 <- mantel_runs(tree = by_tree, traits = traits, unknown = "AB", 
                       sqrtPhylo = TRUE)
  
  return(list(mantel = mtel, mantel_sqrt = mtel1))
})

names(mante_results) <- char_names

# plot
## limits
poly_lim <- range(c(mante_results$Polygyny$mantel$perm.NULL, 
                    mante_results$Polygyny$mantel$r.Mantel))
colo_lim <- range(c(mante_results$Coloniality$mantel$perm.NULL, 
                    mante_results$Coloniality$mantel$r.Mantel))
habi_lim <- range(c(mante_results$Habitat$mantel$perm.NULL, 
                    mante_results$Habitat$mantel$r.Mantel))
diet_lim <- range(c(mante_results$Diet$mantel$perm.NULL, 
                    mante_results$Diet$mantel$r.Mantel))

#pdf("Mantel_test_results/Mantel_test_results.pdf", width = 7, height = 7) # open pdf
jpeg("Mantel_test_results/Mantel_test_results.jpg", width = 16.6, height = 16.6,
     units = "cm", res = 600)

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 0.5, 0.6))

hist(mante_results$Polygyny$mantel$perm.NULL, breaks = 20, xlim = poly_lim, 
     main = "", xlab = "", border = "gray45", col = "gray80")
abline(v = mante_results$Polygyny$mantel$r.Mantel, col = "#20C12F", lwd = 2)
abline(v = mean(mante_results$Polygyny$mantel$perm.NULL), lty = 2, lwd = 1.5)
legend("topright", legend = "a", box.col = "white", bg = "white", inset = 0.001,
       text.font = 2, cex = 1.1)
box(bty = "l")

hist(mante_results$Coloniality$mantel$perm.NULL, breaks = 20, xlim = colo_lim, 
     main = "", xlab = "", ylab = "", border = "gray45", col = "gray80")
abline(v = mante_results$Coloniality$mantel$r.Mantel, col = "#20C12F", lwd = 2)
abline(v = mean(mante_results$Coloniality$mantel$perm.NULL), lty = 2, lwd = 1.5)
legend("topright", legend = "b", box.col = "white", bg = "white", inset = 0.001,
       text.font = 2, cex = 1.1)
box(bty = "l")

hist(mante_results$Habitat$mantel$perm.NULL, breaks = 20, xlim = habi_lim, 
     main = "", xlab = "Mantel r", border = "gray45", col = "gray80")
abline(v = mante_results$Habitat$mantel$r.Mantel, col = "#20C12F", lwd = 2)
abline(v = mean(mante_results$Habitat$mantel$perm.NULL), lty = 2, lwd = 1.5)
legend("topright", legend = "c", box.col = "white", bg = "white", inset = 0.001,
       text.font = 2, cex = 1.1)
box(bty = "l")

hist(mante_results$Diet$mantel$perm.NULL, breaks = 20, xlim = diet_lim, 
     main = "", xlab = "Mantel r", ylab = "", border = "gray45", col = "gray80")
abline(v = mante_results$Diet$mantel$r.Mantel, col = "#20C12F", lwd = 2)
abline(v = mean(mante_results$Diet$mantel$perm.NULL), lty = 2, lwd = 1.5)
legend("topright", legend = "d", box.col = "white", bg = "white", inset = 0.001,
       text.font = 2, cex = 1.1)
box(bty = "l")

legend("right", legend = c("Observed r", "Null mean", "Mantel null"), 
       lty = c(1, 2, NA), pch = c(NA, NA, 22), pt.cex = 2, lwd = c(2, 1.5, NA),
       col = c("#20C12F", "black", "gray45"), pt.bg = c(NA, NA, "gray80"), 
       bty = "n", cex = 0.8, inset = 0.001)

dev.off()

# saving results
dir.create("Mantel_test_results")

save(mante_results, file = "Mantel_test_results/Mantel_results.RData")
# ------------------------------------------------------------------------------
