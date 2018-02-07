# Plot tests

dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

parms<-matrix(c(NA,
                NA),byrow=TRUE,ncol=1,nrow=2)

#RMCriteria::PrometheeI(dados,c(0.3,0.7),c(0,0),parms,FALSE)

#Step 1: Construct the RPrometheeArguments
PromObj <- RPrometheeConstructor(datMat=dados,vecWeights=c(0.3,0.7),vecMaximiz=c(F,T),prefFunction=c(0,0),parms=parms,normalize=FALSE)
res <- RPrometheeI(PromObj)
str(res)

datMatDF <- data.frame(PromObj@datMat)
vecWeightsDF <- data.frame(PromObj@vecWeights)
parmsDF <- data.frame(PromObj@parms)
resDF <- data.frame("PhiPlus" = res@PhiPlus, "PhiMinus" = res@PhiMinus)


phiLabels <- c(rep("PhiPlus", nrow(datMatDF)), rep("PhiMinus", nrow(datMatDF)))
phiNums <- c(resDF[,1], resDF[,2])
alternatives <- c(rep(rownames(datMatDF), 2))
resultsPlot <- data.frame(alternatives, phiLabels, phiNums)
resultsPlot[,2] <- as.factor(resultsPlot[,2])

# Create a dataframe to use as source for both Phi bars in PrometheeI plot
limits <- data.frame(
  class = c("PhiPlus", "PhiPlus", "PhiMinus", "PhiMinus"),
  boundaries = c(0.5, 0.5, 0.5, 0.5),
  pos_neg = c("Pos", "Neg", "Pos", "Neg"))

# Change order of factors
limits$class <- factor(limits$class, levels = c("PhiPlus", "PhiMinus"))
limits$pos_neg <- factor(limits$pos_neg, levels = c("Pos", "Neg"))

# Filter results table to exclude Phi flows
resultsPlot[,2] <- factor(resultsPlot[,2],
                              levels = c("PhiPlus", "PhiMinus"))

# Full bars as in Visual-Promethee.
ggplot(limits) +
  geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
           stat = "identity", width = 0.5) +
  geom_point(data = resultsPlot, aes(x = phiLabels, y = phiNums),
             stat = "identity") +
  geom_line(data = resultsPlot, aes(x = phiLabels, y = phiNums),
            group = resultsPlot[,1], stat = "identity") +
  geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
            label = sprintf("%0.3f",
                            round(resultsPlot$phiNums, digits = 3),
                            position = position_dodge(width = 0.9)),
            hjust = 0, nudge_x = 0.05) +
  scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
  geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
            label = alternatives, hjust = 1, nudge_x = -0.05) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Alternative/Phi")

ggplot(limits) +
  geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
           stat = "identity", width = 0.5) +
  geom_point(data = resultsPlot, aes(x = phiLabels, y = phiNums),
             stat = "identity") +
  geom_line(data = resultsPlot, aes(x = phiLabels, y = phiNums),
            group = resultsPlot[,1], stat = "identity") +
  geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
            label = sprintf("%0.3f",
                            round(resultsPlot$phiNums, digits = 3),
                            position = position_dodge(width = 0.9)),
            hjust = 0, nudge_x = 0.05) +
  scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
  geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
            label = alternatives, hjust = 1, nudge_x = -0.05) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Alternative/Phi")






################################################################################
################################################################################
################################################################################
################################################################################


res <- RPrometheeII(PromObj)

datMatDF <- data.frame(PromObj@datMat)
vecWeightsDF <- data.frame(PromObj@vecWeights)
parmsDF <- data.frame(PromObj@parms)
resDF <- data.frame("Phi" = res@Phi)

phiLabels <- c(rep("Phi", nrow(datMatDF)))
phiNums <- c(resDF[,1])
alternatives <- c(rep(rownames(datMatDF), 2))
resultsPlot <- data.frame(alternatives, phiLabels, phiNums)
resultsPlot[,2] <- as.factor(resultsPlot[,2])


limits <- data.frame(
  class = c("Phi", "Phi"),
  boundaries = c(-1, 1),
  pos_neg = c("Neg", "Pos"))

# Change order of factors
limits$pos_neg <- factor(limits$pos_neg, levels = c("Pos", "Neg"))

resultsPlot[,2] <- factor(resultsPlot[,2], levels = "Phi")

# Full Ranking bar as in Visual-Promethee.
ggplot(limits) +
  geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
           stat = "identity", width = 0.3) +
  geom_point(data = resultsPlot, aes(x = phiLabels, y = phiNums),
             stat = "identity") +
  geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
            label = sprintf("%0.3f",
                            round(resultsPlot$phiNums, digits = 3)),
            hjust = 0, nudge_x = 0.03) +
  scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
  geom_text(data = resultsPlot, aes(x = phiLabels,
                                         y = resultsPlot$phiNums),
            label = resultsPlot$alternatives,
            hjust = 1, nudge_x = -0.03) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Alternative/Phi")


############################################################################
############################################################################
############################################################################
# Walking Weights plot

datMatDF <- data.frame(PromObj@datMat)
vecWeightsDF <- data.frame(PromObj@vecWeights)
parmsDF <- data.frame(PromObj@parms)
resDF <- data.frame("Phi" = res@Phi)

phiLabels <- c(rep("Phi", nrow(datMatDF)))
phiNums <- c(resDF[,1])
alternatives <- c(rep(rownames(datMatDF), 2))
resultsPlot <- data.frame(alternatives, phiLabels, phiNums)
resultsPlot[,2] <- as.factor(resultsPlot[,2])
weightsDF <- setNames(data.frame(c(1:ncol(datMatDF)), vecWeightsDF), c("criterias", "weights"))

plot_a <- ggplot(resultsPlot) +
  geom_bar(aes(x = alternatives, y = phiNums, fill = alternatives),
           stat = "identity") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(aes(x = alternatives, y = phiNums,
                label = sprintf("%0.3f", round(phiNums, digits = 3))),
            vjust = 1, nudge_y = -0.1) +
  labs(x = "Alternatives", y = "Phi")

plot_b <- ggplot(weightsDF) +
  geom_bar(aes(x = as.character(criterias), y = weights), stat = "identity", width = 0.5) +
  geom_text(aes(x = as.character(criterias), y = weights,
                label = sprintf("%0.2f%%", 100*weights),
                vjust = 1)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = "Criterias", y = "Weights")

grid.arrange(plot_a, plot_b, nrow = 2, ncol = 1,
             heights = unit(c(0.7, 0.3), "npc"))


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Promethee III


dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

parms<-matrix(c(NA,
                NA),byrow=TRUE,ncol=1,nrow=2)

PromObj <- RPrometheeConstructor(datMat=dados, vecWeights=c(0.3,0.7), vecMaximiz=c(F,T), prefFunction=c(0,0), parms=parms, normalize=FALSE, alphaVector=c(1,1,1))

res <- RPrometheeIII(PromObj)
str(res)


Plus       <- res@PhiPlus
Minus      <- res@PhiMinus
limInf     <- res@limInf
limSup     <- res@limSup

    # Create dataframes
resDF <- data.frame("PhiPlus" = Plus, "PhiMinus" = Minus, "limInf" = limInf, "limSup" = limSup)

    # Create a dataframe with results from RPrometheeI and arguments
    phiLabels <- c(rep("PhiPlus", nrow(resDF)), rep("PhiMinus", nrow(resDF)))
    phiNums <- c(resDF[,1], resDF[,2])
    errorMin <- c(rep(resDF[,3], 2))
    errorMax <- c(rep(resDF[,4], 2))
    alternatives <- c(as.character(rep(1:nrow(resDF),2)))

    resultsPlot <- data.frame(alternatives, phiLabels, phiNums, errorMin, errorMax)
    resultsPlot[,2] <- as.factor(resultsPlot[,2])


    # Create a dataframe to use as source for the plot
    limits <- data.frame(
      class = c("PhiPlus", "PhiPlus", "PhiMinus", "PhiMinus"),
      boundaries = c(0.5, 0.5, 0.5, 0.5),
      pos_neg = c("Pos", "Neg", "Pos", "Neg"))

    # Change order of factors and levels
    limits$class <- factor(limits$class, levels = c("PhiPlus", "PhiMinus"))
    limits$pos_neg <- factor(limits$pos_neg, levels = c("Pos", "Neg"))
    resultsPlot[,2] <- factor(resultsPlot[,2],
                              levels = c("PhiPlus", "PhiMinus"))

    # Partial bars as in Visual-Promethee
    ggplot(limits) +
      geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
               stat = "identity", width = 0.5) +
      geom_point(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                 stat = "identity") +
      geom_line(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                group = resultsPlot[,1], stat = "identity") +
      geom_errorbar(data = resultsPlot,
                    aes(x = phiLabels,
                        ymin = errorMin, ymax = errorMax, width = 0.2),
                    color = alternatives) +
      geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                label = sprintf("%0.3f",
                                round(resultsPlot$phiNums, digits = 3),
                                position = position_dodge(width = 0.9)),
                hjust = 0, nudge_x = 0.05) +
      scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
      geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                label = alternatives, hjust = 1, nudge_x = -0.05) +
      theme(axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank()) +
      labs(y = "Alternative/Phi")


###########################################

    dados<-matrix(c(5.2,-3.5,
                    4.3,-1.2,
                    6.7,-2.0),byrow = T, ncol=2,nrow=3)

    parms<-matrix(c(NA,
                    NA),byrow=TRUE,ncol=1,nrow=2)

    PromObj <- RPrometheeConstructor(datMat=dados, vecWeights=c(0.3,0.7), vecMaximiz=c(F,T), prefFunction=c(0,0), parms=parms, normalize=FALSE, alphaVector=c(1,1,1))

    res <- RPrometheeIII(PromObj)
    str(res)


    Phi       <- res@Phi
    limInf     <- res@limInf
    limSup     <- res@limSup

    # Create dataframes
    resDF <- data.frame("Phi" = Phi, "limInf" = limInf, "limSup" = limSup)

    phiLabels <- c(rep("Phi", nrow(resDF)))
    phiNums <- c(resDF[,1])
    errorMin <- c(rep(resDF[,2]))
    errorMax <- c(rep(resDF[,3]))
    alternatives <- c(as.character(1:nrow(resDF)))

    resultsPlot <- data.frame(alternatives, phiLabels, phiNums, errorMin, errorMax)
    resultsPlot[,2] <- as.factor(resultsPlot[,2])


    # Create a dataframe to use as source for the plot
    limits <- data.frame(
      class = c("Phi", "Phi"),
      boundaries = c(-1, 1),
      pos_neg = c("Neg", "Pos"))

    # Change order of factors
    limits$pos_neg <- factor(limits$pos_neg, levels = c("Pos", "Neg"))

    resultsPlot[,2] <- factor(resultsPlot[,2], levels = "Phi")

    # Full Ranking bar as in Visual-Promethee
    ggplot(limits) +
      geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
               stat = "identity", width = 0.3) +
      geom_point(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                 stat = "identity") +
      geom_text(data = resultsPlot, aes(x = phiLabels, y = phiNums),
                label = sprintf("%0.3f",
                                round(resultsPlot$phiNums, digits = 3)),
                hjust = 0, nudge_x = 0.03) +
      geom_errorbar(data = resultsPlot, aes(x = phiLabels,
                                            ymin = errorMin,
                                            ymax = errorMax),
                    colour = alternatives, width = 0.07, size = 1) +
      scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
      geom_text(data = resultsPlot, aes(x = phiLabels,
                                        y = resultsPlot$phiNums),
                label = resultsPlot$alternatives,
                hjust = 1, nudge_x = -0.03) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank()) +
      labs(y = "Alternative Phi")




