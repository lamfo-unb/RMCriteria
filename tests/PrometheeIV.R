dados<-matrix(c(5.2,-3.5,
                4.3,-1.2,
                6.7,-2.0),byrow = T, ncol=2,nrow=3)

parms<-matrix(c(1.0,
                1.3),byrow=TRUE,ncol=1,nrow=2)

#RMCriteria::PrometheeIV(dados,c(0.3,0.7),c(0,0),parms,FALSE)

#RMCriteria::PrometheeIV(dados,c(0.3,0.7),c(1,1),parms,FALSE)

#RMCriteria::PrometheeIV(dados,c(0.3,0.7),c(2,2),parms,FALSE)

#RMCriteria::PrometheeIV(dados,c(0.3,0.7),c(3,3),parms,FALSE)

PromObj <- RPrometheeConstructor(datMat=dados,vecWeights=c(0.5,0.5),vecMaximiz=c(F,T),prefFunction=c(0,0),parms=parms,normalize=FALSE)
res <- RPrometheeIV(PromObj)
str(res)

PrometheeIVPlot(res)


Plus     <-   res@PhiPlus
Minus    <-   res@PhiMinus
Index    <-   res@Index

# Create dataframes
resDF <- data.frame("PhiPlus" = Plus, "PhiMinus" = Minus)

# Create a dataframe with results from RPrometheeI and arguments
phiLabels <- c(rep("PhiPlus", nrow(resDF)), rep("PhiMinus", nrow(resDF)))
phiNums <- c(resDF[,1], resDF[,2])
alternatives <- c(as.character(rep(1:nrow(resDF),2)))
resultsPlot <- data.frame(alternatives, phiLabels, phiNums)
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
