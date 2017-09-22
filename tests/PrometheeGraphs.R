rm(list=ls())

library(RMCriteria)
library(ggplot2)
library(dplyr)

# Examplo Ishizaka - Cap 6

dados<-matrix(c(15000,7.5,1,50,
                29000,9.0,4,110,
                38000,8.5,10,90,
                24000,8.0,8,75,
                25500,7.0,8,85),byrow = T, ncol=4,nrow=5)

rownames(dados)<-c("Economic","Sport","Luxury","Touring A","Touring B")
colnames(dados)<-c("Price","Consumption","Comfort","Power")

parms<-matrix(c(2000, 5000,
                0.5,1,
                1,2,
                10,20),byrow=FALSE,ncol=4,nrow=2)

RMCriteria::PrometheeI(dados,c(1,1,1,1),c(5,1,1,1),parms,T)
flows <- data.frame(RMCriteria::PrometheeI(dados,c(1,1,1,1),c(5,1,1,1),parms,T))
colnames(flows) <-  c("PhiPlus", "PhiMinus")

# Create each column for the results table
Phi_labels <- c(rep("PhiPlus", nrow(dados)), rep("PhiMinus", nrow(dados)))
Phi_nums <- c(flows$PhiPlus, flows$PhiMinus)
alternatives <- c(rep(rownames(dados), 2))

# Create the dataframe using columns created before
datatest <- data.frame(alternatives, Phi_labels, Phi_nums)
datatest$Phi_labels <- as.factor(datatest$Phi_labels)

# Filter the results table for only positive flows
datatest_pos <- filter(datatest, datatest$Phi_labels == "PhiPlus")
datatest_pos <- datatest_pos[,-2]

# Filter the results table for only negative flows
datatest_neg <- filter(datatest, datatest$Phi_labels == "PhiMinus")
datatest_neg <- datatest_neg[,-2]

# Create the Promethee II (sum of positive and negative flows)
datatest_sum <- datatest_pos
datatest_sum$Phi_nums <- (datatest_pos$Phi_nums - datatest_neg$Phi_nums)

# Simple bar plot with both flows
ggplot(datatest, aes(x = alternatives, y = Phi_nums, fill = Phi_labels)) +
    geom_bar(stat = "identity", position = position_dodge(), width = 0.9) +
    geom_text(aes(label = sprintf("%0.3f", Phi_nums)),
              vjust = -0.25, position = position_dodge(width = 0.9))

# Simple bar plot with Promethee2 ranking
ggplot(datatest_sum, aes(x = alternatives, y = Phi_nums)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%0.3f", Phi_nums)),
              vjust = -0.25, position = position_dodge(width = 0.9))

# Create a dataframe to use as source for both bars in PrometheeI plot
limits <- data.frame(
    class = c("PhiPlus", "PhiPlus", "PhiMinus", "PhiMinus"),
    boundaries = c(0.5, 0.5, 0.5, 0.5),
    pos_neg = c("Pos", "Neg", "Pos", "Neg"))

# Change order of factors
limits$class <- factor(limits$class, levels = c("PhiPlus", "PhiMinus"))
limits$pos_neg <- factor(limits$pos_neg, levels = c("Pos", "Neg"))

# Full bars as in Visual-Promethee.
# To-do: adjust data labels
ggplot(limits) +
  geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
           stat = "identity", width = 0.5) +
  geom_point(data = datatest, aes(x = Phi_labels, y = Phi_nums),
             stat = "identity") +
  geom_line(data = datatest, aes(x = Phi_labels, y = Phi_nums),
            group = alternatives, stat = "identity") +
  geom_text(data = datatest, aes(x = Phi_labels, y = Phi_nums),
            label = sprintf("%0.3f", round(Phi_nums, digits = 3),
                            position = position_dodge(width = 0.9)),
            hjust = 0, nudge_x = 0.05) +
  scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
  geom_text(data = datatest, aes(x = Phi_labels, y = Phi_nums),
            label = alternatives, hjust = 1, nudge_x = -0.05)
