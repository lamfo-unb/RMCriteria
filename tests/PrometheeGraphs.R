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
flows$Phi <- flows[,1] - flows[,2]
colnames(flows) <-  c("PhiPlus", "PhiMinus", "Phi")

# Create each column for the results table
Phi_labels <- c(rep("PhiPlus", nrow(dados)), rep("PhiMinus", nrow(dados)),
                rep("Phi", nrow(dados)))
Phi_nums <- c(flows[,1], flows[,2], flows[,3])
alternatives <- c(rep(rownames(dados), 3))

# Create the dataframe using columns created before
datatest <- data.frame(alternatives, Phi_labels, Phi_nums)
datatest$Phi_labels <- as.factor(datatest$Phi_labels)


# Create a dataframe to use as source for both Phi bars in PrometheeI plot
limits <- data.frame(
  class = c("PhiPlus", "PhiPlus", "PhiMinus", "PhiMinus"),
  boundaries = c(0.5, 0.5, 0.5, 0.5),
  pos_neg = c("Pos", "Neg", "Pos", "Neg"))

# Change order of factors
limits$class <- factor(limits$class, levels = c("PhiPlus", "PhiMinus"))
limits$pos_neg <- factor(limits$pos_neg, levels = c("Pos", "Neg"))

# Filter results table to exclude Phi flows
prometheeI_temp <- filter(datatest, Phi_labels != "Phi")
prometheeI_temp[,2] <- factor(prometheeI_temp[,2],
                              levels = c("PhiPlus", "PhiMinus"))
prometheeI_temp[,1] <- factor(prometheeI_temp[,1],
                              levels = c("Economic","Sport","Luxury","Touring A","Touring B"))

# Full bars as in Visual-Promethee.
ggplot(limits) +
  geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
           stat = "identity", width = 0.5) +
  geom_point(data = prometheeI_temp, aes(x = Phi_labels, y = Phi_nums),
             stat = "identity") +
  geom_line(data = prometheeI_temp, aes(x = Phi_labels, y = Phi_nums),
            group = prometheeI_temp[,1], stat = "identity") +
  geom_text(data = prometheeI_temp, aes(x = Phi_labels, y = Phi_nums),
            label = sprintf("%0.3f",
                            round(prometheeI_temp$Phi_nums, digits = 3),
                            position = position_dodge(width = 0.9)),
            hjust = 0, nudge_x = 0.05) +
  scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
  geom_text(data = prometheeI_temp, aes(x = Phi_labels, y = Phi_nums),
            label = prometheeI_temp$alternatives, hjust = 1, nudge_x = -0.05)



##################################################################
# Create a dataframe to use as source for bar in PrometheeII plot
limits_II <- data.frame(
  class = c("Phi", "Phi"),
  boundaries = c(-1, 1),
  pos_neg = c("Neg", "Pos"))

# Change order of factors
limits_II$pos_neg <- factor(limits_II$pos_neg, levels = c("Pos", "Neg"))

# Filter results table to exclude Phi flows
prometheeII_temp <- filter(datatest, Phi_labels == "Phi")
prometheeII_temp[,2] <- factor(prometheeII_temp[,2], levels = "Phi")
prometheeII_temp[,1] <- factor(prometheeII_temp[,1],
                              levels = c("Economic","Sport","Luxury","Touring A","Touring B"))


# Full Ranking bar as in Visual-Promethee.
ggplot(limits_II) +
  geom_bar(aes(x = class, y = boundaries, fill = pos_neg),
           stat = "identity", width = 0.5) +
  geom_point(data = prometheeII_temp, aes(x = Phi_labels, y = Phi_nums),
             stat = "identity") +
  geom_text(data = prometheeII_temp, aes(x = Phi_labels, y = Phi_nums),
            label = sprintf("%0.3f",
                            round(prometheeII_temp$Phi_nums, digits = 3)),
            hjust = 0, nudge_x = 0.03) +
  scale_fill_manual(aes(x = class, y = boundaries), values = c("#a1d99b", "#F57170")) +
  geom_text(data = prometheeII_temp, aes(x = Phi_labels,
                                         y = prometheeII_temp$Phi_nums),
            label = prometheeII_temp$alternatives,
            hjust = 1, nudge_x = -0.03)
