  #### 1. Libraries and data ####

library(ape)
library(phytools)
library(tidyverse)
library(tidytree) 
library(data.table)
library(ggtree)

tree <- read.tree("data/mimosa_tree-Vasconcelos2020.tre")
tree_data <- read.csv("data/Mimosa_tree_data-Vasconcelos2020.csv", 
                      header = T, 
                      sep = ",")
pollen_data <- read.csv("data/pollen_data-mimosa.csv")

tree <- drop.tip(tree, c(subset(tree_data, clade == 
                                         "(outgroup)")$cleaned_name))

plot(tree, show.tip.label = F)

## Marking clades

tree_tibble <- as_tibble(tree)
tree_tibble$node.labels <- NA

tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "A"]))] <- "A"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "B"]))] <- "B"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "C"]))] <- "C"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "D"]))] <- "D"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "E"]))] <- "E"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "F"]))] <- "F"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "G"]))] <- "G"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "H"]))] <- "H"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "I"]))] <- "I"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "J"]))] <- "J"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "K"]))] <- "K"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "L"]))] <- "L"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "M"]))] <- "M"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "N"]))] <- "N"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "O"]))] <- "O"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "P"]))] <- "P"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "Q"]))] <- "Q"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "R"]))] <- "R"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "S"]))] <- "S"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "T"]))] <- "T"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "U"]))] <- "U"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "V"]))] <- "V"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "W"]))] <- "W"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "X"]))] <- "X"
tree$node.labels[tree$node == getMRCA(tree, c(tree_data$cleaned_name[tree_data$clade == "Z"]))] <- "Z"

tree_data <- as.treedata(tree)

c(tree$node.labels[!is.na(tree$node.labels)])
c(tree[!is.na(tree$node.labels),"node"])

pdf("tree_clades.pdf")

ggtree(tree_data, branch.length = "none") + 
  geom_nodepoint(aes(subset = node %in% 
                       c(390, 458, 471, 477, 501, 503, 524, 526, 535, 
                         596, 604, 610, 612, 615, 634, 636, 646, 657,
                         662, 663, 693, 748)), 
                 color = "black", size = 5) +
  geom_text(aes(label = node.labels), hjust = 0.5, vjust = 0.5, 
            color = "white",size=3)

dev.off()
