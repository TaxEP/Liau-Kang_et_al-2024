### MIMOSA POLLEN SAMPLING INFORMATIONS
## 04-06-2024

#### Libraries and Reading data ####

library(tidyverse)
library(phytools)
library(TNRS)
library(ape)
library(tidytree) 
library(data.table)
library(ggtree)
library(Claddis)

tree <- read.tree("data/mimosa_tree-Vasconcelos2020.tre")

dat_clade <- read.csv("data/Mimosa_tree_data-Vasconcelos2020.csv", 
                      header = T, sep = ",")

dat_pollen <- read.csv("data/pollen_data-mimosa.csv") %>%
  filter(genus == "Mimosa") %>%
  select(cleaned_name, new_description)



#-----------------------------------------------------------------------------#



#### Preparing data for sampling analysis ####

## formatting tree and clade data to TNRS

dat_clade$cleaned_name <- gsub("_", " ", dat_clade$cleaned_name, fixed = TRUE)

tree$tip.label <- gsub("_", " ", tree$tip.label, fixed = TRUE)

# correcting errors

dat_clade[311, 2] <- "Mimosa montis-carasae" #typo
dat_clade[250, 2] <- "Mimosa invisa" #"accepted name" puts it as Mimosa invisa
dat_clade[107, 4] <- "A" #typo
tree$tip.label[281] <- "Mimosa invisa" #"accepted name" puts it as Mimosa invisa
dat_clade$clade <- sub("U|V|W", "X", dat_clade$clade) #removing X sub-clades
dat_clade$clade <- sub("Z", NA, dat_clade$clade) #removing delimitation, it is
# not from Simon (2011) and groups only 2 species

tree <- drop.tip(tree, "Mimosa hirsutula") #name don't exist

##### Updating and cleaning names #####

dir.create("output/tnrs", showWarnings = F)

#----------------------------#
# LAST RUNNED IN: 04-06-2024 #
#----------------------------#

## getting updated tip labels on TNRS

#resolved_tree_labels <- TNRS(tree$tip.label,
#                           sources = c("wcvp"),
#                           classification = "wfo",
#                           mode = "resolve",
#                           matches = "best",
#                           accuracy = NULL,
#                           skip_internet_check = FALSE)
#write.csv(resolved_tree_labels, "output/tnrs/resolved_tree_labels.csv")

## getting updated barneby data taxa on TNRS

## getting updated clades data on TNRS

#resolved_dat_clade <- TNRS(dat_clade$cleaned_name,
#                            sources = c("wcvp"),
#                            classification = "wfo",
#                            mode = "resolve",
#                            matches = "best",
#                            accuracy = NULL,
#                            skip_internet_check = FALSE)
#write.csv(resolved_dat_clade, "output/tnrs/resolved_dat_clade.csv")

## reading TNRS outputs

resolved_tree_labels <- read.csv("output/tnrs/resolved_tree_labels.csv") %>%
  select(
    Name_submitted, Warnings, Name_matched, Taxonomic_status, Accepted_name)

resolved_dat_clade <- read.csv("output/tnrs/resolved_dat_clade.csv") %>%
  select(
    Name_submitted, Warnings, Name_matched, Taxonomic_status, Accepted_name)

## reviewing outputs

updates_dat_clade <- resolved_dat_clade %>%
  subset(Taxonomic_status != "Accepted" | Accepted_name == "Mimosa")

updates_tree_labels <- resolved_tree_labels %>%
  subset(Taxonomic_status != "Accepted" | Accepted_name == "Mimosa")

## correcting, if necessary

resolved_tree_labels[10, 5] = "Mimosa vernicosa var. ciliata"
resolved_dat_clade[363, 5] = "Mimosa vernicosa var. ciliata"
# "no opinion", but the name is accepted



# removing "(double entry)" and duplicates of cleaned_name

dat_clade$cleaned_name[duplicated(dat_clade$cleaned_name)]

dat_clade <- dat_clade[!duplicated(dat_clade$cleaned_name), ]

## cleaning accepted names and updating

resolved_dat_clade <- resolved_dat_clade %>%
  mutate(Accepted_name = gsub(" subsp\\.| var\\.|-", "", Accepted_name) %>%
           gsub(" ", "_", .))

resolved_tree_labels <- resolved_tree_labels %>%
  mutate(Accepted_name = gsub(" subsp\\.| var\\.|-", "", Accepted_name) %>%
           gsub(" ", "_", .))

dat_clade$accepted_name <- resolved_dat_clade$Accepted_name
tree$tip.label <- resolved_tree_labels$Accepted_name

# removing duplicates of clade data

dat_clade <- dat_clade[!duplicated(dat_clade$accepted_name), ]

## checking if all tree species have clade data

setdiff(tree$tip.label, dat_clade$accepted_name)
# Mimosa eurystegia really don't have clade info... moving on



#-----------------------------------------------------------------------------#



#### Generating analyses data

tree <- keep.tip(tree, tree$tip.label[grep("^Mimosa_", tree$tip.label)])

dat <- as.data.frame(tree$tip.label)
colnames(dat) <- "taxon"

## adding clade information from Vasconcelos 2020 to data

dat <- merge(dat, dat_clade[, c("accepted_name", "clade")], 
             by.x = "taxon", 
             by.y = "accepted_name", all.x = TRUE)

## adding pollen sampling information to data

dat <- dat %>%
  mutate(pollen_data = ifelse(taxon %in% dat_pollen$cleaned_name, "yes", NA))

## saving file for disparity analysis

write.csv(dat, "output/data/spp_clade_data.csv", row.names = F)



#-----------------------------------------------------------------------------#



#### Analyzing sampling ####

## Percentage of tree taxa with pollen data (original)

original_sampling <- dat_pollen %>%
  filter(!new_description %in% "x")

pollen_sampling <- dat %>%
  filter(pollen_data == "yes")

round(sum(original_sampling$cleaned_name %in% pollen_sampling$taxon) / 
        length(tree$tip.label) * 100,2)

## Percentage of tree taxa with pollen data (new description)

round(sum(!is.na(dat$pollen_data)) / length(tree$tip.label) * 100,2)

## Number and name of matrix species not in the phylogeny

length(setdiff(dat_pollen$cleaned_name, dat$taxon))

setdiff(dat_pollen$cleaned_name, dat$taxon)

dir.create("output/data", showWarnings = F)

write.table(setdiff(dat_pollen$cleaned_name, dat$taxon), 
            "output/data/taxa-data_not_tree.csv", 
            col.names = FALSE, 
            row.names = FALSE)

## Percentage of pollen data for each tree clade

clades_percentages <- dat %>%
  group_by(clade) %>%
  summarise(pollen_data = round(mean(!is.na(pollen_data)) * 100,2))

write.csv(clades_percentages,
          "output/data/sampling-clade_percentages.csv",
          row.names = F)

  #### DISPARITY ANALYSIS - MORPHOSPACE AND METRICS - MIMOSA POLLEN ####

dat_pollen <- read.csv("data/pollen_data-mimosa.csv", na.strings = c("NA")) %>%
  filter(genus == "Mimosa")

#### Preparing data for disparity analysis ####

## creating dataframe with important information columns

info_pollen <- dat_pollen %>%
  select(cleaned_name, new_description)

## filtering for characters with missing data minor or equal to 60%

matrix_pollen <- dat_pollen %>%
  select(11:21) %>% #selecting columns with pollen characters
  .[, colMeans(is.na(.)) * 100 <= 40] #filtering

rownames(matrix_pollen) <- dat_pollen$cleaned_name

## making discretization of continuous traits in order to use Claddis package

# we will use the Freedman-Diaconis rule to define an optimum number of breaks

# In statistics, the Freedman–Diaconis rule can be used to select the width of
# the bins to be used in a histogram.
# Freedman-Diaconis rule is designed roughly to minimize the integral of the 
# squared difference between the histogram (i.e., relative frequency density) 
# and the density of the theoretical probability distribution.
# The Freedman-Diaconis rule isn’t axiomatic: If the suggested bin width/number 
# of bins seems too few or too great, judgment is used to scale up or down as 
# needed.

#shorter_diameter_mean

plot(density(matrix_pollen$shorter_diameter_mean, na.rm = TRUE))

hist(matrix_pollen$shorter_diameter_mean,
     breaks = "FD",
     main = "shorter_diameter_mean")

#longer_diameter_mean

plot(density(matrix_pollen$longer_diameter_mean, na.rm = TRUE))

hist(matrix_pollen$longer_diameter_mean,
     breaks = "FD",
     main = "longer_diameter_mean")

#exine_thickness_mean

plot(density(matrix_pollen$exine_thickness_mean, na.rm = TRUE))

hist(matrix_pollen$exine_thickness_mean,
     breaks = "FD",
     main = "exine_thickness_mean")

#running the discretization function

matrix_pollen <- arules::discretizeDF(
  matrix_pollen,
  methods = list
  (
    shorter_diameter_mean = list(
      method = "interval",
      breaks = 14,
      labels = c(0:13)
    ),
    longer_diameter_mean = list(
      method = "interval",
      breaks = 8,
      labels = c(0:7)
    ),
    exine_thickness_mean = list(
      method = "interval",
      breaks = 10,
      labels = c(0:9)
    )
  )
)

#### Disparity analysis - Morphospace ####

# characters ordination

ord <- c(rep("unordered",7))

## building cladistic matrix

cladistic_matrix <- as.matrix(matrix_pollen) %>%
  build_cladistic_matrix(ordering = ord)

## calculating distance matrix and performing principal coordinates analysis

pcoa <- ordinate_cladistic_matrix(
  cladistic_matrix,
  distance_metric = "mord",
  distance_polymorphism_behaviour = "min_difference",
  correction = "cailliez")

## calculating variance explained explained by each principal component

eigen <- pcoa$values$Corr_eig

# PC1 variance explained:

round(eigen[1] / sum(eigen[eigen > 0]) * 100, 2)

# PC2 variance explained:

round(eigen[2] / sum(eigen[eigen > 0]) * 100, 2)

# PC3 variance explained

round(eigen[3] / sum(eigen[eigen > 0]) * 100, 2)

# extracting vectors to plot

info_pollen[, c(3,4)] <- (pcoa$vectors[, c(1,2)])

colnames(info_pollen)[2:4] <- c("Status", "PC1", "PC2")

## preparing "status" column with pollen description info

info_pollen <- info_pollen %>%
  mutate(Status = case_when(
    is.na(Status) ~ "Already described",
    Status == "x" ~ "New description"
  ))

## plotting morphospace

morphospace <- ggplot(info_pollen, mapping = aes(x = PC1, 
                                                 y = PC2,
                                                 color = Status,
                                                 shape = Status)) +
  geom_point() +
  coord_fixed(
    ratio = 1,
    xlim = c(-4.5, 3),
    ylim = c(-4.5, 3),
    expand = TRUE,
    clip = "on"
  ) +
  xlab(paste(
    "PCo", 1, " ", "(", round(pcoa$values$Rel_corr_eig[1] * 100, 2), "%", ")",
    sep = ""
  )) + ylab(paste("PCo", 2, " ", "(", 
                  round(pcoa$values$Rel_corr_eig[2] * 100, 2), "%", ")",
                  sep = ""
  )) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

## saving

dir.create("output/plots", showWarnings = F)

pdf("output/plots/morphospace.pdf")

morphospace

dev.off()



#### Disparity analysis - Disparity Metrics ####

library(dispRity)
library(openxlsx)

## reading clades data file

clades_dat <- read.csv("output/data/spp_clade_data.csv")

## creating lists grouping species of each clade
## (one list for before and other for after new descriptions)

# adding clade information to info_pollen
vectors_clades <- merge(info_pollen, clades_dat[, c("taxon", "clade")], 
             by.x = "cleaned_name", 
             by.y = "taxon", all.x = TRUE)

# removing taxa without clade information
vectors_clade <- na.omit(vectors_clades)

# creating the list of lists for clades

groups_clades <- vectors_clades %>%
  select(cleaned_name, clade)

groups_clades <- split(groups_clades$cleaned_name, groups_clades$clade)

# creating the list of lists for all data (before and after)

group_original <- vectors_clades %>%
  filter(Status == "Already described") %>%
  select(cleaned_name)

group_new <- vectors_clades %>%
  select(cleaned_name)

groups_samplings <- list("Previous data" = c(group_original$cleaned_name), 
                         "Current data" = group_new$cleaned_name)

## calculating sum of variances and statistics

bootstraps <- 1000

set.seed(777)

# clades - SV

groups_clades <- Filter(function(x) { length(x) > 3 }, groups_clades)

sub_clade <- custom.subsets(pcoa$vectors, groups_clades)
bootstrapped_data <-
  boot.matrix(sub_clade, bootstraps = bootstraps, rarefaction = F)
disparity_clade_sv <-
  dispRity(bootstrapped_data, metric = c(sum, variances))
disp_clade_sv <- summary(disparity_clade_sv)

# samplings - SV

sub_sampling <- custom.subsets(pcoa$vectors, groups_samplings)
bootstrapped_data <-
  boot.matrix(sub_sampling, bootstraps = bootstraps, rarefaction = F)
disparity_sampling_sv <-
  dispRity(bootstrapped_data, metric = c(sum, variances))
disp_sampling_sv <- summary(disparity_sampling_sv)

# clades - SR

groups_clades <- Filter(function(x) { length(x) > 3 }, groups_clades)

sub_clade <- custom.subsets(pcoa$vectors, groups_clades)
bootstrapped_data <-
  boot.matrix(sub_clade, bootstraps = bootstraps, rarefaction = F)
disparity_clade_sr <-
  dispRity(bootstrapped_data, metric = c(sum, ranges))
disp_clade_sr <- summary(disparity_clade_sr)

# samplings - SR

sub_sampling <- custom.subsets(pcoa$vectors, groups_samplings)
bootstrapped_data <-
  boot.matrix(sub_sampling, bootstraps = bootstraps, rarefaction = F)
disparity_sampling_sr <-
  dispRity(bootstrapped_data, metric = c(sum, ranges))
disp_sampling_sr <- summary(disparity_sampling_sr)

# saving

library(patchwork)

list_of_datasets_sv <- list("Disparity - Samplings" = disp_sampling_sv, 
                         "Disparity - Clades" = disp_clade_sv)
write.xlsx(list_of_datasets_sv, file = "output/data/sum_of_variances.xlsx")

list_of_datasets_sr <- list("Disparity - Samplings" = disp_sampling_sr, 
                         "Disparity - Clades" = disp_clade_sr)
write.xlsx(list_of_datasets_sr, file = "output/data/sum_of_ranges.xlsx")

pdf("output/plots/metrics.pdf")

par(mfrow=c(2,2), mar=c(2,2,1.5,0.5), xpd=TRUE)

sampling_colors <- c("#4C72B0", "#DD8452")
clades_colors <- c("#4C78A8", "#9ECF99", "#F17C67", 
                   "#B279A2", "#CCB974", "#64B5CD",
                   "#E36C09","#D17FED")

plot(disparity_sampling_sv, main = "Sum of Variances", 
     colors = sampling_colors)
plot(disparity_sampling_sr, main = "Sum of Ranges", yaxt = "n", 
     colors = sampling_colors)
plot(disparity_clade_sv, 
     colors = clades_colors)
plot(disparity_clade_sr, yaxt = "n", 
     colors = clades_colors)

dev.off()

  #### TREE PLOT - CLADES AND POLLEN DATA ####

plot(tree, show.tip.label = F)

## Marking clades

tree_tibble <- as_tibble(tree)
tree_tibble$node.labels <- NA

tree_tibble$node.labels[tree_tibble$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "A"]))] <- "A"
tree_tibble$node.labels[tree_tibble$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "B"]))] <- "B"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "C"]))] <- "C"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "D"]))] <- "D"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "E"]))] <- "E"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "F"]))] <- "F"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "G"]))] <- "G"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "H"]))] <- "H"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "I"]))] <- "I"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "J"]))] <- "J"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "K"]))] <- "K"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "L"]))] <- "L"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "M"]))] <- "M"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "N"]))] <- "N"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "O"]))] <- "O"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "P"]))] <- "P"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "Q"]))] <- "Q"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "R"]))] <- "R"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "S"]))] <- "S"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "T"]))] <- "T"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "U"]))] <- "U"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "V"]))] <- "V"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "W"]))] <- "W"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "X"]))] <- "X"
tree$node.labels[tree$node == getMRCA(tree, c(clades_dat$taxon[clades_dat$clade == "Z"]))] <- "Z"

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
