### MIMOSA POLLEN SAMPLING INFORMATIONS
## 04-06-2024

#### Libraries and Reading data ####

library(tidyverse)
library(phytools)
library(TNRS)

tree <- read.nexus("data/mimosa_tree-Vasconcelos2020.nex")

dat_clade <- read.csv("data/Mimosa_tree_data-Vasconcelos2020.csv", 
                      header = T, sep = ",")

dat_pollen <- read.csv("data/pollen_data-mimosa.csv") %>%
  filter(genus == "Mimosa") %>%
  select(cleaned_name)



#-----------------------------------------------------------------------------#



#### Preparing data for analysis ####

## formatting tree and clade data to TNRS

dat_clade$cleaned_name <- gsub("_", " ", dat_clade$cleaned_name, fixed = TRUE)

tree$tip.label <- gsub("_", " ", tree$tip.label, fixed = TRUE)

# correcting errors

dat_clade[311, 2] <- "Mimosa montis-carasae" #typo
dat_clade[250, 2] <- "Mimosa invisa" #"accepted name" puts it as Mimosa invisa
dat_clade[107, 4] <- "A" #typo
tree$tip.label[261] <- "Mimosa invisa" #"accepted name" puts it as Mimosa invisa

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

## adding clade information from Vasconcelos 2020 to data for analysis

dat <- merge(dat, dat_clade[, c("accepted_name", "clade")], 
              by.x = "taxon", 
              by.y = "accepted_name", all.x = TRUE)

dat <- merge(dat, dat_pollen[, c("cleaned_name", "cleaned_name")], 
             by.x = "taxon", 
             by.y = "cleaned_name", all.x = TRUE)

## inserting the clades for some species

#-----------------------------------------------------------------------------#



#### Analyzing sampling ####

## Percentage of tree taxa with pollen data

round(sum(!is.na(dat$cleaned_name.1)) / length(dat$cleaned_name.1) * 100,2)

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
  summarise(pollen_data = round(mean(!is.na(cleaned_name.1)) * 100,2))

write.csv(clades_percentages,
          "output/data/sampling-clade_percentages.csv",
          row.names = F)

