### DISPARITY ANALYSIS - MORPHOSPACE - MIMOSA POLLEN
## 05-10-2023

#### Libraries and Data ####

library(tidyverse)
library(phytools)
library(Claddis)

dat_pollen <- read.csv("data/pollen_data-mimosa.csv", na.strings = c("", "?")) %>%
  filter(genus == "Mimosa")

#### Preparing data for analysis ####

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

ord <- c("unord", "unord", "unord", "unord", "unord", "unord", "unord")

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

jpeg(
  filename = "output/plots/morphospace.jpeg",
  width = c(174),
  height = c(134),
  units = "mm",
  res = 600
)

morphospace

dev.off()


