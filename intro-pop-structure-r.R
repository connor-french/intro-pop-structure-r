## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----message=FALSE-------------------------------------------------------
#Install packages. Remove any you already have installed
pkgs <-
  c("adegenet",
    "vcfR",
    "dplyr",
    "ggplot2",
    "stringr",
    "rnaturalearth")
#install.packages(pkgs)

#load packages into environment
library(adegenet) #for pca, dapc
library(vcfR) #for reading in data and data conversion
library(dplyr) #for manipulating data
library(stringr) #for some text manipulation
library(ggplot2) #for plotting
library(rnaturalearth) #for our basemap


## ----message=FALSE, results="hide"---------------------------------------
#read in anolis vcf file that we're pulling from github and convert to genlight object
anole_genlight <-
  read.vcfR(
    "https://github.com/ivanprates/2018_Anolis_EcolEvol/blob/master/data/VCFtools_SNMF_punctatus_t70_s10_n46/punctatus_t70_s10_n46_filtered.recode.vcf?raw=true"
  ) %>%
  vcfR2genlight()

#get the population names for your samples. These can be collection site names, a priori hypotheses for the population structure, etc. 
pops <-
  read.csv(
    "https://github.com/ivanprates/2018_Anolis_EcolEvol/raw/master/data/plot_order_punctatus_n46.csv"
  )


## ------------------------------------------------------------------------
#quick summary of genlight to see what info we have
anole_genlight


## ------------------------------------------------------------------------
#centering our data is important. We are retaining 3 principal component axes (nf = 3). This is a somewhat arbitrary decision and you can change the number if you feel like it.
anole_pca <- glPca(anole_genlight, center = TRUE, nf = 3)

#tidy the data for plotting
plot_data <-
  as_tibble(anole_pca$scores, rownames = "individual") %>%
  mutate(population = pops$pop) #make a column for populations. 

#plot PC scores colored according to hypotheses population
pca_plot <-
  ggplot(plot_data, aes(
    x = PC1,
    y = PC2,
    fill = population,
    size = 3
  )) +
  geom_point(shape = 21) +
  scale_fill_viridis_d(direction = -1) + #reverse the color direction to better reflect Prates et al. 2018
  guides(size = FALSE,  #don't want a legend for the size
         fill = guide_legend(override.aes = list(size = 3))) + #the point sizes in the legend are too small, so I'm increasing their size
  theme_bw(base_size = 16) #increase legend font size and choose a more pleasant theme
pca_plot


## ------------------------------------------------------------------------
#find the best number of clusters. you have options here, so pay attention!
#Rule of thumb for me is to use most of the PC axes- k-means clustering likes as much info as it can get. Play around with setting n.pca to different values and see what you get.
num_clust <- find.clusters(anole_genlight, n.pca = 30, n.clust = 3)

#run a dpca, specifying the best number of clusters. you have more options!
#use fewer PCs for plotting (it's somewhat arbitrary- I pick a number that makes the points not too clustered together. In this case, 10 seems alright). Play around with setting n.pca to different values and see what you get.
anole_dapc <- dapc(anole_genlight, num_clust$grp, n.pca = 10, n.da = 3)

#tidy the data for plotting. adegenet has the scatter() function for plotting, but I don't like it.
plot_data_dapc <-
  as_tibble(anole_dapc$ind.coord, rownames = "individual") %>%
  mutate(population = pops$pop,
         group = anole_dapc$grp) #make a column for populations

#plot the data. Can color the points according to your pre-defined populations and the dapc groups to see if it conforms to your hypothesis.
dapc_plot <-
  ggplot(plot_data_dapc, aes(
    x = LD1,
    y = LD2,
    fill = population,
    size = 3
  )) +
  geom_point(shape = 21) +
  scale_fill_viridis_d(direction = -1) + #reverse the color direction to better reflect Prates et al. 2018
  guides(size = FALSE,  #don't want a legend for the size
         fill = guide_legend(override.aes = list(size = 3))) +
  theme_bw(base_size = 16)
dapc_plot



## ------------------------------------------------------------------------
#read in the locality data and merge it with the genetic data
anole_locs <-
  read.csv(
    "https://github.com/ivanprates/2018_Anolis_EcolEvol/raw/master/maps/qmatrix_punctatus.csv"
  ) %>%
  mutate(individual = str_c("punc_", ID)) %>% #Have to change the ID column to reflect our genetic individual ID column
  full_join(plot_data_dapc, by = "individual") #join the two datasets by the "individual" category


#read in a basemap for plotting
sa_map <-
  ne_countries(continent = "south america", returnclass = "sf")

#plot the map
anole_map <- ggplot(data = sa_map) +
  geom_sf() +
  geom_point(
    data = anole_locs,
    aes(
      x = longitude,
      y = latitude,
      fill = population,
      size = 2
    ),
    shape = 21
  ) +
  scale_fill_viridis_d(direction = -1) +
  guides(size = FALSE,  #don't want a legend for the size
         fill = guide_legend(override.aes = list(size = 3))) +
  theme_bw(base_size = 16)
anole_map

