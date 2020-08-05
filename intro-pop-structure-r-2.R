#' ---
#' title: "Exploring population structure in R with adegenet and sNMF"
#' knit: (function(input_file, encoding) {
#'   out_dir <- 'docs';
#'   rmarkdown::render(input_file,
#'  encoding=encoding,
#'  output_file=file.path(dirname(input_file), out_dir, 'index.html'))})
#' author: "Connor French"
#' date: "`r Sys.Date()`"
#' output: 
#'   html_document:
#'     toc: TRUE
#'     toc_float: TRUE
#' ---
#' 
## ---- echo=FALSE---------------------------------------------
# library(downloadthis)
# 
# download_link(
#   link = "https://github.com/connor-french/intro-pop-structure-r/raw/master/docs/pre-workshop-instructions.pdf",
#   button_label = "Pre-workshop instructions",
#   button_type = "danger",
#   has_icon = TRUE,
#   icon = "fa fa-save",
#   self_contained = FALSE
# )

#' 
## ---- echo=FALSE---------------------------------------------
# download_file(
#   path = "/Users/connorfrench/Dropbox/Old_Mac/Science/tutorials/intro-pop-structure-r/intro-pop-structure-r.R",
#   button_label = "R script",
#   button_type = "danger",
#   has_icon = TRUE,
#   icon = "fa fa-save",
#   self_contained = FALSE
# )

#' 
#' 
#' # Introduction
#' 
#' The aim of this tutorial is to give you a head start in exploring structure in your SNP data. Some familiarity of the R language or command line is helpful, but not necessary. I've provided resources at the end for anyone new to R who want to develop their skills further. There are also links to tutorials that explore other population genetics analysis, like estimating genetic variation, genetic divergence between populations, etc.
#' 
#' This tutorial is split into four sections:
#' 
#' **Loading and pre-processing data**
#' 
#' **DAPC for quick population structure inference**
#' 
#' **sNMF for admixture inference**
#' 
#' **Mapping your results**
#' 
#' 
#' We will be exploring the genetic structure of *Anolis punctatus*, an anole species from the Amazon and Atlantic Forest. Our data comes from [Prates et al. 2018, Local adaptation in mainland anole lizards: Integrating population history and genome-enviroment associations](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4650). We want to see if the genetic structure we infer matches that found in *Prates et al.*. 
#' 
#' <img src="https://photos.smugmug.com/HerpGalleries/Peru-Amazon-2013-1/i-58tqn64/1/b4bc26eb/L/Anolis%20punctatus%20Amazon%20Green%20Anole%20MS-L.jpg" alt="*Anolis punctatus*. Source: https://cages.smugmug.com/HerpGalleries/Peru-Amazon-2013-1/i-58tqn64" width="200"/>
#' 
#' ![](images/prates-fig.png)
#' 
#' 
#' 
#' **Workflow**  
#' All code can be run by copying and pasting the text in gray code chunks into your R console. In the top left corner of each chunk there is a button that makes copy-pasting easier. If you would rather run the code from a script, a download button for the R script is provided at the top. As some of the functions in this tutorial are interactive, I don't recommend using an Rmarkdown file to perform these analyses. If you do not have R, you can download it [here](https://cloud.r-project.org/). If you have R and would like to work in the RStudio environment, you can download it [here](https://www.rstudio.com/products/rstudio/download/). For those familiar with R- I am following the [tidyverse](https://www.tidyverse.org/) style of coding, which may look unfamiliar for those used to working in base R. I've provided a link at the bottom (and [here](https://rstudio.cloud/learn/primers)) that can bring you up to speed and go further if you're interested in this reader and user-friendly coding style. Also, feel free to ask questions during the workshop!  
#' 
#' ## 1. Loading and preprocessing data
#' 
#' If you haven't installed the following packages, download [the preworkshop instructions](https://github.com/connor-french/intro-pop-structure-r/raw/master/docs/pre-workshop-instructions.pdf) and follow those steps. After the packages have been installed, they're still not ready to be used. You now need to load them into our R environment. You can do this with the `library()` function.  
#' 
#' **1.1** 
## ----message=FALSE-------------------------------------------
#load packages into environment
library(adegenet) # for dapc
library(vcfR) # for reading in genetic data
library(tidyverse) # for manipulating and plotting data
library(LEA) # For sNMF
library(rnaturalearth) #for mapping

#' 
## ----klippy, echo=FALSE, include=TRUE------------------------
# enable code copying
# klippy::klippy()

#' 
## ----setup, echo=FALSE---------------------------------------
# knitr::opts_chunk$set(collapse = TRUE, tidy = TRUE)

#' 
#' 
#' Now time for some data! The `read.vcfR()` function from the `vcfR` package allows you to read in a vcf file. The `vcfR` package comes with functionality that makes manipulating vcf data more intuitive than other packages, so I use it when I can. Note that the path is to a file hosted on github! The `vcfR2genlight()` function converts the vcf into the genlight format. The genlight format is a compact representation of the genotype matrix that `adegenet` uses to efficiently process the data.  
#' 
#' **1.2**
## ----message=FALSE, results="hide"---------------------------
#read in anolis vcf file that we're pulling from github and convert to genlight object
anole_vcf <-
  read.vcfR(
    "https://github.com/ivanprates/2018_Anolis_EcolEvol/blob/master/data/VCFtools_SNMF_punctatus_t70_s10_n46/punctatus_t70_s10_n46_filtered.recode.vcf?raw=true"
  )
anole_genlight <- vcfR2genlight(anole_vcf)

#' 
#' Lets take a quick look at a summary of our SNP data. Everything looks right, according to *Prates et al.*! There are 3,249 SNPs, with 7.13% missing data for 46 individuals. The SNPs have already been thinned to one per locus, a necessary step to limit linkage among loci.  
#' 
#' **1.3**
## ------------------------------------------------------------
#quick summary of genlight to see what info we have
anole_genlight

#' 
#' Another important component of the analysis is population metadata. These can be locality names, *a priori* hypotheses for population structure (e.g. whether the individuals are on one side of a mountain or not), or other information that you think could help you interpret your results. In this case, you're comparing your inference against *Prates et al.*'s findings. You're going to read in the results of their population structure inference to compare results with.  
#' 
#' Here, they have an individual ID column, the position they plotted the individual in their admixture barplot, and the population name they assigned the individual.  
#' 
#' **1.4**
## ------------------------------------------------------------
pops <-
  read_csv(
    "https://github.com/ivanprates/2018_Anolis_EcolEvol/raw/master/data/plot_order_punctatus_n46.csv"
  )

pops

#' 
#' 
#' ## 2. DAPC for quick population structure inference
#' 
#' Discriminant analysis of principal components (DAPC) is one of many techniques to explore population structure. It distinguishes itself by being applicable to very large datasets and doesn't assume an evolutionary model like [STRUCTURE](https://web.stanford.edu/group/pritchardlab/structure.html) or [ADMIXTURE](http://software.genetics.ucla.edu/admixture/). The output of DAPC is a little different than other programs. DAPC does not assign mixed ancestry- if there is shared ancestry within individuals, cluster colors won't indicate it. DAPC is best used to gain a quick understanding of different levels of structure in your data. **All of these methods are exploratory**, and should be treated as such. Don't fret over ambiguous results too much- use it as an opportunity to explore your data.  
#' 
#' Again, the aim of a DAPC analysis is to find the best number of clusters in your data and maximize the differences between clusters. The `find.clusters()` function first performs a PCA on the data, then performs a clustering algorithm (k-means) that considers a range of potential K values to group the data by. The goodness of fit of each K value is indicated by a BIC score, where the lowest BIC indicates the best fit. The `dapc()` function performs a discriminant analysis with the chosen K value, finding the axis of variation that best discriminates among the groups. This is why you'll typically see stronger grouping in a DAPC analysis versus a PCA.
#' 
#' The `find.clusters()` function will give you visualizations of each step in the process and prompt you to choose the number of PCs and the number of clusters to consider. You should retain as many PCs as your computer processor will allow, as k-means clustering is data-hungry and including more axes improves performance. In this case, with a modest number of individuals and loci, retaining all of the PCs (46) is reasonable- the program runs very quickly.  
#' Picking the number of clusters to retain is a tricky exercise. The "elbow" rule is a good rule of thumb for choosing a K value, but oftentimes an elbow isn't evident in the information criterion plot, or there are multiple elbows. I go more into depth about this in the slide presentation, but I'll emphasize it here- a "best K value" doesn't exist! Pick a few, visualize the results, and interpret the results in the context of your study.  
#' 
#' You'll notice that the best number of clusters isn't immediately obvious for this data set. There is no obvious "knee" and K=3, 4, or 5 could be reasonable. For this exercise, let's choose to retain 3 clusters, the same as *Prates et al*. If you come back to the tutorial after this workshop, explore the different K values and see what patterns pop out!  
#' 
#' **2.1**
## ---- eval=FALSE---------------------------------------------
num_clust <- find.clusters(anole_genlight)

#' 
#' ![](images/pca_1.png)
#' 
#' ![](images/bic.png)
#' 
#' 
## ---- echo=FALSE, include=TRUE-------------------------------
# assign num_clust so the rmarkdown file can run
num_clust <- find.clusters(anole_genlight, n.clust = 3, n.pca = 46)

#' 
#' Now, run the DAPC using the `dapc()` function. You have two options here. First, choose the number of PCs to retain. This choice is important as choosing too many can result in overfitting, while too few can result in overdispersion. A rule of thumb is to look at where the cumulative variance explained by the PCs begins to level off. It's hard to tell with this data set, but 10 PCs looks reasonable. As always, this is an exploratory analysis, so try out a few options with your own data set and see what the results look like.  
#' 
#' Second, choose the number of discriminant functions to retain. For small numbers of clusters, retaining all of them is reasonable. However, as the number of clusters increases (tens of axes), the information contained within each axis decreases, so choosing the first few axes is a better option. With only three clusters, retaining both discriminant functions is reasonable.  
#' 
#' **2.2**
## ---- eval=FALSE---------------------------------------------
anole_dapc <- dapc(anole_genlight, num_clust$grp)

#' 
#' ![](images/pca_2.png)
#' 
#' ![](images/da.png)
#' 
## ---- echo=FALSE, include=FALSE------------------------------
# so rmarkdown can run
anole_dapc <- dapc(anole_genlight, num_clust$grp, n.pca = 10, n.da = 2)

#' 
#' 
#' Now, let's plot the results! First, adegenet's default plotting method. It is customizable (see [this tutorial](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf)) for all of the options), but I like a little more flexibility. For instance, it's not intuitive how to map an *a priori* clustering hypothesis using their method. Since you're comparing your results with *Prates et al.*, you need some of that flexibility.  
#' 
#' **2.3**
## ------------------------------------------------------------
scatter(anole_dapc, posi.da="bottomleft")

#' 
#' 
#' The `anole_dapc` object contains all of the information you need to use a more flexible plotting software, so you are going to do that here with `ggplot2`.  
#' 
#' First, you need to extract and reshape the data to make it amenable to plotting. You need three pieces of information to do this- the individual coordinates across the first two linear discriminant axes (`anole_dapc$ind.coord`), the population assignment of each individual (`pops$pop`), and the cluster assignment (`anole_dapc$grp`).  
#' 
#' **2.4**
## ------------------------------------------------------------
#tidy the data for plotting. 
dapc_data_df <-
  # as_tibble() converts the ind.coord matrix to a special data frame called a tibble. 
  as_tibble(anole_dapc$ind.coord, rownames = "individual") %>%
  # mutate changes or adds columns to a data frame. Here, we're adding the population and group assignment columns to the data frame
  mutate(population = pops$pop,
         group = anole_dapc$grp)

dapc_data_df

#' 
#' 
#' Now that you have an appropriate data frame, you can plot the data! This is using the `ggplot2` package from the `tidyverse`. You're making a scatterplot with the first LD on the x-axis and the second LD on the y-axis, and coloring each point according to the population name given by *Prates et al.*. The rest of the options change the look of the plot. If you're interested in learning more about data wrangling and data visualization in R, check out these excellent [interactive primers](https://rstudio.cloud/learn/primers).  
#' 
#' It looks like our findings agree with *Prates et al*!  
#' 
#' **2.5**
## ------------------------------------------------------------
#plot the data. Can color the points according to your pre-defined populations and the dapc groups to see if it conforms to your hypothesis.
dapc_plot <-
  ggplot(dapc_data_df, aes(
    x = LD1,
    y = LD2,
    fill = population
  )) +
  geom_point(shape = 21, size = 3) +
  #reverse the color direction to better reflect Prates et al. 2018
  scale_fill_viridis_d(direction = -1) + 
  theme_bw(base_size = 16)

dapc_plot

#' 
#' 
#' ## 3. sNMF for admixture inference
#' Now, let's use another approach for inferring population structure and admixture. `sNMF` requires a `.geno` file as input. This format keeps the individuals as columns and their genotypes as rows, where 0 indicates homozygous for the reference allele, 1 indicates heterozygous, 2 indicates homozygous for the alternate allele, and 9 indicates missing data.
#' 
#' ![](images/genotype-matrix.png)
#' 
#' 
#' There are many ways to get a `.geno` file from a vcf or genotype matrix (`LEA` has a function, `vcf2geno` that converts a vcf to a `.geno` file, but it's finicky with certain types of vcf files). Here's a way to convert from the `vcfR` object to the `.geno` format using the `vcfR` package and some data wrangling. This should work with any vcf that's read in using `vcfR`, so for your own data, just replace `anole_vcf` with your vcf!  
#' 
#' **3.1**
## ------------------------------------------------------------
# first need to code the genotypes according to the .geno system
anole_tidy <- anole_vcf %>% 
  extract_gt_tidy() %>% 
  select(-gt_DP, -gt_CATG, -gt_GT_alleles) %>% 
  mutate(gt1 = str_split_fixed(gt_GT, "/", n = 2)[,1],
         gt2 = str_split_fixed(gt_GT, "/", n = 2)[,2],
         geno_code = case_when(
           # homozygous for reference allele = 0
           gt1 == 0 & gt2 == 0 ~ 0,
           # heterozygous = 1
           gt1 == 0 & gt2 == 1 ~ 1,
           gt1 == 1 & gt2 == 0 ~ 1,
           # homozygous for alternate allele = 2
           gt1 == 1 & gt2 == 1 ~ 2,
           # missing data = 9
           gt1 == "" | gt2 == "" ~ 9
         )) %>% 
  select(-gt_GT, -gt1, -gt2)

# now you need to rotate the table so individuals are columns and genotypes are rows
anole_geno <- anole_tidy %>% 
  pivot_wider(names_from = Indiv, values_from = geno_code) %>% 
  select(-Key)

anole_geno

#' 
#' sNMF requires that files be written to your computer, so you're going to write the `.geno` file to a directory of your choice. Replace the path I provide with the path you want to write to. Make sure the file is called "anole_geno.geno".  
#' 
#' **3.2**
## ------------------------------------------------------------
write.table(anole_geno, 
            "~/Desktop/anole_geno.geno",
            col.names = FALSE,
            row.names = FALSE,
            sep = "")

#' 
#' 
#' Now to run an sNMF analysis! First, a quick explanation of each argument:  
#' 
#' * input.file is a path to your `.geno` file  
#' * K is a vector of cluster numbers you want to consider
#' * entropy- set to TRUE to calculate the cross-entropy criterion, which is how you evaluate models  
#' * repetitions- results can vary among runs. You want to conduct multiple runs to make sure your results are stable and/or find the best run  
#' * project- you can add to an existing project ("continue"), remove the current project and start a new one ("new"), or store the results in the current project even if the input file has been modified ("force")  
#' * alpha- the alpha regularization parameter. This penalizes intermediate ancestry proportions. A general rule of thumb for smaller (<10,000 SNPs) datasets is higher alpha (100-1000+) leads to better performance  
#' 
#' sNMF spits out a lot of text to your screen. It's telling you what parameters it is using for each run, where it is writing its data, and how to load the project. For most applications, sNMF has helper functions so you don't directly interact with these files.  
#' 
#' **3.3**
## ---- collapse=FALSE,  results="hide"------------------------
anole_snmf <- snmf(input.file = "~/Desktop/anole_geno.geno",
                   K = 1:10,
                   entropy = TRUE,
                   repetitions = 5,
                   project = "new",
                   alpha = 100
                  )

#' 
#' 
#' sNMF has a `plot()` helper function that makes visualizing the cross-entropy criterion across K values easy. We can see here that K=3 or K=4 are suitable K values. *Prates et al.* performed a more robust exploration of parameter space (alpha values of 1, 10, 100, 1000 with 100 repetitions), and consistently recovered K=3. Let's go forward using K=3, but exploring K=4 or K=5 could recover interesting population structure.
#' 
#' **3.4**
## ------------------------------------------------------------
plot(anole_snmf, cex = 1.2, col = "lightblue", pch = 19)

#' 
#' For K=3, which run has the best fit? In this case, run 4 looks like it has the best fit.  
#' 
#' **3.5**
## ------------------------------------------------------------
# get the cross-entropy of the 5 runs for K = 3
ce <-  cross.entropy(anole_snmf, K = 3)

ce

#' 
#' You can get the index of the fastest run with the `which.min()` function.  This is useful when you're running many different analyses and the best run changes across analyses.  
#' 
#' **3.6**
## ------------------------------------------------------------
# select the run with the lowest cross-entropy for K = 3
best_run <- which.min(ce)

best_run

#' 
#' 
#' Now time for plotting! First, you need the Q-matrix, which contains the ancestry proportions. You're extracting the Q matrix for the best run of K=3. Each column is a population and each row is an individual, so let's add column names to make that more apparent.  
#' 
#' **3.7**
## ------------------------------------------------------------
q_mat <- LEA::Q(anole_snmf, K = 3, run = best_run) 

colnames(q_mat) <- paste0("P", 1:3)

head(q_mat)

#' 
#' 
#' Unfortunately, this format doesn't make plotting easy. `LEA` has a default `barchart()` function, but it isn't very flexible and labeling individuals doesn't work out of the box (even though the "lab" argument in the function makes it seem like it does). To match your barchart to *Prates et al.*'s, you'll have to do some data wrangling.  
#' 
#' First, you need to combine the population metadata from *Prates et al* with the Q matrix.  
#' 
#' **3.8**
## ------------------------------------------------------------
# convert the q matrix to a data frame
q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = pops$ID,
         region = pops$pop,
         order = pops$plot_order)
q_df

#' 
#' Now, you're transforming the data to a "long" format for plotting. The population column names get their own column and the ancestry proportions (q) get their own column.  
#' 
#' **3.9**
## ------------------------------------------------------------
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

q_df_long

#' 
#' 
#' Finally, you're ordering the individuals for plotting. Since we're interested in matching the data to *Prates et al.*'s results, you're going to order the individuals according to their plotting specifications. 
#' 
#' 
#' I also included code that will order the individuals by their ancestry proportions. We won't use it here, but I've provided it for use with your own code.  
#' 
#' **3.10**
## ------------------------------------------------------------
q_df_prates <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_df_prates

#' 
#' **3.10b** (don't run this for the workshop)
## ---- eval=FALSE---------------------------------------------
## q_df_ordered <- q_df_long %>%
##   # assign the population assignment according to the max q value (ancestry proportion) and include the assignment probability of the population assignment
##   group_by(individual) %>%
##   mutate(likely_assignment = pop[which.max(q)],
##          assignment_prob = max(q)) %>%
##   # arrange the data set by the ancestry coefficients
##   arrange(likely_assignment, assignment_prob) %>%
##   # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
##   ungroup() %>%
##   mutate(individual = forcats::fct_inorder(factor(individual)))
## 
## q_df_ordered

#' 
#' 
#' Time to plot! I'm specifying the palette and legend based on *Prates et al.*, who did a great job of making all of their code available on [github](https://github.com/ivanprates/2018_Anolis_EcolEvol). I've also included a hashed-out line that specifies a similar color palette that is applicable in case you want to experiment with different K values and don't want to mess with specifying your color palette over and over.  
#' 
#' **3.11**
## ------------------------------------------------------------

# a custom palette for plotting
q_palette <- c("#fde725", "#35b779", "#440154")


q_df_prates %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette, labels = c("AF", "Eam", "Wam")) +
  #scale_fill_viridis_d() +
  labs(fill = "Region") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
        )


#' 
#' This is very close to *Prates et al.*!
#' 
#' ![](images/prates-fig.png)
#' 
#' ## 4. Mapping your results
#' 
#' Now we'll plot our population structure onto a map to see how genetic diversity is geographically structured. First we have to read in the locality data, which contains latitude and longitude coordinates, then we're merging the locality data and the dapc results by the individual IDs to facilitate plotting. Next, we read in a basemap of South America to plot the locality data on top of. Finally, we plot everything with ggplot2! There are many options and packages to facilitate making maps in R, and I've provided a few links at the bottom that you can explore further.
#' 
#' **4.1**
## ------------------------------------------------------------
#read in the locality data and merge it with the genetic data
anole_locs <-
  read_csv(
    "https://github.com/ivanprates/2018_Anolis_EcolEvol/raw/master/maps/qmatrix_punctatus.csv"
  ) %>%
  mutate(individual = str_c("punc_", ID)) %>% #Have to change the ID column to reflect our genetic individual ID column
  full_join(dapc_data_df, by = "individual") #join the two datasets by the "individual" category


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

#' 
#' 
#' ## Resources
#' 
#' [Introductory R and RStudio](https://rstudio.cloud/learn/primers)
#' 
#' [Drawing maps in R](https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html)
#' 
#' [DAPC tutorial](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf)
#' 
#' [Basic population genetic statistics in R](https://popgen.nescent.org/StartSNP.html)
#' 
#' 
#' 
#' 
