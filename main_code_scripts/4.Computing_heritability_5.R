##############################################################################
# ---- Summarizing GCTA results in a table

load("files_gtca.RData")
heritability <- matrix(nrow = 44, ncol = 14)
rownames(heritability) <- colnames(pheno_file)[c(3:46)]
colnames(heritability) <- c("V(G)","V(G)_SE","V(e)","V(e)_SE","Vp","Vp_SE",
                            "V(G)/Vp", "V(G)/Vp_SE","logL","logL0","LTR",
                            "df","Pval","n")
for (i in 1:44) {
  output <- read.table(
    paste("HELIX_GWAS_FINAL_genotyped_XYchr_without_Y_MT_996_individuals_",
          i,".hsq", sep = ""),
    fill = TRUE,header = T)
  for (j in 1:14) {
      heritability[i,j] <- round(output[1,j+1], digits = 3)
  }
}

# Single urine metabolites
write.table(heritability, 
            file ="summary_heritability_results_by_metabolite.txt",
            sep = "\t", quote = F)


##############################################################################
# ---- Creating circular plot 1

# Loading library
library(tidyverse)
library(viridis)

# Obtaining metabolites for which the likelihood ratio test is
# significant (p-value < 0.05)
sign <- ifelse(heritability[,13] < 0.05, "P-value < 0.05", "P-value > 0.05")

# Introducing line breaks to make metabolite names labels shorter
row.names(heritability)[28] <- "X3.hydroxybutyrate-\n3.aminoisobutyrate"
row.names(heritability)[38] <- "N.methyl-\n2.pyridone-\n5.carboxamide"
row.names(heritability)[30] <- "N.acetyl-\nneuraminic.acid"
row.names(heritability)[34] <- "Proline.\nbetaine"
row.names(heritability)[42] <- "N.methylpicolinic.\nacid"

# Create dataset
data <- data.frame(
  individual=gsub(".urine", "", row.names(heritability), fixed = T),
  group=sign,
  heritability[,c(1,3,5)]
)

row.names(data) <- NULL

# Transform data in a tidy format (long format)
data <- data %>% gather(key = "observation", value="value", -c(1,2)) 

# Setting a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$observation))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType,
                             ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(group, individual)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)

# Getting the name and the y position of each label
label_data <- data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Preparing a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# Preparing a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Making the plot
p <- ggplot(data) +      
  
  # Adding the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=observation), 
           stat="identity", alpha=0.5) +
  scale_fill_viridis(discrete=TRUE) +
  

  ylim(-10,max(label_data$tot, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  coord_polar() +

  # Adding labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=tot+1, label=individual,hjust=hjust),
            color="black", fontface="bold",alpha=0.9, size=6, 
            angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Adding base line information
  geom_segment(data=base_data, aes(x = start, y = -1, xend = end, yend = -1),
               colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -5, label=group), hjust=c(1,0),
            colour = "black", alpha=0.8, size=8, fontface="bold", 
            inherit.aes = FALSE)


# Save at png
ggsave(p, file="circular_plot_heritability_variation.png", width=20, height=20)



##############################################################################
# ---- Creating circular plot 2 (just heritability)

# Changing metabolite names

row.names(heritability)[28] <- "3-hydroxybutyrate-\n3-aminoisobutyrate"
row.names(heritability)[38] <- "N-methyl-2-pyridone-\n5-carboxamide"
row.names(heritability)[39] <- "p-hydroxy-\nphenylacetate"
row.names(heritability)[30] <- "N-acetyl-\nneuraminic acid"
row.names(heritability)[34] <- "Proline\nbetaine"
row.names(heritability)[42] <- "N-methylpicolinic\nacid"

names <- gsub("X", "", row.names(heritability))
names <- gsub(".acid", " acid", names, fixed = T)
names <- gsub(".", "-", names, fixed = T)

# Loading libraries
library(tidyverse)

# Creating dataset

data <- data.frame(
  individual=gsub("-urine", "", names, fixed = T),
  group=sign,
  value = heritability[,c(7)]
)
row.names(data) <- NULL

# Setting a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)))
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Getting the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Preparing a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# Preparing a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Making the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), 
           stat="identity", alpha=0.5) +
  
   # Add a val=100/75/50/25 lines.
  geom_segment(data=grid_data, aes(x = 48, y = 0.80, xend = 1, yend = 0.80),
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = 48, y = 0.60, xend = 1, yend = 0.60),
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = 48, y = 0.40, xend = 1, yend = 0.40),
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = 48, y = 0.20, xend = 1, yend = 0.20),
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = 48, y = 1, xend = 1, yend = 1),
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
   # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),5), y = c(0.20, 0.40, 0.60, 0.80,1), 
           label = c("0.2", "0.4", "0.6", "0.8","1") , color="grey", size=6 ,
           angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), 
           stat="identity", alpha=1) +
  ylim(-0.75,2) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+0.1, label=individual, 
                                 hjust=hjust), 
            color="black", fontface="bold",alpha=0.9, size=6, 
            angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.1, xend = end, 
                                   yend = -0.1), colour = "black", 
               alpha=0.8, size=0.8 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title+1, y = -0.2, label=group),
            hjust=c(1,0), colour = "black", alpha=0.8, size=7, 
            fontface="bold", inherit.aes = FALSE)


p
# Saving at png
ggsave(p, file="circular_plot_just_heritab.png", width=20, height=20)
