#### way from https://www.molecularecologist.com/2015/04/23/procrustes-analyses-in-r/ ####
#Load packages
library(MCMCpack)
library(ggplot2)

#read data in
nov<-read.csv("pca_stats/PC_mix.csv",header=TRUE)
X<-as.matrix(cbind(nov$PC1a,nov$PC2a))
Xstar<-as.matrix(cbind(nov$PC1r,nov$PC2r))
#Procrustes transformation
p<-procrustes(Xstar,X,translation=TRUE,dilation=TRUE)

#X= aciurina PC scores
#Xstar+ rbrush PC scores 
#The procrustes function will perform the Procrustes transformation to align Xstar to X through rotation, translation, and scaling.

# Access the transformed matrix and other results
X.new <- p$X.new
R <- p$R
tt <- p$tt
s <- p$s

####plot the original vs transformed rabbitbrush data####
# Create a data frame for the original data (Xstar) and the transformed data (X.new)
data_original <- data.frame(PC1 = Xstar[, 1], PC2 = Xstar[, 2], Type = "Original")
data_transformed <- data.frame(PC1 = p$X.new[, 1], PC2 = p$X.new[, 2], Type = "Transformed")

# Combine the data frames
combined_data <- rbind(data_original, data_transformed)

# Create a scatter plot
ggplot(combined_data, aes(x = PC1, y = PC2, color = Type)) +
  geom_point() +
  labs(title = "Procrustes Transformation",
       x = "PC1", y = "PC2",
       color = "Data Type") +
  theme_minimal()

####plot the data from X (aciurina) and Xstar (rbrush) before the Procrustes transformation####
# Create data frames for 'X' and 'Xstar'
data_X <- data.frame(PC1 = X[, 1], PC2 = X[, 2], Type = "X",  Population = nov$pop_a)
data_Xstar <- data.frame(PC1 = Xstar[, 1], PC2 = Xstar[, 2], Type = "Xstar", Population = nov$pop_r)

# Combine the data frames
combined_data <- rbind(data_X, data_Xstar)

# Create a scatter plot
ggplot(combined_data, aes(x = PC1, y = PC2, color = Population, shape = Type)) +
  geom_point() +
  labs(title = "Original Data Comparison",
       x = "PC1", y = "PC2",
       color = "Population", shape = "Type") +
  theme_minimal()
####THIS ONE####
#plot the data from X (aciurina) and Xstar (rbrush) after the Procrustes transformatioN
# Create data frames for 'X' and 'Xstar' after transformation
data_X <- data.frame(PC1 = X[, 1], PC2 = X[, 2], Type = "Aciurina", Population = nov$pop_a)
data_X_transformed <- data.frame(PC1 = p$X.new[, 1], PC2 = p$X.new[, 2], Type = "Ericameria", Population = nov$pop_r)

# Combine the data frames
combined_data <- rbind(data_X, data_X_transformed)

# Create a scatter plot
ggplot(combined_data, aes(x = PC1, y = PC2, color = Population, shape = Type)) +
  geom_point() +
  labs(title = "Transformed Data Comparison",
       x = "PC1", y = "PC2",
       color = "Population", shape = "Type") +
  theme_minimal()

# Combine the data frames
combined_data <- rbind(data_X, data_X_transformed)

####making final figure####
legend_colors <- c("#009E73", "#E69F00", "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Create a scatter plot
ggplot(combined_data, aes(x = PC1, y = PC2, color = Population, shape = Type)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = legend_colors) +  # Set color scale using your legend_colors
  labs(x = "PC1", y = "PC2",
       color = "Population", shape = "Type") +
  theme_minimal() +  # Start with the minimal theme
  theme(  # Customize the theme
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove background shading
  )
