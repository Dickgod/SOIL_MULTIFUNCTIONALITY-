
setwd("C:/Users/dickson/Desktop/Mauki/R CLASS/Directory/r4ds/SMF")
#install.packages("remotes")
#remotes::install_version("htmltools", version = "0.5.8")
#install.packages("piecewiseSEM")
#install.packages("htmltools")
#install.packages("GGally")
#devtools::install_github("houyunhuang/ggcor")
#install.packages("remotes")
#remotes::install_github("Hy4m/linkET")
library(GGally)
library(ggplot2)
library(ggcorrplot)
library(igraph)
library(data.table)
library(mice)
library(modEvA)####plotting glm
library(rgbif)
library(tidyverse)
library(tidySEM)
library(lavaan)
library(corrplot)
library(ggplot2)
library(readxl)
library(tibble)
library(readxl)
library(readr)
library(factoextra)
library(ggpubr)
library(MuMIn)
library(semPlot)
library(lme4)
library(lmtest)
library(lmerTest)
library(performance)
library(see)
library(lmtest)
library(car)
library(data.table)
library(FD)
library(readxl)
library(FactoMineR)
library(ggfortify)
library(caret)
library(performance)
library(see)
library(gridExtra)
library(effectsize)
library(sjPlot)
library(semEff)
library(shiny)
library(pscl)
library(lubridate)
library(RColorBrewer)
library(multcompView)
library(readr)
library(viridis)
library(data.table)
library(readxl)
library(readr)
library(dplyr)
library(mice)
library(ggplot2)
library(factoextra)
library(ggpubr)
library(gridExtra)
library(ggrepel)
library(piecewiseSEM)
library(ggcorrplot)
library(rdacca.hp)
library(linkET)
library(csnl)
library(linkET)
library(dplyr)
library(vegan)
library(lattice)
library(permute)
library(ggplot2)
library(FD)
library(vegan)
library(lattice)
library(permute)
###############################################

#Import agregate stability data
Agregates <- read_excel("Agregates.xlsx")
#Import soil data
#SP3_all_data <- read_excel("Soil_properties.xls")
#import GHG data
GHG_Emmision <- read_excel("Emmission.xlsx")
#import decomposition data
Decomposition_rate <- read_excel("Decomposition_rate.xlsx")
#import leaching data
Leaching_data <- read_excel("Leaching_data.xlsx")
#import isotopy data
Isotopy_data <- read_excel("Isotopy_data.xlsx")
#import enzymes data
Enzymes_data <- read_excel("Enzymes_data.xlsx")
##

# Normalize the PlotId column to uppercase for consistent merging
Agregates$PlotID <- toupper(Agregates$PlotID)
#SP3_all_data$PlotID <- toupper(SP3_all_data$PlotID)
GHG_Emmision$PlotID <- toupper(GHG_Emmision$PlotID)
Decomposition_rate$PlotID <- toupper(Decomposition_rate$PlotID)
Leaching_data$PlotID <- toupper(Leaching_data$PlotID)
Isotopy_data$PlotID <- toupper(Isotopy_data$PlotID)
Enzymes_data$PlotID <- toupper(Enzymes_data$PlotID)


# Ensure unique PlotID within each dataset
Agregates <- Agregates %>% distinct(PlotID, .keep_all = TRUE)
#SP3_all_data <- SP3_all_data %>% distinct(PlotID, .keep_all = TRUE)
GHG_Emmision <- GHG_Emmision %>% distinct(PlotID, .keep_all = TRUE)
Decomposition_rate <- Decomposition_rate %>% distinct(PlotID, .keep_all = TRUE)
Leaching_data <- Leaching_data %>% distinct(PlotID, .keep_all = TRUE)
Isotopy_data <- Isotopy_data %>% distinct(PlotID, .keep_all = TRUE)
Enzymes_data <- Enzymes_data %>% distinct(PlotID, .keep_all = TRUE)

# Merge the dataframes on the PlotID column
SMF_merged_data1 <- full_join(Enzymes_data, GHG_Emmision, by = "PlotID")
SMF_merged_data2 <- full_join(SMF_merged_data1, Decomposition_rate, by = "PlotID")
SMF_merged_data3 <- full_join(SMF_merged_data2, Leaching_data, by = "PlotID")
SMF_merged_data4 <- full_join(SMF_merged_data3, Isotopy_data, by = "PlotID")
SMF_merged_data5 <- full_join(SMF_merged_data4, Agregates, by = "PlotID")

#View(SMF_merged_data5)
########################
# Delete rows 61 to 71 from the merged data
SMF_merged_data <- SMF_merged_data5 %>% slice(-c(61:65))

# View the merged data
#View(SMF_merged_data)

# write the merged data to a new CSV file
write.csv(SMF_merged_data, "Merged_SMF_Data.csv", row.names = FALSE)
###################################################################################


# Select columns representing enzymes
# Assuming the enzyme columns are named "Enzyme1", "Enzyme2", etc.
enzymes_variables <- c("Glucosidase", "Phosphotase", "Urease", "Cellobiose") # Replace with actual column names

# Ensure the enzyme columns exist in the dataframe
enzymes_variables <- enzymes_variables[enzymes_variables %in% colnames(SMF_merged_data)]

# Normalize the enzyme data if necessary (e.g., scale to 0-1)
SMF_merged_data[enzymes_variables] <- lapply(SMF_merged_data[enzymes_variables], function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

# Calculate the Enzyme Index as the average of the selected enzyme columns
SMF_merged_data$Enzyme_Index <- rowMeans(SMF_merged_data[enzymes_variables], na.rm = TRUE)

# View the updated data with the Enzyme Index
#View(SMF_merged_data)

# Optionally, write the updated data to a new CSV file
#write.csv(SMF_merged_data, "SMF_Merged_Data_With_Enzyme_Index.csv", row.names = FALSE)


##################################
# Select columns representing
# Assuming the enzyme columns are named "Enzyme1", "Enzyme2", etc.
enzymes_variables <- c("Glucosidase", "Phosphotase", "Urease", "Cellobiose") # Replace with actual column names

# Ensure the enzyme columns exist in the dataframe
enzymes_variables <- enzymes_variables[enzymes_variables %in% colnames(SMF_merged_data)]

# Normalize the enzyme data if necessary (e.g., scale to 0-1)
SMF_merged_data[enzymes_variables] <- lapply(SMF_merged_data[enzymes_variables], function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))

# Calculate the Enzyme Index as the average of the selected enzyme columns
SMF_merged_data$Enzyme_Index <- rowMeans(SMF_merged_data[enzymes_variables], na.rm = TRUE)

# View the updated data with the Enzyme Index
#View(SMF_merged_data)

# Optionally, write the updated data to a new CSV file
#write.csv(SMF_merged_data, "SMF_Merged_Data_With_Enzyme_Index.csv", row.names = FALSE)
#################################################################



####CLIMATE REGULATION AND CARBON STORAGE


# Select relevant columns
Climate_regulation <- SMF_merged_data %>%
  select(FluxCH4,FluxN2O, FluxCO2, SOC)

# Normalize the Carbon Storage directly
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

#Normalize Carbon Storage and other indicators
normalized_SOC <- normalize(Climate_regulation$SOC)
#normalized_Decomposition_rate <- normalize(Climate_regulation$Decomposition_rate/yr)
normalized_FluxCH4 <- normalize(Climate_regulation$FluxCH4)
normalized_FluxN2O <- normalize(Climate_regulation$FluxN2O)
normalized_FluxCO2 <- normalize(Climate_regulation$FluxCO2)

# Invert the fluxes (assuming higher fluxes are worse for climate regulation)
normalized_FluxCH4 <- 1 - normalized_FluxCH4
normalized_FluxN2O <- 1 - normalized_FluxN2O
normalized_FluxCO2 <- 1 - normalized_FluxCO2
Climate_regulation$normalized_FluxCH4<-normalized_FluxCH4
Climate_regulation$normalized_FluxN2O<-normalized_FluxN2O
Climate_regulation$normalized_FluxCO2<-normalized_FluxCO2
Climate_regulation$PlotID<-SMF_merged_data$PlotID
Climate_regulation$Ecot<-SMF_merged_data$Ecot
################
Climate_regulation <- na.omit(Climate_regulation)
anova <- aov(normalized_FluxCH4 ~ Ecot, data = Climate_regulation)
summary(anova)
# Tukey's test and compact letter display
Tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, Tukey)

# Table with the mean, the standard deviation and the letters indications significant differences for each treatment
dt <- group_by(Climate_regulation, Ecot) %>%
  summarise(normalized_FluxCH4_mean=mean(normalized_FluxCH4), se=sd(normalized_FluxCH4) / sqrt(n())) %>%
  arrange(desc(normalized_FluxCH4_mean))
cld <- as.data.frame.list(cld$Ecot)
dt$Tukey <- cld$Letters
#Tukey<-dt$Tukey
print(dt)
#############################

# Define the desired order of Ecot levels
desired_order <- c("Helichrysum", "Erica", "Dist_Podocarpus", "Podocarpus", "Dist_Ocotea", "Ocotea", "Lower_Montane", "Grassland", "Homegarden","Coffee", "Maize", "Savanna")

Climate_regulation$Ecot <- factor(Climate_regulation$Ecot, levels = desired_order)
# Use the custom color palette in scale_color_brewer
custom_palette <- colorRampPalette(c("#92C5DE","#2166AC", "#B2182B"))(length(unique(SMF_merged_data$Ecosystems)))
# Reorder the levels of Ecot factor variable according to the desired order
dt$Ecot <- factor(dt$Ecot, levels = desired_order)

# Plotting the barplot
OSMF<-ggplot(dt, aes(x = Ecot, y = normalized_FluxCH4_mean, fill = Ecot)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = normalized_FluxCH4_mean + se, ymin = normalized_FluxCH4_mean - se),
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label=Tukey, y = normalized_FluxCH4_mean + se+0.03), size = 3, color = "Gray25",
            show.legend = FALSE,
            position = position_dodge(0.9)) +
  ylim(0,0.6)+
  scale_fill_manual(values = custom_palette) +  # Use scale_fill_manual for custom palette
  theme_classic()+
  labs(x = "", y = "Overall Soil Multifunctionality Index", title = "OVERALL SOIL MULTIFUNCTIONALITY") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = "black"),  # X-axis text settings
        axis.text.y = element_text(size = 12, color = "black"),  # Y-axis text settings
        axis.title.y = element_text(size = 12, face = "bold", color = "black"), # Y-axis title settings
        plot.title = element_text(size = 14, face = "bold", color = "black"),   # Title settings
        axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=FALSE)
OSMF

######################################################################################################

# Combine normalized indicators into a data frame
normalized_indicators <- data.frame(
  normalized_SOC,
  #normalized_Decomposition_rate,
  normalized_FluxCH4,
  normalized_FluxN2O,
  normalized_FluxCO2
)

# Calculate the Climate Regulation Index as the average of the normalized indicators
climate_regulation_index <- rowMeans(normalized_indicators, na.rm = TRUE)

# Add the index to the original data
SMF_merged_data$CR_Index <- climate_regulation_index

# View the updated data frame
print(SMF_merged_data)

# View the updated data
View(SMF_merged_data)

# Optionally, write the updated data to a new CSV file
write.csv(data, "SMF_merged_data.csv", row.names = FALSE)
#########################################


####################NUTRIENTS RECYCLING

# Select relevant columns for mineralization
nu_mineralization_indicators <- SMF_merged_data %>%
  select(`15N/14N`,`12C/13C`, Cellobiose, Glucosidase, Urease, Phosphotase,`conc(N-NH4µg/gmresin)`,
         `conc(N-NO3µg/gmresin)` ,
         Decomposition_rate_yr, FluxCO2,  FluxCH4,FluxN2O)

# Normalize the indicators (Min-Max Normalization)
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Normalize mineralization indicators
normalized_Decomposition <- normalize(nu_mineralization_indicators$Decomposition_rate_yr)
normalized_Cisotopy <- normalize(nu_mineralization_indicators$`12C/13C`)
normalized_Nisotopy <- normalize(nu_mineralization_indicators$`15N/14N`)
normalized_Cellobiose <- normalize(nu_mineralization_indicators$Cellobiose)
normalized_Glucosidase <- normalize(nu_mineralization_indicators$Glucosidase)
normalized_Urease <- normalize(nu_mineralization_indicators$Urease)
normalized_Phosphotase <- normalize(nu_mineralization_indicators$Phosphotase)
normalized_Ammonia <- normalize(nu_mineralization_indicators$`conc(N-NH4µg/gmresin)`)
normalized_nitrate <- normalize(nu_mineralization_indicators$`conc(N-NO3µg/gmresin)`)
normalized_FluxCO2 <- normalize(nu_mineralization_indicators$FluxCO2)
normalized_FluxCH4 <- normalize(nu_mineralization_indicators$FluxCH4)
normalized_FluxN2O <- normalize(nu_mineralization_indicators$FluxN2O)
# Combine normalized indicators into a data frame
normalized_nu_mineralization_indicators <- data.frame(
  normalized_Decomposition,
  normalized_Cisotopy,
  normalized_Nisotopy,
  normalized_Cellobiose,
  normalized_Glucosidase ,
  normalized_Urease,
  normalized_Phosphotase,
  normalized_Ammonia,
  normalized_nitrate,
  normalized_FluxCO2,
  normalized_FluxCH4,
  normalized_FluxN2O
  
)

# Assign weights to each indicator (assuming equal weights for this example)
#weights_nu_mineralization_indicators <- rep(1 / ncol(normalized_nu_mineralization_indicators), ncol(normalized_nu_mineralization_indicators))

# Calculate the Mineralization Index
Nu_mineralization_index <- rowMeans(normalized_nu_mineralization_indicators, na.rm = TRUE)

# Add the indices to the original data
SMF_merged_data$NUM_index <- Nu_mineralization_index

# View the updated data
View(SMF_merged_data)

# Optionally, write the updated data to a new CSV file
write.csv(data, "SMF_merged_data.csv", row.names = FALSE)


#############################################

####WATER REGULATION

# Select relevant columns for water regulation
water_regulation_indicators <- SMF_merged_data %>%
  select(Clay, Sand, Silt, Mean_aggregates, `conc(N-NH4µg/gmresin)`, `conc(N-NO3µg/gmresin)`)

# Normalize the indicators (Min-Max Normalization)
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Normalize water regulation indicators
normalized_Clay <- normalize(water_regulation_indicators$Clay)
normalized_Sand <- normalize(water_regulation_indicators$Sand)
normalized_Silt <- normalize(water_regulation_indicators$Silt)
normalized_Mean_aggregates <- normalize(water_regulation_indicators$Mean_aggregates)
normalized_Ammonia <- normalize(water_regulation_indicators$`conc(N-NH4µg/gmresin)`)
normalized_nitrate <- normalize(water_regulation_indicators$`conc(N-NO3µg/gmresin)`)

# Invert the ammonia and nitrate indicators (assuming higher values are negative contributors)
normalized_Ammonia <- 1 - normalized_Ammonia
normalized_nitrate <- 1 - normalized_nitrate
normalized_Sand <- 1 - normalized_Sand

# Combine normalized indicators into a data frame
normalized_water_regulation_indicators <- data.frame(
  normalized_Clay,
  normalized_Sand,
  normalized_Silt,
  normalized_Mean_aggregates,
  normalized_Ammonia,
  normalized_nitrate
)

# Calculate the Water Regulation Index as the average of the normalized indicators
water_regulation_index <- rowMeans(normalized_water_regulation_indicators, na.rm = TRUE)

# Add the index to the original data
SMF_merged_data$WR_Index <- water_regulation_index

# View the updated data
View(SMF_merged_data)

# Optionally, write the updated data to a new CSV file
write.csv(data, "Data_With_Water_Regulation_Index.csv", row.names = FALSE)
#################################


##PRIMARY PRODUCTIVITY### SOIL QUALITY INDEX FROM NEEMA

# Select relevant columns for primary productivity
primary_productivity_indicators <- SMF_merged_data %>%
  select(SON, SOC, SPH, PH, `K_tot(mg/kg)`,Cellobiose, Glucosidase, Urease, Phosphotase, Decomposition_rate/yr
  )

# Normalize the indicators (Min-Max Normalization)
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Normalize primary productivity indicators
normalized_primary_productivity <- primary_productivity_indicators %>%
  mutate(across(everything(), normalize))

# Invert negative contributors for primary productivity (if necessary)
# Extreme pH values are negative contributors, so invert them if necessary
# Example: If pH is too low or too high, create a factor that represents optimal pH contribution
#optimal_pH <- 1 - abs(normalized_primary_productivity$Soil_pH - 7) / 7

# Assign weights to each indicator (assuming equal weights for this example)
#weights_primary_productivity <- rep(1 / ncol(normalized_primary_productivity), ncol(normalized_primary_productivity))


# Calculate the Primary Productivity Index
primary_productivity_index <- rowMeans(normalized_primary_productivity, na.rm = TRUE)

# Add the index to the original data
SMF_merged_data$PP_Index <- primary_productivity_index

# View the updated data
View(SMF_merged_data)

# Optionally, write the updated data to a new CSV file
write.csv(data, "SMF_merged_data.csv", row.names = FALSE)
#######################################################

###BIODIVERSITY CONSERVATION

# Select relevant columns for biodiversity conservation
biodiversity_indicators <- data %>%
  select(Microbial_Biomass_Carbon, Microbial_Biomass_Nitrogen, Microbial_Diversity_Index, Enzyme_Activity1, Enzyme_Activity2, Enzyme_Activity3, Enzyme_Activity4, Enzyme_Activity5, Soil_Organic_Carbon, Aggregate_Stability, Porosity, Soil_Moisture, Nitrogen_Content, Phosphorus_Content, Potassium_Content)

###Adding data from soil earthworms

# Normalize the indicators (Min-Max Normalization)
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Normalize biodiversity indicators
normalized_biodiversity <- biodiversity_indicators %>%
  mutate(across(everything(), normalize))

# Invert negative contributors for biodiversity conservation (if necessary)
# Example: If soil compaction is included, it would be a negative contributor and should be inverted

# Assign weights to each indicator (assuming equal weights for this example)
weights_biodiversity <- rep(1 / ncol(normalized_biodiversity), ncol(normalized_biodiversity))

# Calculate the Biodiversity Conservation Index
biodiversity_conservation_index <- rowSums(as.matrix(normalized_biodiversity) * weights_biodiversity)

# Add the index to the original data
data$Biodiversity_Conservation_Index <- biodiversity_conservation_index

# View the updated data
View(data)

# Optionally, write the updated data to a new CSV file
write.csv(data, "Data_With_Biodiversity_Conservation_Index.csv", row.names = FALSE)
#####################################################################



###############PCA OF THE BUNDLES##

# Select relevant columns for PCA (assuming you have already calculated indices)
pca_data <- SMF_merged_data %>%
  select(`15N/14N`,`12C/13C`, Cellobiose, Glucosidase, Urease, Phosphotase,`conc(N-NH4µg/gmresin)`,
         `conc(N-NO3µg/gmresin)` ,
         Decomposition_rate_yr, FluxCO2,  FluxCH4,FluxN2O,CN, NP, `K(mmol/L)`, `BS%`,Mean_aggregates, AWC.x)


#'conc(N-NO3µg/gmresin)', 'conc(N-NH4µg/gmresin)',
# Perform PCA
for (col in names(pca_data)) {
  pca_data[[col]][is.na(pca_data[[col]]) | is.infinite(pca_data[[col]])] <- median(pca_data[[col]], na.rm = TRUE)
}
#view(pca_data)
# Perform PCA
pca_data<-scale(pca_data)
#view(pca_data)

pca_result <- prcomp(pca_data, scale = TRUE)

#######################
fviz_pca_var(pca_result, col.var = "black")
# Visualize PCA variables for PC3 and PC4
#fviz_pca_var(pca_result,
#axes = c(3, 4),  # Specify to plot PC3 (axis 3) vs PC4 (axis 4)
#col.var = "black")  # Color of the variable arrows

# The plot will display automatically in R

# Contributions of variables to PC1
fviz_contrib(pca_result, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(pca_result, choice = "var", axes = 2, top = 10)
ind <- get_pca_ind(pca_result)
ind
head(ind$coord)

fviz_pca_ind(pca_result, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
##########################
SMF_merged_data$Precipitation<-scale(SMF_merged_data$MAP)
SMF_merged_data$Temperature<-scale(SMF_merged_data$MAT)
SMF_merged_data$Elevation<-scale(SMF_merged_data$Elev)
SMF_merged_data$Land_use_intensity<-scale(SMF_merged_data$LUI)
SMF_merged_data$Slow_fast<-scale(SMF_merged_data$PCA1_RL)
SMF_merged_data$Woody_grassy<-scale(SMF_merged_data$PCA2_RL)
SMF_merged_data$Clay_s<-scale(SMF_merged_data$Clay)
SMF_merged_data$Sand_s<-scale(SMF_merged_data$Sand)
SMF_merged_data$PH_s<-scale(SMF_merged_data$PH)
# Create a data frame with PCA results
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Land_Use_Intensity <- SMF_merged_data$Land_use_intensity
pca_scores$Precipitation <- SMF_merged_data$Precipitation
pca_scores$Temperature <- SMF_merged_data$Temperature
pca_scores$Elevation<-SMF_merged_data$Elevation
pca_scores$Slow_fast<-SMF_merged_data$Slow_fast
pca_scores$Woody_grassy<-SMF_merged_data$Woody_grassy
pca_scores$Ecosystem <- SMF_merged_data$Ecosystems
pca_scores$Clay<-scale(SMF_merged_data$Clay_s)
pca_scores$Sand<-scale(SMF_merged_data$Sand_s)
pca_scores$PH<-scale(SMF_merged_data$PH_s)
###################################################

####USING ENV FIT PACKAGE########################

################
#Envi_data<-SMF_merged_data[,c("MAT","MAP","Land_use_intensity","Elev", "PCA1_RL", "PCA2_RL", "Clay", "Sand", "PH"), drop=FALSE]
#names(Envi_data)<-c("Temperature", "Precipitation","Land use intensity", "Elevation", "Slow_fast", "Woody_grassy", "Clay", "Sand", "Soil_pH")

Envi_data<-SMF_merged_data[,c("Temperature", "Precipitation","Land_use_intensity", "Elevation", "Slow_fast", "Woody_grassy", "Clay_s", "Sand_s", "PH_s"), drop=FALSE]
#Envi_data<-scale(Envi_data)
Envi_data <- as.data.frame(Envi_data)
################
env_fit <- envfit(pca_result, Envi_data, perm = 999, na.rm = TRUE)
# Add environmental variables
env_scores <- as.data.frame(scores(env_fit, display = "vectors"))
env_scores <- cbind(env_scores, rownames(env_scores))
env_scores
write.csv(env_scores, file = 'env_scoresSMF', row.names = FALSE)
colnames(env_scores) <- c("PC1", "PC2", "Variable")
# Add environmental variables
env_scores <- as.data.frame(scores(env_fit, display = "vectors"))
env_scores <- cbind(env_scores, rownames(env_scores))
colnames(env_scores) <- c("PC1", "PC2", "Variable")
# Create the PCA plot with variable vectors
# Define the desired order of Ecot levels
desired_order <- c("Helichrysum", "Erica", "Dist_Podocarpus", "Podocarpus", "Dist_Ocotea", "Ocotea", "Lower_Montane", "Grassland", "Homegarden","Coffee", "Maize", "Savanna")

SMF_merged_data$Ecosystems <- factor(SMF_merged_data$Ecosystems, levels = desired_order)

####"burlywood3", =  "darkolivegreen",= "deepskyblue4"))

# Use the custom color palette in scale_color_brewer
custom_palette <- colorRampPalette(c("deepskyblue4","burlywood3", "darkred"))(length(unique(SMF_merged_data$Ecosystems)))
# Reorder the levels of Ecot factor variable according to the desired order
SMF_merged_data$Ecosystems <- factor(SMF_merged_data$Ecosystems, levels = desired_order)

manual_shapes <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21, 22)
#manual_shapes <- c(0, 1, 2, 4, 5, 6, 0, 1, 2, 4, 5, 6)
p <- autoplot(pca_result, data = SMF_merged_data,shape = "Ecosystems", fill = "Ecosystems",
              loadings = TRUE,
              loadings.label = TRUE,
              loadings.colour = "black",
              loadings.label.size = 6,
              loadings.point.size = 15,
              loadings.label.colour = "black",
              loadings.label.repel = TRUE,
              loadings.size = 6,
              arrow.size= 12) +
  scale_fill_manual(values = custom_palette)+
  scale_shape_manual(values = manual_shapes)+
  theme_classic()+
  labs(title = "PCA of Soil function indices with Environmental Factors",
       x = "Principal Component 1(28.4%)",
       y = "Principal Component 2(13.6%)")+
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14), # Increase legend text size
    legend.key.size = unit(1.0, 'lines') # Increase size of legend key
  ) +
  guides(shape = guide_legend(override.aes = list(size = 3)), # Increase shape size in legend
         fill = guide_legend(override.aes = list(size = 3)))+
  geom_point(aes(shape = Ecosystems, fill = Ecosystems), size = 3)+# Ensure fill size is increased in legend+
  geom_segment(data = env_scores, aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", size = 0.8) + # Increase line thickness
  geom_text(data = env_scores, aes(x = PC1, y = PC2, label = Variable),
            hjust = 0, vjust = 1, color = "Red")+
  geom_point(aes(shape = Ecosystems, fill = Ecosystems), size = 3)
p
ggsave("SMFPCA_Plot.png", plot = p, width = 10, height = 8)
##########################################################################################
p <- autoplot(pca_result, data = SMF_merged_data,shape = "Ecosystems", fill = "Ecosystems",
              loadings = TRUE,
              loadings.label = TRUE,
              loadings.colour = "black",
              loadings.label.size = 6,
              loadings.point.size = 15,
              loadings.label.colour = "black",
              loadings.label.repel = TRUE,
              loadings.size = 6,
              arrow.size= 12) +
  scale_fill_manual(values = custom_palette)+
  scale_shape_manual(values = manual_shapes)+
  theme_classic()+
  labs(title = "PCA of Soil function indices with Environmental Factors",
       x = "Principal Component 1(28.4%)",
       y = "Principal Component 2(13.6%)")+
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14), # Increase legend text size
    legend.key.size = unit(1.0, 'lines') # Increase size of legend key
  )+
  guides(shape = guide_legend(override.aes = list(size = 3)), # Increase shape size in legend
         fill = guide_legend(override.aes = list(size = 3)))+
  geom_point(aes(shape = Ecosystems, fill = Ecosystems), size = 3)
p

##


# Create the PCA plot with ellipses
p <- autoplot(pca_result, data = SMF_merged_data, shape = "Ecosystems", fill = "Ecosystems",
              loadings = TRUE,
              loadings.label = TRUE,
              loadings.colour = "black",
              loadings.label.size = 6,
              loadings.point.size = 15,
              loadings.label.colour = "black",
              loadings.label.repel = TRUE,
              loadings.size = 6,
              arrow.size= 12) +
  scale_fill_manual(values = custom_palette) +
  scale_shape_manual(values = manual_shapes) +
  theme_classic() +
  labs(title = "PCA of Soil Function Indices without Environmental Factors",
       x = "Principal Component 1(28.4%)",
       y = "Principal Component 2(13.6%)") +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.0, 'lines')
  ) +
  guides(shape = guide_legend(override.aes = list(size = 3)),
         fill = guide_legend(override.aes = list(size = 3))) +
  geom_point(aes(shape = Ecosystems, fill = Ecosystems), size = 3) +
  geom_segment(data = env_scores, aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", size = 0.8) +
  geom_text(data = env_scores, aes(x = PC1, y = PC2, label = Variable),
            hjust = 0, vjust = 1, color = "black") +
  stat_ellipse(aes(x = PC1, y = PC2, color = Management_types, fill = Ecosystems),
               type = "t", level = 0.95, alpha = 0.95, show.legend = FALSE)

# Print the plot
p
##########################################################################


# Extract PCA scores from the rda object
pca_scores <- scores(pca_result, display = "sites")

# Extract PC1 and PC2 from the PCA scores
pca_scores_df <- as.data.frame(pca_scores)  # Convert to data frame for easier handling
pca_scores_df$PC1 <- pca_scores_df[, 1]
pca_scores_df$PC2 <- pca_scores_df[, 2]

# Merge PCA scores with your grouping data (e.g., SMF_merged_data)
SMF_merged_data$PCA1 <- pca_scores_df$PC1
SMF_merged_data$PCA2 <- pca_scores_df$PC2

##################################################################################
# Define the desired order of Ecot levels
desired_order <- c("Helichrysum", "Erica", "Dist_Podocarpus", "Podocarpus", "Dist_Ocotea", "Ocotea", "Lower_Montane", "Grassland", "Homegarden","Coffee", "Maize", "Savanna")

SMF_merged_data$Ecosystems <- factor(SMF_merged_data$Ecosystems, levels = desired_order)

####"burlywood3", =  "darkolivegreen",= "deepskyblue4"))

# Use the custom color palette in scale_color_brewer
custom_palette <- colorRampPalette(c("deepskyblue4","burlywood3", "darkred"))(length(unique(SMF_merged_data$Ecosystems)))
# Reorder the levels of Ecot factor variable according to the desired order
SMF_merged_data$Ecosystems <- factor(SMF_merged_data$Ecosystems, levels = desired_order)

manual_shapes <- c(21, 22, 23, 24, 25, 21, 22, 23, 24, 25, 21, 22)

# Create scatter plot of PCA1 vs PCA2 across different ecosystems
pca1 <- ggplot(SMF_merged_data, aes(x = Ecosystems, y = PCA1, fill = Ecosystems, shape = Ecosystems, fill = Ecosystems)) +
  geom_point(size = 3) +  # Add points for each observation
  theme_classic() +  # Use a classic theme for the plot
  labs(
    title = "Distribution of PCA1 Across Ecosystems",
    x = "Ecosystem types",
    y = "Principal Component 1 (PCA1)"
  ) +
  scale_fill_manual(values = custom_palette) +  # Customize colors for ecosystems if needed
  scale_shape_manual(values = manual_shapes) +# Customize shapes for ecosystems if needed
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.0, 'lines')
  )

# Display the plot
pca1

# Print the plot
print(pca1)

##########################################
pca2 <- ggplot(SMF_merged_data, aes(x = Ecosystems, y = PCA2, fill = Ecosystems, shape = Ecosystems , fill = Ecosystems)) +
  geom_point(size = 3) +  # Add points for each observation
  theme_classic() +  # Use a classic theme for the plot
  labs(
    title = "Distribution of PCA2 Across Ecosystems",
    x = "Ecosystem types",
    y = "Principal Component 2 (PCA2)"
  ) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.0, 'lines')
  ) +
  scale_fill_manual(values = custom_palette) +  # Customize colors for ecosystems if needed
  scale_shape_manual(values = manual_shapes)  # Customize shapes for ecosystems if needed

# Print the plot
print(pca2)

####################################################


#####MULTIFUNCTIONALITY
#columns_to_keep <- c("OSMF_index" , "CSMF_index","ASMF_index")
SMF<-pca_data
################################################
####ENVIRONMENTAL DRIVERS
columns_to_keep2 <- c("MAP" , "MAT","LUI", "PCA1_RL", "PCA2_RL")
Drivers<-SMF_merged_data[,columns_to_keep2]
Envi_data<-SMF_merged_data[,c("MAT", "MAP","LUI", "Elev", "PCA1_RL", "PCA2_RL",
                              "Clay", "Sand", "PH"), drop=FALSE]
##########################MANTEL ANALYSIS
mantel <- mantel_test(Envi_data, SMF,
                      spec_select = list(MAT = 1,
                                         MAP = 2,
                                         LUI = 3,
                                         Elev=4,
                                         PCA1_RL=5,
                                         PCA2_RL=6,
                                         Clay=7,
                                         Sand=8,
                                         PH=9)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#> mantel_test() using 'bray' dist method for 'spec'.
#> mantel_test() using 'euclidean' dist method for 'env'.


qcorrplot(correlate(SMF), type = "lower", diag = TRUE) +
  geom_square() +coord_fixed(ratio = 1)+
  geom_mark(only_mark = TRUE,
            sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05, size=3, colour="black")+
  geom_couple(aes(colour = pd, size = rd),
              data = mantel,
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"),
                             order = 2),
         colour = guide_legend(title = "Mantel's p",
                               override.aes = list(size = 2),
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r",order=2))



#############################################################
####PCA
# Get log-likelihood based AIC
evaluate_model <- function(model) {
  # AIC
  aic <- AIC(model)
  # AIC with "dsep" adjustment
  aic_dsep <- AIC(model, AIC.type = "dsep")
  # ANOVA
  anova_result <- anova(model)
  # Fisher's C test
  fisher_c <- fisherC(model)
  # Model coefficients
  coefficients <- coefs(model)
  # Model summary
  model_summary <- summary(model)
  # Model plot
  model_plot <- plot(model)
  
  # Create a list to store the results
  results <- list(
    AIC = aic,
    AIC_dsep = aic_dsep,
    ANOVA = anova_result,
    Fisher_C = fisher_c,
    Coefficients = coefficients,
    Summary = model_summary,
    Plot = model_plot
  )
  # Return the results
  return(results)
}
#################################################
set.seed(12345)
#Correlation_dataset<-merge(x = pca_data, y = Envi_data, by = "common_column", all = FALSE)

#pca_data$PC1<- SMF_merged_data$PCA1
#pca_data$PC2<-SMF_merged_data$PCA2

###########################################

set.seed(1112)
# Identify rows with any NAs
rows_with_na <- apply(SMF_merged_data, 1, function(row) any(is.na(row)))
rows_with_na
# Display the rows with NAs
na_rows <- SMF_merged_data[rows_with_na,  ]
print(na_rows)
SMF_merged_data_SEM<-SMF_merged_data[, -c(81:84)]
# Inspect the structure of your dataset
str(SMF_merged_data_SEM)
pca1_smf1<-psem(
  lm(PCA1~MAT+MAP+Sand+ Clay+ PH+LUI+PCA1_RL+PCA2_RL,data=SMF_merged_data_SEM),
  lm(LUI~MAP+MAT,data =SMF_merged_data_SEM),
  lm(Sand~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(Clay~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PH~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PCA1_RL~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PCA2_RL~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  data=SMF_merged_data_SEM)

summary(pca1_smf1)
fisherC(pca1_smf1)
AIC(pca1_smf1)
summary(pca1_smf1, conserve = T, scale=TRUE)

###################update model################

pca1_smf1<-update(pca1_smf1, Clay %~~%Sand ,  PH %~~%Sand,
                  PCA2_RL %~~% Sand, PH %~~%Clay, PCA2_RL %~~% Clay, PCA1_RL%~~%PH, PCA1_RL %~~% Clay)
summary(pca1_smf1)
fisherC(pca1_smf1)
AIC(pca1_smf1)
evaluate_model(pca1_smf1)

#####CHECKING THE ASSUMPTION OF NORMALITY
A<-lm(PCA1~MAT+MAP+Elev+Sand+ Clay+ LUI+PCA1_RL+PCA2_RL,data=SMF_merged_data)
B<-lm(LUI~MAP+MAT+Elev,family=gaussian,data =SMF_merged_data)
C<-lm(Elev~MAT+MAP+Elev +LUI,data =SMF_merged_data)
C<-lm(Sand~MAT+MAP+Elev +LUI,data =SMF_merged_data)
D<-lm(Clay~MAT+MAP+Elev +LUI,data =SMF_merged_data)
E<-lm(PCA1_RL~MAT+MAP+Elev +LUI,data =SMF_merged_data)
F<-lm(PCA2_RL~MAT+MAP+Elev+LUI,data =SMF_merged_data)

check_heteroscedasticity(F)
check_normality(F)
check_model(A)
#########################################################################
pca1_smf1_coef<-summary(pca1_smf1SEMmodelAGC)$coefficients
pca1_smf1_coef
pca1_smf1_coef <- as.data.frame(pca1_smf1_coef)
write.csv(pca1_smf1_coef, file = "pca1_smf1_coef_coef_table.csv")

########################################################################

##usng rather root traits
set.seed(1332)
pca1_smf2<-psem(
  lm(PCA1~MAT+MAP+Sand+ Clay+ PH+LUI+PCA1_RT+PCA2_RT,data=SMF_merged_data_SEM),
  lm(LUI~MAP+MAT,data =SMF_merged_data_SEM),
  lm(Sand~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(Clay~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PH~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PCA1_RT~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PCA2_RT~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  data=SMF_merged_data_SEM)

summary(pca1_smf2)
fisherC(pca1_smf2)
AIC(pca1_smf2)
summary(pca1_smf2, conserve = T, scale=TRUE)

###################update model################

pca1_smf2<-update(pca1_smf2, Clay %~~%Sand ,  PH %~~%Sand,
                  PCA2_RT %~~% Sand, PH %~~%Clay, PCA1_RL%~~%PH, PCA1_RT %~~% Clay, PCA1_RT %~~% Sand, PCA1_RT %~~%PH)
summary(pca1_smf2)
fisherC(pca1_smf2)
AIC(pca1_smf2)
evaluate_model(pca1_smf2)

#####CHECKING THE ASSUMPTION OF NORMALITY
A<-lm(PCA1~MAT+MAP+Sand+ Clay+ LUI+PCA1_RT+PCA2_RT,data=SMF_merged_data)
B<-lm(LUI~MAP+MAT,family=gaussian,data =SMF_merged_data)
C<-lm(Elev~MAT+MAP+LUI,data =SMF_merged_data)
C<-lm(Sand~MAT+MAP+LUI,data =SMF_merged_data)
D<-lm(Clay~MAT+MAP+LUI,data =SMF_merged_data)
E<-lm(PCA1_RT~MAT+MAP+LUI,data =SMF_merged_data)
F<-lm(PCA2_RT~MAT+MAP+LUI,data =SMF_merged_data)

check_heteroscedasticity(F)
check_normality(F)
check_model(A)
#########################################################################
pca1_smf2_coef<-summary(pca1_smf2)$coefficients
pca1_smf2_coef
pca1_smf2_coef <- as.data.frame(pca1_smf2_coef)
write.csv(pca1_smf2_coef, file = "pca1_smf2_coef_coef_table.csv")


######################################################################
set.seed(11133)
# Identify rows with any NAs
rows_with_na <- apply(SMF_merged_data, 1, function(row) any(is.na(row)))
rows_with_na
# Display the rows with NAs
na_rows <- SMF_merged_data[rows_with_na,  ]
print(na_rows)
SMF_merged_data_SEM<-SMF_merged_data[, -c(81:84)]
# Inspect the structure of your dataset
str(SMF_merged_data_SEM)
set.seed(11133)
pca2_smf3<-psem(
  lm(PCA2~MAT+MAP+Sand+ Clay+ PH+LUI+PCA1_RL+PCA2_RL,data=SMF_merged_data_SEM),
  lm(LUI~MAP+MAT,data =SMF_merged_data_SEM),
  lm(Sand~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(Clay~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PH~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PCA1_RL~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PCA2_RL~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  data=SMF_merged_data_SEM)

summary(pca2_smf3)
fisherC(pca2_smf3)
AIC(pca2_smf3)
summary(pca2_smf3, conserve = T, scale=TRUE)

###################update model################
set.seed(22133)
pca2_smf3<-update(pca2_smf3, Clay %~~%Sand ,  PH %~~%Sand,
                  PCA2_RL %~~% Sand, PH %~~%Clay, PCA2_RL %~~% Clay, PCA1_RL%~~%PH, PCA1_RL %~~% Clay)
summary(pca2_smf3)
fisherC(pca2_smf3)
AIC(pca2_smf3)
evaluate_model(pca2_smf3)

#####CHECKING THE ASSUMPTION OF NORMALITY
A<-lm(PCA2~MAT+MAP+Elev+Sand+ Clay+ LUI+PCA1_RL+PCA2_RL,data=SMF_merged_data)
B<-lm(LUI~MAP+MAT+Elev,family=gaussian,data =SMF_merged_data)
C<-lm(Elev~MAT+MAP+Elev +LUI,data =SMF_merged_data)
C<-lm(Sand~MAT+MAP+Elev +LUI,data =SMF_merged_data)
D<-lm(Clay~MAT+MAP+Elev +LUI,data =SMF_merged_data)
E<-lm(PCA1_RL~MAT+MAP+Elev +LUI,data =SMF_merged_data)
F<-lm(PCA2_RL~MAT+MAP+Elev+LUI,data =SMF_merged_data)

check_heteroscedasticity(A)
check_normality(A)
check_model(A)
#########################################################################
pca2_smf3_coef<-summary(pca2_smf3)$coefficients
pca2_smf3_coef
pca2_smf3_coef <- as.data.frame(pca2_smf3_coef)
write.csv(pca2_smf3_coef, file = "pca2_smf3_coef_coef_table.csv")

########################################################################

##usng rather root traits
set.seed(221443)
pca2_smf4<-psem(
  lm(PCA2~MAT+MAP+Sand+ Clay+ PH+LUI+PCA1_RT+PCA2_RT,data=SMF_merged_data_SEM),
  lm(LUI~MAP+MAT,data =SMF_merged_data_SEM),
  lm(Sand~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(Clay~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PH~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PCA1_RT~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  lm(PCA2_RT~MAT+MAP+LUI,data =SMF_merged_data_SEM),
  data=SMF_merged_data_SEM)

summary(pca2_smf4)
fisherC(pca2_smf4)
AIC(pca2_smf4)
summary(pca2_smf4, conserve = T, scale=TRUE)

###################update model################

pca2_smf4<-update(pca2_smf4, Clay %~~%Sand ,  PH %~~%Sand,
                  PCA2_RT %~~% Sand, PH %~~%Clay, PCA1_RL%~~%PH, PCA1_RT %~~% Clay, PCA1_RT %~~% Sand, PCA1_RT %~~%PH)
summary(pca2_smf4)
fisherC(pca2_smf4)
AIC(pca2_smf4)
evaluate_model(pca2_smf4)

#####CHECKING THE ASSUMPTION OF NORMALITY
A<-lm(PCA2~MAT+MAP+Sand+ Clay+ LUI+PCA1_RT+PCA2_RT,data=SMF_merged_data)
B<-lm(LUI~MAP+MAT,family=gaussian,data =SMF_merged_data)
C<-lm(Elev~MAT+MAP+LUI,data =SMF_merged_data)
C<-lm(Sand~MAT+MAP+LUI,data =SMF_merged_data)
D<-lm(Clay~MAT+MAP+LUI,data =SMF_merged_data)
E<-lm(PCA1_RT~MAT+MAP+LUI,data =SMF_merged_data)
F<-lm(PCA2_RT~MAT+MAP+LUI,data =SMF_merged_data)

check_heteroscedasticity(A)
check_normality(A)
check_model(A)
#########################################################################
pca2_smf4_coef<-summary(pca2_smf4)$coefficients
pca2_smf4_coef
pca2_smf4_coef <- as.data.frame(pca2_smf4_coef)
write.csv(pca2_smf4_coef, file = "pca2_smf4_coef_coef_table.csv")

########################################################################






###############plotting individual SFUNCTIONS INDEX AGAIST VARIABLES


# Plot Nu_mineralization_index against environmental factors with correlation statistics
plot_nu_mineralization_vs_env <- function(data) {
  p1 <- ggplot(data, aes(x = LUI, y = Nu_mineralization_index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 0.05) +
    labs(title = "Nu_mineralization_index vs Land Use Intensity",
         x = "Land Use Intensity",
         y = "Nu_mineralization_index") +
    theme_minimal()
  
  p2 <- ggplot(data, aes(x = MAP, y = Nu_mineralization_index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 700) +
    labs(title = "Nu_mineralization_index vs Precipitation",
         x = "Precipitation",
         y = "Nu_mineralization_index") +
    theme_minimal()
  
  p3 <- ggplot(data, aes(x = MAT, y = Nu_mineralization_index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 3) +
    labs(title = "Nu_mineralization_index vs Temperature",
         x = "Temperature",
         y = "Nu_mineralization_index") +
    theme_minimal()
  
  list(p1, p2, p3)
}

# Plot Water_Regulation_Index against environmental factors with correlation statistics
plot_water_regulation_vs_env <- function(data) {
  p1 <- ggplot(data, aes(x = LUI, y = Water_Regulation_Index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 0.05) +
    labs(title = "Water_Regulation_Index vs Land Use Intensity",
         x = "Land Use Intensity",
         y = "Water_Regulation_Index") +
    theme_minimal()
  
  p2 <- ggplot(data, aes(x = MAP, y = Water_Regulation_Index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 700) +
    labs(title = "Water_Regulation_Index vs Precipitation",
         x = "Precipitation",
         y = "Water_Regulation_Index") +
    theme_minimal()
  
  p3 <- ggplot(data, aes(x = MAT, y = Water_Regulation_Index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 3) +
    labs(title = "Water_Regulation_Index vs Temperature",
         x = "Temperature",
         y = "Water_Regulation_Index") +
    theme_minimal()
  
  list(p1, p2, p3)
}

# Plot Primary_Productivity_Index against environmental factors with correlation statistics
plot_primary_productivity_vs_env <- function(data) {
  p1 <- ggplot(data, aes(x = LUI, y = Primary_Productivity_Index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 0.05) +
    labs(title = "Primary_Productivity_Index vs Land Use Intensity",
         x = "Land Use Intensity",
         y = "Primary_Productivity_Index") +
    theme_minimal()
  
  p2 <- ggplot(data, aes(x = MAP, y = Primary_Productivity_Index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 700) +
    labs(title = "Primary_Productivity_Index vs Precipitation",
         x = "Precipitation",
         y = "Primary_Productivity_Index") +
    theme_minimal()
  
  p3 <- ggplot(data, aes(x = MAT, y = Primary_Productivity_Index)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    stat_cor(method = "pearson", label.x = 3) +
    labs(title = "Primary_Productivity_Index vs Temperature",
         x = "Temperature",
         y = "Primary_Productivity_Index") +
    theme_minimal()
  
  list(p1, p2, p3)
}

# Generate plots
nu_mineralization_plots <- plot_nu_mineralization_vs_env(SMF_merged_data)
water_regulation_plots <- plot_water_regulation_vs_env(SMF_merged_data)
primary_productivity_plots <- plot_primary_productivity_vs_env(SMF_merged_data)

# Display plots
grid.arrange(grobs = nu_mineralization_plots, nrow = 3)
grid.arrange(grobs = water_regulation_plots, nrow= 3)
grid.arrange(grobs = primary_productivity_plots, nrow = 3)

##########################################################################

##CALCULATIING THE SCORES FOR THE ECOSYSTEM MANAGEMENT TYPES

Priority_points = data.table(NCP = c('nu_mineralization', 'water_regulation', 'primary_productivity'),
                             National_parks = c(2, 5, 7),
                             Agricultural_plots=  c(4, 3,8),
                             Savanna_grassland =  c(2,5, 8))

##################################################

Priority_relative = Priority_points[, list(NCP = NCP,
                                           National_parks = National_parks/(sum(National_parks)),
                                           Agricultural_plots = Agricultural_plots/(sum(Agricultural_plots)),
                                           Savanna_grassland = Savanna_grassland/(sum(Savanna_grassland)))]

Priority_relative
#################################################################
# Select only the ecosystem types and calculated indices
Supply_NCP <- SMF_merged_data%>%
  select(Management_types,PlotID, Nu_mineralization_index, Water_Regulation_Index, Primary_Productivity_Index)
Supply_NCP<-as.data.table(Supply_NCP)

# Replace NA values with zero in the entire data table
#Supply_NCP[is.na(Supply_NCP)] <- 0

#Management_type <- list(
# "FER" = "National_parks",
# "HEL" = "National_parks",
# "FPO" = "National_parks",
# "FPD" = "National_parks",
# "FOD" = "National_parks",
#  "FOC" = "National_parks",
# "FLM" = "National_parks",
#  "HOM" = "Agricultural_plots",
# "MAI" = "Agricultural_plots",
# "COF" = "Agricultural_plots",
# "GRA" = "Agricultural",
# "SAV" = "Savanna-grassland"
#)

# Create a new column 'Management_types' based on the mapping
#Supply_NCP[, Management_types := Management_type[as.character(Ecot)]]

# View the updated data table
print(Supply_NCP)

# Now, create the final table with calculations and scaling
Supply_NCP = Supply_NCP[, list(
  Management_types = Management_types,
  PlotID = PlotID,
  Nu_mineralization_index = Nu_mineralization_index,
  Water_Regulation_Index = Water_Regulation_Index,
  Primary_Productivity_Index = Primary_Productivity_Index
)]

Supply_NCP_long = melt.data.table(Supply_NCP, id.vars = c('PlotID', 'Management_types'),
                                  variable.name = 'NCP', value.name = 'NCP_supply')
# Assuming Supply_NCP_long is your data table
# Ensure that columns 'Management_types' and 'NCP' are not lists
Supply_NCP_long <- as.data.table(Supply_NCP_long)
Supply_NCP_long[, Management_types := unlist(Management_types)]
Supply_NCP_long[, NCP := unlist(NCP)]

# Aggregate NCP_supply by Management_types and NCP
Supply_NCP_long_aggregated <- Supply_NCP_long[, .(NCP_supply_aggregated = mean(NCP_supply, na.rm = TRUE)), by = .(Management_types, NCP)]

# View the result
print(Supply_NCP_long_aggregated)

Supply_NCP_long_aggregated = Supply_NCP_long[, list(NCP_supply_aggregated = mean(NCP_supply)), by = list(Management_types, NCP)]

Priority_long = melt.data.table(Priority_relative, id.vars = 'NCP', variable = 'Stakeholder', value.name = 'Priority_score')

Supply_NCP_long_aggregated
#####################################################################

Supply_demand = merge.data.table(Supply_NCP_long_aggregated, Priority_long, by = 'NCP', allow.cartesian=TRUE)

# I create a new column for weighted scores
Supply_demand[, Weighted_scores := NCP_supply_aggregated * Priority_score]

# We can sum weighted scores to calculate MF for each stakeholder in each plot

Multifunctionality_forest_type = Supply_demand[, list(MF = sum(Weighted_scores)), by = list(Stakeholder, Management_types)]
Multifunctionality_forest_type
##############################################################################

ggplot(Multifunctionality_forest_type, aes(y = MF, x = Stakeholder, fill = Forest_type)) +
  geom_col(position = position_dodge() )

#########################################################


#### FIRST CALCULATE OVERALL SOIL MULTIFUNCTIONALITY

# Select relevant columns for minOVERALL SOIL MULTIFUNCTIONALITY

Overall_Indicators <- SMF_merged_data %>%
  select(Enzyme_Index,`conc(N-NH4µg/gmresin)`,
         `conc(N-NO3µg/gmresin)` ,
         Decomposition_rate_yr, total_flux, `BS%`, SOC, Mean_aggregates,NP, CN,
         `K_tot(mg/kg)`, `15N/14N`, `12C/13C`, AWC.x)

# Normalize the indicators (Min-Max Normalization)
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Normalize mineralization indicators
#normalized_Cisotopy <- normalize(nu_mineralization_indicators$`12C/13C`)

normalized_Enzyme_Index <- Overall_Indicators$Enzyme_Index
normalized_Ammonia <- normalize(Overall_Indicators$`conc(N-NH4µg/gmresin)`)
normalized_Nisotopy <- normalize(Overall_Indicators$`15N/14N`)
normalized_Cisotopy <- normalize(Overall_Indicators$`12C/13C`)
normalized_nitrate <- normalize(Overall_Indicators$`conc(N-NO3µg/gmresin)`)
normalized_Decomposition_rate <- normalize(Overall_Indicators$Decomposition_rate_yr)
normalized_Tflux<- normalize(Overall_Indicators$total_flux)
normalized_BS<- normalize(Overall_Indicators$`BS%`)
normalized_Carbon <- normalize(Overall_Indicators$SOC)
normalized_aggregates<- normalize(Overall_Indicators$Mean_aggregates)
normalized_NP<- normalize(Overall_Indicators$NP)
normalized_CN<- normalize(Overall_Indicators$CN)
normalized_Pottasium<- normalize(Overall_Indicators$`K_tot(mg/kg)`)
normalized_AWC<- normalize(Overall_Indicators$AWC.x)



# Invert the ammonia and nitrate indicators (assuming higher values are negative contributors)
normalized_Tflux <- 1 - normalized_Tflux
normalized_nitrate <- 1 - normalized_nitrate
normalized_Ammonia<- 1 - normalized_Ammonia
normalized_aggregates<-1-normalized_aggregates
normalized_Nisotopy<-1-normalized_Nisotopy
normalized_Cisotopy<-1-normalized_Cisotopy
# Combine normalized indicators into a data frame
normalized_Overall_Indicators_indicators <- data.frame(
  normalized_Enzyme_Index,
  normalized_Ammonia,
  normalized_Nisotopy,
  normalized_Cisotopy,
  normalized_nitrate,
  normalized_Decomposition_rate,
  normalized_Tflux,
  normalized_BS,
  normalized_Carbon,
  normalized_aggregates,
  normalized_NP,
  normalized_CN,
  normalized_Pottasium,
  normalized_AWC
  
)
view(normalized_Overall_Indicators_indicators)
# Assign weights to each indicator (assuming equal weights for this example)
#weights_nu_mineralization_indicators <- rep(1 / ncol(normalized_nu_mineralization_indicators), ncol(normalized_nu_mineralization_indicators))

# Calculate the Mineralization Index
Overall_Indicators_SMF_index <- rowMeans(normalized_Overall_Indicators_indicators, na.rm = TRUE)

# Add the indices to the original data
SMF_merged_data$OSMF_index <- Overall_Indicators_SMF_index
plot(Overall_Indicators_SMF_index)
# View the updated data
View(SMF_merged_data)

# Optionally, write the updated data to a new CSV file
write.csv(data, "SMF_merged_data.csv", row.names = FALSE)

#####################################
# analysis of variance
anova <- aov(OSMF_index ~ factor(Ecosystems), data = SMF_merged_data)
summary(anova)
# Tukey's test and compact letter display
Tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, Tukey)

# Table with the mean, the standard deviation and the letters indications significant differences for each treatment
dt <- group_by(SMF_merged_data, Ecosystems) %>%
  summarise(OSMF_index_mean=mean(OSMF_index), se=sd(OSMF_index) / sqrt(n())) %>%
  arrange(desc(OSMF_index_mean))
cld <- as.data.frame.list(cld$`factor(Ecosystems)`)
dt$Tukey <- cld$Letters
Tukey<-dt$Tukey
print(dt)
#############################

# Define the desired order of Ecot levels
desired_order <- c("Helichrysum", "Erica", "Dist_Podocarpus", "Podocarpus", "Dist_Ocotea", "Ocotea", "Lower_Montane", "Grassland", "Homegarden","Coffee", "Maize", "Savanna")

SMF_merged_data$Ecosystems <- factor(SMF_merged_data$Ecosystems, levels = desired_order)
# Use the custom color palette in scale_color_brewer
custom_palette <- colorRampPalette(c("#92C5DE","#2166AC", "#B2182B"))(length(unique(SMF_merged_data$Ecosystems)))
# Reorder the levels of Ecot factor variable according to the desired order
dt$Ecosystems <- factor(dt$Ecosystems, levels = desired_order)

# Plotting the barplot
OSMF<-ggplot(dt, aes(x = Ecosystems, y = OSMF_index_mean, fill = Ecosystems)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = OSMF_index_mean + se, ymin = OSMF_index_mean - se),
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label=Tukey, y = OSMF_index_mean + se+0.03), size = 3, color = "Gray25",
            show.legend = FALSE,
            position = position_dodge(0.9)) +
  ylim(0,0.65)+
  scale_fill_manual(values = custom_palette) +  # Use scale_fill_manual for custom palette
  theme_classic()+
  labs(x = "", y = "Overall Soil Multifunctionality Index", title = "OVERALL SOIL MULTIFUNCTIONALITY") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = "black"),  # X-axis text settings
        axis.text.y = element_text(size = 12, color = "black"),  # Y-axis text settings
        axis.title.y = element_text(size = 12, face = "bold", color = "black"), # Y-axis title settings
        plot.title = element_text(size = 14, face = "bold", color = "black"),   # Title settings
        axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
guides(fill=TRUE)
OSMF

######################################################################################################
# analysis of variance
anova <- aov(OSMF_index ~ factor(Management_types), data = SMF_merged_data)
summary(anova)
# Tukey's test and compact letter display
Tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, Tukey)

# Table with the mean, the standard deviation and the letters indications significant differences for each treatment
dt <- group_by(SMF_merged_data, Management_types) %>%
  summarise(OSMF_index_mean=mean(OSMF_index), se=sd(OSMF_index) / sqrt(n())) %>%
  arrange(desc(OSMF_index_mean))
cld <- as.data.frame.list(cld$`factor(Management_types)`)
dt$Tukey <- cld$Letters
Tukey<-dt$Tukey
print(dt)

##########
# Define the desired order of Ecot levels

#desired_order2 <- c("Agricultural_plots", "National_parks", "Savanna_grassland")

#SMF_merged_data$Management_types <- factor(SMF_merged_data$Management_types, levels = desired_order2)
# Use the custom color palette in scale_color_brewer
#custom_palette2 <- colorRampPalette(c("#92C5DE","#2166AC", "#B2182B"))(length(unique(SMF_merged_data$Management_types)))
# Reorder the levels of Ecot factor variable according to the desired order
#dt$Management_types <- factor(dt$Management_types, levels = desired_order2)

OSMF2<-ggplot(dt, aes(x = Management_types, y = OSMF_index_mean, fill = Management_types)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = OSMF_index_mean + se, ymin = OSMF_index_mean - se),
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label=Tukey, y = OSMF_index_mean + se+0.03), size = 3, color = "Gray25",
            show.legend = FALSE,
            position = position_dodge(0.9)) +
  ylim(0,0.6)+
  scale_fill_manual(values = c("Agricultural_plots" = "#92C5DE",
                               "National_parks" = "#2166AC",
                               "Savanna_grassland" = "#B2182B"))+ # Use scale_fill_manual for custom palette
  theme_classic()+
  labs(x = "", y = "Overall Soil Multifunctionality Index", title = "OVERALL SOIL MULTIFUNCTIONALITY") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = "black"),  # X-axis text settings
        axis.text.y = element_text(size = 12, color = "black"),  # Y-axis text settings
        axis.title.y = element_text(size = 12, face = "bold", color = "black"), # Y-axis title settings
        plot.title = element_text(size = 14, face = "bold", color = "black"),   # Title settings
        axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=FALSE)
OSMF2

############################################################


# Select relevant columns for AGRONOMIC STAKEHOLDERS
Agronomic_Indicators <- SMF_merged_data %>%
  select(Enzyme_Index,
         Decomposition_rate_yr, `BS%`, SOC,NP, CN,
         `K_tot(mg/kg)`, AWC.x)

# Normalize the indicators (Min-Max Normalization)
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Normalize mineralization indicators

normalized_Enzyme_Index <- Agronomic_Indicators$Enzyme_Index
#normalized_Ammonia <- normalize(Agronomic_Indicators$`conc(N-NH4µg/gmresin)`)
#normalized_nitrate <- normalize(Agronomic_Indicators$`conc(N-NO3µg/gmresin)`)
normalized_Decomposition_rate <- normalize(Agronomic_Indicators$Decomposition_rate_yr)
#normalized_Tflux<- normalize(Agronomic_Indicators$total_flux)
normalized_BS<- normalize(Agronomic_Indicators$`BS%`)
normalized_Carbon <- normalize(Agronomic_Indicators$SOC)
#normalized_aggregates<- normalize(Agronomic_Indicators$Mean_aggregates)
normalized_NP<- normalize(Agronomic_Indicators$NP)
normalized_CN<- normalize(Agronomic_Indicators$CN)
normalized_Pottasium<- normalize(Agronomic_Indicators$`K_tot(mg/kg)`)
normalized_AWC<- normalize(Agronomic_Indicators$AWC.x)
# Invert the ammonia and nitrate indicators (assuming higher values are negative contributors)
#normalized_Tflux <- 1 - normalized_Tflux
#normalized_nitrate <- 1 - normalized_nitrate
#normalized_Ammonia<- 1 - normalized_Ammonia
#normalized_aggregates<-1-normalized_aggregates

# Combine normalized indicators into a data frame
normalized_Agronomic_Indicators <- data.frame(
  normalized_Enzyme_Index,
  normalized_Decomposition_rate,
  normalized_BS,
  normalized_Carbon,
  normalized_NP,
  normalized_CN,
  normalized_Pottasium,
  normalized_AWC
  
)
view(normalized_Agronomic_Indicators)
# Assign weights to each indicator (assuming equal weights for this example)
#weights_nu_mineralization_indicators <- rep(1 / ncol(normalized_nu_mineralization_indicators), ncol(normalized_nu_mineralization_indicators))

# Calculate the Mineralization Index
Agronomic_SMF_index <- rowMeans(normalized_Agronomic_Indicators, na.rm = TRUE)

# Add the indices to the original data
SMF_merged_data$ASMF_index <- Agronomic_SMF_index
plot(Agronomic_SMF_index)
# View the updated data
View(SMF_merged_data)

##############################################


# analysis of variance
anova <- aov(ASMF_index ~ factor(Ecosystems), data = SMF_merged_data)
summary(anova)
# Tukey's test and compact letter display
Tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, Tukey)

# Table with the mean, the standard deviation and the letters indications significant differences for each treatment
dt <- group_by(SMF_merged_data, Ecosystems) %>%
  summarise(ASMF_index_mean=mean(ASMF_index), se=sd(ASMF_index) / sqrt(n())) %>%
  arrange(desc(ASMF_index_mean))
cld <- as.data.frame.list(cld$`factor(Ecosystems)`)
dt$Tukey <- cld$Letters
Tukey<-dt$Tukey
print(dt)

##########

ASMF<-ggplot(dt, aes(x = Ecosystems, y = ASMF_index_mean, fill = Ecosystems)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = ASMF_index_mean + se, ymin = ASMF_index_mean - se),
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label=Tukey, y = ASMF_index_mean + se+0.03), size = 3, color = "Gray25",
            show.legend = FALSE,
            position = position_dodge(0.9)) +
  ylim(0,0.55)+
  scale_fill_manual(values = custom_palette) +  # Use scale_fill_manual for custom palette
  theme_classic()+
  labs(x = "", y = "Agronomic Soil Multifunctionality Index", title = "AGRONOMIC SOIL MULTIFUNCTIONALITY") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = "black"),  # X-axis text settings
        axis.text.y = element_text(size = 12, color = "black"),  # Y-axis text settings
        axis.title.y = element_text(size = 12, face = "bold", color = "black"), # Y-axis title settings
        plot.title = element_text(size = 14, face = "bold", color = "black"),   # Title settings
        axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=FALSE)
ASMF


###################################################

# analysis of variance
anova <- aov(ASMF_index ~ factor(Management_types), data = SMF_merged_data)
summary(anova)
# Tukey's test and compact letter display
Tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, Tukey)

# Table with the mean, the standard deviation and the letters indications significant differences for each treatment
dt <- group_by(SMF_merged_data, Management_types) %>%
  summarise(ASMF_index_mean=mean(ASMF_index), se=sd(ASMF_index) / sqrt(n())) %>%
  arrange(desc(ASMF_index_mean))
cld <- as.data.frame.list(cld$`factor(Management_types)`)
dt$Tukey <- cld$Letters
Tukey<-dt$Tukey
print(dt)

##########

ASMF2<-ggplot(dt, aes(x = Management_types, y = ASMF_index_mean, fill = Management_types)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = ASMF_index_mean + se, ymin = ASMF_index_mean - se),
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label=Tukey, y = ASMF_index_mean + se+0.03), size = 3, color = "Gray25",
            show.legend = FALSE,
            position = position_dodge(0.9)) +
  ylim(0,0.45)+
  scale_fill_manual(values = c("Agricultural_plots" = "#92C5DE",
                               "National_parks" = "#2166AC",
                               "Savanna_grassland" = "#B2182B"))+  # Use scale_fill_manual for custom palette
  theme_classic()+
  labs(x = "", y = "Agronomic Soil Multifunctionality Index", title = "AGRONOMIC SOIL MULTIFUNCTIONALITY") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = "black"),  # X-axis text settings
        axis.text.y = element_text(size = 12, color = "black"),  # Y-axis text settings
        axis.title.y = element_text(size = 12, face = "bold", color = "black"), # Y-axis title settings
        plot.title = element_text(size = 14, face = "bold", color = "black"),   # Title settings
        axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=FALSE)
ASMF2


# Select relevant columns for CONSERVATION  STAKEHOLDERS
Conservation_Indicators <- SMF_merged_data %>%
  select(`conc(N-NH4µg/gmresin)`,
         `conc(N-NO3µg/gmresin)` ,
         Decomposition_rate_yr, total_flux, SOC, Mean_aggregates)

# Normalize the indicators (Min-Max Normalization)
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Normalize mineralization indicators

#normalized_Enzyme_Index <- Agronomic_Indicators$Enzyme_Index
normalized_Ammonia <- normalize(Conservation_Indicators$`conc(N-NH4µg/gmresin)`)
normalized_nitrate <- normalize(Conservation_Indicators$`conc(N-NO3µg/gmresin)`)
normalized_Decomposition_rate <- normalize(Conservation_Indicators$Decomposition_rate_yr)
normalized_Tflux<- normalize(Conservation_Indicators$total_flux)
#normalized_BS<- normalize(Agronomic_Indicators$`BS%`)
normalized_Carbon <- normalize(Conservation_Indicators$SOC)
normalized_aggregates<- normalize(Conservation_Indicators$Mean_aggregates)
#normalized_NP<- normalize(Agronomic_Indicators$NP)
#normalized_CN<- normalize(Agronomic_Indicators$CN)
#normalized_Pottasium<- normalize(Agronomic_Indicators$`K_tot(mg/kg)`)

# Invert the ammonia and nitrate indicators (assuming higher values are negative contributors)
normalized_Tflux <- 1 - normalized_Tflux
normalized_nitrate <- 1 - normalized_nitrate
normalized_Ammonia<- 1 - normalized_Ammonia
normalized_aggregates<-1-normalized_aggregates


# Combine normalized indicators into a data frame
normalized_Conservation_Indicators <- data.frame(
  normalized_Decomposition_rate,
  normalized_Carbon,
  normalized_Ammonia,
  normalized_nitrate,
  normalized_Tflux
)

view(normalized_Conservation_Indicators)
# Assign weights to each indicator (assuming equal weights for this example)
#weights_nu_mineralization_indicators <- rep(1 / ncol(normalized_nu_mineralization_indicators), ncol(normalized_nu_mineralization_indicators))

# Calculate the Mineralization Index
Conservation_SMF_index <- rowMeans(normalized_Conservation_Indicators, na.rm = TRUE)

# Add the indices to the original data
SMF_merged_data$CSMF_index <- Conservation_SMF_index
plot(Conservation_SMF_index)
# View the updated data
View(SMF_merged_data)

##############################################


# analysis of variance
anova <- aov(CSMF_index ~ factor(Ecosystems), data = SMF_merged_data)
summary(anova)
# Tukey's test and compact letter display
Tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, Tukey)

# Table with the mean, the standard deviation and the letters indications significant differences for each treatment
dt <- group_by(SMF_merged_data, Ecosystems) %>%
  summarise(CSMF_index_mean=mean(CSMF_index), se=sd(CSMF_index) / sqrt(n())) %>%
  arrange(desc(CSMF_index_mean))
cld <- as.data.frame.list(cld$`factor(Ecosystems)`)
dt$Tukey <- cld$Letters
Tukey<-dt$Tukey
print(dt)

##########

CSMF<-ggplot(dt, aes(x = Ecosystems, y = CSMF_index_mean, fill = Ecosystems)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = CSMF_index_mean + se, ymin = CSMF_index_mean - se),
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label=Tukey, y = CSMF_index_mean + se+0.03), size = 3, color = "Gray25",
            show.legend = FALSE,
            position = position_dodge(0.9)) +
  ylim(0,0.8)+
  scale_fill_manual(values = custom_palette) +  # Use scale_fill_manual for custom palette
  theme_classic()+
  labs(x = "", y = "Conservation Soil Multifunctionality Index", title = "CONSERVATION SOIL MULTIFUNCTIONALITY") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = "black"),  # X-axis text settings
        axis.text.y = element_text(size = 12, color = "black"),  # Y-axis text settings
        axis.title.y = element_text(size = 12, face = "bold", color = "black"), # Y-axis title settings
        plot.title = element_text(size = 14, face = "bold", color = "black"),   # Title settings
        axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=FALSE)
CSMF
############################################################################################
# analysis of variance
anova <- aov(CSMF_index ~ factor(Management_types), data = SMF_merged_data)
summary(anova)
# Tukey's test and compact letter display
Tukey <- TukeyHSD(anova)
cld <- multcompLetters4(anova, Tukey)

# Table with the mean, the standard deviation and the letters indications significant differences for each treatment
dt <- group_by(SMF_merged_data, Management_types) %>%
  summarise(CSMF_index_mean=mean(CSMF_index), se=sd(CSMF_index) / sqrt(n())) %>%
  arrange(desc(CSMF_index_mean))
cld <- as.data.frame.list(cld$`factor(Management_types)`)
dt$Tukey <- cld$Letters
Tukey<-dt$Tukey
print(dt)

##########

CSMF2<-ggplot(dt, aes(x = Management_types, y = CSMF_index_mean, fill = Management_types)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = CSMF_index_mean + se, ymin = CSMF_index_mean - se),
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label=Tukey, y = CSMF_index_mean + se+0.03), size = 3, color = "Gray25",
            show.legend = FALSE,
            position = position_dodge(0.9)) +
  ylim(0,0.8)+
  scale_fill_manual(values = c("Agricultural_plots" = "#92C5DE",
                               "National_parks" = "#2166AC",
                               "Savanna_grassland" = "#B2182B"))+  # Use scale_fill_manual for custom palette
  theme_classic()+
  labs(x = "", y = "Conservation Soil Multifunctionality Index", title = "CONSERVATION SOIL MULTIFUNCTIONALITY") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 12, color = "black"),  # X-axis text settings
        axis.text.y = element_text(size = 12, color = "black"),  # Y-axis text settings
        axis.title.y = element_text(size = 12, face = "bold", color = "black"), # Y-axis title settings
        plot.title = element_text(size = 14, face = "bold", color = "black"),   # Title settings
        axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=FALSE)
CSMF2
############################################################################################


#correlation between SMF AND THE AXIS
# Relate to PCA axis scores using correlation

cor.test(soil_data$mf_agronomist, soil_data$Axis1)
cor.test(soil_data$mf_conservationist, soil_data$Axis2)

# Select five variables for correlation analysis (Replace with your chosen variables)
selected_vars <- SMF_merged_data[, c( "OSMF_index", "ASMF_index", "CSMF_index","PCA1","PCA2")]

# Compute the correlation matrix
correlation_matrix <- cor(selected_vars, use = "complete.obs")
# Set up the output device to save the plot

png(filename = "correlation_plot.png", width = 8, height = 6, units = "in", res = 800)

# Create a correlation plot with correlation coefficients, color, and size variations
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, # Text label color and rotation
         number.cex = 0.8, # Size of the correlation coefficient numbers
         addCoef.col = "black", # Color of the correlation coefficients
         col = colorRampPalette(c("blue", "white", "red"))(200), # Color palette
         diag = FALSE) # Hide the diagonal

# Close the output device to save the plot


smfcorrelation<-corrplot(correlation_matrix, method = "color", type = "upper",
                         tl.col = "black", tl.srt = 45, # Text label color and rotation
                         number.cex = 0.8, # Size of the correlation coefficient numbers
                         addCoef.col = "black", # Color of the correlation coefficients
                         col = colorRampPalette(c("blue", "white", "red"))(200), # Color palette
                         diag = FALSE) #
print(smfcorrelation)
dev.off()
##################################################################

##################STRUCTRUAL EQUATION MODDELLING FOR THE MULTIFUNCTIONALITY####

##usng rather root traits
set.seed(221443)
OSMF_index_SEM<-psem(
  lm(log(OSMF_index)~MAT+MAP+Sand+ Clay+ PH+LUI+PCA1_RL+PCA2_RL+PCA1+PCA2, data=SMF_merged_data),
  lm(LUI~MAP+MAT,data =SMF_merged_data),
  lm(Sand~MAT+MAP+LUI,data =SMF_merged_data),
  lm(Clay~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PH~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA1_RL~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA2_RL~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA1~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA2~MAT+MAP+LUI,data =SMF_merged_data),
  data=SMF_merged_data)

summary(OSMF_index_SEM)
fisherC(OSMF_index_SEM)
AIC(OSMF_index_SEM)
summary(OSMF_index_SEM, conserve = T, scale=TRUE)

###################update model################

OSMF_index_SEM<-update(OSMF_index_SEM,
                       Clay%~~%Sand,
                       PH %~~%Sand,
                       PCA2_RL %~~% Sand,
                       PH %~~%Clay,
                       PCA1_RL %~~% Clay,
                       PCA1_RL %~~% Sand,
                       PCA1_RL %~~%PH,
                       PCA1 %~~% Sand,
                       PCA2_RL %~~% Clay,
                       PCA1 %~~% Clay,
                       PCA1 %~~% PH,
                       PCA1 %~~% PCA1_RL,
                       PCA2 %~~% PCA1_RL,
                       PCA1 %~~% PCA2_RL,
                       PCA2 %~~% PCA2_RL)
summary(OSMF_index_SEM)
fisherC(OSMF_index_SEM)
AIC(OSMF_index_SEM)
evaluate_model(OSMF_index_SEM)

#####CHECKING THE ASSUMPTION OF NORMALITY
A<-lm(log(OSMF_index)~MAT+MAP+Sand+ Clay+ LUI+PCA1_RL+PCA2_RL+PCA1+PCA2, data=SMF_merged_data)
B<-lm(LUI~MAP+MAT,family=gaussian,data =SMF_merged_data)
C<-lm(Elev~MAT+MAP+LUI,data =SMF_merged_data)
C<-lm(Sand~MAT+MAP+LUI,data =SMF_merged_data)
D<-lm(Clay~MAT+MAP+LUI,data =SMF_merged_data)
E<-lm(PCA1_RL~MAT+MAP+LUI,data =SMF_merged_data)
F<-lm(PCA2_RL~MAT+MAP+LUI,data =SMF_merged_data)
E<-lm(PCA1~MAT+MAP+LUI,data =SMF_merged_data)
F<-lm(PCA2~MAT+MAP+LUI,data =SMF_merged_data)

check_heteroscedasticity(E)
check_normality(E)
check_model(A)
#########################################################################
OSMF_index_SEM_coef<-summary(OSMF_index_SEM)$coefficients
OSMF_index_SEM_coef
OSMF_index_SEM_coef <- as.data.frame(OSMF_index_SEM_coef)
write.csv(OSMF_index_SEM_coef, file = "OSMF_index_SEM_coef_coef_table.csv")

################################################################
##usng rather root traits
set.seed(221443)
CSMF_index_SEM<-psem(
  lm(CSMF_index~MAT+MAP+Sand+ Clay+ PH+LUI+PCA1_RL+PCA2_RL+PCA1+PCA2, data=SMF_merged_data),
  lm(LUI~MAP+MAT,data =SMF_merged_data),
  lm(Sand~MAT+MAP+LUI,data =SMF_merged_data),
  lm(Clay~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PH~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA1_RL~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA2_RL~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA1~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA2~MAT+MAP+LUI,data =SMF_merged_data),
  data=SMF_merged_data)

summary(CSMF_index_SEM)
fisherC(CSMF_index_SEM)
AIC(CSMF_index_SEM)
summary(CSMF_index_SEM, conserve = T, scale=TRUE)

###################update model################

CSMF_index_SEM<-update(CSMF_index_SEM,
                       Clay%~~%Sand,
                       PH %~~%Sand,
                       PCA2_RL %~~% Sand,
                       PH %~~%Clay,
                       PCA1_RL %~~% Clay,
                       PCA1_RL %~~% Sand,
                       PCA1_RL %~~%PH,
                       PCA1 %~~% Sand,
                       PCA2_RL %~~% Clay,
                       PCA1 %~~% Clay,
                       PCA1 %~~% PH,
                       PCA1 %~~% Sand,
                       PCA1 %~~% PCA1_RL,
                       PCA2 %~~% PCA1_RL,
                       PCA1 %~~% PCA2_RL,
                       PCA2 %~~% PCA2_RL)

summary(CSMF_index_SEM)
fisherC(CSMF_index_SEM)
AIC(CSMF_index_SEM)
evaluate_model(CSMF_index_SEM)

#####CHECKING THE ASSUMPTION OF NORMALITY
A<-lm(CSMF_index~MAT+MAP+Sand+ Clay+ LUI+PCA1_RL+PCA2_RL+PCA1+PCA2, data=SMF_merged_data)
B<-lm(LUI~MAP+MAT,family=gaussian,data =SMF_merged_data)
C<-lm(Elev~MAT+MAP+LUI,data =SMF_merged_data)
C<-lm(Sand~MAT+MAP+LUI,data =SMF_merged_data)
D<-lm(Clay~MAT+MAP+LUI,data =SMF_merged_data)
E<-lm(PCA1_RT~MAT+MAP+LUI,data =SMF_merged_data)
F<-lm(PCA2_RT~MAT+MAP+LUI,data =SMF_merged_data)

check_heteroscedasticity(A)
check_normality(A)
check_model(A)
#########################################################################
CSMF_index_SEM_coef<-summary(CSMF_index_SEM)$coefficients
CSMF_index_SEM_coef
CSMF_index_SEM_coef <- as.data.frame(CSMF_index_SEM_coef)
write.csv(CSMF_index_SEM_coef, file = "CSMF_index_SEM_coef_coef_table.csv")
#########################################################################



##usng rather root traits
set.seed(225543)
ASMF_index_SEM<-psem(
  lm(ASMF_index~MAT+MAP+Sand+ Clay+ PH+LUI+PCA1_RL+PCA2_RL+PCA1+PCA2, data=SMF_merged_data),
  lm(LUI~MAP+MAT,data =SMF_merged_data),
  lm(Sand~MAT+MAP+LUI,data =SMF_merged_data),
  lm(Clay~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PH~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA1_RL~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA2_RL~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA1~MAT+MAP+LUI,data =SMF_merged_data),
  lm(PCA2~MAT+MAP+LUI,data =SMF_merged_data),
  data=SMF_merged_data)

summary(ASMF_index_SEM)
fisherC(ASMF_index_SEM)
AIC(ASMF_index_SEM)
summary(ASMF_index_SEM, conserve = T, scale=TRUE)

###################update model################

ASMF_index_SEM<-update(ASMF_index_SEM,
                       Clay%~~%Sand,
                       PH %~~%Sand,
                       PCA2_RL %~~% Sand,
                       PH %~~%Clay,
                       PCA1_RL %~~% Clay,
                       PCA1_RL %~~% Sand,
                       PCA1_RL %~~%PH,
                       PCA1 %~~% Sand,
                       PCA2_RL %~~% Clay,
                       PCA1 %~~% Clay,
                       PCA1 %~~% PH,
                       PCA1 %~~% Sand,
                       PCA1 %~~% PCA1_RL,
                       PCA2 %~~% PCA1_RL,
                       PCA1 %~~% PCA2_RL,
                       PCA2 %~~% PCA2_RL)
summary(ASMF_index_SEM)
fisherC(ASMF_index_SEM)
AIC(ASMF_index_SEM)
evaluate_model(ASMF_index_SEM)

#####CHECKING THE ASSUMPTION OF NORMALITY
A<-lm(ASMF_index~MAT+MAP+Sand+ Clay+ LUI+PCA1_RT+PCA2_RT,data=SMF_merged_data)
B<-lm(LUI~MAP+MAT,family=gaussian,data =SMF_merged_data)
C<-lm(Elev~MAT+MAP+LUI,data =SMF_merged_data)
C<-lm(Sand~MAT+MAP+LUI,data =SMF_merged_data)
D<-lm(Clay~MAT+MAP+LUI,data =SMF_merged_data)
E<-lm(PCA1_RT~MAT+MAP+LUI,data =SMF_merged_data)
F<-lm(PCA2_RT~MAT+MAP+LUI,data =SMF_merged_data)

check_heteroscedasticity(A)
check_normality(A)
check_model(A)
#########################################################################
ASMF_index_SEM_coef<-summary(ASMF_index_SEM)$coefficients
ASMF_index_SEM_coef
ASMF_index_SEM_coef <- as.data.frame(ASMF_index_SEM_coef)
write.csv(ASMF_index_SEM_coef, file = "ASMF_index_SEM_coef_coef_table.csv")
###########################################################################


####CORRELATION WITH THE MANTEL TEST
#Sample analysis Kilidataset
# Load necessary libraries
#drivers<-read.csv("Plot characteristics.csv")
#colnames(drivers)
#columns_to_keep <- c("habitat" , "MAT","MAP","LUI", "lu1", "PS1")
#drivers2 <- drivers[, columns_to_keep]

#colnames(drivers2)

# Group by habitat type and calculate summary statistics
##drivers_sum <- drivers2 %>%
#group_by(habitat) %>%
#summarize(
#mean_MAT = mean(MAT,na.rm = TRUE),
#mean_MAP = mean(MAP,na.rm = TRUE),
#mean_LUI = mean(LUI,na.rm = TRUE),
#mean_lu1 = mean(lu1,na.rm = TRUE),
#mean_PS1 = mean(PS1,na.rm = TRUE)
#)

# Ensure NCP contains only numeric columns
#NCP<- NCP[sapply(NCP, is.numeric)]
#####MULTIFUNCTIONALITY
columns_to_keep <- c("OSMF_index" , "CSMF_index","ASMF_index")
SMF<-SMF_merged_data[,columns_to_keep]
################################################
####ENVIRONMENTAL DRIVERS
columns_to_keep2 <- c("MAP" , "MAT","LUI", "PCA1_RL", "PCA2_RL")
Drivers<-SMF_merged_data[,columns_to_keep2]

##########################MANTEL ANALYSIS
mantel <- mantel_test(Drivers, SMF,
                      spec_select = list(MAT = 1,
                                         MAP = 2,
                                         LUI = 3,
                                         PCA1_RL=4,
                                         PCA2_RL=5)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
#> mantel_test() using 'bray' dist method for 'spec'.
#> mantel_test() using 'euclidean' dist method for 'env'.


qcorrplot(correlate(SMF), type = "lower", diag = TRUE) +
  geom_square() +coord_fixed(ratio = 1)+
  geom_mark(only_mark = TRUE,
            sig_level = c(0.05, 0.01, 0.001),
            mark = c("", "", ""),
            sig_thres = 0.05, size=3)+
  geom_couple(aes(colour = pd, size = rd),
              data = mantel,
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"),
                             order = 2),
         colour = guide_legend(title = "Mantel's p",
                               override.aes = list(size = 2),
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r",order=2))

####################################################




##########VARIANCE PARTITIONING
# Install and load the necessary package if not already installed
# install.packages("rdacca.hp")
library(rdacca.hp)

# Assuming your data is ready in dv (response variables) and iv (independent variables)
# No data reduction or selection, so all variables will be used

# Perform the variance analysis without reducing data
## Select five variables for correlation analysis (Replace with your chosen variables)
#OSMF_index_var <- SMF_merged_data[, c( "OSMF_index", "ASMF_index", "CSMF_index","PCA1","PCA2")]
OSMF_index_var <- SMF_merged_data[,"OSMF_index"]
ASMF_index_var <- SMF_merged_data[,"ASMF_index"]
CSMF_index_var <- SMF_merged_data[,"CSMF_index"]
Drivers_vars <- SMF_merged_data[, c("PH","Clay","LUI","MAT","MAP","PCA1","PCA2","Slow_fast","Woody_grassy")]
# Select five variables for correlation analysis (Replace with your chosen variables)
#selected_vars <- SMF_merged_data[, c( "OSMF_index", "ASMF_index", "CSMF_index","PCA1","PCA2")]
result <- rdacca.hp(
  dv = OSMF_index_var,                       # Your full response variables matrix
  iv = Drivers_vars,                       # Your full independent variables (drivers)
  method = "RDA",                # Method: RDA (Redundancy Analysis)
  type = "adjR2",                # Using adjusted R² for variance partitioning
  scale = FALSE,                 # No scaling of the variables (optional)
  add = FALSE,                   # Not adding additional species data
  sqrt.dist = FALSE,             # No sqrt transformation (only needed for dbRDA)
  n.perm = 0,                    # Set permutations to 0 to avoid any randomization
  var.part = TRUE                # Perform variance partitioning, keeping all data
)

# View the result install.packages("VennDiagram")

print(result)
library("VennDiagram")

# Prepare data: Top 3 contributing variables (LUI, PCA1, Woody_grassy)
venn_data <- list(
  LUI = c("A", "B", "C", "D"),           # Example placeholder for the LUI variable
  PCA1 = c("B", "C", "E"),               # Placeholder for PCA1
  Woody_grassy = c("C", "D", "F")        # Placeholder for Woody_grassy
)

# Draw the Venn diagram for 3 sets
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("LUI", "PCA1", "Woody_grassy"),
  filename = NULL,                       # Set to NULL to plot in R instead of saving to a file
  fill = c("red", "green", "blue"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.fontface = "bold"
)

# Plot the Venn diagram
grid.draw(venn.plot)

###############################################################





# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Normalization function (Min-Max Normalization)
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Step 1: Calculate Overall Soil Multifunctionality Index
Overall_Indicators <- SMF_merged_data %>%
  select(Enzyme_Index, `conc(N-NH4µg/gmresin)`,
         `conc(N-NO3µg/gmresin)`,
         Decomposition_rate_yr, total_flux, `BS%`, SOC, Mean_aggregates, NP, CN,
         `K_tot(mg/kg)`, `15N/14N`, `12C/13C`)

# Normalization of indicators
normalized_Enzyme_Index <- normalize(Overall_Indicators$Enzyme_Index)
normalized_Ammonia <- normalize(Overall_Indicators$`conc(N-NH4µg/gmresin)`)
normalized_Nisotopy <- normalize(Overall_Indicators$`15N/14N`)
normalized_Cisotopy <- normalize(Overall_Indicators$`12C/13C`)
normalized_nitrate <- normalize(Overall_Indicators$`conc(N-NO3µg/gmresin)`)
normalized_Decomposition_rate <- normalize(Overall_Indicators$Decomposition_rate_yr)
normalized_Tflux <- normalize(Overall_Indicators$total_flux)
normalized_BS <- normalize(Overall_Indicators$`BS%`)
normalized_Carbon <- normalize(Overall_Indicators$SOC)
normalized_aggregates <- normalize(Overall_Indicators$Mean_aggregates)
normalized_NP <- normalize(Overall_Indicators$NP)
normalized_CN <- normalize(Overall_Indicators$CN)
normalized_Pottasium <- normalize(Overall_Indicators$`K_tot(mg/kg)`)

# Invert indicators where necessary
normalized_Tflux <- 1 - normalized_Tflux
normalized_nitrate <- 1 - normalized_nitrate
normalized_Ammonia <- 1 - normalized_Ammonia
normalized_aggregates <- 1 - normalized_aggregates
normalized_Nisotopy <- 1 - normalized_Nisotopy
normalized_Cisotopy <- 1 - normalized_Cisotopy

# Combine normalized indicators for Overall Soil Multifunctionality
normalized_Overall_Indicators_indicators <- data.frame(
  normalized_Enzyme_Index,
  normalized_Ammonia,
  normalized_Nisotopy,
  normalized_Cisotopy,
  normalized_nitrate,
  normalized_Decomposition_rate,
  normalized_Tflux,
  normalized_BS,
  normalized_Carbon,
  normalized_aggregates,
  normalized_NP,
  normalized_CN,
  normalized_Pottasium
)

# Calculate the Overall Soil Multifunctionality Index (OSMF_index)
Overall_Indicators_SMF_index <- rowMeans(normalized_Overall_Indicators_indicators, na.rm = TRUE)
SMF_merged_data$OSMF_index <- Overall_Indicators_SMF_index
# Calculate the average (mean) of the OSMF_index
#average_OSMF_index <- mean(SMF_merged_data$OSMF_index, na.rm = TRUE)
##average_OSMF_index=0.516
# Print the average index
print(average_OSMF_index)

# Step 2: Calculate Agronomic Stakeholders Multifunctionality Index
Agronomic_Indicators <- SMF_merged_data %>%
  select(Enzyme_Index, Decomposition_rate_yr, `BS%`, SOC, NP, CN, `K_tot(mg/kg)`)

# Normalize Agronomic indicators
normalized_Enzyme_Index <- normalize(Agronomic_Indicators$Enzyme_Index)
normalized_Decomposition_rate <- normalize(Agronomic_Indicators$Decomposition_rate_yr)
normalized_BS <- normalize(Agronomic_Indicators$`BS%`)
normalized_Carbon <- normalize(Agronomic_Indicators$SOC)
normalized_NP <- normalize(Agronomic_Indicators$NP)
normalized_CN <- normalize(Agronomic_Indicators$CN)
normalized_Pottasium <- normalize(Agronomic_Indicators$`K_tot(mg/kg)`)

# Combine normalized indicators for Agronomic Stakeholders
normalized_Agronomic_Indicators <- data.frame(
  normalized_Enzyme_Index,
  normalized_Decomposition_rate,
  normalized_BS,
  normalized_Carbon,
  normalized_NP,
  normalized_CN,
  normalized_Pottasium
)

# Calculate Agronomic Soil Multifunctionality Index (ASMF_index)
Agronomic_SMF_index <- rowMeans(normalized_Agronomic_Indicators, na.rm = TRUE)
SMF_merged_data$ASMF_index <- Agronomic_SMF_index
#average_ASMF_index <- mean(SMF_merged_data$ASMF_index, na.rm = TRUE)
#average_ASMF_index=0.3812


# Step 3: Calculate Conservation Stakeholders Multifunctionality Index
Conservation_Indicators <- SMF_merged_data %>%
  select(`conc(N-NH4µg/gmresin)`, `conc(N-NO3µg/gmresin)`,
         Decomposition_rate_yr, total_flux, SOC, Mean_aggregates)

# Normalize Conservation indicators
normalized_Ammonia <- normalize(Conservation_Indicators$`conc(N-NH4µg/gmresin)`)
normalized_nitrate <- normalize(Conservation_Indicators$`conc(N-NO3µg/gmresin)`)
normalized_Decomposition_rate <- normalize(Conservation_Indicators$Decomposition_rate_yr)
normalized_Tflux <- normalize(Conservation_Indicators$total_flux)
normalized_Carbon <- normalize(Conservation_Indicators$SOC)
normalized_aggregates <- normalize(Conservation_Indicators$Mean_aggregates)

# Invert indicators where higher values are negative contributors
normalized_Tflux <- 1 - normalized_Tflux
normalized_nitrate <- 1 - normalized_nitrate
normalized_Ammonia <- 1 - normalized_Ammonia
normalized_aggregates <- 1 - normalized_aggregates

# Combine normalized indicators for Conservation Stakeholders
normalized_Conservation_Indicators <- data.frame(
  normalized_Ammonia,
  normalized_nitrate,
  normalized_Decomposition_rate,
  normalized_Tflux,
  normalized_Carbon,
  normalized_aggregates
)

# Calculate Conservation Soil Multifunctionality Index (CSMF_index)
Conservation_SMF_index <- rowMeans(normalized_Conservation_Indicators, na.rm = TRUE)
SMF_merged_data$CSMF_index <- Conservation_SMF_index
#average_CSMF_index <- mean(SMF_merged_data$CSMF_index, na.rm = TRUE)

#average_CSMF_index=0.603

# View the updated data frame with all indices
View(SMF_merged_data)
#####################################


# Calculate total index values
Mean_overall_index <- mean(SMF_merged_data$OSMF_index, na.rm = TRUE)
Mean_agronomic_index <- mean(SMF_merged_data$ASMF_index, na.rm = TRUE)
Mean_conservation_index <- mean(SMF_merged_data$CSMF_index, na.rm = TRUE)

# Calculate mean contributions for each index
#mean_contributions_overall <- colSums(normalized_Overall_Indicators_indicators, na.rm = TRUE) / total_overall_index
#mean_contributions_agronomic <- colSums(normalized_Agronomic_Indicators, na.rm = TRUE) / total_agronomic_index
#mean_contributions_conservation <- colSums(normalized_Conservation_Indicators, na.rm = TRUE) / total_conservation_index

########################################
total_overall_index <- sum( normalized_Overall_Indicators_indicators, na.rm = TRUE)
total_agronomic_index <- sum(normalized_Agronomic_Indicators, na.rm = TRUE)
total_conservation_index <- sum(normalized_Conservation_Indicators, na.rm = TRUE)
###################################################

mean_contributions_overall <- colSums(normalized_Overall_Indicators_indicators, na.rm = TRUE)
mean_contributions_agronomic <- colSums(normalized_Agronomic_Indicators, na.rm = TRUE)
mean_contributions_conservation <- colSums(normalized_Conservation_Indicators, na.rm = TRUE)
###################################
mean_contributions_overall <-  (mean_contributions_overall/ total_overall_index)*Mean_overall_index
mean_contributions_agronomic <- (mean_contributions_agronomic / total_agronomic_index)*Mean_agronomic_index
mean_contributions_conservation <- (mean_contributions_conservation/total_conservation_index)* Mean_conservation_index
#############################################################################################


# Combine contributions into a single data frame
contributions_df <- data.frame(
  Variable = c(names(mean_contributions_overall), names(mean_contributions_agronomic), names(mean_contributions_conservation)),
  MeanValue = c(mean_contributions_overall, mean_contributions_agronomic, mean_contributions_conservation),
  Index = c(rep("Overall SMF", length(mean_contributions_overall)),
            rep("Agronomic SMF", length(mean_contributions_agronomic)),
            rep("Conservation SMF", length(mean_contributions_conservation)))
)
view(contributions_df)
###########################################################

# Plot the bar graph
ggplot(contributions_df, aes(x = MeanValue, y = Index, fill = Variable)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Contribution of Variables to Soil Multifunctionality Indices",
    x = "Score",
    y = "Indices"
  ) +
  theme_minimal() +
  theme(legend.position = "right")
###################################################
# Define the color palette with the specified color codes
#custom_colors <- c("#000000", "#FF8000", "#660000", "#FF6666", "#CCFFCC", "#CCE5FF",
#"#006666", "#00CCCC", "#E0E0E0", "#006633", "#990099", "#FFCCCC", "#000066")
custom_colors <- c(
  "#000000",  # Black (darkest)
  "#000066",  # Dark Blue
  "#660000",  # Dark Red
  "#006633",  # Dark Green
  "#006666",  # Teal
  "#FF8000",  # Orange
  "#990099",  # Purple
  "#00CCCC",  # Cyan
  "#CCE5FF",  # Light Blue
  "#CCFFCC",  # Light Green
  "#808080",  # Light Grey
  "#FF6666",  # Light Red
  "#FFCCCC"   # Light Pink (lightest)
)

# Ensure there are enough colors for all variables by repeating the color vector if necessary
num_vars <- length(unique(contributions_df$Variable))
custom_colors <- rep(custom_colors, length.out = num_vars)

# Plot the bar graph with custom colors
ggplot(contributions_df, aes(x = MeanValue, y = Index, fill = Variable)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Contribution of Variables to Soil Multifunctionality Indices",
    x = "Score",
    y = "SMF categories"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = custom_colors)
######################################################

ggplot(contributions_df, aes(x = MeanValue, y = Index, fill = Variable)) +
  geom_bar(stat = "identity", position = "stack", width = 0.3) +  # Adjust the width here (0.7 is an example)
  labs(
    title = "Contribution of Variables to Soil Multifunctionality Indices",
    x = "Score",
    y = "SMF categories"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_fill_manual(values = custom_colors)


##########################################################################################

