# MEDS5420-Final-project

# Project outline   
  Data sets analyzed 
      1. D07 proliferating cells vs D77 neurons 
        -Alternative splicing analyses for 3' and 5' splice sites, mutually exclusive exons, retained introns, and skipped exons
      2. WT D77 neurons vs D77 neurons shRNA _CELF5_ KD neurons
        -Alternative splicing analyses for 3' and 5' splice sites, mutually exclusive exons, retained introns, and skipped exons

      Analyzed data sets have been graphed as volcano plots and heatmaps. Used Dr. Miura's code as a starting point, and Claude and DeepSeek to troubleshoot codes.  

# Load and install relevant packages for analyses and plotting 

```{r}
#install and load relevant packages for analyses
install.packages("ggplot2")
install.packages("here")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("pheatmap")

library("ggplot2")
library("here")
library("ggrepel")
library("dplyr")
library("tidyverse")
library("pheatmap")
```
# Step 1: Make the rMATS file suitable for plotting 

```{r, warning=FALSE, message=FALSE}

#This following code is needed to visualize volcano plots and heatmaps.  
#Load the data (using Dr. Miura's code as a starting point)
#Reads a tab-delimited file containing alternative 5' splice site (A5SS) data from the rMATS tool
#Uses the here() function to specify the file location 
rmats_data <- read.delim (here ("New_folder", "file.MATS.JC.txt"), header = TRUE)

# Inspect the data
#Displays the first few rows of the data
head(rmats_data)
#Shows the structure of the data frame
str(rmats_data)

# Clean up IncLevel columns (split and average)
clean_inc_level <- function(x) {
  if (is.na(x)) {
    return(NA)
  }
  values <- strsplit(as.character(x), ",")[[1]] # Extract values
  numeric_values <- as.numeric(values) # Attempt numeric conversion
  valid_values <- numeric_values[!is.na(numeric_values)] # Remove NAs from conversion
  if (length(valid_values) > 0) {
 return(mean(valid_values))
  } else {
    return(NA)
  }
}

rmats_data$IncLevel1_mean <- sapply(rmats_data$IncLevel1, clean_inc_level) #Creates new columns with the mean inclusion levels for each condition
rmats_data$IncLevel2_mean <- sapply(rmats_data$IncLevel2, clean_inc_level)

# Calculate delta PSI
#Computes the difference in Percent Spliced In (PSI) between the two conditions
rmats_data$deltaPSI <- rmats_data$IncLevel2_mean - rmats_data$IncLevel1_mean

# Convert FDR to -log10(FDR) for volcano plots
#Converts False Discovery Rate (FDR) values to -log10 scale
rmats_data$negLog10FDR <- -log10(rmats_data$FDR)
```   
# Step 2: Make volcano plot
```{r}
library(ggplot2)

volcano_plot <- ggplot(rmats_data, aes(x = deltaPSI, y = negLog10FDR)) +
  geom_point(aes(color = ifelse(FDR < 0.05, "Significant", "Not Significant")), size = 2) + #Colors each point based on whether its FDR value is below 0.05 (significant) or not and sets point size to 2
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) + #colors points based on statistical significance
  labs(title = "Volcano Plot of Differential Splicing",
       x = "Delta PSI (IncLevel2 - IncLevel1)",
       y = "-log10(FDR)") +
  theme_minimal() +
  theme(legend.position="none") + #hides legend
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") #blue line to indicate where the threshold for significance is 

print(volcano_plot)

```

# Step 3: Make volcano plot and label the up and down regulated genes different colors and add corresponding gene labels to the data points 
```{r}
library(ggplot2)
library(ggrepel)  # For non-overlapping text labels

# Define significance threshold as a variable
sig_threshold <- 0.05

# Create color categories based on significance and direction
rmats_data$regulation <- "Not Significant"
rmats_data$regulation[rmats_data$FDR < sig_threshold & rmats_data$deltaPSI > 0] <- "Up-regulated"
rmats_data$regulation[rmats_data$FDR < sig_threshold & rmats_data$deltaPSI < 0] <- "Down-regulated"

# Create the volcano plot with gene labels and direction-based coloring
volcano_plot <- ggplot(rmats_data, aes(x = deltaPSI, y = negLog10FDR)) +
  # Basic points
  geom_point(aes(color = regulation), 
             size = 2, alpha = 0.7) +
  # Add gene name labels with ggrepel to avoid overlapping
  geom_text_repel(
    data = subset(rmats_data, FDR < sig_threshold),  # Only label significant genes
    aes(label = geneSymbol, color = regulation),
    size = 3,
    max.overlaps = 15,  # Limit overlaps for readability
    box.padding = 0.3,
    segment.color = "grey50"
  ) +
scale_color_manual(values = c(
    "Up-regulated" = "#E41A1C",    # Red
    "Down-regulated" = "#377EB8",  # Blue
    "Not Significant" = "gray"
  )) +
  labs(title = "Volcano Plot of Differential Splicing",
       x = "Delta PSI (IncLevel2 - IncLevel1)",
       y = "-log10(FDR)",
       color = "Regulation",
       caption = paste("Significance threshold: FDR <", sig_threshold)) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) + 
  geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "blue")

# Print the plot
print(volcano_plot)
```
# Step 4: Make heatmap
```{r}
install.packages("pheatmap")
library(pheatmap)

# Select significant events based on FDR and include events less than 0.05
significant_events <- subset(rmats_data, FDR < 0.05)

# Create a matrix of PSI values
psi_matrix <- data.frame(
  Sample1 = significant_events$IncLevel1_mean, #mean inclusion levels for sample 1
  Sample2 = significant_events$IncLevel2_mean #mean inclusion levels for sample 2
)

rownames(psi_matrix) <- significant_events$ID #to label splicing events in heatmap

# Generate the heatmap
pheatmap(t(psi_matrix), scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, main = "Heatmap of PSI Values")
```
# Step 5: Make heatmap and annotate the genes
```{r}
library(pheatmap)

# Assuming 'rmats_data' is your data frame as before

# Select significant events based on FDR
significant_events <- subset(rmats_data, FDR < 0.05)

# Create a matrix of PSI values
psi_matrix <- data.frame(
  Sample1 = significant_events$IncLevel1_mean,
  Sample2 = significant_events$IncLevel2_mean
)

rownames(psi_matrix) <- significant_events$ID

# Create annotation data frame for gene symbols
annotation_row <- data.frame(Gene = significant_events$geneSymbol)
rownames(annotation_row) <- significant_events$ID
# Generate the heatmap with annotations
pheatmap(t(psi_matrix),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of PSI Values",
         annotation_col = annotation_row) #add gene symbol annotation
```
# Step 6: Make heatmap of only the top 10 significant genes 
```{r}
library(pheatmap)

# Select top 10 most significant events (by FDR)
top_10_events <- head(rmats_data[order(rmats_data$FDR), ], 10)

# Create a matrix of PSI values for top 10
psi_matrix <- data.frame(
  Sample1 = top_10_events$IncLevel1_mean,
  Sample2 = top_10_events$IncLevel2_mean,
  row.names = top_10_events$ID
)

# Create annotation for gene symbols
annotation_row <- data.frame(
  Gene = top_10_events$geneSymbol,
  row.names = top_10_events$ID
)
# Generate the heatmap with adjusted parameters
pheatmap(t(psi_matrix),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Top 10 Significant Splicing Events",
         annotation_col = annotation_row,
         fontsize_col = 10,  # Larger font for gene labels
         cellwidth = 30,     # Wider cells for better label visibility
         cellheight = 30,     # Taller cells
         angle_col = 45)      # Diagonal labels for readability
```
