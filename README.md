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
# Conclusions for D07 proliferating cells vs D77 neurons
  1. 3' alternate splicing analysis- _SYNJ1_ appears to be one of the top differentially spliced genes between the two samples. Codes for synaptojanin 1 protein which is involved in endocytosis of synaptic vescicles. It is not clear which sample is D77 neurons, but given that this gene encodes a protein important for synpatic vesicular trafficking, I would predict that its expression is higher in the D77 sample compared to D07.
  2. 5' alternate splicing analysis- _AUTS2_, a gene that is implicated in autism spectrum disorders, is differentially spliced between the 2 samples. Since it is important for neuronal migration, I would assume that the D77 sample has higher expression.
  3. Retained intron analysis- _DDX1_ is an RNA helicase that is involved in neuronal differentiation
  4. Mutually exclusive exon analysis- _PPFIA3_ encodes a scaffolding protein that is involved in maintaining the structure of synapses
  5. Skipped econ analysis- top 10 differentially spliced genes do not appear to be different between the 2 samples
Overall, genes that are differentially spliced between the data sets appear to be ones that are involved in neuronal differentiation.

# Conclusions for WT vs _CELF_ shRNA knockdown D77 neurons 
  1. 3' alternate splicing analysis- one of the top differentially spliced genes appears to be a long noncoding RNA
  2. 5' alternate spllicing analysis- _ANK2_ which codes for ankyrin B appears to be differentially spliced between the 2 samples, suggesting it could be a target of _CELF5_
  3. Retained intron analysis- _APRG1_ is another long noncoding RNA that is differentially spliced between the two samples
  4. Mutually exclusive exon analysis- Top genes do not appear to be different between the two samples
  5. Skipped exon analysis- Top genes do not appear to be different between these samples




# Failed codes

```{r}
library(ggplot2)
library(ggrepel) # For smart label placement
library(tidyverse)
library(dplyr)

# Create a copy of the data for labeling
plot_data <- rmats_data %>%
  mutate(
    Significance = ifelse(FDR < 0.05, "Significant", "Not Significant"),
    # Create a column for labels (top 10 most significant genes)
    Label = ifelse(
      FDR < 0.05 & 
      rank(-negLog10FDR) <= 10,  # Label top 10 most significant
      as.character(geneSymbol),   # Use geneSymbol if available, else ID
      ""
    )
  )

# Create the volcano plot with labels
volcano_plot <- ggplot(plot_data, aes(x = deltaPSI, y = negLog10FDR)) +
  geom_point(aes(color = Significance), size = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  
  # Add gene labels (only for significant points)
  geom_text_repel(
    aes(label = Label),
    size = 0.2,
    box.padding = 0.5,   # Adjust spacing
    max.overlaps = 50,   # Increase if needed
    min.segment.length = 0.2,  # Show more lines
    segment.color = "grey50"
  ) +
  
  # Reference lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  # Labels and theme
  labs(
    title = "Volcano Plot of Differential Splicing",
    x = "Delta PSI (IncLevel2 - IncLevel1)",
    y = "-log10(FDR)",
    color = "Significance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

print(volcano_plot)

```

```{r}
library(ggplot2)
library(ggrepel)
library(dplyr)

# Prepare the data with all genes labeled
plot_data <- rmats_data %>%
  mutate(
    Significance = ifelse(FDR < 0.05, "Significant", "Not Significant"),
    Label = ifelse(is.na(geneSymbol) | geneSymbol == "", as.character(ID), geneSymbol)
  )

# Create the volcano plot with all genes labeled
volcano_plot <- ggplot(plot_data, aes(x = deltaPSI, y = negLog10FDR)) +
  geom_point(aes(color = Significance), 
             size = 0.5,  # Slightly larger points
             alpha = 0.6) +
  
  # Label ALL points
  geom_text_repel(
    aes(label = Label),
    size = 1.5,                # Smaller font size
    max.overlaps = 5,        # Show all labels
    min.segment.length = 0,    # Always show connecting lines
    segment.size = 0.1,        # Thinner connecting lines
    segment.alpha = 0.5,       # Semi-transparent lines
    box.padding = 0.1,         # Tighten label padding
    force = 0.5,               # Reduce repelling force
    max.time = 2,              # Allow more computation time
    max.iter = 1e6,            # Increase maximum iterations
    seed = 123                 # For reproducible positioning
  ) +
  
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray70")) +
  
  # Reference lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  
  # Labels and theme
  labs(
    title = "Volcano Plot of Differential Splicing (All Genes Labeled)",
    x = "Delta PSI (IncLevel2 - IncLevel1)",
    y = "-log10(FDR)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 10),
    text = element_text(size = 8)
  )

# For better performance with many points:
options(ggrepel.max.overlaps = Inf)

print(volcano_plot)
```

```{r}
library(pheatmap)

# Select significant events
significant_events <- subset(rmats_data, FDR < 0.05)

# Create PSI matrix
psi_matrix <- data.frame(
  Sample1 = significant_events$IncLevel1_mean,
  Sample2 = significant_events$IncLevel2_mean,
  row.names = significant_events$ID
)

# Calculate optimal plot size based on number of events
num_events <- nrow(psi_matrix)
plot_height <- 8  # Fixed height in inches
plot_width <- max(6, num_events * 0.2)  # Dynamic width (0.2" per event)

# Generate heatmap with adjusted parameters
heatmap <- pheatmap(t(psi_matrix),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of PSI Values",
         fontsize_col = 8,           # Smaller font for columns
         angle_col = 45,             # Diagonal labels
         cellwidth = 5,             # Fixed cell width (adjust as needed)
         treeheight_row = 60,        # Shorter dendrogram for rows
         treeheight_col = 60,        # Shorter dendrogram for columns
         width = plot_width,         # Dynamic width
         height = plot_height)       # Fixed height
```

```{r}
library(pheatmap)
library(grid)

# Define significance threshold
sig_threshold <- 0.05

# Select significant events based on FDR
significant_events <- subset(rmats_data, FDR < sig_threshold)

# If there are too many events, consider taking top N by significance
if(nrow(significant_events) > 50) {
  significant_events <- significant_events[order(significant_events$FDR), ][1:50, ]
  message("Limiting to top 50 most significant events for clarity")
}

# Create a matrix of PSI values
psi_matrix <- data.frame(
  Sample1 = significant_events$IncLevel1_mean,
  Sample2 = significant_events$IncLevel2_mean
)

# Create more informative row names by combining ID with gene symbol
informative_ids <- paste(significant_events$geneSymbol, 
                         substr(significant_events$ID, 1, 10), 
                         sep = "_")
rownames(psi_matrix) <- informative_ids

# Create annotation data frame for gene symbols
annotation_row <- data.frame(Gene = significant_events$geneSymbol)
rownames(annotation_row) <- informative_ids

# Set custom color palette
my_colors <- colorRampPalette(c("#0571b0", "#f7f7f7", "#ca0020"))(100)

# Custom font sizes and rotation for better visibility
custom_fontsize <- 10  # Adjust as needed

# Generate the heatmap with improved label display
pheatmap(t(psi_matrix),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Heatmap of PSI Values",
         annotation_col = annotation_row,
         color = my_colors,
         fontsize = custom_fontsize,
         fontsize_col = custom_fontsize,
         angle_col = 45,  # Rotate column labels for better visibility
         cellwidth = 12,  # Increase cell width to accommodate labels
         cellheight = 25, # Adjust cell height
         treeheight_col = 50,  # More space for dendrogram
         border_color = NA,    # Remove cell borders
         labels_row = c("Sample1", "Sample2"),  # Clear sample labels
         fontsize_row = 12)    # Larger font for row labels

# Alternative approach using display adjustments (uncomment if needed)
# heatmap_result <- pheatmap(t(psi_matrix),
#          scale = "row",
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          main = "Heatmap of PSI Values",
#          annotation_col = annotation_row,
#          fontsize_col = 8,
#          display_numbers = FALSE)
# 
# # Add more space to bottom margin if labels are still cut off
# grid::grid.draw(grid::textGrob("", 
#                             x = 0.5, 
#                             y = 0.1, 
#                             gp = grid::gpar(fontsize = 1)))
```

```{r}

library(ggplot2)
library(ggrepel)  # For smart label placement

# Create the plot with all genes labeled
volcano_plot <- ggplot(rmats_data, aes(x = deltaPSI, y = negLog10FDR)) +
  geom_point(aes(color = ifelse(FDR < 0.05, "Significant", "Not Significant")), 
             size = .2, alpha = 0.7) +
  
  # Add all gene labels (using ggrepel to avoid overlaps)
  geom_text_repel(
    aes(label = ifelse(FDR < 0.05, geneSymbol, "")),  # Use geneSymbol column for labels
    size = 3,                # Adjust font size
    max.overlaps = 25,      # Show all labels regardless of overlaps
    min.segment.length = 0,  # Always show connecting lines
    segment.size = 0.2,      # Thinner connecting lines
    box.padding = 0.5,       # Adjust spacing around labels
    force = 0.5,             # Adjust repelling force
    max.time = 1,            # Limit computation time
    max.iter = 1e5           # Increase iterations for better placement
  ) +
  
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  labs(title = "Volcano Plot of Differential Splicing",
       x = "Delta PSI (IncLevel2 - IncLevel1)",
       y = "-log10(FDR)") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue")

print(volcano_plot)
```

```{r}
library(ggplot2)
library(ggrepel)

# First calculate log2 fold change if not already in your data
# Assuming you have IncLevel1_mean and IncLevel2_mean columns
rmats_data$log2FC <- log2(rmats_data$IncLevel2_mean / rmats_data$IncLevel1_mean)

# Create the volcano plot with proper axes
volcano_plot <- ggplot(rmats_data, aes(x = log2FC, y = negLog10FDR)) +
  geom_point(aes(color = ifelse(FDR < 0.05, "Significant", "Not Significant")), 
             size = .1, alpha = 0.6) +
  
  # Label significant points
  geom_text_repel(
    aes(label = ifelse(FDR < 0.05, geneSymbol, "")),
    size = 3,
    max.overlaps = 25,
    min.segment.length = 0.2,
    segment.size = 0.2,
    box.padding = 0.5,
    force = 0.5
  ) +
  
  # Proper axis labels
  labs(title = "Volcano Plot of Differential Splicing",
       x = "log2(Fold Change)",
       y = "-log10(FDR)",
       color = "Significance") +
  
  # Color scheme
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray60")) +
  
  # Reference lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  
  # Theme adjustments
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank())

print(volcano_plot)
```

```{r}
library(ggplot2)
library(ggrepel)

# Calculate log2 fold change
rmats_data$log2FC <- log2(rmats_data$IncLevel2_mean / rmats_data$IncLevel1_mean)

# Calculate axis limits with 5% padding
x_limits <- c(min(rmats_data$log2FC, na.rm = TRUE) * 1.05, 
              max(rmats_data$log2FC, na.rm = TRUE) * 1.05)
y_limits <- c(0, max(rmats_data$negLog10FDR, na.rm = TRUE) * 1.05)

# Create the volcano plot with expanded axes
volcano_plot <- ggplot(rmats_data, aes(x = log2FC, y = negLog10FDR)) +
  geom_point(aes(color = ifelse(FDR < 0.05, "Significant", "Not Significant")), 
             size = 1,  # Slightly larger than 0.1 for better visibility
             alpha = 0.6) +
  
  geom_text_repel(
    aes(label = ifelse(FDR < 0.05 & abs(log2FC) > 1, geneSymbol, "")),  # Label only strong effects
    size = 3,
    max.overlaps = 25,
    min.segment.length = 0.2,
    segment.size = 0.2,
    box.padding = 0.5,
    force = 0.5
  ) +
  
  labs(title = "Volcano Plot of Differential Splicing",
       x = "log2(Fold Change)",
       y = "-log10(FDR)",
       color = "Significance") +
  
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray60")) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()) +
  
  # Key addition: Expand axes while keeping all data points
  coord_cartesian(xlim = x_limits, ylim = y_limits)

print(volcano_plot)
```

```{r}
library(ggplot2)
library(ggrepel)

# 1. First check for and handle NA values
sum(is.na(rmats_data$IncLevel1_mean))  # Check NAs in numerator
sum(is.na(rmats_data$IncLevel2_mean))  # Check NAs in denominator
sum(is.na(rmats_data$FDR))            # Check NAs in FDR

# 2. Calculate log2FC with NA handling
rmats_data$log2FC <- ifelse(
  rmats_data$IncLevel1_mean == 0 | rmats_data$IncLevel2_mean == 0,
  NA,  # Return NA if either value is 0 to avoid infinite values
  log2(rmats_data$IncLevel2_mean / rmats_data$IncLevel1_mean)
)

# 3. Handle NA values in FDR/negLog10FDR
rmats_data$negLog10FDR <- -log10(rmats_data$FDR)
rmats_data <- rmats_data[!is.na(rmats_data$log2FC) & !is.na(rmats_data$negLog10FDR), ]

# 4. Create plot with NA-safe conditions
volcano_plot <- ggplot(rmats_data, aes(x = log2FC, y = negLog10FDR)) +
  geom_point(
    aes(color = ifelse(!is.na(FDR) & FDR < 0.05, "Significant", "Not Significant")), 
    size = .1,
    alpha = 0.6
  ) +
  
  geom_text_repel(
    aes(label = ifelse(!is.na(FDR) & FDR < 0.05 & !is.na(geneSymbol), geneSymbol, "")),
    size = 1,
    max.overlaps = 10,
    min.segment.length = 0.2
  ) +
  
  labs(title = "Volcano Plot of Differential Splicing",
       x = "log2(Fold Change)", 
       y = "-log10(FDR)") +
  
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray60")) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  
  theme_minimal()

print(volcano_plot)
```

```{r}
library(ggplot2)
library(ggrepel)

# 1. First check for and handle NA values
sum(is.na(rmats_data$IncLevel1_mean))  # Check NAs in numerator
sum(is.na(rmats_data$IncLevel2_mean))  # Check NAs in denominator
sum(is.na(rmats_data$FDR))            # Check NAs in FDR

# 2. Calculate log2FC with NA handling
rmats_data$log2FC <- ifelse(
  rmats_data$IncLevel1_mean == 0 | rmats_data$IncLevel2_mean == 0,
  NA,  # Return NA if either value is 0 to avoid infinite values
  log2(rmats_data$IncLevel2_mean / rmats_data$IncLevel1_mean)
)

# 3. Handle NA values in FDR/negLog10FDR
rmats_data$negLog10FDR <- -log10(rmats_data$FDR)
rmats_data <- rmats_data[!is.na(rmats_data$log2FC) & !is.na(rmats_data$negLog10FDR), ]

# 4. Create plot with NA-safe conditions and set axis limits
volcano_plot <- ggplot(rmats_data, aes(x = log2FC, y = negLog10FDR)) +
  geom_point(
    aes(color = ifelse(!is.na(FDR) & FDR < 0.05, "Significant", "Not Significant")), 
    size = .1,
    alpha = 0.6
  ) +
  
  geom_text_repel(
    aes(label = ifelse(!is.na(FDR) & FDR < 0.05 & !is.na(geneSymbol), geneSymbol, "")),
    size = 1,
    max.overlaps = 10,
    min.segment.length = 0.2
  ) +
  
  labs(title = "Volcano Plot of Differential Splicing",
       x = "log2(Fold Change)", 
       y = "-log10(FDR)") +
  
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray60")) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  
  # Set axis limits here (before theme_minimal)
  coord_cartesian(
    xlim = c(-4, 4),       # log2FC range from -4 to 4
    ylim = c(0, 20)        # -log10(FDR) from 0 to 10
  ) +
  
  theme_minimal()

print(volcano_plot)
```

```{r}
library(ggplot2)
library(ggrepel)  # For non-overlapping text labels

# Define significance threshold as a variable
sig_threshold <- 0.05

# Create the volcano plot with gene labels
volcano_plot <- ggplot(rmats_data, aes(x = deltaPSI, y = negLog10FDR)) +
  # Basic points
  geom_point(aes(color = ifelse(FDR < sig_threshold, "Significant", "Not Significant")), 
             size = 2, alpha = 0.7) +
  # Add gene name labels with ggrepel to avoid overlapping
  geom_text_repel(
    data = subset(rmats_data, FDR < sig_threshold),  # Only label significant genes
    aes(label = gene_name),
    size = 3,
    max.overlaps = 15,  # Limit overlaps for readability
    box.padding = 0.3,
    segment.color = "grey50"
  ) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  labs(title = "Volcano Plot of Differential Splicing",
       x = "Delta PSI (IncLevel2 - IncLevel1)",
       y = "-log10(FDR)",
       caption = paste("Significance threshold: FDR <", sig_threshold)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) + 
  geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "blue")

# Print the plot
print(volcano_plot)
```

```{r}
library(ggplot2)
library(ggrepel)
library(dplyr)

# Step 1: Prepare the data for the volcano plot
# --------------------------------------
# Define significance threshold
p_threshold <- 0.05

# Calculate log2 fold change and -log10(p-value)
# Assuming you have PSI values for two conditions
volcano_data <- rmats_data %>%
  mutate(
    # Calculate log2 fold change from PSI values
    # Adding small value to avoid log(0)
    log2FC = log2((IncLevel2_mean + 0.01) / (IncLevel1_mean + 0.01)),
    
    # Calculate -log10(p-value)
    negLog10P = -log10(PValue),
    
    # Label significant points
    significance = case_when(
      PValue < p_threshold & log2FC > 1 ~ "Up-regulated",
      PValue < p_threshold & log2FC < -1 ~ "Down-regulated",
      TRUE ~ "Not Significant"
    )
  )

# Step 2: Generate the volcano plot
# --------------------------------------
# Create a more attractive color palette
my_colors <- c(
  "Up-regulated" = "#D53E4F",  # Bright red
  "Down-regulated" = "#3288BD", # Bright blue
  "Not Significant" = "#BBBBBB"  # Light gray
)

# Create the volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = negLog10P)) +
  # Add a horizontal line for p-value threshold
  geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "#666666") +
  
  # Add vertical lines for fold change thresholds
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#666666") +
  
  # Add points
  geom_point(aes(color = significance), 
             size = 2.5, 
             alpha = 0.7) +
  
  # Add gene labels for significant genes
  geom_text_repel(
    data = subset(volcano_data, PValue < p_threshold & abs(log2FC) > 1),
    aes(label = geneSymbol, color = significance),
    size = 3,
    max.overlaps = 25,
    box.padding = 0.35,
    segment.color = "grey50",
    segment.alpha = 0.6
  ) +
  
  # Set colors
  scale_color_manual(values = my_colors) +
  
  # Set plot labels
  labs(
    title = "Volcano Plot: Differential Splicing Analysis",
    subtitle = paste("p-value threshold:", p_threshold, ", log2FC threshold: ±1"),
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"(p-value)"),
    color = "Regulation"
  ) +
  
  # Set theme
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(fill = NA, color = "gray80")
  )

# Print the plot
print(volcano_plot)

# Step 3: Summary statistics (optional)
# --------------------------------------
# Count genes in each category for informative labeling
up_regulated <- sum(volcano_data$significance == "Up-regulated", na.rm = TRUE)
down_regulated <- sum(volcano_data$significance == "Down-regulated", na.rm = TRUE)
not_significant <- sum(volcano_data$significance == "Not Significant", na.rm = TRUE)

cat(paste0("Summary Statistics:\n",
          "Up-regulated events: ", up_regulated, "\n",
          "Down-regulated events: ", down_regulated, "\n",
          "Not significant events: ", not_significant, "\n",
          "Total events analyzed: ", nrow(volcano_data), "\n"))

# Step 4: Enhanced version with faceting by event type (optional)
# --------------------------------------
# If your rmats_data has an "event_type" column (e.g., SE, A3SS, A5SS, RI, MXE)
# Uncomment code below to create faceted volcano plots by event type

# Check if event_type column exists
if("event_type" %in% colnames(rmats_data)) {
  # Create faceted volcano plot
  faceted_volcano <- ggplot(volcano_data, aes(x = log2FC, y = negLog10P)) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "#666666") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#666666") +
    geom_point(aes(color = significance), size = 1.8, alpha = 0.7) +
    scale_color_manual(values = my_colors) +
    labs(
      title = "Volcano Plot by Splicing Event Type",
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(p-value)"),
      color = "Regulation"
    ) +
    facet_wrap(~event_type, scales = "free") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  # Print the faceted plot
  print(faceted_volcano)
}
```
