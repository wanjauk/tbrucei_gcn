library(tidyverse)
library(janitor)
library(ggplot2)
library(ggforce)
library(here)

# read in the data to create the pie chart.
gene_counts <- read.csv(here("data","intermediate","module_gene_counts.csv"), sep = ",")

# preprocess the gene counts to add sector percentage for each module colour
gene_counts$module_number <- as.character(gene_counts$module_number)
gene_counts <- gene_counts %>% select(module_colour, everything())
gene_counts <- gene_counts %>% mutate(perc = round(gene_count/colSums(gene_counts[,"gene_count", drop = FALSE])*100, 1))


# colours to fill each sector for the corresponding module colours
fill_colours <- c("Black" = "#000000",
                  "Blue" = "#0000ff",
                  "Brown" = "#964b00",
                  "Cyan" = "#00FFFF",
                  "Darkgreen" = "#013220",
                  "Darkgrey" = "#a9a9a9",
                  "Darkorange" = "#ff8c00",
                  "Darkred" = "#8b0000",
                  "Darkturquoise"= "#00ced1",
                  "Green" = "#00FF00",
                  "Greenyellow" = "#adff2f",
                  "Grey" = "#808080",
                  "Grey60" = "#999999",
                  "Lightcyan" = "#e0ffff",
                  "Lightgreen" = "#90ee90",
                  "Lightyellow" = "#ffffe0",
                  "Magenta" = "#ff00ff",
                  "Midnightblue" = "#191970",
                  "Orange" = "#ffa500",
                  "Pink" = "#ffc0cb",
                  "Purple" = "#6a0dad",
                  "Red" = "#ff0000",
                  "Royalblue" = "#4169e1",
                  "Salmon" = "#ff8c69",
                  "Tan" = "#d2b48c",
                  "Turquoise" = "#30d5c8",
                  "White" = "#ffffff",
                  "Yellow" = "#FFFF00")

gene_counts <- gene_counts %>%
  mutate(end = 2 * pi * cumsum(gene_count)/sum(gene_count),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))

# create a label with module colour and gene count
gene_counts$label <- paste(gene_counts$module_colour, paste(gene_counts$gene_count), sep=": ")

# plot the pie chart
png(file = here("results","figures","gene_counts_pie_chart.png"), type = "cairo", res =300, units = 'in',
    width = 10, height = 6, pointsize = 7)
ggplot(gene_counts) +
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = module_colour), color = "transparent") +
  geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = label,
                hjust = hjust, vjust = vjust)) +
  # coord_fixed() +
  scale_x_continuous(limits = c(-1.6, 1.4),  # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1.2, 1.2),      # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  labs(fill="Module Colour",x=NULL,y=NULL,title="",caption="") +
  scale_fill_manual(values= fill_colours) +
  theme_void()
dev.off()
