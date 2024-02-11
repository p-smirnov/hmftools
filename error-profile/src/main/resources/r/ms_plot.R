library(dplyr)
library(tidyr)
library(ggplot2)

# Parse and check inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 2)
{
  print("Requires arguments 1=SampleId, 2=Directory")
  stop()
}

sampleId <- args[1]
outDir <- args[2]

cupDataFile <- paste0(outDir, sampleId, '.cup.data.csv')

if (!file.exists(cupDataFile))
{
  print(sprintf('Missing CUP sample data file: %s', cupDataFile))
  stop()
}

df <- read.table("~/hartwig/repeat_profiles/ACTN01020270R.ms_table.tsv.gz", header = TRUE, sep = "\t",
                    check.names = FALSE)

# fill in missing rows
df <- df %>%
  merge(expand.grid(unit = unique(df$unit), numUnits = c(4:20)), all=TRUE) %>%
  replace(is.na(.), 0)

# https://stackoverflow.com/questions/58515587/heat-map-with-additional-values-in-r
# change to long form
df <- df %>% 
  pivot_longer(
    cols = starts_with('count'), 
    names_to = "jitter",
    values_to = "count")

df$jitter <- as.character(gsub("count", "j", df$jitter))
df$rate <- pmax(df$count / df$readCount * 100, 0.001)

df$rate <- replace(df$rate, is.na(df$rate), 0)

# make sure we have a row for read count for every unit, numUnits pair
df <- df %>%
  group_by(unit, numUnits) %>%
  summarise(jitter="readCount", rate=first(readCount)) %>%
  rbind(df)

# need to give it the sorting of the jitter values
#df$jitter <- factor(df$jitter, levels=c(sprintf("j%+d", c(-10:10)), "readCount"))

units <- c("A/T", "C/G", "AT/TA", "AG/GA/CT/TC", "AC/CA/GT/TG", "CG/GC", "3bp repeat", "4bp repeat", "5bp repeat")
#units <- head(units, 1)

# Create a list to store individual plots
plot_list <- list()

for (msUnit in units) {
  print(msUnit)
  jitter_df <- df %>%
    filter(unit == msUnit, numUnits %in% 4:20) %>%
    select(-c(unit))
  
  jitter_df$jitter <- factor(jitter_df$jitter, levels=c(sprintf("j%+d", c(-10:10)), "readCount"))
  
  if(FALSE) {
    # to fix: https://stackoverflow.com/questions/23478497/ggplot2-y-axis-order-changes-after-subsetting
  # Create a plot for each unit
  p <- ggplot(jitter_df, aes(x = jitter, y = reorder(numUnits, desc(numUnits)))) +
    geom_tile(aes(fill = rate), color = "white", subset(jitter_df, jitter != "readCount" & readCount > 0)) +
    scale_fill_gradient(low = "white", high = "deepskyblue", trans = "log", limits=c(0.001,100), guide="none") +
    geom_text(aes(label = sprintf("%.2f", rate)), size=3, subset(jitter_df, jitter != "readCount" & readCount > 0)) +
    geom_text(aes(label = format(rate, big.mark=",", scientific=FALSE)), size=3, subset(jitter_df, jitter == "readCount")) +
    theme_minimal() +
    labs(title = paste("unit =", msUnit), x="jitter", y="num units")
  }
  
  # Create a plot for each unit
  p <- ggplot(jitter_df, aes(x = jitter, y = reorder(numUnits, desc(numUnits)))) +
    geom_tile(aes(fill = rate), color = "white") +
    scale_fill_gradient(low = "white", high = "deepskyblue", trans = "log", limits=c(0.001,100), guide="none") +
    geom_text(aes(label = sprintf("%.2f", rate)), size=3, subset(jitter_df, jitter != "readCount" & readCount > 0)) +
    geom_text(aes(label = format(rate, big.mark=",", scientific=FALSE)), size=3, subset(jitter_df, jitter == "readCount")) +
    theme_minimal() +
    labs(title = paste("unit =", msUnit), x="jitter", y="num units")

  
  # Add the plot to the list
  plot_list[[msUnit]] <- p
}

# Arrange the plots in a 3x3 grid
gridExtra::grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)
