# Run with Rscript

library(dplyr)
library(ggplot2)
library(ggpubr)
theme_set(theme_pubclean())

table<-read.table("predictions_ratio.txt", header=TRUE)
pdf("genomad_vir2_pred.pdf", width = 5, height = 5, pointsize = 8)
ggboxplot(table, x = "pred_type", y = "total_ratio", color = "#02618a", add = "jitter", ylab="Total predictions ratio (geNomad / VIRify2)", xlab='') + stat_compare_means(aes(group = pred_type))
dev.off()

table<-read.table("intersection_ratio.txt", header=TRUE, sep='\t')
pdf("genomad_vir2_intersec.pdf", width = 5, height = 5, pointsize = 10)
ggboxplot(table, x = "Tool", y = "intersec_ratio", color = "Tool", palette = c("#528c1f", "#8c1f52"), add = "jitter", ylab="Intersection ratio (intersection / total predictions)", xlab='', facet.by = "Type", short.panel.labs = FALSE) + stat_compare_means(label = "p.format", label.x = 1.5)
dev.off()

summary_stats <- table %>%
  group_by(Tool, Type) %>%
  summarize(
    average = mean(intersec_ratio),
    std_deviation = sd(intersec_ratio),
    min_value = min(intersec_ratio),
    max_value = max(intersec_ratio),
    median_value = median(intersec_ratio),
    iqr_value = IQR(intersec_ratio),
    .groups = 'drop'  # Ungroup the result
  )
write.table(summary_stats, file = 'intersec_ratio_stats.txt', sep = "\t", quote = FALSE, row.names = FALSE)

table<-read.table("concat_gev2_venn_counts.txt", header=TRUE)

viral_data = subset(table, prediction_type == "viral_genome")
pdf("viral_len.pdf", width = 5, height = 5, pointsize = 10)
ggboxplot(viral_data, x = "label", y = "length", color = "label", palette = c("#528c1f", "#8c1f52"), ylab="Viral genomes prediction len (bp)", xlab='', facet.by = "tool", short.panel.labs = FALSE) + stat_compare_means(label = "p.format", label.x = 1.5) + scale_y_log10()
dev.off()


prophages_data = subset(table, prediction_type == "prophage")
pdf("prophages_len.pdf", width = 5, height = 5, pointsize = 10)
ggboxplot(prophages_data, x = "label", y = "length", color = "label", palette = c("#528c1f", "#8c1f52"), ylab="Prophage prediction len (bp)", xlab='', facet.by = "tool", short.panel.labs = FALSE) + stat_compare_means(label = "p.format", label.x = 1.5) + scale_y_log10()
dev.off()


### To generate stacked barplots
table<-read.table("corr_concat_gev2_venn_counts_taxo.txt", header=TRUE)

# For virify:
virify_data = subset(table, tool == "Virify")

filtered_data <- virify_data %>%
  filter(taxonomy != "unclassified")

total_counts <- filtered_data %>%
  group_by(prediction_type, tool, label) %>%
  summarize(total_count = sum(counts), .groups = 'drop') %>%
  ungroup()

virify_data_with_total <- filtered_data %>%
  left_join(total_counts, by = c("prediction_type", "tool", "label"))

virify_data_with_total <- virify_data_with_total %>%
  mutate(relative_abundance = (counts / total_count) * 100)

sum_relative_abundance <- virify_data_with_total %>%
  group_by(taxonomy) %>%
  summarize(total_abundance = sum(relative_abundance)) %>%
  ungroup() %>%
  arrange(desc(total_abundance))

top_20_taxonomy <- sum_relative_abundance %>%
  head(20) %>%
  pull(taxonomy)

virify_data_top_20 <- virify_data_with_total %>%
  filter(taxonomy %in% top_20_taxonomy)

# To sort the labels:
virify_data_top_20$taxonomy <- factor(
  virify_data_top_20$taxonomy,
  levels = rev(top_20_taxonomy)  # Reverse order to plot most abundant on top
)


my_palette <- c("#000000", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6")

pdf("virify_taxo.pdf", width = 7, height = 5, pointsize = 5)
ggplot(virify_data_top_20, aes(x = label, y = relative_abundance, fill = taxonomy)) +
  geom_col() +
  facet_grid(. ~ prediction_type, labeller = labeller(prediction_type = c("prophage" = "Prophage", "viral_genome" = "Viral Genome"))) +
  scale_fill_manual(values = my_palette) +
  labs(fill = "", y = "Relative abundance (%)", x = "")
dev.off()


### To write the table without filtering out unclass
sum_relative_abundance <- virify_data_with_total %>%
  group_by(taxonomy) %>%
  summarize(total_abundance = sum(relative_abundance))


sorted_taxonomy_levels <- sum_relative_abundance %>%
  arrange(desc(total_abundance)) %>%
  pull(taxonomy)

virify_data_with_total$taxonomy <- factor(
  virify_data_with_total$taxonomy,
  levels = sorted_taxonomy_levels
)

write.table(virify_data_with_total, file = 'virify_data_with_total.txt', sep = "\t", quote = FALSE, row.names = FALSE)

pdf("virify_taxo.pdf", width = 15, height = 15, pointsize = 8)
ggplot(virify_data_with_total, aes(x = label, y = relative_abundance, fill = taxonomy)) + geom_col() + facet_grid(. ~ prediction_type)
dev.off()

# For genomad:
table<-read.table("corr_concat_gev2_venn_counts_taxo.txt", header=TRUE)

genomad_data = subset(table, tool == "geNomad")

filtered_data <- genomad_data %>%
  filter(taxonomy != "Unclassified")

total_counts <- filtered_data %>%
  group_by(prediction_type, tool, label) %>%
  summarize(total_count = sum(counts), .groups = 'drop') %>%
  ungroup()

genomad_data_with_total <- filtered_data %>%
  left_join(total_counts, by = c("prediction_type", "tool", "label"))

genomad_data_with_total <- genomad_data_with_total %>%
  mutate(relative_abundance = (counts / total_count) * 100)

sum_relative_abundance <- genomad_data_with_total %>%
  group_by(taxonomy) %>%
  summarize(total_abundance = sum(relative_abundance)) %>%
  ungroup() %>%
  arrange(desc(total_abundance))

top_20_taxonomy <- sum_relative_abundance %>%
  head(20) %>%
  pull(taxonomy)

genomad_data_top_20 <- genomad_data_with_total %>%
  filter(taxonomy %in% top_20_taxonomy)

# To sort the labels:
genomad_data_top_20$taxonomy <- factor(
  genomad_data_top_20$taxonomy,
  levels = rev(top_20_taxonomy)  # Reverse order to plot most abundant on top
)


my_palette <- c("#000000", "#6A3D9A", "#FFFF99", "#B15928", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6")

pdf("genomad_taxo.pdf", width = 7, height = 5, pointsize = 5)
ggplot(genomad_data_top_20, aes(x = label, y = relative_abundance, fill = taxonomy)) +
  geom_col() +
  facet_grid(. ~ prediction_type, labeller = labeller(prediction_type = c("prophage" = "Prophage", "viral_genome" = "Viral Genome"))) +
  scale_fill_manual(values = my_palette) +
  labs(fill = "", y = "Relative abundance (%)", x = "")
dev.off()


# Taxonomy barplots

table<-read.table("corr_concat_gev2_venn_counts_taxo.txt", header=TRUE)
genomad_data = subset(table, tool == "geNomad")

total_counts <- genomad_data %>%
  group_by(prediction_type, tool, label) %>%
  summarize(total_count = sum(counts), .groups = 'drop') %>%
  ungroup()

genomad_data_with_total <- genomad_data %>%
  left_join(total_counts, by = c("prediction_type", "tool", "label"))

genomad_data_with_total <- genomad_data_with_total %>%
  mutate(relative_abundance = (counts / total_count) * 100)

# To sort the labels:
sum_relative_abundance <- genomad_data_with_total %>%
  group_by(taxonomy) %>%
  summarize(total_abundance = sum(relative_abundance))

sorted_taxonomy_levels <- sum_relative_abundance %>%
  arrange(desc(total_abundance)) %>%
  pull(taxonomy)

genomad_data_with_total$taxonomy <- factor(
  genomad_data_with_total$taxonomy,
  levels = sorted_taxonomy_levels
)

write.table(genomad_data_with_total, file = 'genomad_data_with_total.txt', sep = "\t", quote = FALSE, row.names = FALSE)

pdf("genomad_taxo.pdf", width = 15, height = 15, pointsize = 8)
ggplot(genomad_data_with_total, aes(x = label, y = relative_abundance, fill = taxonomy)) + geom_col() + facet_grid(. ~ prediction_type)
dev.off()



