
# Script Summary
## Plotting mutations

library(tidyverse)
library(here)
library(patchwork)
theme_set(
  theme_classic(base_size = 16) +
    theme(plot.subtitle = element_text(size = 14))
)

# read mutation data -output from variant_filtering script
mutations <- readxl::read_excel(
  here("output/true_mutation_calls_tag_amp_seq_my_name.xlsx"),
  col_type = "text"
) %>% 
  mutate( #extract the variant and gene name - adapt depending on the format of variant description
    variant_gene = str_remove(VariantDescription, ":.*"),
    variant_name = str_extract(VariantDescription, "p\\.[A-Z]\\d+[A-Z]"),
    variant_name = str_remove(variant_name, fixed("p."))
  ) %>% 
  select( # reordering columns
    sample_id := Sample,
    variant_gene, variant_name, 
    variant_id := Variant,
    variant_description := VariantDescription,
    everything()
  )

# number of total patients (regardless of true mutations)
n_patients <- readxl::read_excel(
  here("output/mutation_calls_tag_amp_seq_my_name.xlsx"),
  col_type = "text"
) %>% 
  pull(Sample) %>%
  n_distinct()

total_mutated_patients <- n_distinct(mutations$sample_id)

n_real_mutations <- nrow(mutations)


# create output dir

dir.create("output/mutations-plots", showWarnings = FALSE)

## gene level
gene_level_summary <- mutations %>% 
  group_by(variant_gene) %>% 
  summarise(
    n_mutated_patients = n_distinct(sample_id),
    p_mutated_patients = n_mutated_patients / n_patients
  )
p1 <- gene_level_summary %>% 
  ggplot(aes(fct_reorder(variant_gene, -p_mutated_patients), p_mutated_patients)) +
  geom_col(width = .85) +
  geom_text(
    aes(
      label = paste0(
        round(p_mutated_patients*100, 1), "%\n(",
        n_mutated_patients, " patients)"
      ),
      y = p_mutated_patients + .125
    )
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, .1)) +
  labs(
    x = NULL,
    y = NULL,
    subtitle = paste0(
      "Fraction of patients with each mutated gene\n",
      "\t", total_mutated_patients, " of ", n_patients,
      " patients had at least one mutation (",
      round(100 * total_mutated_patients/n_patients, 1),
      "%)"
    )
  )

ggsave(
  here("output/mutations-plots/gene-frequency.png"),
  p1, width = 5.5, height = 4.5, dpi = 600
)

# protein level
variant_level_summary <- mutations %>% 
  group_by(variant_gene, variant_name) %>% 
  summarise(
    n_mutated_patients = n_distinct(sample_id),
    p_mutated_patients = n_mutated_patients / n_patients
  ) %>% 
  ungroup()
p2 <- variant_level_summary %>% 
  ggplot(aes(fct_reorder(variant_name, p_mutated_patients), p_mutated_patients)) +
  geom_col(width = .85) +
  facet_wrap(
    ~paste0(
      variant_gene, "\n",
      column_to_rownames(gene_level_summary, "variant_gene")[variant_gene, "n_mutated_patients"],
      " mutated patients (",
      round(
        100 * column_to_rownames(gene_level_summary, "variant_gene")[variant_gene, "p_mutated_patients"],
        1
      ),
      "%)"
    ),
    scales = "free_y"
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL,
    subtitle = "Fraction of patients with each mutation"
  ) +
  coord_flip()

ggsave(
  here("output/mutations-plots/variant-frequency.png"),
  p2, width = 8, height = 11, dpi = 600
)

# plot 

plot_gene <- function(.gene) {
  .variants_summary <- variant_level_summary %>% 
    filter(variant_gene == .gene)
  .gene_n_mutated <- gene_level_summary %>% 
    filter(variant_gene == .gene) %>% 
    pull(p_mutated_patients)
  .gene_p_mutated <- gene_level_summary %>% 
    filter(variant_gene == .gene) %>% 
    pull(p_mutated_patients)
  
  # add not mutated 
  .plot_df <- .variants_summary %>% 
    add_row(
      variant_gene = .gene,
      variant_name = "Not mutated",
      n_mutated_patients = n_patients - .gene_n_mutated,
      p_mutated_patients = 1 - .gene_p_mutated
    ) %>%
    mutate(
      variant_gene = paste0(
        variant_gene,
        "\n(",
        round(.gene_p_mutated * 100, 1),
        "%)"
      ),
      variant_name = paste0(
        variant_name,
        " (",
        round(100 * p_mutated_patients, 1),
        "%)"
      ),
      variant_name = fct_reorder(variant_name, p_mutated_patients)
    )
  
  vars <- unique(.plot_df[["variant_name"]])
  n_vars <- length(vars)
  colors <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(8, "Dark2")
  )(n_vars)
  names(colors) <- vars
  colors[str_detect(names(colors), "Not mutated")] <- "gray80"
  
  .plot_df %>% 
    ggplot(aes(variant_gene, n_mutated_patients, fill = variant_name)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(10)
    ) +
    scale_fill_manual(
      values = colors
    ) +
    guides(
      fill = guide_legend(ncol = 1, reverse = TRUE)
    ) +
    coord_cartesian(ylim = c(0, 825)) +
    theme_classic(base_size = 14) +
    labs(
      x = NULL, y = "Number of patients",
      fill = paste0(.gene, " variant")
    ) +
    theme(legend.text = element_text(size = 8)) +
    guides(fill = guide_legend(ncol = 1 + as.numeric(n_vars > 10)))
}

gene_plots <- list()
for (.gene in unique(gene_level_summary$variant_gene)) {
  gene_plots[[.gene]] <- plot_gene(.gene)
}

ggsave(
  here("output/mutations-plots/stacked-barplot-by-variant-and-gene.png"),
  patchwork::wrap_plots(plotlist = gene_plots, nrow = 1),
  width = 16, height = 6
)