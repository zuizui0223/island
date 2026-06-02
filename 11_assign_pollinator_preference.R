############################################################
## Pollinator preference analysis (Dirichlet posterior)
##  - Exclude pollination_main == "unknown"
##  - Apply to:
##      (A) flower_color_simple (categorical)
##      (B) flower_shape_simple (categorical)
############################################################

## -------------------------
## 0) Packages
## -------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(MCMCpack)   # rdirichlet
})

set.seed(1)

## -------------------------
## 1) Load data
## -------------------------
infile <- file.path(Sys.getenv("HOME"),
                    "Desktop/list_with_traits_island_distance_bom_with_PDI_reclassified_isolated.csv")
dat <- fread(infile)

## -------------------------
## 2) Safety check
## -------------------------
need_cols <- c("pollination_main",
               "flower_color_simple",
               "flower_shape_simple")
miss <- setdiff(need_cols, names(dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))

## -------------------------
## 3) Exclude unknown pollination
## -------------------------
dat2 <- dat %>%
  filter(!is.na(pollination_main),
         pollination_main != "unknown")

## Fix pollination order (no unknown)
poll_levels <- c("bees","birds","bumblebees","butterflies",
                 "flies","mixed","moths","self","wind")
dat2 <- dat2 %>%
  mutate(pollination_main = factor(pollination_main, levels = poll_levels))

## -------------------------
## 4) Category definitions
## -------------------------
## (A) Color categories
color_levels <- c(
  "blue_purple",
  "green_brown_inconspicuous",
  "red_pink",
  "white",
  "yellow_orange"
)

## (B) Shape categories（←必要に応じて微調整してOK）
shape_levels <- c(
  "actinomorphic_open",
  "brush_puff",
  "composite_head",
  "reduced_wind",  "tubular",
  "zygomorphic"
)

## Keep only defined categories
dat2 <- dat2 %>%
  filter(flower_color_simple %in% color_levels,
         flower_shape_simple %in% shape_levels) %>%
  mutate(
    flower_color_simple = factor(flower_color_simple, levels = color_levels),
    flower_shape_simple = factor(flower_shape_simple, levels = shape_levels)
  )

## -------------------------
## 5) Dirichlet preference function
## -------------------------
dirichlet_preference <- function(counts, cats, ndraw = 50000, alpha0 = 1){
  k <- length(counts)
  u <- rep(1/k, k)     # uniform baseline
  
  post <- rdirichlet(n = ndraw, alpha = counts + alpha0)
  
  tibble(
    category = cats,
    prob_greater_than_uniform =
      colMeans(post > matrix(u, nrow(post), k, byrow = TRUE)),
    post_median = apply(post, 2, median),
    post_q025   = apply(post, 2, quantile, probs = 0.025),
    post_q975   = apply(post, 2, quantile, probs = 0.975)
  )
}

## -------------------------
## 6A) COLOR analysis
## -------------------------
tab_color <- dat2 %>%
  count(pollination_main, flower_color_simple, name="n") %>%
  complete(pollination_main, flower_color_simple, fill=list(n=0)) %>%
  group_by(pollination_main) %>%
  mutate(
    N = sum(n),
    prop = n / N
  ) %>%
  ungroup()

pref_color <- tab_color %>%
  arrange(pollination_main, flower_color_simple) %>%
  group_by(pollination_main) %>%
  summarise(
    pref = list(
      dirichlet_preference(
        counts = n,
        cats   = as.character(flower_color_simple)
      )
    ),
    .groups="drop"
  ) %>%
  unnest(pref) %>%
  rename(flower_color_simple = category)

## -------------------------
## 6B) SHAPE analysis
## -------------------------
tab_shape <- dat2 %>%
  count(pollination_main, flower_shape_simple, name="n") %>%
  complete(pollination_main, flower_shape_simple, fill=list(n=0)) %>%
  group_by(pollination_main) %>%
  mutate(
    N = sum(n),
    prop = n / N
  ) %>%
  ungroup()

pref_shape <- tab_shape %>%
  arrange(pollination_main, flower_shape_simple) %>%
  group_by(pollination_main) %>%
  summarise(
    pref = list(
      dirichlet_preference(
        counts = n,
        cats   = as.character(flower_shape_simple)
      )
    ),
    .groups="drop"
  ) %>%
  unnest(pref) %>%
  rename(flower_shape_simple = category)

## -------------------------
## 7) Plots
## -------------------------
## (A) Color heatmap
p_color_heat <- ggplot(pref_color,
                       aes(x = flower_color_simple,
                           y = pollination_main,
                           fill = prob_greater_than_uniform)) +
  geom_tile(color="grey85", linewidth=0.3) +
  scale_fill_gradient2(
    low="#2c7bb6", mid="white", high="#d7191c",
    midpoint=0.5, limits=c(0,1),
    name="P(pref > uniform)"
  ) +
  labs(title="Bayesian evidence for pollinator COLOR preference",
       x="Flower color category", y="Pollination mode") +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

## (B) Shape heatmap
p_shape_heat <- ggplot(pref_shape,
                       aes(x = flower_shape_simple,
                           y = pollination_main,
                           fill = prob_greater_than_uniform)) +
  geom_tile(color="grey85", linewidth=0.3) +
  scale_fill_gradient2(
    low="#2c7bb6", mid="white", high="#d7191c",
    midpoint=0.5, limits=c(0,1),
    name="P(pref > uniform)"
  ) +
  labs(title="Bayesian evidence for pollinator SHAPE preference",
       x="Flower shape category", y="Pollination mode") +
  theme_classic(base_size=14) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

## -------------------------
## 8) Save outputs
## -------------------------
out_dir <- file.path(Sys.getenv("HOME"), "Desktop/pollinator_preference_no_unknown")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

fwrite(as.data.table(tab_color), file.path(out_dir, "tab_color_counts_prop.csv"))
fwrite(as.data.table(pref_color), file.path(out_dir, "pref_color_dirichlet.csv"))

fwrite(as.data.table(tab_shape), file.path(out_dir, "tab_shape_counts_prop.csv"))
fwrite(as.data.table(pref_shape), file.path(out_dir, "pref_shape_dirichlet.csv"))

ggsave(file.path(out_dir, "heatmap_color_preference.pdf"),
       p_color_heat, width=10, height=6)
ggsave(file.path(out_dir, "heatmap_shape_preference.pdf"),
       p_shape_heat, width=10, height=6)

cat("\n✅ Finished. Results saved in:\n", out_dir, "\n")