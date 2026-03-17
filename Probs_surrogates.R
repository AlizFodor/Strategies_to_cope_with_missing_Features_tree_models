
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(patchwork)

# --- user chosen parameters 

m <- 1  # m = number of features at a node
K <- 5 # number of stored surrogates per node
k <- 20 # tree size
p_values    <- seq(0, 0.9, by = 0.01)

f_fractions <- c(0.05, 0.2, 0.4, 0.6, 0.8) # fractions of missing features

# ----------------------------------------

# probability that the primary split is missing at a node
q_node <- function(p, m = 1) 1 - (1 - p)^m

# probabiltiy that at least one surrogate feature is observed given K candidates
s_surrogate_obs <- function(p, K = 5) 1 - p^K

# number of nodes visized as function of training feature fraction
k_approxi <- function(f_frac, maxdepth, scale = 5, f_ref_max = max(f_fractions)) {
  round(maxdepth * scale * (f_frac / f_ref_max))}

# primary / surrogate / fallback probabilities at a tree path level
theory_probs <- function(p, m = 1, K = 5, k = 10) {
  q <- q_node(p, m)
  r <- 1 - (1 - q)^k        
  s <- s_surrogate_obs(p, K)
  
  primary   <- (1 - q)^k
  surrogate <- r * s
  fallback  <- r * (1 - s)
  
  tibble(primary = primary, surrogate = surrogate, fallback = fallback)
}


df_p <- tibble(p = p_values) %>%
  rowwise() %>%
  mutate(out = list(theory_probs(p, m = m, K = K, k = k))) %>%
  unnest(out) %>%
  ungroup() %>%
  pivot_longer(cols = c(primary, surrogate, fallback),
               names_to = "route", values_to = "prob") %>%
  mutate(route = factor(route, levels = c("primary", "surrogate", "fallback")))


plot_p <- ggplot(df_p, aes(x = p, y = prob, color = route)) +
  geom_line(linewidth = 1.2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_viridis_d(
    option = "D", end = 0.9,
    labels = c(primary = "Primary", surrogate = "Surrogate", fallback = "Fallback")
  ) +
  labs(
    title = "Expected routing probabilities with increasing feature missingness",
    subtitle = paste0("m = ", m, " number of features per split, K = ", K, " surrogates per node, k = ", k, " visited nodes"),
    x = "Probability a feature is missing (p)",
    y = "Probability",
    color = "Routing"
  ) +
  theme_minimal(base_size = 27) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

plot_p


ggsave(
  filename = "expected_probs4.png",
  plot = plot_p,
  path = "Plots/",
  width = 4800,
  height = 3000,
  units = "px",
  dpi = 300
)


















