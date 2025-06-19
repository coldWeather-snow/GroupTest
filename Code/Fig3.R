library(ggplot2)
library(dplyr)
library(tidyr)

eff_labels <- tribble(
  ~M, ~q, ~rho, ~efficiency,
  150, 0, 0, 1,
  150, 0, 0.001, 0.998,
  150, 0, 0.01, 0.978,
  150, 0.8, 0, 1,
  150, 0.8, 0.001, 0.979,
  150, 0.8, 0.01, 0.86,
  61, 0, 0, 1,
  61, 0, 0.001, 0.998,
  61, 0, 0.01, 0.981,
  61, 0.8, 0, 1,
  61, 0.8, 0.001, 0.964,
  61, 0.8, 0.01, 0.791
)

# Add formatted factors for merging with support_long
eff_labels <- eff_labels %>%
  mutate(
    q_f = factor(q, labels = c("q = 0", "q = 0.8")),
    rho_f = as.factor(rho),
    M_f = as.factor(M)
  )

# Manual entry of support points (expand as needed)
support_data <- tribble(
  ~M,  ~q,     ~rho,      ~support_points,
  150, 0,      0,         c(1, 17, 150),
  150, 0,      0.001,     c(1, 17, 148, 149, 150),
  150, 0,      0.01,      c(1, 16, 17, 18, 141,142,143,144,145,146,147,148,149,150),
  150, 0.8,    0,         c(1, 8, 78),
  150, 0.8,    0.001,     c(1, 7, 8, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86),
  150, 0.8,    0.01,      c(1,6,7,8,9,10,11,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106),
  61,  0,      0,         c(1, 16, 61),
  61,  0,      0.001,     c(1, 16, 61),
  61,  0,      0.01,      c(1, 15, 16, 61),
  61,  0.8,    0,         c(1, 7, 61),
  61,  0.8,    0.001,     c(1, 7, 8, 59, 60, 61),
  61,  0.8,    0.01,      c(1,6,7,8,9,10,53,54,55,56,57,58,59,60,61)
)

# Expand the list columns to long format
support_long <- support_data %>%
  unnest(cols = support_points)

# For pretty facets
support_long$q_f <- factor(support_long$q, labels=c("q = 0", "q = 0.8"))
support_long$rho_f <- as.factor(support_long$rho)
support_long$M_f <- as.factor(support_long$M)

library(ggplot2)
library(dplyr)

# Find max support per facet for positioning text
eff_labels <- eff_labels %>%
  left_join(support_long %>%
              group_by(M_f, q_f) %>%
              summarize(xmax = max(support_points, na.rm = TRUE)),
            by = c("M_f", "q_f"))

library(dplyr)



# Get the maximum rho value for positioning the label
eff_labels_pos <- eff_labels %>%
  group_by(M_f, q_f) %>%
  mutate(y_eff = max(as.numeric(as.character(rho_f))) + 0.005)  # adjust offset as needed

# Define a named vector for q_f
q_labels <- c("q = 0" = "q: 0", "q = 0.8" = "q: 0.8")
M_labels <- c("61" = "M: 61", "150" = "M: 150")

ggplot(support_long, aes(x = support_points, y = rho_f, color = rho_f)) +
  geom_jitter(height = 0.1, size = 3, width = 0, alpha = 0.8) +
  facet_grid(M_f ~ q_f, scales = "free_x", 
             labeller = labeller(q_f = q_labels, M_f = M_labels)) +
  # (rest of your plotting code unchanged)
  geom_text(
    data = eff_labels,
    aes(x = xmax + 5, y = rho_f, label = paste("Eff:", efficiency)),
    color = "black", size = 3.5, hjust = 0, fontface = "bold", inherit.aes = FALSE
  ) +
  labs(
    # title = "Support Points for Robust E-optimal Designs",
    x = "Support Point",
    y = expression(rho),
    color = expression(rho)
  ) +
  scale_color_manual(values = c("red", "green3", "blue")) +
  theme_bw(base_size = 14) +
  theme(strip.text = element_text(face = "bold")) + 
  xlim(0, 200) + 
  theme(legend.position = "none")
