# =============================================================================
# TCMDATA Logo Generator - Fresh & Clean Design
# =============================================================================

library(hexSticker)
library(ggplot2)
library(showtext)

# Clean academic font
font_add_google("Montserrat", "montserrat")
showtext_auto()

# -----------------------------------------------------------------------------
# Fresh Herb Leaf Design - Light & Clean
# -----------------------------------------------------------------------------

# Elegant leaf shape
leaf_data <- data.frame(
  x = c(0, 0.25, 0.4, 0.35, 0.2, 0, -0.2, -0.35, -0.4, -0.25),
  y = c(0.9, 0.65, 0.25, -0.15, -0.45, -0.6, -0.45, -0.15, 0.25, 0.65)
)

# Network nodes on the leaf
nodes_data <- data.frame(
  x = c(0, 0.12, -0.12, 0.08, -0.08, 0),
  y = c(0.45, 0.15, 0.15, -0.15, -0.15, -0.35),
  size = c(3.5, 2.5, 2.5, 2, 2, 1.8)
)

# Connecting lines
edges_data <- data.frame(
  x =    c(0,     0,     0.12,  -0.12, 0.12,  -0.12),
  xend = c(0.12,  -0.12, 0.08,  -0.08, 0,     0),
  y =    c(0.45,  0.45,  0.15,  0.15,  0.15,  0.15),
  yend = c(0.15,  0.15,  -0.15, -0.15, -0.35, -0.35)
)

p_herb <- ggplot() +
  # Leaf shape - soft mint green
  geom_polygon(data = leaf_data, aes(x, y),
               fill = "#7FBFB3", color = "#5A9E91", 
               alpha = 0.9, linewidth = 0.8) +
  # Center vein
  geom_segment(aes(x = 0, xend = 0, y = 0.8, yend = -0.5),
               color = "#5A9E91", linewidth = 0.5, alpha = 0.4) +
  # Network edges - subtle sage
  geom_segment(data = edges_data, 
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "#3D6B5E", linewidth = 0.9, alpha = 0.85) +
  # Network nodes - clean white
  geom_point(data = nodes_data, aes(x, y, size = size),
             color = "#FFFFFF", alpha = 0.95) +
  scale_size_identity() +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(legend.position = "none")

# Generate final logo
sticker(
  subplot = p_herb,
  package = "TCMDATA",
  p_size = 22,
  p_y = 1.48,
  p_color = "#2D5A4E",       # Dark teal text
  p_family = "montserrat",
  p_fontface = "bold",
  s_x = 1,
  s_y = 0.72,
  s_width = 1.1,
  s_height = 0.95,
  h_fill = "#E8F5F1",        # Very light mint background
  h_color = "#7FBFB3",       # Soft green border
  h_size = 1.8,
  spotlight = FALSE,
  filename = "man/figures/logo.png",
  dpi = 300
)

message("✓ Logo saved to: man/figures/logo.png")
