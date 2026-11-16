library(ggplot2)
library(grid)

GN_TITLE_SIZE    <- 18
GN_SUBTITLE_SIZE <- 14
GN_BASE_SIZE     <- 16
GN_SUBTITLE_COL  <- "grey35"
GN_STRIP_FILL    <- "#f0f0f0"
GN_GRID_COL      <- "grey85"

X_BETA <- "Beta (95% CI)"

BRK_BETA <- seq(-1, 1, by = 0.1)
BRK_R2   <- seq(-1, 1, by = 0.005)

theme_gn <- function(base       = GN_BASE_SIZE,
                     strip      = GN_BASE_SIZE,
                     legend     = "bottom",
                     grid       = c("x", "y", "both", "none"),
                     base_theme = theme_bw) {
  grid <- match.arg(grid)
  gline <- element_line(linetype = "dotted", colour = GN_GRID_COL)

  base_theme(base_size = base) +
    theme(
      plot.title       = element_text(size = GN_TITLE_SIZE, face = "bold"),
      plot.subtitle    = element_text(size = GN_SUBTITLE_SIZE, face = "plain",
                                      colour = GN_SUBTITLE_COL),
      plot.caption     = element_text(size = GN_SUBTITLE_SIZE, colour = GN_SUBTITLE_COL),

      axis.title       = element_text(size = base),
      axis.text        = element_text(size = base, colour = "black"),

      strip.text       = element_text(size = strip, face = "bold",
                                      margin = margin(t = 8, b = 8)),
      strip.background = element_rect(fill = GN_STRIP_FILL, colour = NA),

      legend.position  = legend,
      legend.direction = if (legend %in% c("bottom", "top")) "horizontal" else "vertical",
      legend.title     = element_text(size = base, face = "bold"),
      legend.text      = element_text(size = base, face = "plain"),

      panel.grid.minor   = element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "both")) gline else element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "both")) gline else element_blank(),

      panel.spacing    = unit(1.2, "lines"),
      plot.margin      = margin(t = 15, r = 25, b = 10, l = 15)
    )
}

theme_gn_corr <- function(base = GN_BASE_SIZE, strip = GN_BASE_SIZE, legend = "right") {
  theme_gn(base = base, strip = strip, legend = legend,
           grid = "none", base_theme = theme_minimal)
}

MAT_COL     <- "#8C4FAE"
PAT_COL     <- "#2E9E5D"
OVR_COL     <- "#DEB221"
CHD_POP_COL <- "#F4A87A"
CHD_DIR_COL <- "#D85A30"

SERIES4_COL <- c("Child PGS (Population)"         = CHD_POP_COL,
                 "Child PGS (Direct)"             = CHD_DIR_COL,
                 "Maternal PGS (Genetic Nurture)" = MAT_COL,
                 "Paternal PGS (Genetic Nurture)" = PAT_COL)

DOM_COL <- c("Cognitive Abilities"                = "#4C72B0",
             "Education Achievement & Attainment" = "#DD8452",
             "Height & BMI"                       = "#55A868",
             "Psychopathology"                    = "#C44E52")

DOM_SHORT <- c("Cognitive Abilities"                = "Cognitive",
               "Education Achievement & Attainment" = "Education",
               "Height & BMI"                       = "Height & BMI",
               "Psychopathology"                    = "Psychopathology")

lighten_gn <- function(hex, f = 0.45) {
  rgb_v <- grDevices::col2rgb(hex)[, 1]
  grDevices::rgb(t(rgb_v + (255 - rgb_v) * f), maxColorValue = 255)
}

sample_ramp <- function(col) {
  c(All_Trios     = lighten_gn(col, 0.60),
    All_Unrelated = lighten_gn(col, 0.30),
    Whole_Sample  = col)
}

SEX_COL <- c("Female Trios" = "#EC4C8C", "Male Trios" = "#3B76EF")
ZYG_COL <- c("MZ Trios"     = "#E1A11A", "DZ Trios"   = "#1F78B4")
SIB_COL     <- "#6FA8DC"

save_gn <- function(p, fname, width, height, dir = outFileStem) {
  ggsave(file.path(sub("/+$", "", dir), paste0(fname, ".pdf")),
         plot = p, device = "pdf",
         width = width, height = height, units = "in", bg = "white")
  invisible(p)
}

fdr_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < .001 ~ "***",
    p < .01  ~ "**",
    p < .05  ~ "*",
    TRUE     ~ ""
  )
}

cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
