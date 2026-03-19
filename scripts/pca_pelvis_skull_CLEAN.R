# open pelvis csizes
csize_pelvis <- readRDS("/your/working/directory/csize_pelvises.rds")
# pelvis data
pelvis_all <- readRDS( "/your/working/directory/pelvises_all.rds")
# skull data
skull_all<- readRDS("/your/working/directory/skull_all.rds")
# skull csize
csize_skulls <- readRDS( "/your/working/directory/csize_skulls.rds" )
# csize neurocranium
csize_neuro <- readRDS("/your/working/directory/csize_neuro.rds")
# classifiers
classifiers_all <- readRDS("/your/working/directory/classifiers_all.rds")


# Mean-centering by Generation nested within Strain ----------------------

regr1 <- procD.lm(pelvis_all$coords ~ classifiers_all$Generation)
pelvis_all_resid <- residualize(regr1$residuals, pelvis_all)
pelvis_all <- gpagen(pelvis_all_resid$coords)

regr1 <- procD.lm(skull_all$coords ~ classifiers_all$Generation)
skull_all_resid <- residualize(regr1$residuals, skull_all)
skull_all <- gpagen(skull_all_resid$coords)


# PCA ---------------------------------------------------------------------

pca_pelvis <- gm.prcomp(pelvis_all$coords)
pca_skull  <- gm.prcomp(skull_all$coords)

plot(pca_pelvis)
plot(pca_skull)


# Dataframes for plotting -------------------------------------------------

df_pel <- data.frame(
  X = pca_pelvis$x[,1],
  Y = pca_pelvis$x[,2],
  Generation = classifiers_all$Generation,
  Strain = classifiers_all$Strain,
  Sex = classifiers_all$sex
)

df_sku <- data.frame(
  X = pca_skull$x[,1],
  Y = pca_skull$x[,2],
  Generation = classifiers_all$Generation,
  Strain = classifiers_all$Strain,
  Sex = classifiers_all$sex
)


# Plot settings -----------------------------------------------------------

custom_colors <- c(
  Parental = "blue",
  F1 = "gold1",
  F2 = "red",
  B1 = "olivedrab3",
  B2 = "purple",
  B3 = "darkgreen"
)

ellipse_fill_alpha <- 0.10
ellipse_line_alpha <- 0.40


# Function to generate PCA plots -----------------------------------------

plot_pca <- function(df){
  
  ggplot(df, aes(x = X, y = Y, color = Generation)) +
    
    stat_ellipse(
      aes(fill = Generation),
      geom = "polygon",
      alpha = ellipse_fill_alpha,
      color = NA,
      level = 0.99
    ) +
    
    stat_ellipse(
      geom = "path",
      alpha = ellipse_line_alpha,
      linewidth = 0.6,
      level = 0.99
    ) +
    
    geom_point(size = 1.5, alpha = 0.4) +
    
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    
    theme_classic(base_size = 14) +
    theme(legend.position = "none")
}


# Full datasets -----------------------------------------------------------

pel_all <- plot_pca(df_pel)
sku_all <- plot_pca(df_sku)

pel_all
sku_all

# Females -----------------------------------------------------------------

df_pel_f <- df_pel %>% dplyr::filter(Sex == "F")
df_sku_f <- df_sku %>% dplyr::filter(Sex == "F")

fem_pel <- plot_pca(df_pel_f)
fem_sku <- plot_pca(df_sku_f)

fem_pel
fem_sku

# Males -------------------------------------------------------------------

df_pel_m <- df_pel %>% dplyr::filter(Sex == "M")
df_sku_m <- df_sku %>% dplyr::filter(Sex == "M")

mas_pel <- plot_pca(df_pel_m)
mas_sku <- plot_pca(df_sku_m)

mas_pel
mas_sku

# Combine plots -----------------------------------------------------------

grid <- plot_grid(
  pel_all, sku_all,
  fem_pel, fem_sku,
  mas_pel, mas_sku,
  ncol = 2
)

grid

ggsave(
  plot = grid,
  filename = "figs/pca_pelvis_skull_expanded_conf_interval.pdf",
  height = 18,
  width = 20,
  units = "cm"
)

