library(geomorph)
library(dplyr)
library(ggplot2)

# ============================================================
# OPEN FILES
# ============================================================
csize_pelvises           <- readRDS("/your/working/directory/csize_pelvises2.rds")
pelvis_all_csized        <- readRDS("/your/working/directory/pelvis_all_csized.rds")
skull_all_csized         <- readRDS("/your/working/directory/skull_all_csized.rds")
csize_skulls             <- readRDS("/your/working/directory/csize_skulls336.rds")
csize_neuro              <- readRDS("/your/working/directory/csize_neuro.rds")
classifiers_all_csized   <- readRDS("/your/working/directory/classifiers_all_csized_comp_w.rds")

# ============================================================
# SHARED SETTINGS
# ============================================================
custom_colors <- c(
  "Parental" = "blue",
  "F1"       = "gold1",
  "F2"       = "red",
  "B1"       = "olivedrab3",
  "B2"       = "purple",
  "B3"       = "darkgreen"
)
gen_levels_recode <- c("F0", "F1", "F2", "B1", "B2", "B3")
ellipse_fill_alpha <- 0.1
ellipse_line_alpha <- 0.40

# helper: compute correlation + CI + MSE per generation
compute_cor_mse <- function(df_pls, classifiers) {
  df <- data.frame(
    X          = df_pls[, 1],
    Y          = df_pls[, 2],
    Generation = as.factor(classifiers$Generation)
  )
  
  cor_df <- df %>%
    group_by(Generation) %>%
    summarize(
      n      = n(),
      r      = cor(X, Y, method="pearson"),
      z      = atanh(r),
      se     = 1 / sqrt(n - 3),
      r_low  = tanh(z - 1.96*se),
      r_high = tanh(z + 1.96*se),
      .groups = "drop"
    ) %>%
    mutate(Generation = recode(Generation, Parental="F0"),
           Generation = factor(Generation, levels=gen_levels_recode)) %>%
    rename(correlation = r)
  
  mse_df <- df %>%
    group_by(Generation) %>%
    do({ data.frame(MSE = mean(resid(lm(Y ~ X, data=.))^2)) }) %>%
    ungroup() %>%
    mutate(Generation = recode(Generation, Parental="F0"),
           Generation = factor(Generation, levels=gen_levels_recode))
  
  left_join(cor_df, mse_df, by="Generation")
}

# helper: PLS scatter plot with correlation lines
plot_pls <- function(pls_obj, classifiers, title="") {
  df <- data.frame(
    X          = pls_obj$XScores[, 1],
    Y          = pls_obj$YScores[, 1],
    Generation = as.factor(classifiers$Generation)
  )
  slopes <- df %>%
    group_by(Generation) %>%
    summarize(r=cor(X, Y), sdX=sd(X), sdY=sd(Y), .groups="drop") %>%
    mutate(slope = r * (sdY / sdX))
  
  ggplot(df, aes(x=X, y=Y, color=Generation)) +
    stat_ellipse(aes(fill=Generation), geom="polygon",
                 alpha=ellipse_fill_alpha, color=NA, level=0.99) +
    stat_ellipse(geom="path", alpha=ellipse_line_alpha,
                 linewidth=0.6, level=0.99) +
    geom_point(size=3, alpha=0.4) +
    geom_abline(data=slopes,
                aes(slope=r, intercept=0, color=Generation),
                linewidth=1.2) +
    scale_color_manual(values=custom_colors) +
    scale_fill_manual(values=custom_colors) +
    theme_classic(base_size=14) +
    theme(axis.text.x=element_text(size=21),
          axis.text.y=element_text(size=21),
          legend.position="none") +
    labs(x="Pelvis PLS1 Scores", y="Skull PLS1 Scores", title=title)
}

# helper: correlation + MSE summary plot
plot_cor_mse <- function(plot_data) {
  ggplot(plot_data, aes(x=Generation)) +
    geom_col(aes(y=correlation), fill="grey", width=0.5) +
    geom_errorbar(aes(ymin=r_low, ymax=r_high),
                  width=0.15, color="grey22", linewidth=0.4) +
    geom_point(aes(y=correlation), color="grey22", size=2) +
    geom_line(aes(y=(MSE / max(MSE)) * max(correlation), group=1),
              color="red", linewidth=0.4) +
    geom_point(aes(y=(MSE / max(MSE)) * max(correlation)),
               color="red", size=2) +
    scale_y_continuous(
      name     = "Correlation (Pearson r)",
      sec.axis = sec_axis(~. / max(plot_data$correlation) * max(plot_data$MSE),
                          name="Mean Square Error")
    ) +
    labs(x="Generation") +
    theme_minimal(base_size=14) +
    theme(axis.title.y.left  = element_text(color="black"),
          axis.title.y.right = element_text(color="black"),
          legend.position    = "none")
}

# helper: morph + heatmap
plot_morph_heatmap <- function(pls_obj, side=c("pelvis","skull"),
                               gpa_resid, atlas_lms, atlas_mesh,
                               scale_factor, zoom=0.6) {
  side <- match.arg(side)
  sdev <- sd(if(side=="pelvis") pls_obj$XScores[,1] else pls_obj$YScores[,1])
  vec  <- t(matrix(if(side=="pelvis") pls_obj$left.pls.vectors[,1]
                   else pls_obj$right.pls.vectors[,1], nrow=3))
  
  pos_morph <- gpa_resid$consensus + vec * sdev * scale_factor
  neg_morph <- gpa_resid$consensus - vec * sdev * scale_factor
  
  sup_pos  <- pPsup(atlas_lms, pos_morph)
  mesh_pos <- tps3d(atlas_mesh, atlas_lms, sup_pos$Mp2)
  
  sup_neg  <- pPsup(atlas_lms, neg_morph)
  mesh_neg <- tps3d(atlas_mesh, atlas_lms, sup_neg$Mp2)
  
  open3d(zoom=zoom); title3d("POSITIVE")
  shade3d(mesh_pos, override=F, color="grey", alpha=1, specular=1)
  rglwidget()
  
  open3d(zoom=zoom); title3d("NEGATIVE")
  shade3d(mesh_neg, override=F, color="grey", alpha=1, specular=1)
  rglwidget()
  
  dist.vec <- meshDist(mesh_neg, mesh_pos, plot=FALSE)$dists
  to   <-  max(abs(dist.vec))
  from <- -to
  open3d(zoom=0.9)
  par3d("windowRect"=c(0,100,600,800))
  meshDist(mesh_pos, distvec=-dist.vec,
           rampcolors=c("blue","white","firebrick2"),
           from=from, to=to, steps=50,
           titleplot="Closest Point Distance (mm)",
           xaxt=5, plot=T, center=T, symmetric=T)
  rglwidget()
}

# ============================================================
# ALL SPECIMENS ####
# ============================================================
regr_pelvis_all <- procD.lm(pelvis_all_csized$coords ~
                              classifiers_all_csized$sex +
                              classifiers_all_csized$Strain +
                              classifiers_all_csized$parous.non_parous.y,
                            SS.type="III")
summary(regr_pelvis_all)
pelvis_all_resid_csized <- residualize(regr_pelvis_all$residuals, pelvis_all_csized)

regr_skull_all <- procD.lm(skull_all_csized$coords ~
                             classifiers_all_csized$sex +
                             classifiers_all_csized$Strain +
                             classifiers_all_csized$parous.non_parous.y,
                           SS.type="III")
summary(regr_skull_all)
skull_all_resid_csized <- residualize(regr_skull_all$residuals, skull_all_csized)

pls_all <- two.b.pls(pelvis_all_resid_csized$coords,
                     skull_all_resid_csized$coords,
                     iter=999, print.progress=TRUE)
summary(pls_all)
pls_all$svd$d[1]^2 / sum(pls_all$svd$d^2)

# plot
all_pls <- plot_pls(pls_all, classifiers_all_csized)
all_pls
ggsave(all_pls, filename="./figs/all_pls_corr2.pdf", height=20, width=20, unit="cm")

# correlation + MSE summary
plot_data_all <- compute_cor_mse(
  data.frame(pls_all$XScores[,1], pls_all$YScores[,1]),
  classifiers_all_csized)
plot_cor_mse(plot_data_all)

# morphs
pelvis_gpa_all <- gpagen(pelvis_all_resid_csized$coords, print.progress=FALSE)
skull_gpa_all  <- gpagen(skull_all_resid_csized$coords, print.progress=FALSE)

plot_morph_heatmap(pls_all, "pelvis", pelvis_gpa_all,
                   pelvis_atlas_lms, pelvis_atlas_mesh, scale_factor=3)
plot_morph_heatmap(pls_all, "skull", skull_gpa_all,
                   atlas_skull_lm[1:336,], skull_mesh, scale_factor=12)

# ============================================================
# NON-PAROUS FEMALES ####
# ============================================================
idx_Fnp <- which(classifiers_all_csized$sex == "F" &
                   classifiers_all_csized$parous.non_parous.y == "N")

pelvis_Fnp_gpa      <- gpagen(pelvis_all_csized$coords[,,idx_Fnp], print.progress=FALSE)
skull_Fnp_gpa       <- gpagen(skull_all_csized$coords[,,idx_Fnp], print.progress=FALSE)
classifiers_Fnp     <- classifiers_all_csized[idx_Fnp, ]

regr_pelvis_Fnp <- procD.lm(pelvis_Fnp_gpa$coords ~
                              classifiers_Fnp$Strain +
                              pelvis_Fnp_gpa$Csize, SS.type="III")
summary(regr_pelvis_Fnp)
pelvis_resid_Fnp <- residualize(regr_pelvis_Fnp$residuals, pelvis_Fnp_gpa)

regr_skull_Fnp <- procD.lm(skull_Fnp_gpa$coords ~
                             classifiers_Fnp$Strain +
                             pelvis_Fnp_gpa$Csize, SS.type="III")
summary(regr_skull_Fnp)
skull_resid_Fnp <- residualize(regr_skull_Fnp$residuals, skull_Fnp_gpa)

pls_Fnp <- two.b.pls(pelvis_resid_Fnp$coords, skull_resid_Fnp$coords,
                     iter=999, print.progress=TRUE)
summary(pls_Fnp)
pls_Fnp$svd$d[1]^2 / sum(pls_Fnp$svd$d^2)

# plot
fnp_pls <- plot_pls(pls_Fnp, classifiers_Fnp)
fnp_pls
ggsave(fnp_pls, filename="./figs/fnp_pls.pdf", height=20, width=20, unit="cm")

# correlation + MSE summary
plot_data_Fnp <- compute_cor_mse(
  data.frame(pls_Fnp$XScores[,1], pls_Fnp$YScores[,1]),
  classifiers_Fnp)
plot_cor_mse(plot_data_Fnp)

# morphs
pelvis_gpa_Fnp <- gpagen(pelvis_resid_Fnp$coords, print.progress=FALSE)
skull_gpa_Fnp  <- gpagen(skull_resid_Fnp$coords, print.progress=FALSE)

plot_morph_heatmap(pls_Fnp, "pelvis", pelvis_gpa_Fnp,
                   pelvis_atlas_lms, pelvis_atlas_mesh, scale_factor=2)
plot_morph_heatmap(pls_Fnp, "skull", skull_gpa_Fnp,
                   atlas_skull_lm[1:336,], skull_mesh, scale_factor=7)

# ============================================================
# MALES ####
# ============================================================
idx_M <- which(classifiers_all_csized$sex == "M" &
                 !is.na(classifiers_all_csized$weight))

pelvis_M_gpa    <- gpagen(pelvis_all_csized$coords[,,idx_M], print.progress=FALSE)
skull_M_gpa     <- gpagen(skull_all_csized$coords[,,idx_M], print.progress=FALSE)
classifiers_M   <- classifiers_all_csized[idx_M, ]

regr_pelvis_M <- procD.lm(pelvis_M_gpa$coords ~
                            classifiers_M$Strain +
                            classifiers_M$weight, SS.type="III")
summary(regr_pelvis_M)
pelvis_resid_M <- residualize(regr_pelvis_M$residuals, pelvis_M_gpa)

regr_skull_M <- procD.lm(skull_M_gpa$coords ~
                           classifiers_M$Strain +
                           classifiers_M$weight, SS.type="III")
summary(regr_skull_M)
skull_resid_M <- residualize(regr_skull_M$residuals, skull_M_gpa)

pls_M <- two.b.pls(pelvis_resid_M$coords, skull_resid_M$coords,
                   iter=999, print.progress=TRUE)
summary(pls_M)
pls_M$svd$d[1]^2 / sum(pls_M$svd$d^2)

# plot
mas_pls <- plot_pls(pls_M, classifiers_M)
mas_pls
ggsave(mas_pls, filename="./figs/mal_pls.pdf", height=20, width=20, unit="cm")

# correlation + MSE summary
plot_data_M <- compute_cor_mse(
  data.frame(pls_M$XScores[,1], pls_M$YScores[,1]),
  classifiers_M)
plot_cor_mse(plot_data_M)

# morphs
pelvis_gpa_M <- gpagen(pelvis_resid_M$coords, print.progress=FALSE)
skull_gpa_M  <- gpagen(skull_resid_M$coords, print.progress=FALSE)

plot_morph_heatmap(pls_M, "pelvis", pelvis_gpa_M,
                   pelvis_atlas_lms, pelvis_atlas_mesh, scale_factor=4)
plot_morph_heatmap(pls_M, "skull", skull_gpa_M,
                   atlas_skull_lm[1:336,], skull_mesh, scale_factor=10)