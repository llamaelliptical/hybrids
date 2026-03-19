library(shapes)
library(geomorph)
library(ggplot2)
library(dplyr)

# ============================================================
# OPEN FILES
# ============================================================
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

# ============================================================
# SHARED SETTINGS
# ============================================================
custom_colors <- c(
  "F0" = "blue", "F1" = "gold1", "F2" = "red",
  "B1" = "olivedrab3", "B2" = "purple", "B3" = "darkgreen"
)
sex_shapes         <- c("M" = 19, "F" = 17)
ellipse_fill_alpha <- 0.08
ellipse_line_alpha <- 0.40

# ============================================================
# HELPERS
# ============================================================

# morph a mesh given a target shape
make_morph <- function(target, atlas_lms, atlas_mesh) {
  sup <- pPsup(atlas_lms, target)
  tps3d(atlas_mesh, atlas_lms, sup$Mp2)
}

# show a mesh in rgl
show_mesh <- function(mesh, color="grey", zoom=0.6) {
  open3d(zoom=zoom)
  shade3d(mesh, override=FALSE, color=color, alpha=1, specular=1)
  rglwidget()
}

# heatmap between two meshes
show_heatmap <- function(mesh_ref, mesh_comp, distvec_sign=1, zoom=0.9,
                         steps=50, title="Closest Point Distance (mm)") {
  dist.vec <- meshDist(mesh_comp, mesh_ref, plot=FALSE)$dists
  to   <-  max(abs(dist.vec))
  from <- -to
  open3d(zoom=zoom)
  par3d("windowRect"=c(0, 100, 600, 800))
  meshDist(mesh_ref,
           distvec    = distvec_sign * dist.vec,
           rampcolors = c("blue", "white", "firebrick2"),
           from=from, to=to, steps=steps,
           titleplot=title, xaxt=5,
           plot=TRUE, center=TRUE, symmetric=TRUE)
  rglwidget()
}

# sex dimorphism morphs + heatmap for one generation
sexdim_morphs <- function(female_consensus, male_consensus,
                          atlas_lms, atlas_mesh,
                          scale=0.5, heatmap_steps=30) {
  vector      <- female_consensus - male_consensus
  superfemale <- female_consensus + scale * vector
  supermale   <- male_consensus   - scale * vector
  
  mesh_f <- make_morph(superfemale, atlas_lms, atlas_mesh)
  mesh_m <- make_morph(supermale,   atlas_lms, atlas_mesh)
  
  show_mesh(mesh_f, color="grey")
  show_mesh(mesh_m, color="grey")
  show_heatmap(mesh_f, mesh_m, distvec_sign=1, steps=heatmap_steps)
}

# ============================================================
# SEXUAL DIMORPHISM: all specimens
# ============================================================
keep                     <- complete.cases(classifiers_all$weight)
coords_clean             <- pelvis_all$coords[,,keep]
classifiers_weight_clean <- classifiers_all[keep, ]

sexdiff <- procD.lm(
  coords_clean ~ classifiers_weight_clean$sex +
    classifiers_weight_clean$Strain +
    classifiers_weight_clean$weight,
  SS.type="III"
)
summary(sexdiff)

# ============================================================
# PCA: all pelvises by generation
# ============================================================
regr_strain     <- procD.lm(pelvis_all$coords ~ classifiers_all$Strain)
pelvis_all_resid <- residualize(regr_strain$residuals, pelvis_all)

pca_all <- gm.prcomp(pelvis_all_resid$coords)

df_pca <- data.frame(
  PC1        = pca_all$x[,1],
  PC2        = pca_all$x[,2],
  PC3        = pca_all$x[,3],
  Generation = factor(classifiers_all$Generation,
                      levels=c("Parental","F1","F2","B1","B2","B3"),
                      labels=c("F0","F1","F2","B1","B2","B3")),
  Sex        = factor(classifiers_all$sex, levels=c("F","M"))
)

centroids <- df_pca %>%
  group_by(Generation, Sex) %>%
  summarise(mean_PC1=mean(PC1), mean_PC2=mean(PC2), .groups="drop")

pca_plot <- ggplot(df_pca, aes(x=PC1, y=PC2, color=Generation)) +
  stat_ellipse(aes(fill=Generation), geom="polygon",
               alpha=ellipse_fill_alpha, color=NA, level=0.99) +
  stat_ellipse(geom="path", alpha=ellipse_line_alpha,
               linewidth=0.6, level=0.99) +
  geom_point(aes(shape=Sex), size=2.5, alpha=0.3) +
  geom_line(data=centroids,
            aes(x=mean_PC1, y=mean_PC2, group=Generation),
            linewidth=0.8) +
  geom_point(data=centroids,
             aes(x=mean_PC1, y=mean_PC2, shape=Sex),
             size=4, stroke=1.2) +
  scale_color_manual(values=custom_colors) +
  scale_fill_manual(values=custom_colors) +
  scale_shape_manual(values=sex_shapes) +
  coord_equal() +
  xlim(-0.09, 0.10) + ylim(-0.05, 0.05) +
  theme_classic(base_size=22) +
  theme(legend.position="none")
pca_plot

# ============================================================
# PROCRUSTES DISTANCE OF SEXUAL DIMORPHISM BY GENERATION
# ============================================================
gens       <- levels(as.factor(classifiers_all$Generation))
vectors_sex <- matrix(NA, nrow=length(gens), ncol=1,
                      dimnames=list(gens, "procdist"))

for(i in seq_along(gens)){
  idx_g  <- classifiers_all$Generation == gens[i]
  cons_F <- gpagen(pelvis_all$coords[,,idx_g & classifiers_all$sex=="F"],
                   print.progress=FALSE)$consensus
  cons_M <- gpagen(pelvis_all$coords[,,idx_g & classifiers_all$sex=="M"],
                   print.progress=FALSE)$consensus
  vectors_sex[i,1] <- procdist(cons_F, cons_M, type="full", reflect=FALSE)
}
print(vectors_sex)

# permutation test of sex dimorphism within generation
procdist_perms <- function(coords, groups, sex, nperm=1000) {
  gens    <- levels(as.factor(groups))
  results <- list()
  for(g in gens){
    idx    <- which(groups == g)
    cons_F <- gpagen(coords[,,idx[sex[idx]=="F"]], print.progress=FALSE)$consensus
    cons_M <- gpagen(coords[,,idx[sex[idx]=="M"]], print.progress=FALSE)$consensus
    obs    <- procdist(cons_F, cons_M, type="full", reflect=FALSE)
    perms  <- numeric(nperm)
    for(p in 1:nperm){
      ps     <- sample(sex[idx])
      cons_Fp <- gpagen(coords[,,idx[ps=="F"]], print.progress=FALSE)$consensus
      cons_Mp <- gpagen(coords[,,idx[ps=="M"]], print.progress=FALSE)$consensus
      perms[p] <- procdist(cons_Fp, cons_Mp, type="full", reflect=FALSE)
    }
    results[[g]] <- list(obs=obs, perm=perms)
  }
  results
}

res_list <- procdist_perms(pelvis_all$coords, classifiers_all$Generation,
                           classifiers_all$sex, nperm=500)
pvals <- sapply(res_list, function(x) mean(x$perm >= x$obs))
print(pvals)

# pairwise comparison between generations
pairwise_procdist <- function(results) {
  gens  <- names(results)
  combs <- combn(gens, 2, simplify=FALSE)
  out   <- data.frame(Gen1=character(), Gen2=character(),
                      Diff=numeric(), pvalue=numeric())
  for(pair in combs){
    g1       <- pair[1]; g2 <- pair[2]
    obs_diff  <- results[[g1]]$obs - results[[g2]]$obs
    perm_diff <- results[[g1]]$perm - results[[g2]]$perm
    pval      <- mean(abs(perm_diff) >= abs(obs_diff))
    out       <- rbind(out, data.frame(Gen1=g1, Gen2=g2,
                                       Diff=obs_diff, pvalue=pval))
  }
  out
}
pairwise_results <- pairwise_procdist(res_list)
print(pairwise_results)

# plot Procrustes distances across generations
df_sex <- data.frame(
  Generation = factor(rownames(vectors_sex),
                      levels=c("Parental","F1","F2","B1","B2","B3"),
                      labels=c("F0","F1","F2","B1","B2","B3")),
  Procdist = vectors_sex[,1]
)

sexdimo <- ggplot(df_sex, aes(x=Generation, y=Procdist, group=1, color=Generation)) +
  geom_line(linewidth=2, color="steelblue") +
  geom_point(size=6) +
  scale_color_manual(values=custom_colors) +
  theme_classic(base_size=16) +
  labs(x="Generation", y="Procrustes Distance") +
  theme(legend.position="none")
sexdimo

# ============================================================
# PCA MORPHS: PC1 and PC2
# ============================================================
pc_morph_heatmap <- function(pc_vec, sdev_val, consensus,
                             atlas_lms, atlas_mesh,
                             scale=4, pdf_path=NULL) {
  pos_morph <- consensus + pc_vec * sdev_val * scale
  neg_morph <- consensus - pc_vec * sdev_val * scale
  
  mesh_pos <- make_morph(pos_morph, atlas_lms, atlas_mesh)
  mesh_neg <- make_morph(neg_morph, atlas_lms, atlas_mesh)
  
  show_mesh(mesh_pos)
  show_mesh(mesh_neg)
  
  dist.vec <- meshDist(mesh_neg, mesh_pos, plot=FALSE)$dists
  to <- max(abs(dist.vec)); from <- -to
  open3d(zoom=0.9)
  par3d("windowRect"=c(0,100,600,800))
  meshDist(mesh_pos, distvec=dist.vec,
           rampcolors=c("blue","white","firebrick2"),
           from=from, to=to, steps=50,
           plot=TRUE, center=TRUE, symmetric=TRUE, xaxt=0, titleplot="")
  rglwidget()
  
  if(!is.null(pdf_path)){
    pdf(pdf_path, width=1.5, height=4)
    meshDist(mesh_neg, distvec=-dist.vec,
             rampcolors=c("blue","white","firebrick2"),
             from=from, to=to, steps=50,
             plot=TRUE, center=TRUE, symmetric=TRUE, xaxt=0, titleplot="")
    dev.off()
  }
}

PC1_vec <- t(matrix(pca_all$rotation[,1], nrow=3))
PC2_vec <- t(matrix(pca_all$rotation[,2], nrow=3))

pc_morph_heatmap(PC1_vec, pca_all$sdev[1], pelvis_all$consensus,
                 pelvis_atlas_lms, pelvis_atlas_mesh,
                 pdf_path="./figs/meshdist_PC1.pdf")

pc_morph_heatmap(PC2_vec, pca_all$sdev[2], pelvis_all$consensus,
                 pelvis_atlas_lms, pelvis_atlas_mesh,
                 pdf_path="./figs/meshdist_PC2.pdf")

# ============================================================
# SEX DIMORPHISM MORPHS BY GENERATION
# ============================================================

# generation index lists
gen_idx <- list(
  F0 = c(CAST, WSB, CZE, SPRET),
  F1 = c(CASTXCZE, CASTXSPRET, CASTXWSB, WSBXCZE, WSBXSPRET),
  F2 = c(CASTXCZE_F2, CASTXWSB_F2),
  B1 = c(CASTXCZE_CAST, CASTXCZE_CZE, CASTXWSB_CAST, CASTXWSB_WSB, WSBXCZE_CZE, WSBXCZE_WSB),
  B2 = c(CAXW__CAXW_CA, CAXW_CAXW_W),
  B3 = c(CA_CAXW_CA, W_CAXW_CA)
)

# B2 uses distvec_sign = -1, others = 1 (kept as in original)
heatmap_sign <- c(F0=-1, F1=1, F2=1, B1=1, B2=-1, B3=-1)

for(gen in names(gen_idx)){
  idx         <- gen_idx[[gen]]
  classif_gen <- classifiers_all[idx, ]
  gpa_gen     <- gpagen(pelvis_all$coords[,,idx], print.progress=FALSE)
  
  female_cons <- gpagen(gpa_gen$coords[,,classif_gen$sex=="F"],
                        print.progress=FALSE)$consensus
  male_cons   <- gpagen(gpa_gen$coords[,,classif_gen$sex=="M"],
                        print.progress=FALSE)$consensus
  
  cat("\n===", gen, "===\n")
  sexdim_morphs(female_cons, male_cons,
                pelvis_atlas_lms, pelvis_atlas_mesh,
                scale=0.5, heatmap_steps=30)
}

# ============================================================
# F0: PCA + PARITY MORPHS
# ============================================================
F0_idx              <- gen_idx$F0
classifiers_F0      <- classifiers_all[F0_idx, ]
F0_pelvis           <- gpagen(pelvis_all$coords[,,F0_idx], print.progress=FALSE)
F0_pelvis_f         <- gpagen(F0_pelvis$coords[,,classifiers_F0$sex=="F"],
                              print.progress=FALSE)

# regress out strain + csize for F0 PCA
regrp           <- procD.lm(F0_pelvis$coords ~ classifiers_F0$Strain + F0_pelvis$Csize,
                            SS.type="III")
summary(regrp)
F0_resid        <- residualize(regrp$residuals, F0_pelvis)
F0_resid_gpa    <- gpagen(F0_resid$coords, print.progress=FALSE)
sexdim_pca      <- gm.prcomp(F0_resid_gpa$coords)

pca_df_F0 <- data.frame(
  PC1    = sexdim_pca$x[,1],
  PC2    = sexdim_pca$x[,2],
  Sex    = factor(classifiers_F0$sex, levels=c("M","F")),
  Strain = factor(classifiers_F0$Strain)
)

ggplot(pca_df_F0, aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(group=Sex, fill=Sex), geom="polygon", alpha=0.1, color=NA) +
  stat_ellipse(aes(group=Sex, color=Sex), linetype="solid",
               linewidth=0.85, alpha=0.3, show.legend=FALSE) +
  geom_point(aes(color=Sex, shape=Strain), size=3, alpha=0.85) +
  scale_color_manual(values=c("M"="blue", "F"="hotpink")) +
  scale_fill_manual(values=c("M"="blue", "F"="hotpink")) +
  scale_shape_manual(values=c("CAST"=19, "WSB"=17, "CZECHI"=13, "SPRET"=21)) +
  coord_cartesian(xlim=c(-0.10,0.10), ylim=c(-0.05,0.05)) +
  labs(x="PC1", y="PC2", color="Sex", fill="Sex", shape="Strain") +
  theme_minimal(base_size=14) +
  theme(legend.position="right", panel.grid=element_blank(),
        axis.line=element_line())

# parity effect (F0 females only)
classifiers_F0_f <- classifiers_all[F0_idx[classifiers_F0$sex=="F"], ]

parity <- procD.lm(
  F0_pelvis_f$coords ~ classifiers_F0_f$parous.non_parous.y +
    classifiers_F0_f$Strain + F0_pelvis_f$Csize,
  SS.type="III"
)
summary(parity)

size_vector         <- matrix(parity$coefficients[2,], ncol=3, byrow=TRUE)
size_vector_centered <- scale(size_vector, center=TRUE, scale=FALSE)
size_vector_std     <- size_vector_centered / sqrt(sum(size_vector_centered^2))

nonparous_morph <- F0_pelvis_f$consensus - 0.02 * size_vector_std
parous_morph    <- F0_pelvis_f$consensus + 0.02 * size_vector_std

show_mesh(make_morph(nonparous_morph, pelvis_atlas_lms, pelvis_atlas_mesh))
show_mesh(make_morph(parous_morph,    pelvis_atlas_lms, pelvis_atlas_mesh))