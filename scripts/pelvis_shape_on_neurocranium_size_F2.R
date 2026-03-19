# open pelvises, skulls/neurocranium and csizes
csize_pelvises <- readRDS("/your/working/directory/csize_pelvises2.rds")
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


# BOAS coordinates #####
# females non-parous F2 only ####
idx_F2_fnp <- classifiers_all$Generation == "F2" &
  classifiers_all$sex == "F" &
  classifiers_all$parous.non_parous.y == "N"

classifiers_F2_fnp <- classifiers_all[idx_F2_fnp, ]

# pelvic centroid size (F2 females non-parous)
csize_pelvis_F2_fnp <- csize_pelvises[idx_F2_fnp]
names(csize_pelvis_F2_fnp) <- classifiers_F2_fnp$Spec.ID

# GPA on BOAS coords (temporary size removal)
pelvis_F2_fnp_gpa <- gpagen(
  pelvis_all$coords[, , idx_F2_fnp]
)

match(
  dimnames(pelvis_F2_fnp_gpa$coords)[[3]],
  classifiers_F2_fnp$Spec.ID
)

# regress out strain from pelvis shape
regrp <- procD.lm(
  pelvis_F2_fnp_gpa$coords ~ classifiers_F2_fnp$Strain + classifiers_F2_fnp$weight ,
  SS.type = "III"
)
summary(regrp)

pelvis_F2_regr_fnp <- residualize(
  regrp$residuals,
  pelvis_F2_fnp_gpa
)

# re-add centroid size → BOAS residuals
pelvis_F2_regr_fnp_boas <- array(
  NA,
  dim = dim(pelvis_F2_regr_fnp$coords)
)

for (i in seq_len(dim(pelvis_F2_regr_fnp$coords)[3])) {
  pelvis_F2_regr_fnp_boas[, , i] <-
    pelvis_F2_regr_fnp$coords[, , i] * csize_pelvis_F2_fnp[i]
}

dimnames(pelvis_F2_regr_fnp_boas)[[3]] <-
  dimnames(pelvis_F2_regr_fnp$coords)[[3]]

# neurocranium size (F2 females non-parous)
csize_neuro_F2_fnp <- csize_neuro[idx_F2_fnp]

# regress out strain from neuro size
regrn <- lm(
  csize_neuro_F2_fnp ~ classifiers_F2_fnp$Strain + classifiers_F2_fnp$weight 
)
summary(regrn)

# pelvis BOAS residuals ~ neuro size
regr <- procD.lm(
  pelvis_F2_regr_fnp_boas ~ regrn$residuals
)
summary(regr)


# plot the regression for females non parous F2 
strain <- as.factor(classifiers_F2_fnp$Strain)

levels(strain) <- c(
  "CAST" = "purple",
  "CZECHI" = "lightsalmon",
  "SPRET" = "khaki4", 
  "WSB" = "darkseagreen"
)


# from package geomorph and function procD.lm (follow instructions)
plot(
  regr,
  type = "regression",
  predictor = regrn$residuals,
  reg.type = "RegScore",
  pch = 19,
  cex= 1.8,
  col = as.character(strain),
  xlab = "Neurocranial size residuals",
  ylab = "Pelvis shape regression score",
  main = "Pelvis form on neuro size (F2 females, non-parous), 
  corrected for body weight"
)

abline(regrn$residuals, scores, col="black", lwd=3)
# abline(regrn$residuals[51:72], scores[51:72], col="darkseagreen", lwd=3)
# abline(regrn$residuals[16:31], scores[16:31], col="lightsalmon", lwd=3)
# abline(regrn$residuals[32:50], scores[32:50], col="khaki4", lwd=3)

# Project pelvic shapes onto size_vector to get scores
scores <- two.d.array(pelvis_F2_regr_fnp_boas) %*% as.vector(size_vector)

# Recompute size_vector as you already did
size_vector <- matrix(regr$coefficients[2,], ncol=3, byrow=TRUE)

# Project pelvic shapes onto size_vector to get scores
scores <- two.d.array(pelvis_F2_regr_fnp_boas) %*% as.vector(size_vector)

#calculate std dev of vector
sdev= sd(scores[,1])

# morphs
small_head<- pelvis_all$consensus - 0.2*sdev*size_vector
big_head<- pelvis_all$consensus + 0.2*sdev*size_vector

# morph
# Center the two matrices at the origin 
sup.pos <- pPsup(pelvis_atlas_lms, small_head)
# Map the mean face onto the target
sup.ilium <- tps3d(pelvis_atlas_mesh, pelvis_atlas_lms, sup.pos$Mp2) 

sup2 <- pPsup(pelvis_atlas_lms, big_head)
# Map the mean face onto the target
sup2.ilium <- tps3d(pelvis_atlas_mesh, pelvis_atlas_lms, sup2$Mp2) 


open3d(zoom = 0.6)
shade3d(sup.ilium, override=F, color="grey", alpha=1, specular=1)
rglwidget()

open3d(zoom=0.6)
shade3d(sup2.ilium, override=F, color="grey", alpha=1, specular=1)
rglwidget()


# plot heatmap with only extreme values 
dist.vec <- meshDist(sup2.ilium, sup.ilium, plot = FALSE)$dists
to   <- max(abs(dist.vec))
from <- -to
open3d(zoom = 0.9)
par3d("windowRect"= c(0, 100, 600, 800))
meshDist(sup.ilium, distvec = dist.vec, rampcolors = c("blue", "white", "firebrick2"), from = from,
         to = to, steps = 50, plot = TRUE, center = TRUE, symmetric = TRUE, xaxt = 0,          # <-- suppress axis ticks
         titleplot = ""
)
legend3d( "right",  legend = c( sprintf("%.2f", from), "0.00", sprintf("%.2f", to) ),
          col = c("blue", "white", "firebrick2"), pch = 15, cex = 3
)
rglwidget()


regrsize<- lm (csize_pelvis_F2_fnp ~ classifiers_F2_fnp$weight)
summary(regrsize)

regrsize<- lm (csize_neuro_F2_fnp ~ classifiers_F2_fnp$weight)
summary(regrsize)












