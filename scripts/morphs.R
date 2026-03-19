library(Rvcg)
library(rgl)
library(remotes)
library(morpho.tools.GM)

# ============================================================
## PELVIS ATLAS
# ============================================================
setwd("/your/working/directory/")

pelvis_atlas_tag  <- morpho.tools.GM::tag2lm("/your/working/directory/pelvis_average_lms.tag")
pelvis_atlas_lms  <- as.matrix(pelvis_atlas_tag)
pelvis_atlas_mesh <- Morpho::file2mesh("/your/working/directory/pelvis_mesh.ply")
pelvis_atlas_mesh <- vcgQEdecim(pelvis_atlas_mesh, percent=0.5)

rownames(pelvis_atlas_lms) <- 1:213

real.lms <- c(1,2,3,4,5,15,16,17,18,45,46,47,74,75,92,103,104,125,
              126,144,145,127,128,162,173,174,195,196,197,198)
semi.lms <- c(1:213)[-real.lms]

# plot pelvis atlas
open3d()
shade3d(pelvis_atlas_mesh, col="grey", specular="black", alpha=1)
spheres3d(pelvis_atlas_lms[real.lms, 1],
          pelvis_atlas_lms[real.lms, 2],
          pelvis_atlas_lms[real.lms, 3],
          col="royalblue", radius=0.1)
spheres3d(pelvis_atlas_lms[semi.lms, 1],
          pelvis_atlas_lms[semi.lms, 2],
          pelvis_atlas_lms[semi.lms, 3],
          col="gold", radius=0.1)
text3d(pelvis_atlas_lms[,1],
       pelvis_atlas_lms[,2],
       pelvis_atlas_lms[,3],
       rownames(pelvis_atlas_lms),
       col="black", cex=1.2, pos=3)
rglwidget()

# ============================================================
# SKULL ATLAS
# ============================================================
skull_mesh     <- geomorph::read.ply("/your/working/directory/Global_Adult_ONLY_Skull_Atlas_lowres.ply")
atlas_skull_lm <- read.table("/your/working/directory/Adult_Cranium_Atlas_Landmarks.tag",
                             skip=5, sep=" ", header=FALSE)[, 2:4]

LM_type_skull    <- suppressWarnings(
  read.table("/your/working/directory/Adult_Cranium_Atlas_Landmarks.tag",
             skip=5, sep=" ", header=FALSE))[, 8]

vec_LM_skull    <- which(LM_type_skull == "LANDMARK")
vec_curve_skull <- which(LM_type_skull == "curve_semilandmark")
vec_surf_skull  <- c(which(LM_type_skull == "surface_semilandmarks"),
                     which(LM_type_skull == "surface_semilandmarks;"))

# remove surface semilandmarks
atlas_skull_lm <- as.matrix(atlas_skull_lm[-vec_surf_skull, ])

skull_fixed.lm  <- vec_LM_skull
skull_curves.lm <- vec_curve_skull

# plot full skull landmark scheme
open3d()
shade3d(skull_mesh, override=F, color="grey", alpha=1, specular="black")
spheres3d(atlas_skull_lm[skull_fixed.lm, 1],
          atlas_skull_lm[skull_fixed.lm, 2],
          atlas_skull_lm[skull_fixed.lm, 3],
          col="royalblue", radius=0.3)
spheres3d(atlas_skull_lm[skull_curves.lm, 1],
          atlas_skull_lm[skull_curves.lm, 2],
          atlas_skull_lm[skull_curves.lm, 3],
          col="gold", radius=0.3)
rglwidget()

# ============================================================
# SKULL ATLAS - NEUROCRANIUM SUBSET
# ============================================================
neurocranium <- c(1:11, 21, 22, 53:58, 64, 65, 71:74, 79, 89:93, 95:115,
                  144:155, 166:170, 181:191, 204:207, 234:237,
                  254:263, 266:275, 296:336)

neurocranium_fixed  <- neurocranium[neurocranium %in% vec_LM_skull]
neurocranium_curves <- neurocranium[neurocranium %in% vec_curve_skull]

# plot neurocranium subset
open3d()
shade3d(skull_mesh, override=F, color="grey", alpha=1, specular="black")
spheres3d(atlas_skull_lm[neurocranium_fixed, 1],
          atlas_skull_lm[neurocranium_fixed, 2],
          atlas_skull_lm[neurocranium_fixed, 3],
          col="royalblue", radius=0.3)
spheres3d(atlas_skull_lm[neurocranium_curves, 1],
          atlas_skull_lm[neurocranium_curves, 2],
          atlas_skull_lm[neurocranium_curves, 3],
          col="gold", radius=0.3)
rglwidget()




