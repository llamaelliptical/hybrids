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


library(dplyr)
library(ggplot2)
library(ggrepel)


# ------------------------------------------------------------
# Add centroid sizes to classifiers_all
# ------------------------------------------------------------

# pelvis centroid size
classifiers_all$csize_pelvis <- csize_pelvis[
  match(classifiers_all$Spec.ID, names(csize_pelvis))
]

# neurocranium centroid size
classifiers_all$csize_neuro <- csize_neuro[
  match(classifiers_all$Spec.ID, names(csize_neuro))
]
# -------------------------------------------------------------------
# Split datasets
# -------------------------------------------------------------------

classifiers_all_f <- classifiers_all %>% filter(sex == "F")
classifiers_all_m <- classifiers_all %>% filter(sex == "M")

# -------------------------------------------------------------------
# Means by generation
# -------------------------------------------------------------------

csize_bygen <- classifiers_all %>%
  group_by(Generation) %>%
  summarise(csize = mean(csize_pelvis))

csizeneuro_bygen <- classifiers_all %>%
  group_by(Generation) %>%
  summarise(csize = mean(csize_neuro))

csize_bygen_F <- classifiers_all_f %>%
  group_by(Generation) %>%
  summarise(csize = mean(csize_pelvis))

csizeneuro_bygen_F <- classifiers_all_f %>%
  group_by(Generation) %>%
  summarise(csize = mean(csize_neuro))

csize_bygen_m <- classifiers_all_m %>%
  group_by(Generation) %>%
  summarise(csize = mean(csize_pelvis))

csizeneuro_bygen_m <- classifiers_all_m %>%
  group_by(Generation) %>%
  summarise(csize = mean(csize_neuro))

# -------------------------------------------------------------------
# Means by strain
# -------------------------------------------------------------------

csize_bystr <- classifiers_all %>%
  group_by(Strain) %>%
  summarise(csize = mean(csize_pelvis))

csizeneuro_bystr <- classifiers_all %>%
  group_by(Strain) %>%
  summarise(csize = mean(csize_neuro))

csize_bystr_F <- classifiers_all_f %>%
  group_by(Strain) %>%
  summarise(csize = mean(csize_pelvis))

csizeneuro_bystr_F <- classifiers_all_f %>%
  group_by(Strain) %>%
  summarise(csize = mean(csize_neuro))

csize_bystr_m <- classifiers_all_m %>%
  group_by(Strain) %>%
  summarise(csize = mean(csize_pelvis))

csizeneuro_bystr_m <- classifiers_all_m %>%
  group_by(Strain) %>%
  summarise(csize = mean(csize_neuro))

# -------------------------------------------------------------------
# Strain numbering
# -------------------------------------------------------------------

df <- data.frame(
  Strain = c("CAST","CZECHI","SPRET","WSB",
             "WSBXCZE","WSBXSPRET","CASTXCZE","CASTXWSB","CASTXSPRET",
             "CASTXCZE_F2","CASTXWSB_F2",
             "WSBXCZE_CZE","WSBXCZE_WSB",
             "CASTXCZE_CAST","CASTXCZE_CZE","CASTXWSB_CAST","CASTXWSB_WSB",
             "CAXW_CAXW_W","CAXW__CAXW_CA",
             "CA_CAXW_CA","W_CAXW_CA"),
  Number = 1:21
)

# -------------------------------------------------------------------
# Colors
# -------------------------------------------------------------------

custom_colors <- c(
  "Parental"="blue",
  "F0"="blue",
  "F1"="gold1",
  "F2"="red",
  "B1"="olivedrab3",
  "B2"="purple",
  "B3"="darkgreen"
)

# -------------------------------------------------------------------
# Helper function to build strain plotting data
# -------------------------------------------------------------------

build_strain_df <- function(csize_bystr, csize_bygen, classifiers){
  
  f0 <- csize_bygen$csize[csize_bygen$Generation=="Parental"]
  
  strain_order <- unique(classifiers$Strain)
  
  csize_bystr <- csize_bystr %>%
    slice(match(strain_order, Strain))
  
  strain_size <- data.frame(
    Strain = csize_bystr$Strain,
    Csize_plot = csize_bystr$csize / f0 * 100,
    numbers = df$Number[match(csize_bystr$Strain, df$Strain)],
    Type = "Strain"
  ) %>%
    mutate(
      Strain = classifiers$Strain[match(Strain, classifiers$Strain)],
      Generation = classifiers$Generation[match(Strain, classifiers$Strain)]
    ) %>%
    mutate(
      Generation = factor(
        Generation,
        levels=c("Parental","F1","F2","B1","B2","B3"),
        labels=c("F0","F1","F2","B1","B2","B3")
      )
    )
  
  return(strain_size)
}

# -------------------------------------------------------------------
# Build plotting datasets
# -------------------------------------------------------------------

strain_all_pel <- build_strain_df(csize_bystr, csize_bygen, classifiers_all)
strain_all_sku <- build_strain_df(csizeneuro_bystr, csizeneuro_bygen, classifiers_all)

strain_fem_pel <- build_strain_df(csize_bystr_F, csize_bygen_F, classifiers_all_f)
strain_fem_sku <- build_strain_df(csizeneuro_bystr_F, csizeneuro_bygen_F, classifiers_all_f)

strain_mas_pel <- build_strain_df(csize_bystr_m, csize_bygen_m, classifiers_all_m)
strain_mas_sku <- build_strain_df(csizeneuro_bystr_m, csizeneuro_bygen_m, classifiers_all_m)

# -------------------------------------------------------------------
# Plot function
# -------------------------------------------------------------------

make_plot <- function(data){
  
  ggplot(data, aes(x=Generation,y=Csize_plot,color=Generation)) +
    geom_boxplot(
      aes(group=Generation, fill=Generation),
      width=0.5,
      outlier.shape=NA,
      alpha=0.15
    ) +
    geom_point(
      aes(shape=Type),
      size=3,
      alpha=0.5,
      position=position_jitter(width=0.15)
    ) +
    geom_text_repel(
      aes(label=numbers),
      position=position_jitter(width=0.15),
      size=4,
      segment.color="grey50",
      show.legend=FALSE
    ) +
    scale_color_manual(values=custom_colors) +
    scale_fill_manual(values=custom_colors) +
    scale_shape_manual(values=c("Strain"=16)) +
    ylim(93,112) +
    theme_classic() +
    theme(
      axis.text.x=element_text(size=0),
      axis.text.y=element_text(size=9),
      legend.position="none",
      text=element_text(size=0)
    )
}

# -------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------

all_pel <- make_plot(strain_all_pel)
all_sku <- make_plot(strain_all_sku)

fem_pel <- make_plot(strain_fem_pel)
fem_sku <- make_plot(strain_fem_sku)

mas_pel <- make_plot(strain_mas_pel)
mas_sku <- make_plot(strain_mas_sku)

# -------------------------------------------------------------------
# Statistical tests
# -------------------------------------------------------------------

all_pvalue <- pairwise.t.test(
  classifiers_all$csize_pelvis,
  classifiers_all$Generation,
  p.adjust.method="bonferroni"
)

F_pvalue <- pairwise.t.test(
  classifiers_all_f$csize_pelvis,
  classifiers_all_f$Generation,
  p.adjust.method="bonferroni"
)

M_pvalue <- pairwise.t.test(
  classifiers_all_m$csize_pelvis,
  classifiers_all_m$Generation,
  p.adjust.method="bonferroni"
)

all_pvalue_neuro <- pairwise.t.test(
  classifiers_all$csize_neuro,
  classifiers_all$Generation,
  p.adjust.method="bonferroni"
)

F_pvalue_neuro <- pairwise.t.test(
  classifiers_all_f$csize_neuro,
  classifiers_all_f$Generation,
  p.adjust.method="bonferroni"
)

M_pvalue_neuro <- pairwise.t.test(
  classifiers_all_m$csize_neuro,
  classifiers_all_m$Generation,
  p.adjust.method="bonferroni"
)


