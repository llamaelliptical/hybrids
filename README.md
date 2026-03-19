# hybrids

This README file refers to the repository for the paper "Cephalo-pelvic covariation and sexual dimorphism are disrupted in hybrid mice: implications for the human obstetrical dilemma" by
XXX 


## Repository structure

data/
    raw morphometric datasets (.rds) for all pelvises and skulls (including neurocranium landmarks) used in the analyses. 

scripts/
    R scripts used for analyses and figure generation. The scripts function_heatmap_pPsup.R is used to geenrato heatmaps and morphs. The function residualize.R was written by the author.
    
atlas/
    mesh files (.ply) and landmark sets for morphing the pelvis shapes and skull shapes onto a reference mesh. The script morpho.R provides a guide to visualize the pelvis and skull landmark schemes in 3D.
