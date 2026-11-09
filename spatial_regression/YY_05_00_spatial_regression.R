library(Rcpp)
library(RcppEigen)
library(sf)
library(spatstat)
library(tidyverse)
library(GWmodel)
library(geojsonio)
library(spatialRF)
library(dplyr)
library(corrplot)
library(RColorBrewer)

setwd("/home/cenk/Withnell_spice/")

df <- read.csv("results/xenium_breast_5k_metadata_with_distances.csv") %>%
    select(cell_type, cell_id, "array_row", "array_col", starts_with("min_dist"), "EMT_score")

pca <- read.csv("scevan_breast_5k/breast_pca_df.csv", row.names = 1)

# Subset to first 34 PCs
pca <- pca[, 1:34]
# Add cell_id back to pca_df
pca$cell_id <- rownames(pca)

merged <- merge(df, pca, by = "cell_id") %>% filter(cell_type == "Tumour cell")
x_pca_df <- merged %>% select(-cell_id)

min_max_scale <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

x_pca_df <- x_pca_df %>% mutate(across(starts_with("min"), ~ -min_max_scale(.x)))

predictor_var_names <- names(x_pca_df) %>%
    .[grepl("^PC", .) | grepl("^min_dist_to", .)]

options(repr.plot.width = 12, repr.plot.height = 10)

corr_mat <- cor(x_pca_df[, predictor_var_names])

corrplot(corr_mat, type = "upper", order = "hclust", col = brewer.pal(n = 8, name = "RdBu"))

formula_str <- paste("EMT_score", "~", paste(predictor_var_names, collapse = " + "))
formula <- as.formula(formula_str)

sf_data <- st_as_sf(x_pca_df, coords = c("array_col", "array_row"), crs = NA)
bw <- bw.gwr(formula, approach = "AIC", adaptive = TRUE, data = sf_data)
save(bw, file = "results/breast_spatial_regression_bw.RData")
load("results/breast_spatial_regression_bw.RData")
m_gwr <- gwr.basic(formula, adaptive = TRUE, data = sf_data, bw = bw)

save(m_gwr, file = "results/breast_spatial_regression_gwr.RData")

# predictor_var_names <- spatialRF::auto_cor(
#     x = x_pca_df[, predictor_var_names], cor.threshold = 0.90
# ) %>%
#     spatialRF::auto_vif(vif.threshold = 10)

predictor_var_names <- colnames(x_pca_df)[grepl("min_", colnames(x_pca_df))]
# list of bws0 based on the number of variables plus 1
bws0_list <- rep(100, length(predictor_var_names) + 1)

gw_ms <- gwr.multiscale(formula,
    data = sf_data,
    adaptive = TRUE,
    max.iterations = 1000,
    criterion = "CVR",
    kernel = "bisquare",
    bws0 = bws0_list,
    verbose = TRUE,
    predictor.centered = rep(TRUE, 4)
)

# Extract the coefficient and SE columns for 'min_distance_to_T.cells.CD8.'
coefficients <- gw_ms$SDF$min_distance_to_T.cells.CD8. # Coefficients
t_values <- gw_ms$SDF$min_distance_to_T.cells.CD8._TV # Standard errors
# print size of coefficients
length(coefficients)
length(t_values)
df <- gw_ms$GW.diagnostic$edf - gw_ms$GW.diagnostic$enp
p_values_gw_ms <- 2 * (1 - pt(abs(t_values), df = df))

# Coefficient, standard error, t-value, and p-value for 'min_distance_to_'
result <- data.frame(
    Coefficient = coefficients,
    TValue = t_values,
    PValue = p_values_gw_ms
)
# filter result for p-values < 0.05
result_filt <- result %>% filter(PValue < 0.05)
# print shape of result_filt
dim(result_filt)
write.csv(result, "results/msgwr_tcell4_coef.csv")


df <- m_gwr$GW.diagnostic$edf - m_gwr$GW.diagnostic$enp

# do this for m_gwr
coefficients <- m_gwr$SDF$min_dist_to_Activated.CD8..T.cell # Coefficients
t_values <- m_gwr$SDF$min_dist_to_Activated.CD8..T.cell._TV # Standard errors
p_values <- 2 * (1 - pt(abs(t_values), df = df))
result <- data.frame(
    Location = 1:length(coefficients),
    Coefficient = coefficients,
    TValue = t_values,
    PValue = p_values
)
# save
# filter
result_filt <- result %>% filter(PValue < 0.05)
# print shape of result_filt
dim(result_filt)

write.csv(result, "results/gwr_tcell4_coef.csv")
