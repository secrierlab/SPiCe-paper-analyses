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

xenium <- read.csv("results/xenium_lung_5k_metadata_with_distances.csv")
pca_df <- read.csv("scevan_lung_5k/lung_pca_df.csv", row.names = 1)

# Subset to first 10 PCs
pca_df <- pca_df[, 1:10]
# Add cell_id back to pca_df
pca_df$cell_id <- rownames(pca_df)
merged <- merge(xenium, pca_df, by = "cell_id")
# filter for tumour cells
xenium <- merged %>%
    filter(cell_type == "Tumour cell")

xenium[1:3, ]

xenium <- xenium %>%
    mutate(
        across(
            starts_with("min"),
            ~ -((. - min(.)) / (max(.) - min(.)))
        ),
        across(
            starts_with("PC"),
            ~ (. - min(.)) / (max(.) - min(.))
        )
    )

sf_data <- st_as_sf(
    xenium,
    coords = c("x_centroid", "y_centroid"),
    crs = NA
)

predictor_var_names <- names(xenium) %>%
    .[grepl("^PC", .) | grepl("^min_dist_to", .)]

options(repr.plot.width = 12, repr.plot.height = 10)
corr_mat <- cor(xenium[, predictor_var_names])

corrplot(corr_mat, type = "upper", order = "hclust", col = brewer.pal(n = 8, name = "RdBu"))

predictor_var_names <- spatialRF::auto_cor(
    x = xenium[, predictor_var_names], cor.threshold = 0.90
) %>%
    spatialRF::auto_vif(vif.threshold = 10)

dependent_variable <- "EMT_score"
formula_str <- paste(dependent_variable, "~", paste(predictor_var_names, collapse = " + "))
formula <- as.formula(formula_str)

sf_data <- st_as_sf(
    xenium,
    coords = c("array_row", "array_col"),
    crs = NA
)

bw <- bw.gwr(formula,
    approach = "AIC",
    adaptive = TRUE,
    data = sf_data
)

m_gwr <- gwr.basic(formula,
    adaptive = TRUE,
    data = sf_data,
    bw = bw
)

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


df <- m.gwr$GW.diagnostic$edf - m.gwr$GW.diagnostic$enp

# do this for m.gwr
coefficients <- m.gwr$SDF$min_distance_to_T.cells.CD8. # Coefficients
t_values <- m.gwr$SDF$min_distance_to_T.cells.CD8._TV # Standard errors
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
