L15 <- read_rds("/Users/Mathias/Documents/PhD/Xenium epilepsy/FINAL_TO_USE_EPILEPSY/L15_Seurat_obj.rds")
xenium_l15 <- L15@assays$Xenium
sum(xenium_l15@layers$counts@x)

blank_l15 <- L15@assays$BlankCodeword
sum(blank_l15@counts@x)
sum(L15@assays$ControlCodeword@counts@x)

sum(L15@assays$ControlProbe@counts@x)
options(future.globals.maxSize = 8000 * 1024 ^ 2)

cropped.coords <- Crop(L15[["L15"]], x = c(7000, 7250), y = c(6500, 7000), coords = "plot")
L15[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(L15[["zoom"]]) <- "segmentation"

ImageDimPlot(L15, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("GRIA1", "GRIA2"))

rm(c(p2,p3)

   