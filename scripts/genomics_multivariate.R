
library(FactoMineR)
library(factoextra)
library(ggforce)
library(ggrepel)

path_output="output/geno/"


# Data import ####

af <- read.table("data/data_geno/allGenotypes_AlleleFreq.csv", header = T,
                 stringsAsFactors = F) %>%
  filter_all(all_vars(. != -1)) # Removes the row if any -1


depth_files <- paste0("data/data_geno/depth/",
                      list.files("data/data_geno/depth/"))

depth <- map_df(depth_files, 
                ~ bind_cols(strain_name = gsub(pattern = "(data/data_geno/depth/)|(_compre.+$)",
                                               replacement = "", .), 
                            read.table(., header = T))) %>%
  mutate(strain_name = convnames(strain_name, eqt_convnames, "dp", "std"))


stats_map_files <- paste0("data/data_geno/stats/",
                          list.files("data/data_geno/stats/"))

stats_map <- map(stats_map_files, readLines) %>%
  map(~ grep("(reads mapped:)|(reads unmapped:)", ., value = T)) %>%
  map_df(~ setNames(., c("mapped","unmapped"))) %>%
  mutate_all(~ sub(pattern = "^[^0-9]+", replacement = "", perl = T, .) %>%
               as.numeric()) %>%
  mutate(strain_name = gsub(pattern = "(data/data_geno/stats/)|(.stat$)",
                            replacement = "", stats_map_files)) %>%
  mutate(strain_name = convnames(strain_name, eqt_convnames, "dp", "std")) %>%
  mutate(prop_mapped = mapped/(mapped + unmapped))

stats_trim_files <- paste0("data/data_geno/00_trimming/",
                           list.files("data/data_geno/00_trimming/"))

stats_trim <- map(stats_trim_files, 
                  ~ list.append(strain_name = gsub("(^(.+)/)|(\\.json)","",.),
                                rjson::fromJSON(file = .))) %>%
  map(~ list(strain_name = .x$strain_name,
           total = .x$summary$before_filtering$total_reads,
           accepted = .x$summary$after_filtering$total_reads,
           too_short = .x$filtering_result$too_short_reads,
           too_long = .x$filtering_result$too_long_reads,
           low_quality = .x$filtering_result$low_quality_reads,
           too_many_N = .x$filtering_result$too_many_N_reads)) %>%
  bind_rows() %>%
  mutate(acc_prct = accepted/total,
         strain_name = convnames(strain_name, eqt_convnames, "dp", "std"))



# Analyses ####
## ACP Alleles frequencies ####

af_mat <- t(af[, names(af)[which(grepl("^af.", names(af), perl = T))]])
rownames(af_mat) <- convnames(rownames(af_mat), data = eqt_convnames,
                              from = "af", to = "std")


af_ALL_pca <- PCA(af_mat[grep("^(LB|MB|AM|PDC|Hir|Ins|Bio)", rownames(af_mat)),], graph = F)
af_LB_pca <- PCA(af_mat[grep("^(LB)", rownames(af_mat)),], graph = F)
af_MB_pca <- PCA(af_mat[grep("^(MB|Hir)", rownames(af_mat)),], graph = F)
af_AM_pca <- PCA(af_mat[grep("^(AM|Hir)", rownames(af_mat)),], graph = F)
af_AM_MB_pca <- PCA(af_mat[grep("^(AM|MB|Hir)", rownames(af_mat)),], graph = F)
af_PDC_pca <- PCA(af_mat[grep("^(PDC|Ins)", rownames(af_mat)),], graph = F)


pca_output <- function(pca, subset) {
  
  pca_ind <- pca$ind$coord %>% 
    as.data.frame() %>% 
    rownames_to_column("strain_name") %>% 
    mutate(baker_id = strains$baker_id[match(strain_name, strains$strain_name)],
           backslopping = strains$backslopping[match(strain_name, 
                                                     strains$strain_name)],
           flour_id = strains$flour_id[match(strain_name, 
                                             strains$strain_name)],
           flour_type = flours$mill_type[match(flour_id, 
                                               flours$flour_id)])
  
  plot_d12 <- ggplot(pca_ind) +
    aes(x = Dim.1, y = Dim.2, color = backslopping, label = strain_name) +
    geom_point() +
    geom_text() +
    geom_line(aes(group = paste(flour_id, backslopping, baker_id))) +
    labs(x = paste0("Dim1 (", round(pca$eig[1,2], 1), "%)"),
         y = paste0("Dim2 (", round(pca$eig[2,2], 1), "%)")) +
    theme_bw()
  
  plot_d23 <- ggplot(pca_ind) +
    aes(x = Dim.3, y = Dim.2, color = backslopping, label = strain_name) +
    geom_point() +
    geom_text() +
    geom_line(aes(group = paste(flour_id, backslopping, baker_id))) +
    labs(x = paste0("Dim3 (", round(pca$eig[3,2], 1), "%)"),
         y = paste0("Dim2 (", round(pca$eig[2,2], 1), "%)")) +
    theme_bw()
  
  plot_d34 <- ggplot(pca_ind) +
    aes(x = Dim.3, y = Dim.4, color = backslopping, label = strain_name) +
    geom_point() +
    geom_text() +
    geom_line(aes(group = paste(flour_id, backslopping, baker_id))) +
    labs(x = paste0("Dim3 (", round(pca$eig[3,2], 1), "%)"),
         y = paste0("Dim4 (", round(pca$eig[4,2], 1), "%)")) +
    theme_bw()
  
  plot_eig <- fviz_eig(pca)
  
  pdf(file = paste0(path_output, "SNP_pca_", subset, ".pdf"), 
      width = 14.14286, height = 10)
  plot(plot_d12)
  plot(plot_d23)
  plot(plot_d34)
  plot(plot_eig)
  dev.off()
  
}

pca_output(af_ALL_pca, "ALL")
pca_output(af_MB_pca, "MB")
pca_output(af_LB_pca, "LB")
pca_output(af_AM_pca, "AM")
pca_output(af_PDC_pca, "PDC")
pca_output(af_AM_MB_pca, "AM_MB")


## Across the genome ####
### Allele frequencies & Depth ####

af_long <- af %>%
  as.data.frame() %>%
  pivot_longer(cols = starts_with("af"), names_to = "strain_name", 
               values_to = "af") %>%
  mutate(strain_name = convnames(strain_name, eqt_convnames, "af", "std"))

af_long_nest <- af_long %>% 
  mutate(chr = eqt_chrom$std[match(seq.id, eqt_chrom$vcf)]) %>%
  filter(af != 0) %>%
  filter(af != 1) %>%
  group_by(strain_name) %>%
  nest()

depth_nest <- depth %>% 
  group_by(strain_name) %>% 
  nest()

af_dp_nest <- full_join(af_long_nest, depth_nest, 
                        by = "strain_name", 
                        suffix = c("_af", "_dp")) %>% 
  rowwise() %>% 
  mutate(af_plot = list(ggplot(data_af) + 
                          aes(x = cum.loc, y = af) + 
                          geom_violin(data = data_af[!(data_af$af %in% 0:1),], 
                                      mapping = aes(group = chr), alpha = 0.3,
                                      fill = "grey" , scale = "width",
                                      color = "transparent", adjust = 0.5) +
                          geom_point(alpha = 0.7, size = 0.1, shape = ".") + 
                          geom_vline(xintercept = chrlen$cumlength,
                                     linetype = 2, size = 0.2) + 
                          # facet_zoom(xlim = range(chrlen$cumlength[chrlen$chr %in% c("mit","A","B","C")]), zoom.data = T) +
                          geom_text_repel(data = chrlen, 
                                          aes(x = cumlength - length/2, 
                                              y = 0.02, label = chr), 
                                          size = 3, 
                                          force = 0,
                                          nudge_x = c(rep(0,16), 
                                                      -5e5,-1e5,2e5,5e5),
                                          nudge_y = c(rep(0,16), 
                                                      rep(5e-2,4))) + 
                          theme_minimal()),
         
         dp_plot = list(ggplot(data_dp) +
                          aes(x = location, y = depth) +
                          geom_point(alpha = 0.5, size = 0.1, shape = ".") +
                          geom_vline(xintercept = chrlen$cumlength,
                                     linetype = 2, size = 0.2) +
                          # facet_zoom(xlim = range(chrlen$cumlength[chrlen$chr %in% c("mit","A","B","C")]), zoom.data = T) +
                          geom_text_repel(data = chrlen, 
                                          aes(x = cumlength - length/2, 
                                              y = 0.02, label = chr), 
                                          size = 3, force = 0,
                                          nudge_x = c(rep(0,16), 
                                                      -5e5,-1e5,2e5,5e5),
                                          nudge_y = c(rep(0,16), 
                                                      rep(.25,4))) +
                          scale_y_log10(limits = c(0, 300)) +
                          theme_minimal()),
         title = list(ggdraw() + 
                        draw_label(strain_name, fontface='bold',
                                   hjust = 0)),
         trim = list(ggdraw() + 
                       draw_label(paste0("Prct good reads = ",
                                         round(stats_trim$acc_prct[stats_trim$strain_name == strain_name] * 100, 1),
                                         " %."),
                                  hjust = 0)),
         mapped = list(ggdraw() + 
                         draw_label(paste0("Prct mapped reads = ",
                                           round(stats_map$prop_mapped[stats_map$strain_name == strain_name] * 100, 1),
                                           " %."),
                                    hjust = 0)), 
         plot = list(plot_grid(title, trim, mapped, af_plot, dp_plot, 
                               rel_heights = c(0.05, 0.05, 0.05, 1, 1), 
                               labels = c("", "", "","Allele frequency", "Depth"),
                               ncol = 1, align = "v"))) %>%
  select(-af_plot, -dp_plot)


for (i in 1:nrow(af_dp_nest)) {
  png(filename = paste0(path_output, "af_dp_plot_", 
                        af_dp_nest$strain_name[[i]], ".png"), 
      width = 21, height = 29.7, units = 'cm', res = 200)
  
  print(af_dp_nest$plot[[i]])
  
  dev.off()
}



### D de Tajima ####

library(PopGenome)




# SANDBOX ####
## PNGs to PDF ####

library(magick)

af_plots <- paste0("output/", list.files(path = "output/", pattern = "af_plot_(.+)png"))
image_write(image_read(af_plots), format = "pdf", "output/af_plots.pdf")

depth_plots <- paste0("output/", list.files(path = "output/", pattern = "depth_plot_(.+)png"))
image_write(image_read(depth_plots)[[1]], format = "pdf", "output/depth_plots.pdf")


## AF plot ridges Sandbox (à mettre dans un mutate précédé d'un rowwise) ####

af_plot_ridges = list(ggplot(data_af) +
                        aes(x = af, y = chr) +
                        ggridges::geom_density_ridges(scale = 7) +
                        theme_minimal())

## Allele frequencies

map2(af_long_plot$strain_name, 
     af_long_plot$plot, 
     ~ ggsave(.y + ggtitle(.x, ) + theme(title = element_text(size = 50)), 
              filename = paste0("output/af_plot_", .x, ".png"),
              width = 29.7, height = 21, bg = "white"))



## Depth ####


rowwise() %>% 
  mutate(depth_plot = list(ggplot(data) +
                             aes(x = location, y = depth) +
                             geom_point(alpha = 0.5) +
                             geom_vline(xintercept = chrlen$cumlength,
                                        linetype = 2) +
                             geom_text_repel(data = chrlen, 
                                             aes(x = cumlength - length/2, 
                                                 y = 0.02, label = chr), 
                                             size = 10, box.padding = 0, 
                                             force = 1) +
                             scale_y_log10() +
                             theme_minimal()),
         depth_plot = list(depth_plot + 
                             ggtitle(strain_name) + 
                             theme(title = element_text(size = 50))))

map(depth_plot$strain_name, 
    ~ ggsave(.x,
             filename = paste0("output/depth_plot_", 
                               .x, ".png"),
             width = 29.7, height = 21, bg = "white"))

