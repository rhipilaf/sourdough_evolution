
library(FactoMineR)
library(factoextra)
library(ggforce)
library(ggrepel)
library(cividis)

path_output="output/geno/"


# Data import ####

## Allelic freq data ####

if(file.exists("data/data_geno/af.rds")) {
  
  af <- readRDS("data/data_geno/af.rds")
  af_names <- names(af)
  
} else {
  
  af_names <- flatten_chr(read.table("data/data_geno/af_allGenotypes_AlleleFreq_clean.csv", 
                                     nrows = 1, stringsAsFactors = F, sep = " "))
  
  af <- read.table("data/data_geno/af_allGenotypes_AlleleFreq_clean.csv", header = T,
                   stringsAsFactors = F, sep = ",") %>%
    filter_all(all_vars(. != -1)) # Removes the row if any -1
  
  saveRDS(af, file = "data/data_geno/af.rds")
  
}


## Depth data ####

if(file.exists("data/data_geno/06_depth/dp.rds")) {
  
  dp <- readRDS("data/data_geno/dp.rds")
  
} else {
  
  depth_files <- paste0("data/data_geno/depth/",
                        list.files("data/data_geno/depth/", 
                                   pattern = "compressed"))
  
  dp <- map_df(depth_files, 
               ~ bind_cols(strain_name = gsub(pattern = "(data/data_geno/depth/)|(_compressed.+$)",
                                              replacement = "", .), 
                           read.table(., header = T, sep = " "))) %>%
    mutate(strain_name = convnames(strain_name, eqt_convnames, "sample", "strain_name")) %>%
    filter(!is.na(strain_name))
  
  saveRDS(dp, file = "data/data_geno/dp.rds")
  
}


## Stats about reads mapping ####

stats_map_files <- paste0("data/data_geno/stats/",
                          list.files("data/data_geno/stats/"))

stats_map <- map(stats_map_files, readLines) %>%
  map(~ grep("(reads mapped:)|(reads unmapped:)", ., value = T)) %>%
  map_df(~ setNames(., c("mapped","unmapped"))) %>%
  mutate_all(~ sub(pattern = "^[^0-9]+", replacement = "", perl = T, .) %>%
               as.numeric()) %>%
  mutate(strain_name = gsub(pattern = "(data/data_geno/stats/)|(.stat$)",
                            replacement = "", stats_map_files)) %>%
  mutate(strain_name = convnames(strain_name, eqt_convnames, "sample", "strain_name")) %>%
  mutate(prop_mapped = mapped/(mapped + unmapped))


## Stats about reads trimming ####

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
         strain_name = convnames(strain_name, eqt_convnames, "sample", "strain_name"))


# Analyses ####
## ACP Alleles frequencies ####

af_mat <- t(af[, names(af)[which(grepl("^af.", names(af), perl = T))]])
# rownames(af_mat) <- convnames(rownames(af_mat), data = eqt_convnames, from = "af", to = "std")

af_TOUT_pca <- PCA(af_mat, graph = F)
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

pca_output(af_TOUT_pca, "TOUT")
pca_output(af_ALL_pca, "ALL")
pca_output(af_MB_pca, "MB")
pca_output(af_LB_pca, "LB")
pca_output(af_AM_pca, "AM")
pca_output(af_PDC_pca, "PDC")
pca_output(af_AM_MB_pca, "AM_MB")

fviz_pca_ind(af_TOUT_pca)
ggsave(filename = "SNP_pca_TOUT.pdf", 
       width = 8, height = 8)


## Across the genome ####

af_long <- af %>%
  pivot_longer(cols = starts_with("af"), 
               names_to = "strain_name", values_to = "af") %>%
  # melt(., measure.vars = names(.)[grepl("^af", names(.))], 
  # value.name = "af", variable.name = "strain_name") %>%
  mutate(strain_name = convnames(strain_name, eqt_convnames, "af", "strain_name"),
         chr = eqt_chrom$std[match(seq.id, eqt_chrom$vcf)])

af_long_nest_strain <- af_long %>% 
  # filter(af != 0) %>%
  # filter(af != 1) %>%
  group_by(strain_name) %>%
  nest()

af_long_nest_chr <- af_long %>% 
  # filter(af != 0) %>%
  # filter(af != 1) %>%
  group_by(strain_name, chr) %>%
  nest()


### Ploidy estimation from allelic frequencies ####

af_strain_dens <- af_long %>%
  filter(af > 0, af < 1) %>%
  group_by(strain_name) %>%
  summarise(dens = list(hist(af, breaks = seq(0,1,0.01), plot = F))) %>%
  mutate(af = map(dens, ~ .x$mids),
         dens = map(dens, ~ .x$density)) %>%
  unnest(c(af, dens)) %>%
  mutate(dens = ifelse(dens > 8, 8, dens))

## Ordering labels with the first PCA axis (works quite well)
af_strain_dens_order <- af_strain_dens %>% 
  pivot_wider(names_from = "af", values_from = "dens") %>%
  column_to_rownames("strain_name") %>%
  as.matrix() %>%
  prcomp(., scale. = T)
af_strain_dens_order <- af_strain_dens_order$x[,"PC2"]

## Ordering labels with hierarchical clustering
mat <- af_strain_dens %>%
  pivot_wider(names_from = "af", values_from = "dens") %>% 
  column_to_rownames("strain_name") %>%
  as.matrix()
nmd <- dist(mat, method = "euclidean") %>%
  ecodist::nmds(., mindim=1, maxdim=1)
tmp <- setNames(as.vector(nmd$conf[[4]]), row.names(mat))

## Ordering labels according to the variance (works less well)
af_strain_dens_order <- af_long %>% 
  group_by(strain_name) %>%
  summarise(var = var(af)) %>%
  column_to_rownames("strain_name") %>%
  as.matrix()
af_strain_dens_order <- setNames(af_strain_dens_order[,1], 
                                 rownames(af_strain_dens_order))

af_strain_dens_plot_evoexp <- af_strain_dens %>%
  filter(strain_name %in% strains_scer_all$strain_name[strains_scer_all$evoexp]) %>%
  arrange(strain_name) %>%
  ggplot() +
  aes(x = strain_name, y = af) +
  geom_tile(aes(fill = dens), color = 'black', show.legend = F) +
  scale_fill_cividis() +
  labs(x = "", y = "Fréquence allélique", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust = 1, size = 7, 
                                   color = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        text = element_text(family = "serif"))


ladder <- c(N1 = "CBS_1387", N2 = "614", N3 = "DBVPG_6861", N4 = "A-18", N5 = "CH13")
af_strain_dens_plot_ladder <- af_strain_dens %>%
  filter(strain_name %in% ladder) %>%
  mutate(strain_name = names(ladder)[match(strain_name, ladder)],
         strain_name = gsub("N","N = ", strain_name)) %>%
  ggplot() +
  aes(x = strain_name, 
      y = af) +
  geom_tile(aes(fill = dens), color = "black") +
  scale_fill_cividis() +
  labs(x = "ECHELLE", y = "", fill = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust = 1, size = 9,
                                   color = "black"),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none",
        text = element_text(family = "serif"),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, 
                                                    b = 0, l = 0)))

af_strain_dens_plot <- plot_grid(af_strain_dens_plot_evoexp, 
                                 af_strain_dens_plot_ladder,
                                 align = "h", nrow = 1, 
                                 axis = "rlbt", 
                                 rel_widths = c(1,0.3))


myggsave(filename = "output/geno/af_strain_dens_plot_all",
         width = 10, height = 4, dpi = 600, device = "pdf")
myggsave(filename = "output/geno/af_strain_dens_plot_all",
         width = 6, height = 3, dpi = 600, device = "png")
myggsave(filename = "output/geno/af_strain_dens_plot_all",
         width = 12, height = 6, dpi = 600, device = "svg")



af_strain_modes <- af_long_nest_strain %>%
  mutate(modes = lapply(data, function(x) get_modes(x$af, 
                                                    adjust = 1, 
                                                    signifi = 2, 0, 1))) %>%
  select(strain_name, modes) %>%
  unnest(modes) %>%
  pivot_longer(cols = starts_with("m"), names_to = "mode") %>%
  filter(!is.na(value))

af_strain_modes_nb <- af_strain_modes %>%
  group_by(strain_name, chr) %>%
  summarise(nb_modes = n())


af_chr_modes <- af_long_nest_chr %>%
  mutate(modes = lapply(data, function(x) get_modes(x$af, 
                                                    adjust = 1, 
                                                    signifi = 2, 0, 1))) %>%
  select(strain_name, chr, modes) %>%
  unnest(modes) %>%
  pivot_longer(cols = starts_with("m"), names_to = "mode") %>%
  filter(!is.na(value)) %>%
  mutate(value = as.numeric(value))

af_chr_modes_nb <- af_chr_modes %>%
  group_by(strain_name, chr) %>%
  summarise(nb_modes = n())

density(x = af_chr_modes$value, adjust = 0.1, na.rm = T) %>% plot()

### Allele frequencies & Depth ####

dp_nest <- dp %>% 
  # filter(depth > 0) %>%
  group_by(strain_name) %>% 
  nest()

af_dp_nest <- left_join(af_long_nest_strain, dp_nest, 
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
                          scale_y_log10(limits = c(1, 1000)) +
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
  cat(af_dp_nest$strain_name[[i]], "\n")
  png(filename = paste0(path_output, "af_dp_plot_", 
                        make.filenames(af_dp_nest$strain_name[[i]]), ".png"), 
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

dp_plots <- paste0("output/", list.files(path = "output/", pattern = "dp_plot_(.+)png"))
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


# Population structure

library(PopGenome)

strns <- c("AM-GT5-T5-8.fasta", "AM-maison-T1-22.fasta", "AM-maison-T5-2.fasta", 
           "AM-PG9-T5-6.fasta", "AM-PJ8-T5-1.fasta", "AM-PR12-T5-28.fasta", 
           "bioreal.fasta", "Hirondelle.fasta", "instant.fasta", "LB-GP4-T1-3.fasta", 
           "LB-GP4-T32-1.fasta", "LB-MAISON-T1-1.fasta", "LB-MAISON-T32-1.fasta", 
           "LB-PD7-T1-1.fasta", "LB-PD7-T32-1.fasta", "LB-RG15-T1-1.fasta", 
           "LB-RG15-T32-1.fasta", "LB-RT17-T1-1.fasta", "LB-RT17-T32-1.fasta", 
           "LB-SP22-T32-1.fasta", "LB-SP22-T5-2.fasta", "LB-SR24-T1-9.fasta", 
           "LB-SR24-T32-1.fasta", "LB-Temoin-T1-1.fasta", "LB-Temoin-T32-1.fasta", 
           "MB-GD1-T32-9.fasta", "MB-GR6-T32-10.fasta", "MB-mais-T32-11.fasta", 
           "MB-maison-T1-2.fasta", "MB-PJ8-T32-12.fasta", "MB-RD13-T32-13.fasta", 
           "MB-RR18-T32-14.fasta", "MB-ST23-T1-1.fasta", "MB-ST23-T32-15.fasta", 
           "PDC-mais-T32-9.fasta", "PDC-SG21-T1-4.fasta", "SRR5678552.fasta", 
           "SRR5678580.fasta")

T1 <- strns[grepl("T1", strns)]
T5 <- strns[grepl("T5", strns)]
T32 <- strns[grepl("T32", strns)]

genomes <- readData(path = "data/data_geno/af_LaurianeGenotypes_fastas", 
                    populations = list(T1,T5,T32))

genomes <- F_ST.stats(genomes)

genomes <- diversity.stats(genomes)

