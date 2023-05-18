#This script is to calculate the trajectories for a single source of data ('origin'). What we do is read in the data, subset to the origin of choice (T47D or tumor), and then identify trajectories separately in each of the three treatment conditions. You need to check the lineages for each of the subsets, and confirm that it is properly going from G0 to G2/M cells. If more than one trajectory is identified, pick the lineage culminating in G2/M.
#setup - make big SCE for coordinates
library(Seurat)
library(grid)
library(slingshot)
library(dplyr)
library(scater)

ti <- function(df, pathname) {
        
        #first we make a SCE object of all the data in order to plot the curves atop the PHATE
        big = subset(df, df$Origin==orig)
        cols = c('pRB_over_RB', 'Ki67', 'pRB', 'RB', 'CDK2', 'CDK4', 'cycD1', 'cycE', 'Cdt1', 'E2F1', 'DNA', 'cycA', 'cycB1', 'p21')
        raw_df_big = subset(big, select=cols)
        
        phate_graph = subset(big, select=c(PHATE_1, PHATE_2))
        sce_big <- SingleCellExperiment(t(raw_df_big))
        colData(sce_big)$phase = big$phase
        
        reducedDims(sce_big) <- list(PHATE=as.matrix(phate_graph))
        #we dont use the trajectories found here, this is just for plotting
        
        #this function finds trajectories 
        trj <- function(data, sce_want) {
                #subset to the features, not using the metadata
                raw_df = subset(data, select=cols)
                #phate coordinates - we use first two dimensions
                phate = subset(data, select=c(PHATE_1, PHATE_2))
                #create SCE object
                sce <- SingleCellExperiment(assays = as.matrix(t(raw_df)))
                colData(sce)$phase = data$phase
                
                reducedDims(sce) <- list(PHATE=as.matrix(phate))
                #run Slingshot, only giving starting cluster of G0
                sce <- slingshot(
                        sce,
                        reducedDim = 'PHATE',
                        clusterLabels = 'phase',
                        start.clus = 'G0', allow.breaks=F
                )
                slo <- SlingshotDataSet(sce)
                if (sce_want==T) {return(sce)} #if you want SCE object or Slingshot object to be returned
                else {return(slo)} 
        }
        
        #first we just run to get the curve coordinates, so we want a Slingshot object
        untrt = subset(df, df$Origin==orig & df$well==0)
        slo_untrt = trj(untrt, sce_want=F)
        
        ten = subset(df, df$Origin==orig & df$well==10)
        slo_ten = trj(ten, sce_want=F)
        
        hundo = subset(df, df$Origin==orig & df$well==100)
        slo_hundo = trj(hundo, sce_want=F)
        
        
        ##plot lineages together - untrt as dashed line, trt as solid
        curve_untrt = slingCurves(slo_untrt, as.df = TRUE) %>%
                dplyr::rename("X" = "PHATE_1", "Y" = "PHATE_2")
        curve_ten = slingCurves(slo_ten, as.df = TRUE) %>%
                dplyr::rename("X" = "PHATE_1", "Y" = "PHATE_2")
        curve_hundo = slingCurves(slo_hundo, as.df = TRUE) %>%
                dplyr::rename("X" = "PHATE_1", "Y" = "PHATE_2")
        
        #subset to lineages we want - look at them to confirm
        slo_untrt@lineages
        slo_ten@lineages
        slo_hundo@lineages
        
        curve_untrt = subset(curve_untrt, curve_untrt$Lineage == 1)
        curve_ten = subset(curve_ten, curve_ten$Lineage == 1)
        curve_hundo = subset(curve_hundo, curve_hundo$Lineage == 1)
        
        #plot the trajectories atop the phate coordinates, colored by cell cycle phase
        #pdf(paste0(pathname,'_phate_trajectory.pdf'), width = 10, height = 10)
        plotReducedDim(sce_big, dimred = "PHATE", colour_by = 'phase', point_alpha=.7) +
                geom_path(data = curve_untrt %>% arrange(Order),
                          aes(group = Lineage), linewidth = 3, col='black') +
                geom_path(data = curve_ten %>% arrange(Order),
                          aes(group = Lineage), linewidth = 3, col='#A554F5') +
                geom_path(data = curve_hundo %>% arrange(Order),
                          aes(group = Lineage), linewidth = 3, col='#E3822D') +
                theme(axis.text.x=element_blank(), #remove x axis labels
                      axis.ticks.x=element_blank(), #remove x axis ticks
                      axis.text.y=element_blank(),  #remove y axis labels
                      axis.ticks.y=element_blank(),  #remove y axis ticks
                      axis.line.y=element_blank(),
                      axis.line.x=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank()
                ) + theme(legend.position="none") + theme(plot.title = element_text(size = 40)) +
                theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values = c("G0" = '#5CAD92', 'G1' = '#594997', 'G2M'='#E7739A', 'S'='#0099CC'))
        ggsave(paste0(pathname,'_phate_trajectory.pdf'), width = 10, height = 10)
        #dev.off()
        
        #code to plot one feature with trajectories from each treatment condition
        graphone <- function(dt, var) {
                thm = theme(
                        axis.title.x = element_blank(),
                        axis.line = element_line(colour = "grey"),
                        panel.background = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.text.y = element_text(size = 18),
                        axis.title.y = element_text(size = 22),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank())
                g <- (ggplot(subset(dt, dt$well==0), aes(x = pseudo_align, y = get(var))) +
                              geom_smooth(subset(dt, dt$well==0),mapping=aes(x=pseudo_align, y = get(var)), linetype = 'solid', color='black', se = F) +
                              geom_smooth(subset(dt, dt$well==10),mapping=aes(x=pseudo_align, y = get(var)), linetype = 'solid', color= '#A554F5',se = F) +
                              geom_smooth(subset(dt, dt$well==100),mapping=aes(x=pseudo_align, y = get(var)), linetype = 'solid', color='#E3822D', se = F) +
                              theme(legend.position="none") + thm + xlim(0, max(join$pseudo_align)))
                
        }
        
        
        source("TrAGEDy_functions_tmz.R")
        untrt = subset(df, df$Origin==orig & df$well==0)
        slo_untrt2 = trj(untrt, sce_want=T)
        untrt = subset(untrt, !is.na(slo_untrt2$slingPseudotime_1))
        slo_untrt2 = slo_untrt2[,!is.na(slo_untrt2$slingPseudotime_1)]
        slo_untrt2@colData@rownames = untrt$X
        
        ten = subset(df, df$Origin==orig & df$well==10)
        slo_ten2 = trj(ten, sce_want=T)
        ten = subset(ten, !is.na(slo_ten2$slingPseudotime_1))
        slo_ten2 = slo_ten2[,!is.na(slo_ten2$slingPseudotime_1)]
        slo_ten2@colData@rownames = ten$X
        
        hundo = subset(df, df$Origin==orig & df$well==100)
        slo_hundo2 = trj(hundo, sce_want=T)
        hundo = subset(hundo, !is.na(slo_hundo2$slingPseudotime_1))
        slo_hundo2 = slo_hundo2[,!is.na(slo_hundo2$slingPseudotime_1)]
        slo_hundo2@colData@rownames = hundo$X 
        
        
        all = rbind(untrt, ten, hundo)
        
        ten_hundo = rbind(ten, hundo)
        slo_ten_hundo = trj(ten_hundo, sce_want=T)
        slo_ten_hundo@colData@rownames = ten_hundo$X
        #first we align the treatment conditions to each other
        pdf(paste0(pathname,'_trt.pdf'))
        slo_ten_hundo$slingPseudotime_1 = align(slo_ten2, slo_hundo2, paste0(pathname, '_alignment_trt.pdf'))
        dev.off()
        #then we align the untreated to the treatment conditions
        all$pseudo_align = align(slo_untrt2, slo_ten_hundo, paste0(pathname, '_alignment_all.pdf'))
        all = subset(all, select=-c(PHATE_1, PHATE_2, prb_ratio, 
                                    Origin, RB_status))
        write.csv(all, paste0(pathname,'_aligned.csv'))
        
        join <- all
        join$dummy_g0 = 0
        join$dummy_g1 = 0
        join$dummy_s = 0
        join$dummy_g2 = 0
        join$dummy_g0 = ifelse(join$phase=='G0', 1, join$dummy_g0)
        join$dummy_g1 = ifelse(join$phase=='G1', 1, join$dummy_g1)
        join$dummy_s = ifelse(join$phase=='S', 1, join$dummy_s)
        join$dummy_g2 = ifelse(join$phase=='G2M', 1, join$dummy_g2)
        
        g_cdk2 <- ggplotGrob(graphone(join, 'CDK2') + ylab('CDK2')  )
        
        g_cdk4 <- ggplotGrob(graphone(join, 'CDK4') + ylab('CDK4')  )
        
        g_cdt1 <- ggplotGrob(graphone(join, 'Cdt1') + ylab('Cdt1')  )
        
        g_dna <- ggplotGrob(graphone(join, 'DNA') + ylab('DNA')  )
        
        g_e2f1 <- ggplotGrob(graphone(join, 'E2F1') + ylab('E2F1')  )
        
        g_ki67 <- ggplotGrob(graphone(join, 'Ki67') + ylab('Ki67')  )
        
        g_rb <-ggplotGrob(graphone(join, 'RB') + ylab('RB')  )
        
        g_prb <- ggplotGrob(graphone(join, 'pRB') + ylab('pRB')  )
        
        g_prb_rb <- ggplotGrob(graphone(join, 'pRB_over_RB') + ylab('pRB_over_RB')  )
        
        g_p21 <- ggplotGrob(graphone(join, 'p21') + ylab('p21')  )
        
        g_cycA <- ggplotGrob(graphone(join, 'cycA') + ylab('cycA')  )
        
        g_cycB1 <- ggplotGrob(graphone(join, 'cycB1') + ylab('cycB1')  )
        
        g_cycD1<- ggplotGrob(graphone(join, 'cycD1') + ylab('cycD1')  )
        
        g_cycE <- ggplotGrob(graphone(join, 'cycE') + ylab('cycE'))
        
        # Combine the plots   
        g = cbind(rbind(g_prb_rb, g_cycE, size = "last"), 
                  rbind(g_ki67, g_cdt1, size = "last"),
                  rbind(g_prb, g_e2f1, size='last'),
                  rbind(g_rb, g_dna, size='last'),
                  rbind(g_cdk2, g_cycA, size='last'),
                  rbind(g_cdk4, g_cycB1, size='last'),
                  rbind(g_cycD1, g_p21, size='last'),
                  size = "first")
        
        # draw it
        grid.newpage()
        pdf(paste0(pathname, "_allfeats.pdf"), width = 20, height = 5) # Open a new pdf file
        grid.draw(g)
        dev.off()

        if (file.exists('Rplots.pdf')) {file.remove('Rplots.pdf')}
}

if (!dir.exists("output")){
        dir.create("output")
}
path = 'output'
df = read.csv("data/sketched_integrated_df.csv")
orig = 'Tumor' #enter origin of choice 
df_tum = subset(df, Origin==orig)
pathname = paste0(path, '/', orig)
ti(df_tum, pathname)

orig = 'T47D' #enter origin of choice 
df_t47 = subset(df, Origin==orig)
pathname = paste0(path, '/', orig)
ti(df_t47, pathname)

