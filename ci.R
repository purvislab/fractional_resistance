library(ggplot2)
library(ggridges)
library(tidyverse)
library(ggbeeswarm)
library(grid)
library(gridExtra)

violin_cis <- function(dat, pathname, pathname_ci) {
        
        col_order = c('pRB_over_RB', 'Ki67', 'pRB', 'RB', 'CDK2', 'CDK4', 'cycD1', 'cycE', 'Cdt1', 'E2F1', 'DNA', 'cycA', 'cycB1', 'p21')
        #get rid of metadata columns
        well_sub <- function(df, wl) {
                df = subset(df, well==wl)
                df = subset(df, select=col_order)
                
                return(df)
        }
        
        dat_0 = well_sub(dat, 0)
        dat_10 = well_sub(dat, 10)
        dat_100 = well_sub(dat, 100)

        tt_dat_0vs10_low = rep(NA,ncol(dat_0))
        tt_dat_0vs10_high = rep(NA,ncol(dat_0))
        
        tt_dat_0vs100_low = rep(NA, ncol(dat_0))
        tt_dat_0vs100_high = rep(NA, ncol(dat_0))
        
        for (i in (1:ncol(dat_0))) {
                
                zero_10 = t.test(dat_10[,i], dat_0[,i],
                                 alternative='two.sided',
                                 var.equal=TRUE)
                
                zero_100 = t.test(dat_100[,i], dat_0[,i],
                                  alternative='two.sided',
                                  var.equal=TRUE)
                tt_dat_0vs10_low[i] = zero_10$conf.int[1]
                tt_dat_0vs10_high[i] = zero_10$conf.int[2]  
                
                tt_dat_0vs100_low[i] = zero_100$conf.int[1]
                tt_dat_0vs100_high[i] = zero_100$conf.int[2]
        }
        dat_colnames <- col_order #feature list
        
        gf = data.frame(features=dat_colnames, lower=tt_dat_0vs10_low, upper=tt_dat_0vs10_high,
                        mid=(tt_dat_0vs10_high + tt_dat_0vs10_low)/2)
        
        gf$features = factor(gf$features, levels=gf$features)
        
        gf$group = '0vs10'
        
        d2 <- ggplot(gf, aes(features, mid)) + 
                geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
                geom_hline(yintercept=0, linetype="dashed", color = "red") +
                ggtitle("95% CI between 0 vs. 10") +
                ylab("95% CI for Effect Size of Feature Difference from two-sample t-test") + xlab("Feature")
        
        gf2 = data.frame(features=dat_colnames, lower=tt_dat_0vs100_low, upper=tt_dat_0vs100_high,
                         mid=(tt_dat_0vs100_high + tt_dat_0vs100_low)/2)
        gf2$features = factor(gf2$features, levels=gf2$features)
        gf2$group = '0vs100'
        
        d3 <- ggplot(gf2, aes(features, mid)) + geom_point(size=2) + 
                geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
                geom_hline(yintercept=0, linetype="dashed", color = "red") +
                ylab("95% CI for Effect Size from two-sample t-test") + xlab("Feature")
        
        
        gf_tot = rbind(gf, gf2)
        ci <- ggplot(data = gf_tot,
                     aes(
                             x = group,
                             y = mid,
                             ymin = lower,
                             ymax = upper
                     )) +
                geom_pointrange(aes(col = group), size=2) +
                geom_hline(yintercept = 0,
                           linetype = "dashed",
                           color = "black") +
                geom_errorbar(aes(ymin = lower, ymax = upper, col = group),
                              width = 0.2,
                              cex = 3) +
                scale_color_manual(values = c("#A554F5", "#E3822D"), name="95% CI",
                                   breaks=c("0vs10", "0vs100"),
                                   labels=c("0 vs. 10 nM", "0 vs. 100 nM")) + 
                theme_classic() + scale_x_discrete(position = "top") + 
                
                facet_wrap(
                        ~ features,
                        strip.position = "top",
                        nrow = 1,
                        scales = "free_x"
                ) +
                theme(strip.background = element_blank(), strip.text.x = element_blank()) +
                theme(
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(colour = "black")
                ) +
                theme(axis.text.x = element_blank(),                                                                                     axis.ticks.x = element_blank()) +
                labs(x = NULL, y = NULL) + theme(text = element_text(size = 54)) + theme(legend.position = 'none')
        
        
        dat_sub = dat[, c(col_order, 'prb_ratio', 'well')]
        
        dat_colnames = col_order
        
        #vertically cap the data so it's not too stretched out
        for (col in col_order) {
                
                perc = 3.8
                dat_sub[,col][dat_sub[,col]>=perc] <- perc
                
        }
        
        dat_long <- gather(dat_sub, key="measure", value="value", all_of(dat_colnames))
        dat_long_high = subset(dat_long, prb_ratio==1)
        dat_long_high$well = factor(dat_long_high$well, levels=c("0","10","100"), labels=c("0","10","100"))
        dat_long_high$measure2 <- factor(dat_long_high$measure, levels=col_order)
        
        
        violins <-
                 ggplot(dat_long_high,
                        aes(
                                x = well,
                                y = value,
                                color = as.factor(well)
                        )) +
                 geom_quasirandom(cex = 1) + scale_color_brewer(palette = "Set1") +
                 facet_grid(cols=vars(measure2), scale = 'free_y') + theme(strip.text = element_text(size = 8)) + theme_classic() + labs(x = NULL, y = NULL) + #scale_x_reordered() +
                 theme(strip.background = element_blank(), strip.text.x = element_blank()) +
                 theme(legend.position = "none") + theme(axis.text.x = element_blank(),                                                                                     axis.ticks.x = element_blank()) +
                 theme(text = element_text(size = 54))


         g1 <- ggplotGrob(violins)
         g2 <- ggplotGrob(ci)
         maxWidth = grid::unit.pmax(g1$widths[2:5], g2$widths[2:5])
         g1$widths[2:5] <- as.list(maxWidth)
         g2$widths[2:5] <- as.list(maxWidth)
         png(pathname, width=16, height=8, units="in", res=1200)
         grid.arrange(g1, g2, ncol=1)
         dev.off()
        
        # #save CI values
         write_csv(gf_tot, pathname_ci)
         if (file.exists('Rplots.pdf')) {file.remove('Rplots.pdf')}
}

set.seed(1453)
if (!dir.exists("output")){
        dir.create("output")
}
select_df <- function(df) {
        df_prb_high = subset(df, prb_ratio==1)
        well_temp = as.factor(df_prb_high$well)
        df_prb_high$well = well_temp
        return(df_prb_high)
}
path = 'output/'
dat = read.csv('data/sketched_integrated_df.csv')
dat = select_df(dat)

orig = 'T47D'
dat_t47 = subset(dat, Origin == orig)
pathname_t47 = paste0(path, orig,"_varequal.png")
pathname_ci_t47 = paste0(path, orig,"_CIs.csv")

orig = 'Tumor'
dat_tum = subset(dat, Origin == orig)
pathname_tum = paste0(path, orig,"_varequal.png")
pathname_ci_tum = paste0(path, orig,"_CIs.csv")

print('T47D')
violin_cis(dat_t47, pathname_t47, pathname_ci_t47)
print('Tumor')
violin_cis(dat_tum, pathname_tum, pathname_ci_tum)

#replicate
t47_rep = read.csv('data/sketched_rep_df.csv')
t47_rep = select_df(t47_rep)
orig = 'T47D_rep'
pathname_rep = paste0(path, orig,"_varequal.png")
pathname_ci_rep = paste0(path, orig,"_CIs.csv")
print('T47D Replicate')
violin_cis(t47_rep, pathname_rep, pathname_ci_rep)