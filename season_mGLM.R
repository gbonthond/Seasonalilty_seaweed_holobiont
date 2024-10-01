##### PACKAGES #####
library("mvabund")  # mulivariate GLMs (mGLMs)
library("ggplot2")  # plots
library("iNEXT")    # diversity indices
library("mobr")       # calculating evenness
library("DHARMa")     # for diagnostic plots when running GLMMs
library("car")        # generalized linear mixed models (GLMMs) additional package
library("DescTools")  # calculating pseudo R2
library("effects")    # plotting
library("emmeans")    # package to get means and confidence intervals for plotting; post-hoc test

##################


##### FUNCTIONS #####
permanent.core <- function(model, smry, aov = NULL, p = 0.01, coef = 0, pooled = F, p.adjust = "none"){
  dif_alga <- data.frame("taxon"         = colnames(model$y),
                         "abundance"     = colSums(model$y),
                         "occurrence"    = colSums(model$y != 0)/nrow(model$y),
                         "p_val_2"       = p.adjust(smry$uni.p[, paste("sample_type", levels(model$data$sample_type)[2], sep = "")], method = p.adjust),
                         "p_val_3"       = p.adjust(smry$uni.p[, paste("sample_type", levels(model$data$sample_type)[3], sep = "")], method = p.adjust),
                         "p_fisher"      = p.adjust(pchisq(-2*(log(smry$uni.p[, paste("sample_type", levels(model$data$sample_type)[2], sep = "")]) + log(smry$uni.p[, paste("sample_type", levels(model$data$sample_type)[3], sep = "")])), df = 4, lower.tail = F), method = p.adjust),
                         "coefficient_2" = smry$est[paste("sample_type", levels(model$data$sample_type)[2], sep = ""), -1],
                         "coefficient_3" = smry$est[paste("sample_type", levels(model$data$sample_type)[3], sep=""), -1],
                         "std_error_2"   = smry$est.stderr[paste("sample_type", levels(model$data$sample_type)[2], sep = ""), ],
                         "std_error_3"   = smry$est.stderr[paste("sample_type", levels(model$data$sample_type)[3], sep = ""), ])
  
  if(p == 0.01){ ci <- 2.807 } else { ci <- 1.96 }
  
  ### calculate lower and upper CI limits
  dif_alga$lower_2 <- dif_alga$coefficient_2 - ci * dif_alga$std_error_2
  dif_alga$upper_2 <- dif_alga$coefficient_2 + ci * dif_alga$std_error_2
  dif_alga$lower_3 <- dif_alga$coefficient_3 - ci * dif_alga$std_error_3
  dif_alga$upper_3 <- dif_alga$coefficient_3 + ci * dif_alga$std_error_3

  ### pooling stder
  dif_alga$averaged_coefficient <- (dif_alga$coefficient_2 + dif_alga$coefficient_3)/2
  dif_alga$pooled_std_error     <- sqrt((dif_alga$std_error_2^2 + dif_alga$std_error_3^2)/2)
  dif_alga$pooled_lower <- dif_alga$averaged_coefficient - ci * dif_alga$pooled_std_error
  dif_alga$pooled_upper <- dif_alga$averaged_coefficient + ci * dif_alga$pooled_std_error
  
  # subsetting
  if(is.null(aov) == T){
  dif_alga_core <- dif_alga[dif_alga$p_val_2 < p &
                            dif_alga$p_val_3 < p &
                            dif_alga$coefficient_2 < -coef &
                            dif_alga$coefficient_3 < -coef &
                            dif_alga$pooled_lower < 0 & dif_alga$pooled_upper < 0, ]
  } else {
  dif_alga_core <- dif_alga[dif_alga$p_val_2 < p &
                            dif_alga$p_val_3 < p &
                            dif_alga$coefficient_2 < -coef &
                            dif_alga$coefficient_3 < -coef &
                            dif_alga$pooled_lower < 0 & dif_alga$pooled_upper < 0, ]
  }
  
  dif_alga_stacked <- data.frame(
    "taxon" = rep(rownames(dif_alga_core), 2),
    "abundance" = rep(dif_alga_core$abundance, 2),
    "coefficient" = stack(dif_alga_core, c("coefficient_2", "coefficient_3"))[, 1],
    "sample_type" = stack(dif_alga_core, c("coefficient_2", "coefficient_3"))[, 2],
    "lower" = stack(dif_alga_core, select = c("lower_2", "lower_3"))[, 1],
    "upper" = stack(dif_alga_core, select = c("upper_2", "upper_3"))[, 1],
    "std_error" = stack(dif_alga_core, c("std_error_2", "std_error_3"))[, 1]
  )
  
  dif_alga_stacked <- dif_alga_stacked[order(dif_alga_stacked$abundance, decreasing = T), ]
 
#  dif_alga_pooled <- data.frame(
#    "taxon" = rownames(dif_alga_core),
#    "abundance" = dif_alga_core$abundance,
#    "coefficient" = dif_alga_core$averaged_coefficient,
#    "lower" = dif_alga_core$pooled_lower,
#    "upper" = dif_alga_core$pooled_upper,
#    "p_wat" = dif_alga_core$p_val_water,
#    "p_sed" = dif_alga_core$p_val_sediment,
#    "p_fisher" = pchisq(-2*(log(dif_alga_core$p_val_water) + log(dif_alga_core$p_val_sediment)), df = 4, lower.tail = F)
#  )
  
  dif_alga_core <- dif_alga_core[order(dif_alga_core$abundance, decreasing = T), ]
  
  colnames(dif_alga_core) <- c("taxon",
                               "abundance",
                               "occurrence",
                               paste("p_val", levels(model$data$sample_type)[2], sep = "_"),
                               paste("p_val", levels(model$data$sample_type)[3], sep = "_"),
                               "p_val_fisher",
                               paste("coefficient", levels(model$data$sample_type)[2], sep = "_"),
                               paste("coefficient", levels(model$data$sample_type)[3], sep = "_"),
                               paste("std_error", levels(model$data$sample_type)[2], sep = "_"),
                               paste("std_error", levels(model$data$sample_type)[3], sep = "_"),
                               paste("lower", levels(model$data$sample_type)[2], sep = "_"),
                               paste("lower", levels(model$data$sample_type)[3], sep = "_"),
                               paste("upper", levels(model$data$sample_type)[2], sep = "_"),
                               paste("upper", levels(model$data$sample_type)[3], sep = "_"),
                               "averaged_coefficient",
                               "pooled_std_error",
                               "pooled_lower",
                               "pooled_upper")
  
  if(pooled == T){
    return(dif_alga_core)
  }  else  {
    return(dif_alga_stacked)
  }
}
seasonal.core  <- function(model, smry, aov = NULL, p = 0.01, coef = 0, p.adjust = "none"){
  if(p == 0.01){ ci <- 2.807 } else { ci <- 1.96 }
  if(is.null(aov) == T){
    dif_seas <- data.frame("taxon"  = rownames(smry$uni.p),
                           "abundance"  = colSums(model$y),
                           "p_val_s" = p.adjust(smry$uni.p[, "seasonwinter"], method = p.adjust),
                           "coefficient"= smry$est["seasonwinter", -1],
                           "std_error"  = smry$est.stderr["seasonwinter", ])
    ### calculate lower and upper CI limits
    dif_seas$lower <- dif_seas$coefficient - ci * dif_seas$std_error
    dif_seas$upper <- dif_seas$coefficient + ci * dif_seas$std_error
    
    dif_seas_core <- dif_seas[dif_seas$p_val_s < p &
                              abs(dif_seas$coefficient) > coef &
                            ((dif_seas$lower < 0 & dif_seas$upper < 0) |
                             (dif_seas$lower > 0 & dif_seas$upper > 0)), ]
    dif_seas_core$core <- ifelse(dif_seas_core$coefficient < 0, "summer", "winter")
    dif_seas_core <- dif_seas_core[order(dif_seas_core$abundance, decreasing = T), ]
    return(dif_seas_core)
    
    } else {
    dif_seas <- data.frame("taxon"  = rownames(smry$uni.p),
                           "abundance"  = colSums(model$y),
                           "p_val"   = p.adjust(aov$uni.p["season", ], method = p.adjust),
                           "p_val_s" = p.adjust(smry$uni.p[, "seasonwinter"], method = p.adjust),
                           "coefficient"= smry$est["seasonwinter", -1],
                           "std_error"  = smry$est.stderr["seasonwinter", ])
    ### calculate lower and upper 95% CI limits
    dif_seas$lower <- dif_seas$coefficient - ci * dif_seas$std_error
    dif_seas$upper <- dif_seas$coefficient + ci * dif_seas$std_error
    
    dif_seas_core <- dif_seas[dif_seas$p_val < p &
                                dif_seas$p_val_s < p &
                                abs(dif_seas$coefficient) > coef &
                              ((dif_seas$lower < 0 & dif_seas$upper < 0) |
                               (dif_seas$lower > 0 & dif_seas$upper > 0)), ]
    dif_seas_core$core <- ifelse(dif_seas_core$coefficient < 0, "summer", "winter")
    dif_seas_core <- dif_seas_core[order(dif_seas_core$abundance, decreasing = T), ]
    return(dif_seas_core)
    }
  

  
}
inv.logit <- function(x) {exp(x)/(1+exp(x))}
###################


##### DATA PREPARATION #####

### variables
permanent_otu <- read.csv(file = "C:/chantal/analysis/season_otu.csv", row.names = 1);dim(permanent_otu)
permanent_var <- read.csv(file = "C:/chantal/analysis/season_var.csv", row.names = 2, stringsAsFactors = T)[, -1];dim(permanent_var)
permanent_var$seq_depth <- rowSums(permanent_otu)
rownames(permanent_var) -> rownames(permanent_otu)

permanent_ko <- read.csv(file="C:/chantal/analysis/16S_ko.csv", header = T, row.names = 1);dim(permanent_ko)
all.equal(rownames(permanent_otu), rownames(permanent_ko))

### taxonomic info 
permanent_tax <- read.csv(file = "C:/chantal/processing/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fit.optifit_mcc.0.03.cons.taxonomy", sep = "\t");dim(permanent_tax)
rownames(permanent_tax) <- permanent_tax$OTU 
permanent_tax$Taxonomy <- gsub(")(", "(", gsub(")_", "", gsub("_(", "", permanent_tax[, 3], fixed = TRUE), fixed = TRUE), fixed = TRUE)
permanent_tax$Taxonomy <- gsub("Group(DSEG","Group_DSEG", permanent_tax[, 3], fixed = TRUE)
permanent_tax$Taxonomy <- gsub("SAR324_clade(Marine_group_B","SAR324_clade_Marine_group_B", permanent_tax[, 3], fixed = TRUE)
permanent_tax$Taxonomy <- gsub("OM60(NOR5clade","OM60_NOR5_clade", permanent_tax[, 3],fixed = TRUE)
permanent_tax$Taxonomy <- gsub("[(]",";", gsub("[)]","", permanent_tax[, 3]))
permanent_tax <- data.frame(permanent_tax$OTU, do.call('rbind', strsplit(as.character(permanent_tax$Taxonomy),';',fixed = TRUE)),row.names = 1)[, c(1, 3, 5, 7, 9, 11)];dim(permanent_tax)
permanent_tax <- permanent_tax[colnames(permanent_otu), ];dim(permanent_tax)
colnames(permanent_tax) <- c("domain", "phylum", "class", "order", "family", "genus")

# verify if variables are correctly interpreted 
summary(permanent_var)

# double check if data frames match
unique(row.names(permanent_otu) == row.names(permanent_var)) 

# (optional) reduce to some proportion most abundant taxa to save computation time
permanent_otu_sub <- permanent_otu[, cumsum(colSums(permanent_otu))/sum(colSums(permanent_otu)) < 0.95 & 
                               colSums(permanent_otu != 0)/nrow(permanent_otu) > 0.25];dim(permanent_otu_sub)
permanent_tax_sub <- permanent_tax[colnames(permanent_otu_sub), ];dim(permanent_tax_sub)

# aggregating to higher taxonomic ranks
permanent_gen <- aggregate(t(permanent_otu), by = permanent_tax[, c(1:6)], FUN = "sum");dim(permanent_gen)
rownames(permanent_gen) <- paste(permanent_gen$phylum, permanent_gen$genus, 1:nrow(permanent_gen), sep = "_")
permanent_fam <- aggregate(t(permanent_otu), by = permanent_tax[, c(1:5)], FUN = "sum");dim(permanent_fam)
rownames(permanent_fam) <- paste(permanent_fam$phylum, permanent_fam$fam, 1:nrow(permanent_fam), sep = "_")
permanent_ord <- aggregate(t(permanent_otu), by = permanent_tax[, c(1:4)], FUN = "sum");dim(permanent_ord)
rownames(permanent_ord) <- paste(permanent_ord$phylum, permanent_ord$order, 1:nrow(permanent_ord), sep = "_")
permanent_cls <- aggregate(t(permanent_otu), by = permanent_tax[, c(1:3)], FUN = "sum");dim(permanent_cls)
rownames(permanent_cls) <- paste(permanent_cls$phylum, permanent_cls$class, 1:nrow(permanent_cls), sep = "_")
permanent_phl <- aggregate(t(permanent_otu), by = permanent_tax[, c(1:2)], FUN = "sum");dim(permanent_phl)
rownames(permanent_phl) <- paste(permanent_phl$phylum, permanent_phl$phylum, 1:nrow(permanent_phl), sep = "_")

### subset data for "seasonal core" to sample_type == "alga" and t6, t1 (winter) and t3, t4 (summer) only
season_var <- droplevels(permanent_var[permanent_var$sample_type == "alga" & permanent_var$timepoint != "t2" & permanent_var$timepoint != "t5", ]);dim(season_var)
season_otu <- permanent_otu_sub[rownames(season_var), colSums(permanent_otu_sub[rownames(permanent_var), ]) > 0];dim(season_otu)
all.equal(row.names(season_otu), row.names(season_var)) 

season_gen <- cbind(permanent_gen[, 1:6], permanent_gen[, rownames(season_var)]);dim(season_gen)
season_fam <- cbind(permanent_fam[, 1:5], permanent_fam[, rownames(season_var)]);dim(season_fam)
season_ord <- cbind(permanent_ord[, 1:4], permanent_ord[, rownames(season_var)]);dim(season_ord)
season_cls <- cbind(permanent_cls[, 1:3], permanent_cls[, rownames(season_var)]);dim(season_cls)
season_phl <- cbind(permanent_phl[, 1:2], permanent_phl[, rownames(season_var)]);dim(season_phl)

### make a new variable season with two levels (summer and winter)
season_var$season <- with(season_var, ifelse(timepoint == "t6", "winter", ifelse(timepoint == "t1", "winter", "summer")))

##########################


##### FIT mGLMs permanent CORES #####
permanent_var$sample_type <- relevel(permanent_var$sample_type, ref = "alga")
permanent_alg_otu_mglm <- manyglm(mvabund(permanent_otu_sub) ~ offset(log(seq_depth)) + sample_type + timepoint + year_timepoint, family = "negative.binomial", data = permanent_var)
levels(permanent_alg_otu_mglm$data$sample_type)

permanent_var$sample_type <- relevel(permanent_var$sample_type, ref = "sediment")
permanent_sed_otu_mglm <- manyglm(mvabund(permanent_otu_sub) ~ offset(log(seq_depth)) + sample_type + timepoint + year_timepoint, family = "negative.binomial", data = permanent_var)
levels(permanent_sed_otu_mglm$data$sample_type)

permanent_var$sample_type <- relevel(permanent_var$sample_type, ref = "water")
permanent_wat_otu_mglm <- manyglm(mvabund(permanent_otu_sub) ~ offset(log(seq_depth)) + sample_type + timepoint + year_timepoint, family = "negative.binomial", data = permanent_var)
permanent_var$sample_type <- relevel(permanent_var$sample_type, ref = "alga")

levels(permanent_wat_otu_mglm$data$sample_type)

### save and load the output to load model instantly next time
#save(permanent_alg_otu_mglm, file = "C:/chantal/analysis/Rdata/permanent_alg_otu_mglm.Rdata") # save model to call it next time directly
#save(permanent_sed_otu_mglm, file = "C:/chantal/analysis/Rdata/permanent_sed_otu_mglm.Rdata") # save model to call it next time directly
#save(permanent_wat_otu_mglm, file = "C:/chantal/analysis/Rdata/permanent_wat_otu_mglm.Rdata") # save model to call it next time directly
load(file = "C:/chantal/analysis/Rdata/permanent_alg_otu_mglm.Rdata") # load saved model
load(file = "C:/chantal/analysis/Rdata/permanent_sed_otu_mglm.Rdata") # load saved model
load(file = "C:/chantal/analysis/Rdata/permanent_wat_otu_mglm.Rdata") # load saved model

#permanent_alg_otu_mglm_s <- summary.manyglm(permanent_alg_otu_mglm, nBoot = 999, block = permanent_alg_otu_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#permanent_sed_otu_mglm_s <- summary.manyglm(permanent_sed_otu_mglm, nBoot = 999, block = permanent_sed_otu_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#permanent_wat_otu_mglm_s <- summary.manyglm(permanent_wat_otu_mglm, nBoot = 999, block = permanent_wat_otu_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#save(permanent_alg_otu_mglm_s, file = "C:/chantal/analysis/Rdata/permanent_alg_otu_mglm_s.Rdata")
#save(permanent_sed_otu_mglm_s, file = "C:/chantal/analysis/Rdata/permanent_sed_otu_mglm_s.Rdata")
#save(permanent_wat_otu_mglm_s, file = "C:/chantal/analysis/Rdata/permanent_wat_otu_mglm_s.Rdata")

load(file = "C:/chantal/analysis/Rdata/permanent_alg_otu_mglm_s.Rdata")
load(file = "C:/chantal/analysis/Rdata/permanent_sed_otu_mglm_s.Rdata")
load(file = "C:/chantal/analysis/Rdata/permanent_wat_otu_mglm_s.Rdata")

#################################


##### GET permanent CORE OTUs ##### 

permanent_alg_core <- permanent.core(permanent_alg_otu_mglm, permanent_alg_otu_mglm_s, p = 0.01, coef = 0, pooled = F, p.adjust = "fdr");dim(permanent_alg_core)
permanent_sed_core <- permanent.core(permanent_sed_otu_mglm, permanent_sed_otu_mglm_s, p = 0.01, coef = 0, pooled = F, p.adjust = "fdr");dim(permanent_sed_core)
permanent_wat_core <- permanent.core(permanent_wat_otu_mglm, permanent_wat_otu_mglm_s, p = 0.01, coef = 0, pooled = F, p.adjust = "fdr");dim(permanent_wat_core)
permanent_alg_core_otu <- permanent.core(permanent_alg_otu_mglm, permanent_alg_otu_mglm_s, p = 0.01, coef = 0, pooled = T, p.adjust = "fdr");dim(permanent_alg_core_otu)
permanent_sed_core_otu <- permanent.core(permanent_sed_otu_mglm, permanent_sed_otu_mglm_s, p = 0.01, coef = 0, pooled = T, p.adjust = "fdr");dim(permanent_sed_core_otu)
permanent_wat_core_otu <- permanent.core(permanent_wat_otu_mglm, permanent_wat_otu_mglm_s, p = 0.01, coef = 0, pooled = T, p.adjust = "fdr");dim(permanent_wat_core_otu)

# plot coefficients and corresponding confidence intervals
{

# asthetics 
{
theme <- theme(legend.title = element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               legend.position = "top",
               panel.background = element_rect(fill = "transparent"),
               panel.border = element_rect(fill = "transparent",colour = "black",linewidth = 1),
               axis.ticks.length = unit(1, "mm"),
               axis.ticks.x = element_line(colour="#333333", linewidth = 1),
               axis.ticks.y = element_line(colour="#333333", linewidth = 1))
}

  gg_permanent_alg_otu_sw <-  ggplot(permanent_alg_core,
                                     aes(x = reorder(taxon, coefficient, col = sample_type), y = coefficient)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper, col = sample_type), alpha = .5, width = 0, linewidth = 2) +
    geom_point(aes(y = coefficient, col = sample_type), size = 2) +
    coord_flip() +
    geom_hline(yintercept=0, linetype="dashed", color = "red",linewidth = 0.75) + theme +
    ggtitle("fold change with respect to alga + 99%CIs");gg_permanent_alg_otu_sw
  
  gg_permanent_sed_otu_sw <-  ggplot(permanent_sed_core[1:100, ],
                                     aes(x = reorder(taxon, coefficient, col = sample_type), y = coefficient)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper, col = sample_type), alpha = .5, width = 0, linewidth = 2) +
    geom_point(aes(y = coefficient, col = sample_type), size = 2) +
    coord_flip() +
    geom_hline(yintercept=0, linetype="dashed", color = "red",linewidth = 0.75) + theme +
    ggtitle("fold change with respect to alga + 99%CIs");gg_permanent_sed_otu_sw
  
  gg_permanent_wat_otu_sw <-  ggplot(permanent_wat_core,
                                     aes(x = reorder(taxon, coefficient, col = sample_type), y = coefficient)) + 
    geom_errorbar(aes(ymin = lower, ymax = upper, col = sample_type), alpha = .5, width = 0, linewidth = 2) +
    geom_point(aes(y = coefficient, col = sample_type), size = 2) +
    coord_flip() +
    geom_hline(yintercept=0, linetype="dashed", color = "red",linewidth = 0.75) + theme +
    ggtitle("fold change with respect to alga + 99%CIs");gg_permanent_wat_otu_sw
  
  gg_permanent_alg_otu <- ggplot(permanent_alg_core_otu[1:25, ],
                                 aes(x = reorder(taxon, -averaged_coefficient, y = -averaged_coefficient))) +
    geom_errorbar(aes(ymin = -pooled_lower, ymax = -pooled_upper), col = "darkgreen", width = 0, size = 1.25) +
    geom_point(aes(y = -averaged_coefficient), size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red",size = 0.75) + theme +
    ggtitle("permanent algal core - OTUs");gg_permanent_alg_otu

  gg_permanent_sed_otu <- ggplot(permanent_sed_core_otu[1:25, ],
                                 aes(x = reorder(taxon, -averaged_coefficient, y = -averaged_coefficient))) +
    geom_errorbar(aes(ymin = -pooled_lower, ymax = -pooled_upper), col = "orange", width = 0, size = 1.25) +
    geom_point(aes(y = -averaged_coefficient), size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red",size = 0.75) + theme +
    ggtitle("permanent sediment core - OTUs");gg_permanent_sed_otu
  
  gg_permanent_wat_otu <- ggplot(permanent_wat_core_otu[1:25, ],
                                 aes(x = reorder(taxon, -averaged_coefficient, y = -averaged_coefficient))) +
    geom_errorbar(aes(ymin = -pooled_lower, ymax = -pooled_upper), col = "darkcyan", width = 0, size = 1.25) +
    geom_point(aes(y = -averaged_coefficient), size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red",size = 0.75) + theme +
    ggtitle("permanent water core - OTUs");gg_permanent_wat_otu
  }
################################


##### FIT mGLMs SEASONAL CORES ##### 

# 1) fitting
#season_otu_mglm <- manyglm(mvabund(season_otu)               ~ offset(log(seq_depth)) + season + year_timepoint, family = "negative.binomial", data = season_var)
#season_gen_mglm <- manyglm(mvabund(t(season_gen[, -c(1:6)])) ~ offset(log(seq_depth)) + season + year_timepoint, family = "negative.binomial", data = season_var)
#season_fam_mglm <- manyglm(mvabund(t(season_fam[, -c(1:5)])) ~ offset(log(seq_depth)) + season + year_timepoint, family = "negative.binomial", data = season_var)
#season_ord_mglm <- manyglm(mvabund(t(season_ord[, -c(1:4)])) ~ offset(log(seq_depth)) + season + year_timepoint, family = "negative.binomial", data = season_var)
#season_cls_mglm <- manyglm(mvabund(t(season_cls[, -c(1:3)])) ~ offset(log(seq_depth)) + season + year_timepoint, family = "negative.binomial", data = season_var)
#season_phl_mglm <- manyglm(mvabund(t(season_phl[, -c(1:2)])) ~ offset(log(seq_depth)) + season + year_timepoint, family = "negative.binomial", data = season_var)
#save(season_otu_mglm, file = "C:/chantal/analysis/Rdata/season_otu_mglm.Rdata") # save model to call it next time directly
#save(season_gen_mglm, file = "C:/chantal/analysis/Rdata/season_gen_mglm.Rdata") # save model to call it next time directly
#save(season_fam_mglm, file = "C:/chantal/analysis/Rdata/season_fam_mglm.Rdata") # save model to call it next time directly
#save(season_ord_mglm, file = "C:/chantal/analysis/Rdata/season_ord_mglm.Rdata") # save model to call it next time directly
#save(season_cls_mglm, file = "C:/chantal/analysis/Rdata/season_cls_mglm.Rdata") # save model to call it next time directly
#save(season_phl_mglm, file = "C:/chantal/analysis/Rdata/season_phl_mglm.Rdata") # save model to call it next time directly
load(file = "C:/chantal/analysis/Rdata/season_otu_mglm.Rdata") # load saved model
load(file = "C:/chantal/analysis/Rdata/season_gen_mglm.Rdata") # load saved model
load(file = "C:/chantal/analysis/Rdata/season_fam_mglm.Rdata") # load saved model
load(file = "C:/chantal/analysis/Rdata/season_ord_mglm.Rdata") # load saved model
load(file = "C:/chantal/analysis/Rdata/season_cls_mglm.Rdata") # load saved model
load(file = "C:/chantal/analysis/Rdata/season_phl_mglm.Rdata") # load saved model

### summary.manyglm() 
#season_otu_mglm_s <- summary.manyglm(season_otu_mglm, nBoot = 999, block = season_otu_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#season_gen_mglm_s <- summary.manyglm(season_gen_mglm, nBoot = 999, block = season_gen_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#season_fam_mglm_s <- summary.manyglm(season_fam_mglm, nBoot = 999, block = season_fam_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#season_ord_mglm_s <- summary.manyglm(season_ord_mglm, nBoot = 999, block = season_ord_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#season_cls_mglm_s <- summary.manyglm(season_cls_mglm, nBoot = 999, block = season_cls_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#season_phl_mglm_s <- summary.manyglm(season_phl_mglm, nBoot = 999, block = season_phl_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#save(season_otu_mglm_s, file = "C:/chantal/analysis/Rdata/season_otu_mglm_s.Rdata")
#save(season_gen_mglm_s, file = "C:/chantal/analysis/Rdata/season_gen_mglm_s.Rdata")
#save(season_fam_mglm_s, file = "C:/chantal/analysis/Rdata/season_fam_mglm_s.Rdata")
#save(season_ord_mglm_s, file = "C:/chantal/analysis/Rdata/season_ord_mglm_s.Rdata")
#save(season_cls_mglm_s, file = "C:/chantal/analysis/Rdata/season_cls_mglm_s.Rdata")
#save(season_phl_mglm_s, file = "C:/chantal/analysis/Rdata/season_phl_mglm_s.Rdata")
load(file = "C:/chantal/analysis/Rdata/season_otu_mglm_s.Rdata")
load(file = "C:/chantal/analysis/Rdata/season_gen_mglm_s.Rdata")
load(file = "C:/chantal/analysis/Rdata/season_fam_mglm_s.Rdata")
load(file = "C:/chantal/analysis/Rdata/season_ord_mglm_s.Rdata")
load(file = "C:/chantal/analysis/Rdata/season_cls_mglm_s.Rdata")
load(file = "C:/chantal/analysis/Rdata/season_phl_mglm_s.Rdata")

#################################


##### GET SEASONAL CORE OTUs ##### 
season_core_otu <- seasonal.core(season_otu_mglm, season_otu_mglm_s, p = 0.01, coef = 0, p.adjust = "fdr");dim(season_core_otu)

# exclude permanent sediment and water cores
season_core_otu <- season_core_otu[rownames(season_core_otu) %in% rownames(permanent_sed_core_otu) == F &
                                   rownames(season_core_otu) %in% rownames(permanent_wat_core_otu) == F, ];dim(season_core_otu)

# plot coefficients and corresponding confidence intervals
# aesthetics 
{
  theme <- theme(legend.title = element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 #legend.position = "top",
                 panel.background = element_rect(fill = "transparent"),
                 panel.border = element_rect(fill = "transparent",colour = "black", linewidth = 1),
                 axis.ticks.length = unit(1, "mm"),
                 axis.ticks.x = element_line(colour = "#333333", linewidth = 1),
                 axis.ticks.y = element_line(colour = "#333333", linewidth = 1))
}


gg_season_otu_sw <-  ggplot(season_core_otu[1:50, ],
                            aes(x = reorder(taxon, -coefficient), y = coefficient)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, col = core), alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = coefficient), size = 2) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + theme +
  ggtitle("fold change with respect to season + 99%CIs");gg_season_otu_sw

################################


##### FIGURES PAPER #####

permanent_alg_core_otu$label <- paste(permanent_alg_core_otu$taxon,
                                      permanent_tax_sub[permanent_alg_core_otu$taxon, "genus"])

gg_permanent_otu_sw <- ggplot(permanent_alg_core_otu[1:25, ],
                         aes(x = reorder(label, -averaged_coefficient), y = -averaged_coefficient)) +
  geom_errorbar(aes(ymin = -pooled_lower, ymax = -pooled_upper), col = "darkgreen", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = -averaged_coefficient), size = 3, stroke = 2, shape = 21, colour = "black", fill = "white") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + theme +
  ggtitle("permanent");gg_permanent_otu_sw

season_core_otu$label <- paste(season_core_otu$taxon, permanent_tax_sub[season_core_otu$taxon, "genus"])
  
gg_summer_otu_sw <-  ggplot(season_core_otu[season_core_otu$core == "summer", ][1:25, ],
                            aes(x = reorder(label, -coefficient), y = -coefficient)) + 
  geom_errorbar(aes(ymin = -lower, ymax = -upper), col = "darkred", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = -coefficient), size = 3, stroke = 2, shape = 21, colour = "black", fill = "white") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + theme +
  ggtitle("Summer");gg_summer_otu_sw

gg_winter_otu_sw <-  ggplot(season_core_otu[season_core_otu$core == "winter", ][1:25, ],
                            aes(x = reorder(label, coefficient), y = coefficient)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), col = "darkblue", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = coefficient), size = 3, stroke = 2, shape = 21, colour = "black", fill = "white") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + theme +
  ggtitle("winter");gg_winter_otu_sw
#######################


##### TABLE S6 ####

# core table excluding sediment and water samples
core_table <- permanent_tax
core_table$otu <- rownames(core_table)
core_table[, "abundance"] <- colSums(permanent_otu[, rownames(core_table)]/rowSums(permanent_otu[, rownames(core_table)]))/nrow(permanent_otu)
core_table[, "occurrence"] <- colSums(permanent_otu[permanent_var$sample_type == "alga", rownames(core_table)] != 0)/nrow(permanent_otu[permanent_var$sample_type == "alga", rownames(core_table)])
#with(core_table, plot(log10(abundance), occurrence))

core_table[, "occurrence_permanent"] <- ifelse((colSums(permanent_otu[permanent_var$sample_type == "alga", rownames(core_table)] != 0)/nrow(permanent_otu[permanent_var$sample_type == "alga", rownames(core_table)]) == 1) == TRUE, "epi", "")
core_table[, "occurrence_season"]    <- ifelse(!rownames(core_table) %in% colnames(permanent_otu_sub),"",
                                              ifelse((colSums(season_otu[season_var$sample_type == "alga" & (season_var$timepoint == "t6" | season_var$timepoint == "t1"), 
                                                                  colnames(permanent_otu_sub) %in% rownames(core_table)] != 0)/nrow(season_otu[season_var$sample_type == "alga" & 
                                                                                                                (season_var$timepoint == "t6" | season_var$timepoint == "t1"), 
                                                                                                                colnames(permanent_otu_sub) %in% rownames(core_table)]) == 1) == TRUE, "winter", 
                                              ifelse((colSums(season_otu[season_var$sample_type == "alga" & (season_var$timepoint == "t3" | season_var$timepoint == "t4"), 
                                                                         colnames(permanent_otu_sub) %in% rownames(core_table)] != 0)/nrow(season_otu[season_var$sample_type == "alga" & (season_var$timepoint == "t3" | season_var$timepoint == "t4"), 
                                                                                                                       colnames(permanent_otu_sub) %in% rownames(core_table)]) == 1) == TRUE, "summer", "")))

for(i in 1:nrow(core_table)) {
  core_table[i, "permanent"] <- ifelse(rownames(permanent_tax[i, ]) %in% permanent_alg_core_otu$taxon, "epi",
                                ifelse(rownames(permanent_tax[i, ]) %in% permanent_sed_core_otu$taxon, "sediment",
                                ifelse(rownames(permanent_tax[i, ]) %in% permanent_wat_core_otu$taxon, "water", "")))
  core_table[i, "perm_coef"] <- ifelse(rownames(permanent_tax[i, ]) %in% permanent_alg_core_otu$taxon, permanent_alg_core_otu[permanent_alg_core_otu$taxon == rownames(permanent_tax[i, ]), "averaged_coefficient"],
                                ifelse(rownames(permanent_tax[i, ]) %in% permanent_sed_core_otu$taxon, permanent_sed_core_otu[permanent_sed_core_otu$taxon == rownames(permanent_tax[i, ]), "averaged_coefficient"], 
                                ifelse(rownames(permanent_tax[i, ]) %in% permanent_wat_core_otu$taxon, permanent_wat_core_otu[permanent_wat_core_otu$taxon == rownames(permanent_tax[i, ]), "averaged_coefficient"], "")))
  core_table[i, "p_fisher"] <-  ifelse(rownames(permanent_tax[i, ]) %in% permanent_alg_core_otu$taxon, permanent_alg_core_otu[permanent_alg_core_otu$taxon == rownames(permanent_tax[i, ]), "p_val_fisher"],
                                ifelse(rownames(permanent_tax[i, ]) %in% permanent_sed_core_otu$taxon, permanent_sed_core_otu[permanent_sed_core_otu$taxon == rownames(permanent_tax[i, ]), "p_val_fisher"], 
                                ifelse(rownames(permanent_tax[i, ]) %in% permanent_wat_core_otu$taxon, permanent_wat_core_otu[permanent_wat_core_otu$taxon == rownames(permanent_tax[i, ]), "p_val_fisher"], "")))
  }

for(i in 1:nrow(core_table)) {
  core_table[i, "season"] <- ifelse(rownames(permanent_tax[i, ]) %in% season_core_otu$taxon, season_core_otu[rownames(permanent_tax[i, ]), "core"], "")
  core_table[i, "seas_coef"] <- ifelse(rownames(permanent_tax[i, ]) %in% season_core_otu$taxon, season_core_otu[rownames(permanent_tax[i, ]), "coefficient"], "")
  core_table[i, "seas_p"] <- ifelse(rownames(permanent_tax[i, ]) %in% season_core_otu$taxon, season_core_otu[rownames(permanent_tax[i, ]), "p_val_s"], "")
}

View(core_table)

write.csv(core_table[order(core_table$occurrence, decreasing = T), ], file = "C:/chantal/analysis/core_table.csv")


##################


#### DIVERSITY ####
permanent_var$otu_chao <- ChaoRichness(t(permanent_otu))[, 2]
permanent_var$otu_PIE <- calc_PIE(permanent_otu, ENS = F)
permanent_var$ko_chao <- ChaoRichness(t(permanent_ko))[, 2]
permanent_var$ko_PIE <- calc_PIE(permanent_ko, ENS = F)

alga_var <- droplevels(permanent_var[permanent_var$sample_type == "alga", ]);dim(alga_var)

glm_otu_chao   <- glm(otu_chao ~ log(seq_depth)+ timepoint * year_timepoint * pop, data = alga_var, family = gaussian(link= "log"), na.action = na.fail)
plot(simulateResiduals(fittedModel = glm_otu_chao , n = 1000))
Anova(glm_otu_chao)
glm_otu_chao_effect <- data.frame(emmeans(glm_otu_chao, specs = "timepoint"));glm_otu_chao_effect
glm_otu_chao_pp <- emmeans(glm_otu_chao, pairwise ~ timepoint)
glm_otu_chao_pseudo_R2 <- 1-(glm_otu_chao$deviance/glm_otu_chao$null.deviance); glm_otu_chao_pseudo_R2


glm_otu_PIE   <- glm(logit(otu_PIE) ~ log(seq_depth)+ timepoint * year_timepoint * pop, data = alga_var, family = gaussian(link= "identity"), na.action = na.fail)
plot(simulateResiduals(fittedModel = glm_otu_PIE , n = 1000))
Anova(glm_otu_PIE)
glm_otu_PIE_effect <- data.frame(emmeans(glm_otu_PIE, specs = "timepoint"));glm_otu_PIE_effect
glm_otu_PIE_pp <- emmeans(glm_otu_PIE, pairwise ~ timepoint)
glm_otu_PIE_pseudo_R2 <- 1-(glm_otu_PIE$deviance/glm_otu_PIE$null.deviance); glm_otu_PIE_pseudo_R2


glm_ko_chao   <- glm(ko_chao ~ log(seq_depth)+ timepoint * year_timepoint * pop, data = alga_var, family = gaussian(link= "log"), na.action = na.fail)
plot(simulateResiduals(fittedModel = glm_ko_chao , n = 1000))
Anova(glm_ko_chao)
glm_ko_chao_effect <- data.frame(emmeans(glm_ko_chao, specs = "timepoint"));glm_ko_chao_effect
glm_ko_chao_pp <- emmeans(glm_ko_chao, pairwise ~ timepoint)
glm_otu_chao_pseudo_R2 <- 1-(glm_ko_chao$deviance/glm_ko_chao$null.deviance); glm_otu_chao_pseudo_R2


glm_ko_PIE   <- glm(logit(ko_PIE) ~ log(seq_depth)+ timepoint * year_timepoint * pop, data = alga_var, family = gaussian(link= "identity"), na.action = na.fail)
plot(simulateResiduals(fittedModel = glm_ko_PIE , n = 1000))
Anova(glm_ko_PIE, test.statistic = "LR")
glm_ko_PIE_effect <- data.frame(emmeans(glm_ko_PIE, specs = "timepoint"));glm_ko_PIE_effect
glm_ko_PIE_pp <- emmeans(glm_ko_PIE, pairwise ~ timepoint)
glm_otu_chao_pseudo_R2 <- 1-(glm_ko_PIE$deviance/glm_ko_PIE$null.deviance); glm_otu_chao_pseudo_R2


### figures
ggplot_glm_otu_chao <- ggplot(glm_otu_chao_effect, aes(x = timepoint, y = exp(emmean))) + 
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), width = 0, size = 1 , col = 1) +
  geom_point(size=2, stroke=1, col=1) + 
  ylim(c(0,4100)) +
  theme(title = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent",colour = "black",size = 1),
        axis.ticks.length = unit(1, "mm"),
        axis.ticks = element_line(colour="#333333",size = 1),
        plot.margin = unit(c(1,1,1,1), "mm"))

ggplot_glm_otu_PIE <- ggplot(glm_otu_PIE_effect, aes(x = timepoint, y = inv.logit(emmean))) + 
  geom_errorbar(aes(ymin = inv.logit(lower.CL), ymax = inv.logit(upper.CL)), width = 0, size = 1 , col = 1) +
  geom_point(size=2, stroke=1, col=1) + 
  #ylim(c(0,4100)) +
  theme(title = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent",colour = "black",size = 1),
        axis.ticks.length = unit(1, "mm"),
        axis.ticks = element_line(colour="#333333",size = 1),
        plot.margin = unit(c(1,1,1,1), "mm"))

ggplot_glm_ko_chao <- ggplot(glm_ko_chao_effect, aes(x = timepoint, y = exp(emmean))) + 
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), width = 0, size = 1 , col = 1) +
  geom_point(size=2, stroke=1, col=1) + 
  theme(title = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent",colour = "black",size = 1),
        axis.ticks.length = unit(1, "mm"),
        axis.ticks = element_line(colour="#333333",size = 1),
        plot.margin = unit(c(1,1,1,1), "mm"))

ggplot_glm_ko_PIE <- ggplot(glm_ko_PIE_effect, aes(x = timepoint, y = inv.logit(emmean))) + 
  geom_errorbar(aes(ymin = inv.logit(lower.CL), ymax = inv.logit(upper.CL)), width = 0, size = 1 , col = 1) +
  geom_point(size=2, stroke=1, col=1) + 
  theme(title = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent",colour = "black",size = 1),
        axis.ticks.length = unit(1, "mm"),
        axis.ticks = element_line(colour="#333333",size = 1),
        plot.margin = unit(c(1,1,1,1), "mm"))


library("egg")
ggarrange(ggplot_glm_otu_chao,
          ggplot_glm_otu_PIE,
          ggplot_glm_ko_chao,
          ggplot_glm_ko_PIE,
          ncol=2,widths = c(1,1),
          labels = c("A","B","C","D"),
          label.args = list(gp = grid::gpar(font = 1, cex = 2), x=unit(1,"line")))
ggsave(file="C:/chantal/figs/diversity.pdf",device="pdf",width = 85,height = 100,units="mm",
       ggarrange(ggplot_glm_otu_chao,
                 ggplot_glm_otu_PIE,
                 ggplot_glm_ko_chao,
                 ggplot_glm_ko_PIE,
                 ncol=2,widths = c(1,1),
                 labels = c("A","B","C","D"),
                 label.args = list(gp = grid::gpar(font = 1, cex = 2), x=unit(1,"line"))))


write.csv(
  cbind(
  data.frame(glm_otu_chao_pp$contrasts)[, c(1, 2, 5, 6)],
  data.frame(glm_otu_PIE_pp$contrasts)[, c(2, 5, 6)],
  data.frame(glm_ko_chao_pp$contrasts)[, c(2, 5, 6)],
  data.frame(glm_ko_PIE_pp$contrasts)[, c(2, 5, 6)]),
  file = "C:/chantal/analysis/post-hoc.csv")


############################



##### BARPLOT #####
test <- stacked.bars(taxon_level = "family",
             var_list = c("year_timepoint", "timepoint", "pop", "SampleID"),
             number = 30,
             input_community = season_otu,
             input_tax = season_tax[colnames(season_otu), ],
             input_variable = season_var,
             tax_disp  = c("phylum", "family"),
             X_levels = NULL, trans, bars = T) 

test + theme(axis.text.x = element_blank(), axis.text.y = element_text(lineheight = 7), axis.title = element_blank(), title = element_blank(), legend.title = element_blank(), panel.spacing.y = unit(c(0, 0), "mm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "darkgrey"), panel.border = element_rect(fill = "transparent", colour = "black",linewidth = 1), axis.ticks.length.y = unit(1, "mm"), axis.ticks.y = element_line(colour = "#333333", linewidth = .1), axis.ticks.length.x = unit(0, "mm"), plot.margin = unit(c(1, 1, 1, 1), "mm")) +
  scale_fill_manual(values = c("#da814e","#b15928",
                                         "#d2ebf2","#c0dce6","#aecddb","#98bccd","#85acc1","#729cb4","#5e8ca7","#4b7d9c","#386d8f","#145078",
                                         "#c6dfaf","#7ec06e","#33a02c",
                                         "#ffe88d","#ffde5c","#ffcc00",
                                         "#dccce4","#a384bf","#6a3d9a",
                                         "#cbe5d0","#94c1a6","#247852",
                                         "#FDBF6F","#ff7f00",
                                         "#fbcecd","#fb9a99","#e31a1c","#800080"))
#################


##### SUMMARY ######
seas_otu <- season_tax[colnames(season_otu), ]
seas_otu$abundance <- colSums(season_otu[, rownames(seas_otu)])
seas_otu$prop <- seas_otu$abundance/sum(seas_otu$abundance)
seas_otu$perc <- seas_otu$prop*100

seas_gen <- aggregate(seas_otu[, 7:9], by =      seas_otu[, 1:6], FUN = "sum")
seas_fam <- aggregate(seas_otu[, 7:9], by =      seas_otu[, 1:5], FUN = "sum")
seas_ord <- aggregate(seas_otu[, 7:9], by =      seas_otu[, 1:4], FUN = "sum")
seas_cls <- aggregate(seas_otu[, 7:9], by =      seas_otu[, 1:3], FUN = "sum")
seas_phl <- aggregate(seas_otu[, 7:9], by =      seas_otu[, 1:2], FUN = "sum")
seas_dom <- aggregate(seas_otu[, 7:9], by = list(seas_otu[, 1  ]), FUN = "sum")

View(seas_fam)
View(seas_phl)

