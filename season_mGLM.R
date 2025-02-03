##### PACKAGES #####
library("suncalc")  # get day lengths
library("mobr")     # calculating diversity 
library("car")      # generalized linear models additional package
library("effects")  # plotting
library("emmeans")  # means and confidence intervals and post-hoc test
library("ggplot2")  # plots
library("egg")
library("mvabund")  # mulivariate GLMs (mGLMs)
library("vegan")  

##################


##### FUNCTIONS #####
permanent.core <- function(model, smry, aov = NULL, p = 0.01, coef = 0, pooled = F, p.adjust = "none"){
  dif_alga <- data.frame("taxon"         = colnames(model$y),
                         "abundance"     = colSums(model$y),
                         "occupancy"     = colSums(model$y != 0)/nrow(model$y),
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
 
  dif_alga_core <- dif_alga_core[order(dif_alga_core$abundance, decreasing = T), ]
  
  colnames(dif_alga_core) <- c("taxon",
                               "abundance",
                               "occupancy",
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
occupancy <- function(community, taxonomy) {
  com <- community[, rownames(taxonomy)]
  com[com != 0] <- 1
  occupancy <- colSums(com)/nrow(com)
}
abundance <- function(community, taxonomy, proportional = T) {
  com <- community[, rownames(taxonomy)]
  if(proportional == F){
    abundance <- colSums(com)}
  else
    abundance <- colSums(com)/sum(colSums(com))
}
rarefaction <- function(x, sample, replicate) {
  library("progress")
  library("vegan")
  pb <- progress_bar$new(total = replicate)
  rar_list <- list()
  for (i in 1:replicate) {
    rar_list[[i]] <- rrarefy(x, sample)
    pb$tick()
    Sys.sleep(1 / 100)
  }
  rar_table <- Reduce(`+`, rar_list)/replicate
  rar_table <- data.frame(rar_table[, colSums(rar_table) > 0])
  rm(rar_list)
  invisible(rar_table)
}
###################


##### DATA PREPARATION #####
permanent_otu <- read.csv(file = "https://raw.githubusercontent.com/gbonthond/Seasonalilty_seaweed_holobiont/refs/heads/main/seasonality_otu.csv", header = T, row.names = 1);dim(permanent_otu)
#permanent_rar <- rarefaction(permanent_otu, sample = 1000, replicate = 100);dim(permanent_rar)
#write.csv(permanent_rar, file = "permanent_rar.csv");dim(permanent_rar)
permanent_rar <- read.csv(file = "https://raw.githubusercontent.com/gbonthond/Seasonalilty_seaweed_holobiont/refs/heads/main/permanent_rar.csv", header = T, row.names = 1);dim(permanent_rar)
permanent_var <- read.csv(file = "https://raw.githubusercontent.com/gbonthond/Seasonalilty_seaweed_holobiont/refs/heads/main/seasonality_var.csv", header = T, row.names = 1, stringsAsFactors = T);dim(permanent_var)
permanent_var$day_length <- with(getSunlightTimes(date = seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day"), lat = 54.3233, lon = 10.1228, keep = c("sunrise", "sunset")), as.numeric(difftime(sunset, sunrise, units = "hours")))[permanent_var$day]
permanent_var <- permanent_var[rownames(permanent_rar), ]

permanent_tax <- read.csv(file = "https://raw.githubusercontent.com/gbonthond/Seasonalilty_seaweed_holobiont/refs/heads/main/seasonality_tax.csv", header = T, row.names = 1);dim(permanent_tax)
permanent_ko  <- round(read.csv(file = "https://raw.githubusercontent.com/gbonthond/Seasonalilty_seaweed_holobiont/refs/heads/main/seasonality_ko.csv", header = T, row.names = 1));dim(permanent_ko)
permanent_ko  <- permanent_ko[, colSums(permanent_ko) > 0];dim(permanent_ko)
permanent_ko  <- permanent_ko[, order(colSums(permanent_ko), decreasing = T)];dim(permanent_ko)
#permanent_ko_rar <- rarefaction(permanent_ko, sample = 1000, replicate = 100);dim(permanent_ko_rar)
#write.csv(permanent_rar, file = "C:/Users/Bonthond/Documents/GitHub/Seasonalilty_seaweed_holobiont/permanent_ko_rar.csv");dim(permanent_rar)
permanent_rar <- read.csv(file = "https://raw.githubusercontent.com/gbonthond/Seasonalilty_seaweed_holobiont/refs/heads/main/permanent_ko_rar.csv", header = T, row.names = 1);dim(permanent_ko_rar)

# Reduced to 95% cumulative most abundant OTUs of at least 0.25 prevelance
permanent_tax_sub <- permanent_tax[cumsum(permanent_tax$abundance) < .95 & permanent_tax$occupancy > .25, ];dim(permanent_tax_sub)
permanent_otu_sub <- permanent_otu[, rownames(permanent_tax_sub)];dim(permanent_otu_sub)
permanent_ko_sub <- permanent_ko[, cumsum(colSums(permanent_ko))/sum(colSums(permanent_ko)) < 0.5];dim(permanent_ko_sub)

### subset data for "seasonal core" to sample_type == "alga" and t6, t1 (winter) and t3, t4 (summer) only
season_var    <- droplevels(permanent_var[permanent_var$sample_type == "alga" & permanent_var$timepoint != "t2" & permanent_var$timepoint != "t5", ]);dim(season_var)
season_otu    <- permanent_otu_sub[rownames(season_var), colSums(permanent_otu_sub[rownames(permanent_var), ]) > 0];dim(season_otu)
season_ko     <- permanent_ko_sub[rownames(season_var), ];dim(season_ko)

### create new variable season with two levels (summer and winter)
season_var$season <- with(season_var, ifelse(timepoint == "t6", "winter", ifelse(timepoint == "t1", "winter", "summer")))
permanent_var$otu_PIE  <- calc_PIE(permanent_otu)
permanent_var$otu_rich <- calc_comm_div(permanent_otu, index = "S")[1:262, "value"]
permanent_var$ko_PIE   <- calc_PIE(round(permanent_ko))
permanent_var$ko_rich  <- calc_comm_div(permanent_ko, index = "S")[1:262, "value"]
alga_var <- droplevels(permanent_var[permanent_var$sample_type == "alga", ]);dim(alga_var)
alga_rar <- droplevels(permanent_rar[permanent_var$sample_type == "alga", ])
alga_rar <- alga_rar[,colSums(alga_rar) >0 ];dim(alga_rar)
alga_ko <- droplevels(permanent_ko[permanent_var$sample_type == "alga", ]);dim(alga_ko)
alga_ko <- alga_ko[, colSums(alga_ko) > 0 ];dim(alga_ko)
alga_ko_rar <- droplevels(permanent_ko_rar[permanent_var$sample_type == "alga", ]);
alga_ko_rar <- alga_ko_rar[, colSums(alga_ko_rar) >0 ];dim(alga_ko_rar)

##########################


#### DIVERSITY ####
glm_otu_rich   <- glm(otu_rich ~ log(seq_depth) + timepoint * year * pop, data = alga_var, 
                      family = gaussian(link = "identity"), na.action = na.fail)
qqPlot(resid(glm_otu_rich), pch = 19, col.lines = 1)
Anova(glm_otu_rich)
glm_otu_rich_pp <- emmeans(glm_otu_rich, pairwise ~ timepoint);glm_otu_rich_pp
glm_otu_rich_effect_general <- data.frame(effect(glm_otu_rich, term = "timepoint"))
glm_otu_rich_effect_pop <- data.frame(effect(glm_otu_rich, term = "timepoint:pop"))
glm_otu_rich_pseudo_R2 <- 1-(glm_otu_rich$deviance/glm_otu_rich$null.deviance);glm_otu_rich_pseudo_R2

# plots
plot_otu_rich <- ggplot(glm_otu_rich_effect_pop, aes(x = timepoint, y = fit)) +
  geom_point(size = 5, stroke = 2, aes(fill = pop), shape = 21, position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = pop), width = 0, size = 1, color = "black",position = position_dodge(width = .5)) +
  geom_line(data = glm_otu_rich_effect_general, aes(x = timepoint, y = fit, group = 1), size = 1, color = "black", linetype = "dashed") + 
  geom_errorbar(data = glm_otu_rich_effect_general, aes(ymin = lower, ymax = upper), width = 0, size = 3, color = "black") +  
  geom_point(data = glm_otu_rich_effect_general, size = 10, stroke = .1, aes(x = timepoint, y = fit, color = timepoint)) +
  ylim(0, 3000) + 
  ylab("OTUs") +
  scale_color_manual(values = c("#0670a4","#84af2b","#ad2b03","#e38434","#fbc005","#77c6e3")) +
  scale_fill_manual(values = c("white", "black")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", colour = "black", linewidth = .1),
        axis.ticks.length = unit(.1, "mm"),
        axis.ticks.x = element_line(colour = "#333333", linewidth = .1),
        axis.ticks.y = element_line(colour = "#333333", linewidth = .1));plot_otu_rich

glm_otu_PIE   <- glm(logit(otu_PIE) ~ timepoint * year * pop, data = alga_var, 
                      family = gaussian(link = "identity"), na.action = na.fail)
qqPlot(resid(glm_otu_PIE), pch = 19, col.lines = 1)
Anova(glm_otu_PIE)
glm_otu_PIE_pp <- emmeans(glm_otu_PIE, pairwise ~ timepoint);glm_otu_PIE_pp
glm_otu_PIE_effect_general <- data.frame(effect(glm_otu_PIE, term = "timepoint"))
glm_otu_PIE_effect_pop <- data.frame(effect(glm_otu_PIE, term = "timepoint:pop"))
glm_otu_PIE_pseudo_R2 <- 1-(glm_otu_PIE$deviance/glm_otu_PIE$null.deviance);glm_otu_PIE_pseudo_R2

# plots
plot_otu_PIE <- ggplot(glm_otu_PIE_effect_pop, aes(x = timepoint, y = fit)) +
  geom_point(size = 5, stroke = 2, aes(fill = pop), shape = 21, position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = pop), width = 0, size = 1, color = "black",position = position_dodge(width = .5)) +
  geom_line(data = glm_otu_PIE_effect_general, aes(x = timepoint, y = fit, group = 1), size = 1, color = "black", linetype = "dashed") + 
  geom_errorbar(data = glm_otu_PIE_effect_general, aes(ymin = lower, ymax = upper), width = 0, size = 3, color = "black",position = position_dodge(width = 0.5)) +  
  geom_point(data = glm_otu_PIE_effect_general, size = 10, stroke = .1, aes(x = timepoint, y = fit, color = timepoint)) +
  ylim(2, 5) + 
  ylab("OTU logit PIE") +
  scale_color_manual(values = c("#0670a4","#84af2b","#ad2b03","#e38434","#fbc005","#77c6e3")) +
  scale_fill_manual(values = c("white", "black")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", colour = "black", linewidth = .1),
        axis.ticks.length = unit(.1, "mm"),
        axis.ticks.x = element_line(colour = "#333333", linewidth = .1),
        axis.ticks.y = element_line(colour = "#333333", linewidth = .1));plot_otu_PIE

glm_ko_rich   <- glm(ko_rich ~ log(seq_depth) + timepoint * year * pop, data = alga_var, 
                      family = gaussian(link = "identity"), na.action = na.fail)
qqPlot(resid(glm_ko_rich), pch = 19, col.lines = 1)
Anova(glm_ko_rich)
glm_ko_rich_pp <- emmeans(glm_ko_rich, pairwise ~ timepoint);glm_ko_rich_pp
glm_ko_rich_effect_general <- data.frame(effect(glm_ko_rich, term = "timepoint"))
glm_ko_rich_effect_pop <- data.frame(effect(glm_ko_rich, term = "timepoint:pop"))
glm_ko_rich_pseudo_R2 <- 1-(glm_ko_rich$deviance/glm_ko_rich$null.deviance);glm_ko_rich_pseudo_R2

# plots
plot_ko_rich <- ggplot(glm_ko_rich_effect_pop, aes(x = timepoint, y = fit)) +
  geom_point(size = 5, stroke = 2, aes(fill = pop), shape = 21, position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = pop), width = 0, size = 1, color = "black",position = position_dodge(width = .5)) +
  geom_line(data = glm_ko_rich_effect_general, aes(x = timepoint, y = fit, group = 1), size = 1, color = "black", linetype = "dashed") + 
  geom_errorbar(data = glm_ko_rich_effect_general, aes(ymin = lower, ymax = upper), width = 0, size = 3, color = "black",position = position_dodge(width = 0.5)) +  
  geom_point(data = glm_ko_rich_effect_general, size = 10, stroke = .1, aes(x = timepoint, y = fit, color = timepoint)) +
  ylim(6000, 6800) + 
  ylab("KOs") +
  scale_color_manual(values = c("#0670a4","#84af2b","#ad2b03","#e38434","#fbc005","#77c6e3")) +
  scale_fill_manual(values = c("white", "black")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", colour = "black", linewidth = .1),
        axis.ticks.length = unit(.1, "mm"),
        axis.ticks.x = element_line(colour = "#333333", linewidth = .1),
        axis.ticks.y = element_line(colour = "#333333", linewidth = .1));plot_ko_rich

glm_ko_PIE   <- glm(logit(ko_PIE) ~ timepoint * year * pop, data = alga_var, 
                      family = gaussian(link = "identity"), na.action = na.fail)
qqPlot(resid(glm_ko_PIE), pch = 19, col.lines = 1)
Anova(glm_ko_PIE)
glm_ko_PIE_pp <- emmeans(glm_ko_PIE, pairwise ~ timepoint);glm_ko_PIE_pp
glm_ko_PIE_effect_general <- data.frame(effect(glm_ko_PIE, term = "timepoint"))
glm_ko_PIE_effect_pop <- data.frame(effect(glm_ko_PIE, term = "timepoint:pop"))
glm_ko_PIE_pseudo_R2 <- 1-(glm_ko_PIE$deviance/glm_ko_PIE$null.deviance);glm_ko_PIE_pseudo_R2

# plots
plot_ko_PIE <- ggplot(glm_ko_PIE_effect_pop, aes(x = timepoint, y = fit)) +
  geom_point(size = 5, stroke = 2, aes(fill = pop), shape = 21, position = position_dodge(width = .5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = pop), width = 0, size = 1, color = "black",position = position_dodge(width = .5)) +
  geom_line(data = glm_ko_PIE_effect_general, aes(x = timepoint, y = fit, group = 1), size = 1, color = "black", linetype = "dashed") + 
  geom_errorbar(data = glm_ko_PIE_effect_general, aes(ymin = lower, ymax = upper), width = 0, size = 3, color = "black",position = position_dodge(width = 0.5)) +  
  geom_point(data = glm_ko_PIE_effect_general, size = 10, stroke = .1, aes(x = timepoint, y = fit, color = timepoint)) +
  ylim(7.2, 7.45) + 
  ylab("KO logit PIE") +
  scale_color_manual(values = c("#0670a4","#84af2b","#ad2b03","#e38434","#fbc005","#77c6e3")) +
  scale_fill_manual(values = c("white", "black")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", colour = "black", linewidth = .1),
        axis.ticks.length = unit(.1, "mm"),
        axis.ticks.x = element_line(colour = "#333333", linewidth = .1),
        axis.ticks.y = element_line(colour = "#333333", linewidth = .1));plot_ko_PIE

ggarrange(plot_otu_rich, plot_otu_PIE, plot_ko_rich, plot_ko_PIE)

##################


##### mGLMs permanent CORES ####

### for each level in sample_type an mGLMs as reference level
#permanent_var$sample_type <- relevel(permanent_var$sample_type, ref = "alga")
#permanent_alg_otu_mglm <- manyglm(mvabund(permanent_otu_sub) ~ offset(log(seq_depth)) + sample_type + timepoint + year, family = "negative.binomial", data = permanent_var)
#save(permanent_alg_otu_mglm, file = "permanent_alg_otu_mglm.Rdata") # load saved model
load(file = "permanent_alg_otu_mglm.Rdata") # load saved model

#permanent_var$sample_type <- relevel(permanent_var$sample_type, ref = "sediment")
#permanent_sed_otu_mglm <- manyglm(mvabund(permanent_otu_sub) ~ offset(log(seq_depth)) + sample_type + timepoint + year, family = "negative.binomial", data = permanent_var)
#save(permanent_sed_otu_mglm, file = "permanent_sed_otu_mglm.Rdata") # load saved model
load(file = "permanent_sed_otu_mglm.Rdata") # load saved model

#permanent_var$sample_type <- relevel(permanent_var$sample_type, ref = "water")
#permanent_wat_otu_mglm <- manyglm(mvabund(permanent_otu_sub) ~ offset(log(seq_depth)) + sample_type + timepoint + year, family = "negative.binomial", data = permanent_var)
#save(permanent_wat_otu_mglm, file = "permanent_wat_otu_mglm.Rdata") # load saved model
load(file = "permanent_wat_otu_mglm.Rdata") # load saved model

#permanent_alg_otu_mglm_s <- summary.manyglm(permanent_alg_otu_mglm, nBoot = 999, block = permanent_alg_otu_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#save(permanent_alg_otu_mglm_s, file = "permanent_alg_otu_mglm_s.Rdata") # load saved model
load(file = "permanent_alg_otu_mglm_s.Rdata")

#permanent_sed_otu_mglm_s <- summary.manyglm(permanent_sed_otu_mglm, nBoot = 999, block = permanent_sed_otu_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#save(permanent_sed_otu_mglm_s, file = "permanent_sed_otu_mglm_s.Rdata") # load saved model
load(file = "permanent_sed_otu_mglm_s.Rdata")

#permanent_wat_otu_mglm_s <- summary.manyglm(permanent_wat_otu_mglm, nBoot = 999, block = permanent_wat_otu_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#save(permanent_wat_otu_mglm_s, file = "permanent_wat_otu_mglm_s.Rdata") # load saved model
load(file = "permanent_wat_otu_mglm_s.Rdata")

permanent_alg_core <- permanent.core(permanent_alg_otu_mglm, permanent_alg_otu_mglm_s, p = 0.01, coef = 0, pooled = F, p.adjust = "fdr");dim(permanent_alg_core)
permanent_sed_core <- permanent.core(permanent_sed_otu_mglm, permanent_sed_otu_mglm_s, p = 0.01, coef = 0, pooled = F, p.adjust = "fdr");dim(permanent_sed_core)
permanent_wat_core <- permanent.core(permanent_wat_otu_mglm, permanent_wat_otu_mglm_s, p = 0.01, coef = 0, pooled = F, p.adjust = "fdr");dim(permanent_wat_core)
permanent_alg_core_otu <- permanent.core(permanent_alg_otu_mglm, permanent_alg_otu_mglm_s, p = 0.01, coef = 0, pooled = T, p.adjust = "fdr");dim(permanent_alg_core_otu)
permanent_sed_core_otu <- permanent.core(permanent_sed_otu_mglm, permanent_sed_otu_mglm_s, p = 0.01, coef = 0, pooled = T, p.adjust = "fdr");dim(permanent_sed_core_otu)
permanent_wat_core_otu <- permanent.core(permanent_wat_otu_mglm, permanent_wat_otu_mglm_s, p = 0.01, coef = 0, pooled = T, p.adjust = "fdr");dim(permanent_wat_core_otu)

# plot coefficients and corresponding confidence intervals
{

  # aesthetics 
  {
    theme <- theme(legend.title = element_blank(),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "transparent"),
                   panel.border = element_rect(fill = "transparent",colour = "black",linewidth = 1),
                   axis.ticks.length = unit(1, "mm"),
                   axis.ticks.x = element_line(colour="#333333", linewidth = 1),
                   axis.ticks.y = element_line(colour="#333333", linewidth = 1))
    }
  
  gg_permanent_alg_otu <- ggplot(permanent_alg_core_otu[1:25, ],
                                 aes(x = reorder(taxon, -averaged_coefficient, y = -averaged_coefficient))) +
    geom_errorbar(aes(ymin = -pooled_lower, ymax = -pooled_upper), col = "darkgreen", width = 0, size = 1.25) +
    geom_point(aes(y = -averaged_coefficient), size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red",size = 0.75) + theme +
    ggtitle("permanent algal core - fold change + 99%CIs");gg_permanent_alg_otu

  gg_permanent_sed_otu <- ggplot(permanent_sed_core_otu[1:25, ],
                                 aes(x = reorder(taxon, -averaged_coefficient, y = -averaged_coefficient))) +
    geom_errorbar(aes(ymin = -pooled_lower, ymax = -pooled_upper), col = "orange", width = 0, size = 1.25) +
    geom_point(aes(y = -averaged_coefficient), size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red",size = 0.75) + theme +
    ggtitle("permanent sediment core - fold change + 99%CIs");gg_permanent_sed_otu
  
  gg_permanent_wat_otu <- ggplot(permanent_wat_core_otu[1:25, ],
                                 aes(x = reorder(taxon, -averaged_coefficient, y = -averaged_coefficient))) +
    geom_errorbar(aes(ymin = -pooled_lower, ymax = -pooled_upper), col = "darkcyan", width = 0, size = 1.25) +
    geom_point(aes(y = -averaged_coefficient), size = 2) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red",size = 0.75) + theme +
    ggtitle("permanent water core - fold change + 99%CIs");gg_permanent_wat_otu
  }
###############################


##### mGLMs SEASONAL CORES #####
#season_otu_mglm <- manyglm(mvabund(season_otu) ~ offset(log(seq_depth)) + season + year, family = "negative.binomial", data = season_var)
#save(season_otu_mglm, file = "season_otu_mglm.Rdata") # save model to call it next time directly
load(file = "season_otu_mglm.Rdata") # load saved model
#season_otu_mglm_s <- summary.manyglm(season_otu_mglm, nBoot = 999, block = season_otu_mglm$data$pop, test = "wald", resamp = "case", show.time = "all", p.uni = "unadjusted")
#save(season_otu_mglm_s, file = "season_otu_mglm_s.Rdata")
load(file = "season_otu_mglm_s.Rdata")

### GET SEASONAL CORE OTUs 
season_core_otu <- seasonal.core(season_otu_mglm, season_otu_mglm_s, p = 0.01, coef = 0, p.adjust = "fdr");dim(season_core_otu)

# exclude permanent sediment and water cores
season_core_otu <- season_core_otu[rownames(season_core_otu) %in% rownames(permanent_sed_core_otu) == F &
                                   rownames(season_core_otu) %in% rownames(permanent_wat_core_otu) == F, ];dim(season_core_otu)

# plot coefficients and corresponding confidence intervals
# aesthetics 
gg_season_otu_sw <-  ggplot(season_core_otu[1:50, ],
                            aes(x = reorder(taxon, -coefficient), y = coefficient)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, col = core), alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = coefficient), size = 2) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + theme +
  ggtitle("fold change with respect to season + 99%CIs");gg_season_otu_sw

##############################


##### CORES - FIGURE 5D-F & TABLE S6 #####
permanent_alg_core_otu$label <- paste(permanent_alg_core_otu$taxon,
                                      permanent_tax_sub[permanent_alg_core_otu$taxon, "genus"])

### figure 5D
gg_permanent_otu_sw <- ggplot(permanent_alg_core_otu[1:25, ],
                         aes(x = reorder(label, -averaged_coefficient), y = -averaged_coefficient)) +
  geom_errorbar(aes(ymin = -pooled_lower, ymax = -pooled_upper), col = "darkgreen", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = -averaged_coefficient), size = 3, stroke = 2, shape = 21, colour = "black", fill = "white") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + theme +
  ggtitle("permanent");gg_permanent_otu_sw

season_core_otu$label <- paste(season_core_otu$taxon, permanent_tax_sub[season_core_otu$taxon, "genus"])
  
### figure 5E
gg_summer_otu_sw <-  ggplot(season_core_otu[season_core_otu$core == "summer", ][1:25, ],
                            aes(x = reorder(label, -coefficient), y = -coefficient)) + 
  geom_errorbar(aes(ymin = -lower, ymax = -upper), col = "darkred", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = -coefficient), size = 3, stroke = 2, shape = 21, colour = "black", fill = "white") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + theme +
  ggtitle("Summer");gg_summer_otu_sw

### figure 5F
gg_winter_otu_sw <-  ggplot(season_core_otu[season_core_otu$core == "winter", ][1:25, ],
                            aes(x = reorder(label, coefficient), y = coefficient)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), col = "darkblue", alpha = .5, width = 0, linewidth = 2) +
  geom_point(aes(y = coefficient), size = 3, stroke = 2, shape = 21, colour = "black", fill = "white") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.75) + theme +
  ggtitle("winter");gg_winter_otu_sw

### mGLM input TABLE S6
core_table <- permanent_tax
core_table$otu <- rownames(core_table)
core_table[, "occupancy_permanent"] <- ifelse((colSums(permanent_otu[permanent_var$sample_type == "alga", rownames(core_table)] != 0)/nrow(permanent_otu[permanent_var$sample_type == "alga", rownames(core_table)]) == 1) == TRUE, "epi", "")
core_table[, "occupancy_season"]    <- ifelse(!rownames(core_table) %in% colnames(permanent_otu_sub),"",
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


########################################


##### NMDS #####

# nmds OTU
#alga_nmds1 <- metaMDS(alga_rar, distance = "bray", trymax = 200, try = 100, autotransform = F)
#alga_nmds2 <- metaMDS(alga_rar, distance = "bray", trymax = 200, try = 100, autotransform = F, previous.best = alga_nmds1)
#save(alga_nmds2, file= "alga_nmds2.Rdata")
load("alga_nmds2.Rdata")
alga_nmds2

alga_env <- envfit(alga_nmds2, alga_var[, c("pH", "salinity", "temperature", "day_length")], permutations = 9999, na.rm = TRUE)

# Basic NMDS plot
custom_colors <- c("#036aaa", "#85af31", "#ad2c03", "#ec8935", "#fcc601", "#79c7db")

# Plot with custom colors for timepoint
plot(alga_nmds2$points, 
     col = custom_colors[as.factor(alga_var$timepoint)],
     pch = ifelse(alga_var$pop == "nor" & alga_var$year == "year1", 22, 
           ifelse(alga_var$pop == "nor" & alga_var$year == "year2", 2, 
           ifelse(alga_var$pop == "nor" & alga_var$year == "year3", 1, 
           ifelse(alga_var$pop == "hei" & alga_var$year == "year1", 15, 
           ifelse(alga_var$pop == "hei" & alga_var$year == "year2", 17, 19))))),
     cex = 1.5, xlab = "NMDS1", ylab = "NMDS2", main = "OTU")

legend("topright", legend = c("Year 1 (nor)", "Year 2 (nor)", "Year 3 (nor)", "Year 1 (hei)", "Year 2 (hei)", "Year 3 (hei)"), pch = c(22, 2, 1, 15, 17, 19), pt.cex = 1.5, title = "Year and Pop", cex = 0.8, bty = "n")
legend("bottomleft", legend = c("t1", "t2", "t3", "t4", "t5", "t6"), fill = custom_colors[as.factor(c("t1", "t2", "t3", "t4", "t5", "t6"))], title = "season", cex = 1, bty = "n")  # Remove border around legend

# envfit vectors
env_scores <- scores(alga_env, display = "vectors") * .5
arrows(1.5, 0.5, 1.5 + env_scores[, "NMDS1"], 0.5 + env_scores[, "NMDS2"], col = "black", lwd = 2, length = 0.1)
text(1.5 + env_scores[, "NMDS1"], 0.5 + env_scores[, "NMDS2"], labels = rownames(env_scores), col = "black", cex = 1.2, pos = 4)

#add ellipse
{
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t1 nor", col = "#036aaa", lwd = 2, kind = "se", conf = 0.95)  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t2 nor", col = "#85af31", lwd = 2, kind = "se", conf = 0.95)  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t3 nor", col = "#ad2c03", lwd = 2, kind = "se", conf = 0.95)  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t4 nor", col = "#ec8935", lwd = 2, kind = "se", conf = 0.95)  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t5 nor", col = "#fcc601", lwd = 2, kind = "se", conf = 0.95)  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t6 nor", col = "#79c7db", lwd = 2, kind = "se", conf = 0.95)  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t1 hei", col = "#036aaa", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t2 hei", col = "#85af31", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t3 hei", col = "#ad2c03", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t4 hei", col = "#ec8935", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t5 hei", col = "#fcc601", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_nmds2, groups = paste(alga_var$timepoint, alga_var$pop), show.groups = "t6 hei", col = "#79c7db", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
}

#trajectories
group_names <- unique(paste(alga_var$year, alga_var$timepoint, alga_var$pop))
centroids <- sapply(group_names, function(group) {
  rows <- which(paste(alga_var$year, alga_var$timepoint, alga_var$pop) == group)
  colMeans(alga_nmds2$points[rows, ])}, simplify = TRUE)
centroid_matrix <- t(centroids)
rownames(centroid_matrix) <- group_names

trajectory_order_nor <- c("year1 t1 nor", "year1 t2 nor", "year1 t3 nor", "year1 t4 nor", "year1 t5 nor", "year1 t6 nor", "year2 t1 nor", "year2 t2 nor", "year2 t3 nor", "year2 t4 nor", "year2 t5 nor", "year2 t6 nor", "year3 t1 nor", "year3 t2 nor", "year3 t3 nor", "year3 t4 nor", "year3 t5 nor", "year3 t6 nor")
trajectory_order_hei <- c("year1 t1 hei", "year1 t2 hei", "year1 t3 hei", "year1 t4 hei", "year1 t5 hei", "year1 t6 hei", "year2 t1 hei", "year2 t2 hei", "year2 t3 hei", "year2 t4 hei", "year2 t5 hei", "year2 t6 hei", "year3 t1 hei", "year3 t2 hei", "year3 t3 hei", "year3 t4 hei", "year3 t5 hei", "year3 t6 hei")

trajectory_centroids_nor <- centroid_matrix[trajectory_order_nor, ]
lines(trajectory_centroids_nor[1:6, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 4)
lines(trajectory_centroids_nor[6:12, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 2)
lines(trajectory_centroids_nor[12:18, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 1)

trajectory_centroids_hei <- centroid_matrix[trajectory_order_hei, ]
lines(trajectory_centroids_hei[1:6, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 4)
lines(trajectory_centroids_hei[6:12, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 2)
lines(trajectory_centroids_hei[12:18, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 1)


# nmds KO
#alga_ko_nmds1 <- metaMDS(alga_ko_rar, distance = "bray", trymax = 200, try = 100, autotransform = F)
#alga_ko_nmds2 <- metaMDS(alga_ko_rar, distance = "bray", trymax = 200, try = 100, autotransform = F, previous.best = alga_ko_nmds1)
#save(alga_ko_nmds2, file= "alga_ko_nmds2.Rdata")
load("alga_ko_nmds2.Rdata")
alga_ko_nmds2
alga_ko_env <- envfit(alga_ko_nmds2, alga_var[, c("pH", "salinity", "temperature", "day_length")], permutations = 9999, na.rm = TRUE)

custom_colors <- c("#036aaa", "#85af31", "#ad2c03", "#ec8935", "#fcc601", "#79c7db")

# Plot with custom colors for timepoint
plot(alga_ko_nmds2$points, 
     col = custom_colors[as.factor(alga_var$timepoint)], 
     pch = ifelse(alga_var$pop == "nor" & alga_var$year == "year1", 22, 
           ifelse(alga_var$pop == "nor" & alga_var$year == "year2", 2, 
           ifelse(alga_var$pop == "nor" & alga_var$year == "year3", 1, 
           ifelse(alga_var$pop == "hei" & alga_var$year == "year1", 15, 
           ifelse(alga_var$pop == "hei" & alga_var$year == "year2", 17, 19))))), 
   #  xlim = c(-.2, .2), ylim = c(-.2, .15),
     cex = 1.5, xlab = "NMDS1", ylab = "NMDS2", main = "OTU")

legend("bottomright", legend = c("Year 1 (nor)", "Year 2 (nor)", "Year 3 (nor)", "Year 1 (hei)", "Year 2 (hei)", "Year 3 (hei)"), pch = c(22, 2, 1, 15, 17, 19), pt.cex = 1.5, title = "Year and Pop", cex = 0.8, bty = "n")
legend("bottomleft",  legend = c("t1", "t2", "t3", "t4", "t5", "t6"), fill = custom_colors[as.factor(c("t1", "t2", "t3", "t4", "t5", "t6"))], title = "season", cex = 1, bty = "n")

# envfit vectors
env_ko_scores <- scores(alga_ko_env, display = "vectors") * .25  # Automatically scales arrows appropriately
arrows(-0.3, -0.2, env_ko_scores[, "NMDS1"] - 0.3, env_ko_scores[, "NMDS2"] - 0.2, col = "black", lwd = 2, length = 0.1)
text(env_ko_scores[, "NMDS1"] - 0.3, env_ko_scores[, "NMDS2"] - 0.2, labels = rownames(env_ko_scores), col = "black", cex = 1.2, pos = 4)


#add ellipses
{
  ordiellipse(alga_ko_nmds2, groups = alga_var$timepoint, show.groups = "t1", col = "#036aaa", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_ko_nmds2, groups = alga_var$timepoint, show.groups = "t2", col = "#85af31", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_ko_nmds2, groups = alga_var$timepoint, show.groups = "t3", col = "#ad2c03", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_ko_nmds2, groups = alga_var$timepoint, show.groups = "t4", col = "#ec8935", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_ko_nmds2, groups = alga_var$timepoint, show.groups = "t5", col = "#fcc601", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  ordiellipse(alga_ko_nmds2, groups = alga_var$timepoint, show.groups = "t6", col = "#79c7db", lwd = 2, kind = "se", conf = 0.95, draw = "polygon")  
  }

#trajectories
group_names_ko <- unique(paste(alga_var$year, alga_var$timepoint))
centroids_ko <- sapply(group_names_ko, function(group) {
  rows <- which(paste(alga_var$year, alga_var$timepoint) == group)
  colMeans(alga_ko_nmds2$points[rows, ])
}, simplify = TRUE)
centroid_matrix_ko <- t(centroids_ko)
rownames(centroid_matrix_ko) <- group_names_ko

trajectory_order_ko <- c("year1 t1", "year1 t2", "year1 t3", "year1 t4", "year1 t5", "year1 t6", "year2 t1", "year2 t2", "year2 t3", "year2 t4", "year2 t5", "year2 t6", "year3 t1", "year3 t2", "year3 t3", "year3 t4", "year3 t5", "year3 t6")

trajectory_centroids_ko <- centroid_matrix_ko[trajectory_order_ko, ]
lines(trajectory_centroids_ko[1:6, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 4)
lines(trajectory_centroids_ko[6:12, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 2)
lines(trajectory_centroids_ko[12:18, ], col = "black", lwd = 2, type = "o", pch = 1, lty = 1)

##############

