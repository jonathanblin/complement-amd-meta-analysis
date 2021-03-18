
library(meta)
library(metafor)
library(dplyr)
library(dmetar)
library(pryr)

data <- read.table("data/tidy_data/meta-analysis.csv", header = TRUE, sep = ",")

data$study <- as.factor(data$study)

amd <- metacont(N_amd,
                mean_amd,
                SD_amd,
                N_ctrl,
                mean_ctrl,
                SD_ctrl,
                data = data,
                studlab=paste(study),
                comb.fixed = FALSE,
                comb.random = TRUE,
                method.tau = "REML",
                hakn = TRUE,
                prediction = TRUE,
                sm = "SMD")

ga <- metacont(N_GA,
               mean_GA,
               SD_GA,
               N_ctrl,
               mean_ctrl,
               SD_ctrl,
               data = data,
               studlab=paste(study),
               comb.fixed = FALSE,
               comb.random = TRUE,
               method.tau = "REML",
               hakn = TRUE,
               prediction = TRUE,
               sm = "SMD")

ga_sens <- metacont(N_GA,
               mean_GA_sens,
               SD_GA_sens,
               N_ctrl,
               mean_ctrl,
               SD_ctrl,
               data = data,
               studlab=paste(study),
               comb.fixed = FALSE,
               comb.random = TRUE,
               method.tau = "REML",
               hakn = TRUE,
               prediction = TRUE,
               sm = "SMD")

# I did above to calculate Hedge's g (SMD) for each effect size and combined as below.

data2 <- as.data.frame(cbind(data, amd$TE, amd$seTE, ga$TE, ga$seTE, ga_sens$TE, ga_sens$seTE))
data2$study <- as.factor(data2$study)
data2$factor <- as.factor(data2$factor)
data2$factor_cat <- as.factor(data2$factor_cat)
data2$cat <- as.factor(data2$cat)

data3 <- mutate(data2, amd_V = amd$seTE^2) %>%
        mutate(data2, ga_V = ga$seTE^2) %>%
        mutate(data2, ga_sens_V = ga_sens$seTE^2)

# FOREST PLOT FOR ACTIVATED PRODUCTS IN NONEXUDATIVE AMD VS CONTROL

data3_amd_activation <- filter(data3, data3$cat == "amd", data3$factor_cat == "activation")

amd_activation <- rma.mv(yi = data3_amd_activation$`amd$TE`,
                         V = amd_V,
                         random = ~ 1 | study/factor,
                         tdist = TRUE,
                         data = data3_amd_activation,
                         method = "REML",
                         slab = paste(factor))
summary(amd_activation)

sav <- confint(amd_activation)
W <- diag(1/data3_amd_activation$amd_V)
X <- model.matrix(amd_activation)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * amd_activation$sigma2 / (amd_activation$sigma2 + (amd_activation$k-amd_activation$p)/sum(diag(P)))  ### study and factor I^2
100 * sav[[1]]$random[1,2:3] / (sav[[1]]$random[1,2:3] + (amd_activation$k-amd_activation$p)/sum(diag(P))) ### CI for study-level I^2
100 * sav[[2]]$random[1,2:3] / (sav[[2]]$random[1,2:3] + (amd_activation$k-amd_activation$p)/sum(diag(P))) ### CI for factor-level I^2

plot.amd_activation %<a-%{
forest(amd_activation,
       addcred = TRUE,
       showweights = FALSE,
       ilab = cbind(data3_amd_activation$N_ctrl, data3_amd_activation$N_amd),
       ilab.xpos = c(-1.45,-.75),
       header = c("Study (Factor)", "(95% CI)"),
       xlim = c(-6,4),
       xlab = "Standardized Mean Difference (SMD)",
       pch = 15,
       psize = 1.5,
       mlab = "",
       col = "blue",
       annosym=c(" (", " to ", ")"))
text(c(-1.45,-0.75), 7.75, c("Ctrl", "AMD"), font = 2, cex = 0.8)
text(-1, 8.5, "Sample Size", font = 2)
text(1.95, 8, "SMD", font = 2)
text(c(-6, -5.4), -1, pos = 4, c(expression(I[study]^2), "= 41% (95% CI: 0-98%)"))
text(-5.75, 9.1, "A", font = 2, cex = 1.3)
     }

# FOREST PLOT FOR NEG REGULATORS IN NONEXUDATIVE AMD VS CONTROL

data3_amd_neg_reg <- filter(data3, data3$cat == "amd", data3$factor_cat == "neg_reg")

amd_neg_reg <- rma.mv(yi = data3_amd_neg_reg$`amd$TE`,
                      V = data3_amd_neg_reg$amd_V,
                         random = ~ 1 | study/factor,
                         tdist = TRUE,
                         data = data3_amd_neg_reg,
                         method = "REML",
                         slab = paste(factor))
summary(amd_neg_reg)

sav <- confint(amd_neg_reg)
W <- diag(1/data3_amd_neg_reg$amd_V)
X <- model.matrix(amd_neg_reg)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * amd_neg_reg$sigma2 / (amd_neg_reg$sigma2 + (amd_neg_reg$k-amd_neg_reg$p)/sum(diag(P)))  ### study and factor I^2
100 * sav[[1]]$random[1,2:3] / (sav[[1]]$random[1,2:3] + (amd_neg_reg$k-amd_neg_reg$p)/sum(diag(P))) ### CI for study-level I^2
100 * sav[[2]]$random[1,2:3] / (sav[[2]]$random[1,2:3] + (amd_neg_reg$k-amd_neg_reg$p)/sum(diag(P))) ### CI for factor-level I^2

plot.amd_neg_reg %<a-%{
forest(amd_neg_reg,
       addcred = TRUE,
       showweights = FALSE,
       ilab = cbind(data3_amd_neg_reg$N_ctrl, data3_amd_neg_reg$N_amd),
       ilab.xpos = c(-2.5,-2),
       header = c("Study (Factor)", "(95% CI)"),
       xlim = c(-5,2.5),
       xlab = "Standardized Mean Difference (SMD)",
       pch = 15,
       psize = 1.5,
       mlab = "",
       col = "blue",
       annosym=c(" (", " to ", ")"))
text(c(-2.5,-2), 6.75, c("Ctrl", "AMD"), font = 2, cex = 0.8)
text(-2.25, 7.5, "Sample Size", font = 2)
text(0.80, 7.05, "SMD", font = 2)
text(c(-5, -4.55), -1, pos = 4, c(expression(I[study]^2), "= 54% (95% CI: 0-98%)"))
text(-4.8, 8.1, "A", font = 2, cex = 1.3)
        }

# FOREST PLOT FOR ACTIVATED PRODUCTS IN GA VS CONTROL

data3_ga_activation <- filter(data3, data3$cat == "ga", data3$factor_cat == "activation")

ga_activation <- rma.mv(yi = data3_ga_activation$`ga$TE`,
                      V = data3_ga_activation$ga_V,
                      random = ~ 1 | study/factor,
                      tdist = TRUE,
                      data = data3_ga_activation,
                      method = "REML",
                      slab = paste(factor))
summary(ga_activation)

sav <- confint(ga_activation)
W <- diag(1/data3_ga_activation$ga_V)
X <- model.matrix(ga_activation)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * ga_activation$sigma2 / (ga_activation$sigma2 + (ga_activation$k-ga_activation$p)/sum(diag(P)))  ### study and factor I^2
100 * sav[[1]]$random[1,2:3] / (sav[[1]]$random[1,2:3] + (ga_activation$k-ga_activation$p)/sum(diag(P))) ### CI for study-level I^2
100 * sav[[2]]$random[1,2:3] / (sav[[2]]$random[1,2:3] + (ga_activation$k-ga_activation$p)/sum(diag(P))) ### CI for factor-level I^2

plot.ga_activation %<a-%{
forest(ga_activation,
          addcred = TRUE,
          showweights = FALSE,
          ilab = cbind(data3_ga_activation$N_ctrl, data3_ga_activation$N_GA),
          ilab.xpos = c(-1.5,-1),
          header = c("Study (Factor)", "(95% CI)"),
          xlim = c(-4,3),
          xlab = "Standardized Mean Difference (SMD)",
          pch = 15,
          psize = 1.5,
          mlab = "",
          col = "blue",
          annosym=c(" (", " to ", ")"))
text(c(-1.5,-1), 8.75, c("Ctrl", "GA"), font = 2, cex = 0.8)
text(-1.25, 9.5, "Sample Size", font = 2)
text(1.5, 9.05, "SMD", font = 2)
text(c(-4, -3.6), -1, pos = 4, c(expression(I[study]^2), "= 47% (95% CI: 0-99%)"))
text(-3.8, 10.1, "B", font = 2, cex = 1.3)
        }

# FOREST PLOT FOR NEG REGULATORS IN GA VS CONTROL

data3_ga_neg_reg <- filter(data3, data3$cat == "ga", data3$factor_cat == "neg_reg")

ga_neg_reg <- rma.mv(yi = data3_ga_neg_reg$`ga$TE`,
                        V = data3_ga_neg_reg$ga_V,
                        random = ~ 1 | study/factor,
                        tdist = TRUE,
                        data = data3_ga_neg_reg,
                        method = "REML",
                        slab = paste(factor))
summary(ga_neg_reg)

sav <- confint(ga_neg_reg)
W <- diag(1/data3_ga_neg_reg$ga_V)
X <- model.matrix(ga_neg_reg)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * ga_neg_reg$sigma2 / (ga_neg_reg$sigma2 + (ga_neg_reg$k-ga_neg_reg$p)/sum(diag(P)))  ### study and factor I^2
100 * sav[[1]]$random[1,2:3] / (sav[[1]]$random[1,2:3] + (ga_neg_reg$k-ga_neg_reg$p)/sum(diag(P))) ### CI for study-level I^2
100 * sav[[2]]$random[1,2:3] / (sav[[2]]$random[1,2:3] + (ga_neg_reg$k-ga_neg_reg$p)/sum(diag(P))) ### CI for factor-level I^2

plot.ga_neg_reg %<a-%{
forest(ga_neg_reg,
       addcred = TRUE,
       showweights = FALSE,
       ilab = cbind(data3_ga_neg_reg$N_ctrl, data3_ga_neg_reg$N_GA),
       ilab.xpos = c(-1.5,-1),
       header = c("Study (Factor)", "(95% CI)"),
       xlim = c(-4,3),
       xlab = "Standardized Mean Difference (SMD)",
       pch = 15,
       psize = 1.5,
       mlab = "",
       col = "blue",
       annosym=c(" (", " to ", ")"))
text(c(-1.5,-1), 6.75, c("Ctrl", "GA"), font = 2, cex = 0.8)
text(-1.25, 7.5, "Sample Size", font = 2)
text(1.5, 7.05, "SMD", font = 2)
text(c(-4, -3.6), -1, pos = 4, c(expression(I[study]^2), "= 0% (95% CI: 0-93%)"))
text(-3.8, 8.1, "B", font = 2, cex = 1.3)
        }

# Making figures with multiple forest plots
figure2 %<a-%{
split.screen(c(1,2))
screen(1)
plot.amd_activation

screen(2)
plot.ga_activation

close.screen(all=TRUE)
}

pdf("output/figures/figure-2.pdf", height = 5, width = 14)
figure2
dev.off()

figure3 %<a-%{
        split.screen(c(1,2))
        screen(1)
        plot.amd_neg_reg

        screen(2)
        plot.ga_neg_reg

        close.screen(all=TRUE)
}

pdf("output/figures/figure-3.pdf", height = 4.5, width = 14)
figure3
dev.off()


# Sensitivity analysis of 10th/90th percentile as max/min

data3_ga_sens_activation <- filter(data3, data3$cat == "ga", data3$factor_cat == "activation")

ga_sens_activation <- rma.mv(yi = data3_ga_sens_activation$`ga_sens$TE`,
                         V = data3_ga_sens_activation$ga_sens_V,
                         random = ~ 1 | study/factor,
                         tdist = TRUE,
                         data = data3_ga_sens_activation,
                         method = "REML",
                         slab = paste(factor))
summary(ga_sens_activation)
forest(ga_sens_activation)

data3_ga_sens_neg_reg <- filter(data3, data3$cat == "ga", data3$factor_cat == "neg_reg")

ga_sens_neg_reg <- rma.mv(yi = data3_ga_sens_neg_reg$`ga_sens$TE`,
                             V = data3_ga_sens_neg_reg$ga_sens_V,
                             random = ~ 1 | study/factor,
                             tdist = TRUE,
                             data = data3_ga_sens_neg_reg,
                             method = "REML",
                             slab = paste(factor))
summary(ga_sens_neg_reg)
forest(ga_sens_neg_reg)
