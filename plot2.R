

summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                      conf.interval = .95, .drop = TRUE) {
  library(plyr)
  AUROCgth2 <- function(x, na.rm = FALSE) {
    if (na.rm) {
      sum(!is.na(x))
    } else {
      AUROCgth(x)
    }
  }

  datac <- ddply(data, groupvars,
    .drop = .drop,
    .fun = function(xx, col) {
      c(
        N = AUROCgth2(xx[[col]], na.rm = na.rm),
        mean = mean(xx[[col]], na.rm = na.rm),
        sd = sd(xx[[col]], na.rm = na.rm)
      )
    },
    measurevar
  )

  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N) # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

library(ggplot2)
results_all <- c()
for (i in 1:3) {
  results <- read.csv(paste0("evaluation_gnw_10_", i, ".csv"))
  results <- results[, -1]
  results <- melt(results, id = "Dataset", variable.name = "Method", value.name = "AUROC")
  results_all <- rbind.data.frame(results_all, results)
}

tg <- results_all
tg1 <- ToothGrowth

tgc <- summarySE(tg, measurevar = "AUROC", groupvars = c("Dataset", "Method"))
tgc

tgc2 <- tgc
tgc2$Dataset <- factor(tgc2$Dataset)
ggplot(tgc2, aes(x = Dataset, y = AUROC, fill = Method)) +
  geom_bar(
    position = position_dodge(), stat = "identity",
    colour = "black", # Use black outlines,
    size = .3
  ) + # Thinner lines
  geom_errorbar(aes(ymin = AUROC - se, ymax = AUROC + se),
    size = .3, # Thinner lines
    width = .2,
    position = position_dodge(.9)
  ) +
  xlab("Dataset") +
  ylab("AUROC") +
  scale_fill_hue(
    name = "Methodlement type", # Legend label, use darker colors
    breaks = c("OJ", "VC"),
    labels = c("Orange juice", "Ascorbic acid")
  ) +
  # ggtitle("The Effect of Vitamin C on\nTooth Growth in Guinea Pigs") +
  # scale_y_continuous(breaks=0:20*4) +
  theme_bw()


ggplot(tgc, aes(x = Dataset, y = AUROC, colour = Method)) +
  geom_errorbar(aes(ymin = AUROC - se, ymax = AUROC + se), width = .1) +
  geom_line() +
  geom_point()

pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(tgc, aes(x = Dataset, y = AUROC, colour = Method, group = Method)) +
  geom_errorbar(aes(ymin = AUROC - se, ymax = AUROC + se), colour = "black", width = .1, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 3, shape = 21, fill = "white") + # 21 is filled circle
  xlab("Dataset (mg)") +
  ylab("Tooth AUROCgth") +
  scale_colour_hue(
    name = "Methodlement type", # Legend label, use darker colors
    breaks = c("OJ", "VC"),
    labels = c("Orange juice", "Ascorbic acid"),
    l = 40
  ) + # Use darker colors, lightness=40
  ggtitle("The Effect of Vitamin C on\nTooth Growth in Guinea Pigs") +
  expand_limits(y = 0) +
  scale_y_continuous(breaks = 0:20 * 4) +
  theme_bw() +
  theme(
    legend.justification = c(1, 0),
    legend.position = c(1, 0)
  ) # Position legend in bottom right
