# ====================================
# EPISKILLS: probability distributions
# ====================================

library(data.table) # for reading, subsetting, joining data etc
library(RColorBrewer) # for chart colours

# normal distribution with mean 178 and sd 7.6
# --------------------------------------------

png('height_of_men.png', height = 5, width = 5, units = 'in', res = 300)
hist(rnorm(33e6, mean = 178, sd = 7.6), main = 'Height of men in the UK', breaks = seq(13, 23, 0.1) * 10, xlab = 'Height in cm', ylab = 'Frequency')
dev.off()

# standard normal distribution (mean = 0, sd = 1)
# -----------------------------------------------

x <- -500:500 / 100
plot(x, dnorm(x), type = 'l', ylab = 'Density')

png('standard_normal.png', height = 6, width = 6.5, units = 'in', res = 300)
plot(x, dnorm(x), type = 'l', ylab = 'Density', xlab = 'x / sd', main = 'mean = 0, sd = 1 (the "standard normal" distribution)')
abline(v = -4:4, lty = 3, col = 'red')
x1 <- x[x >= -1.96 & x <= 1.96]
polygon(x = c(x1, rev(x1)), y = c(rep(0, length(x1)), rev(dnorm(x1))), col = 'lightblue')
text(0, 0.1, '95% of the distribution')
dev.off()

# distribution quiz
# -----------------

png('distribution_quiz.png', height = 8, width = 8, units = 'in', res = 300)
par(mfrow = c(2, 2))
plot(0:20, dbinom(0:20, 100, 0.05), type = 'l', ylab = 'Density', xlab = 'x', main = 'A\nTrials = 100, probability of success = 5%')
plot(seq(-3, 3, 0.1), dnorm(seq(-3, 3, 0.1)), type = 'l', ylab = 'Density', xlab = 'x', main = 'B\nmean = 0, sd = 1')
plot(10:20, dunif(10:20, 10, 20), type = 'l', xlab = 'x', ylab = 'Density', main = "C\nMin = 10, Max = 20")
plot(0:20, dpois(0:20, lambda = 4.5), type = 'l', ylab = 'Density', xlab = 'x', main = 'D\nmean = 4.5')
dev.off()

# theoretical poisson distribution
# --------------------------------

plot(1:20, dpois(1:20, 5), type = 'l', ylab = 'Density', xlab = 'x')

png('poisson.png', height = 6, width = 6, units = 'in', res = 300)
plot(1, type = 'n', xlim = c(0, 20), ylim = c(0, 0.4), ylab = 'Density', xlab = 'Count', main = 'The Poisson distribution')
mapply(lines, x = list(0:20), y = lapply(1:8, function (x) dpois(0:20, x)), col = 1:8)
ys <- seq(0.2, 0.35, length.out = 8)
segments(10, ys, x1 = 12.5, col = 1:8)
text(13, ys, 1:8)
text(10, 0.38, 'Mean (aka lambda)', adj = 0)
dev.off()

# Deaths per LSOA in 2018
# -----------------------

# compare empirical and theoretical distribution

# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/adhocs/15370deathsbylowerlayersuperoutputarealsoaenglandandwalesmidyearperiods1julyto30june2010to2021
d <- fread("https://raw.githubusercontent.com/danlewer/episkills/refs/heads/main/lsoa_deaths.csv")
d <- d[year == 2018, -'year']
d <- d[, lapply(.SD, sum), by = c('la', 'lsoa')]
d$deaths <- rowSums(d[, -c('la', 'lsoa')])
png('deaths_per_LSOA1.png', height = 6, width = 6, units = 'in', res = 300)
hist(d$deaths, breaks = 0:140, xlab = 'Count', ylab = 'Number of LSOAs', ylim = c(0, 3500), xlim = c(0, 80), main = 'Deaths per LSOA in 2018', border = 'white')
abline(v = mean(d$deaths))
lines(0:140, dpois(0:140, 16) * sum(table(d$deaths)), col = 'red')
dev.off()

# add population and stratify by 65+

# https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareamidyearpopulationestimatesnationalstatistics
p <- fread("https://raw.githubusercontent.com/danlewer/episkills/refs/heads/main/pop_mid_2015.csv")
p$p65 <- p$`F65 and over` + p$`M65 and over`
d$d65 <- d$`Persons 65 to 69` + d$`Persons 70 to 74` + d$`Persons 75 to 79` + d$`Persons 80 to 84` + d$`Persons over 85`
d <- p[, c('lsoa', 'p65')][d, on = 'lsoa']
d <- d[, c('lsoa', 'p65', 'd65')]
d[, rt := d65 / p65 * 100]
png('deaths_per_LSOA2.png', height = 6, width = 6, units = 'in', res = 300)
hist(d$rt, breaks = 100, xlab = 'Rate per 100 person-years', ylab = 'Count of LSOAs', border = 'white', xlim = c(0, 20), main = 'Deaths at age 65+ per 100 person-years in 2018')
abline(v = mean(d$rt, na.rm = T))
lines(0:140, dpois(0:140, 4.55) * sum(table(d$rt)), col = 'red')
dev.off()

# poisson regression to calculate average count per LSOA and rate per 100 person-years

d[, rt := NULL]
d <- d[substr(lsoa, 0, 1) == 'E']

glm(d65 ~ 1, data = d, family = 'poisson')
exp(2.593) # mean count of deaths

glm(d65 ~ 1 + offset(log(p65)), data = d, family = 'poisson')
exp(-3.094) * 100 # rate of deaths per 100 person years

# adding an independent variable (IMD quintile)

imd <- fread("https://raw.githubusercontent.com/danlewer/episkills/refs/heads/main/imd2019.csv")
d <- imd[d, on = 'lsoa']
d[, imd5 := factor(imd10, 1:10, rep(1:5, each = 2))]
d[, imd5 := factor(imd5, 1:5, 5:1)]
d[, imd5 := relevel(imd5, ref = 5)]
d[, imd10 := NULL]

m <- glm(d65 ~ imd5 + offset(log(p65)), data = d, family = 'poisson')
summary(m)

exp(0.421236)

# try using linear model instead of poisson model

m <- lm(d65 / p65 * 100 ~ imd5, data = d)
round(cbind(coef(m), confint(m)), 3)

# residuals look wrong (not normally distributed)

hist(m$residuals, xlim = c(-7, 15), breaks = 200, xlab = 'Residual', main = NA, border = 'grey35')
abline(v = 0)

# compare poisson model

m2 <- glm(d65 ~ imd5 + offset(log(p65)), data = d, family = 'poisson')
round(exp(cbind(coef(m2), confint(m2))), 3)

# compare rate, rate ratio, and rate difference calculated by hand
# poisson model looks right, linear model looks wrong

actual <- d[!is.na(d65) & !is.na(p65), .(rate = sum(d65) / sum(p65) * 100), imd5][order(rate)]
actual[, ratio := rate / rate[1]]
actual[, difference := rate - rate[1]]
actual[, rate := round(rate, 3)]
actual[, ratio := round(ratio, 3)]
actual[, difference := round(difference, 3)]
actual

# range of normal distributions with different parameters
# -------------------------------------------------------

png('normal.png', height = 6, width = 6, units = 'in', res = 300)
x <- seq(-3, 3, 0.01)
inputs <- expand.grid(sd = 1:2, mean = -1:1)
y <- mapply(dnorm, x = list(x), mean = inputs$mean, sd = inputs$sd, SIMPLIFY = F)
plot(1, type = 'n', xlim = c(-3, 3), ylim = c(0, 0.4), xlab = 'x', ylab = 'Density', main = 'The normal distribution')
abline(v = 0, lty = 3, lwd = 0.5)
mapply(lines, x = list(x), y = y, col = brewer.pal(6, 'Paired'))
mapply(text, x = rep(-1:1, 2), y = sapply(y, max), labels = LETTERS[1:6], font = 2, cex = 2)
dev.off()

# demonstrate that CI = u +/- 1.96 * sqrt(u/n)
# --------------------------------------------

# example data

x <- rpois(10, 90)
u <- mean(x)

# basic formula

ci <-  1.96 * sqrt(u / 10)
c(u, u - ci, u + ci)

# monte carlo simulation

B <- sapply(x, function (z) rpois(1e6, z))
quantile(rowMeans(B), c(0.5, 0.025, 0.975))

# poisson regression

m <- glm(d65 ~ 1 + offset(log(p65)), data = d, family = 'poisson')
exp(c(coef(m), confint(m))) * 100

# non parametric bootstrap

B <- matrix(sample(x, length(x) * 1e5, replace = T), ncol = 1e5)
quantile(colMeans(B), c(0.5, 0.025, 0.975))
