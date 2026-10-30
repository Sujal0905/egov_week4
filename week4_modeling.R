# ==============================================================
# Week 4 Task - Statistical Analysis and Predictive Modeling in E-Governance
# Question : Does national income explain e-government performance, and can
#            we predict a country's EGDI (and "Very High" status) from
#            GDP per capita, region and population?
# Data     : UN EGDI 2024 (transformed in Week 2) + World Bank GDP and
#            population (GitHub mirrors of World Bank data, latest year 2019-2023)
# ==============================================================
library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
set.seed(2024)                         # reproducible folds and train/test split
dir.create("figures", showWarnings = FALSE)

# ---------- 1. LOAD AND PREPARE DATA ----------
egdi <- read_csv("data/egdi_2024_transformed.csv", show_col_types = FALSE)
gdp  <- read_csv("data/wb_gdp_raw.csv", show_col_types = FALSE) %>%
  rename(wb_name = `Country Name`, iso3 = `Country Code`, year = Year, gdp_usd = Value)
pop  <- read_csv("data/wb_pop_raw.csv", show_col_types = FALSE) %>%
  select(iso3 = `Country Code`, year = Year, population = Value)
name_map <- read_csv("data/wb_name_map.csv", show_col_types = FALSE)   # EGDI name -> World Bank name

# Keep the latest GDP value per country, but only if it is recent (2019-2023)
indicators <- gdp %>%
  filter(year >= 2019, year <= 2023, !is.na(gdp_usd)) %>%
  group_by(iso3) %>% slice_max(year, n = 1) %>% ungroup() %>%
  left_join(pop, by = c("iso3", "year")) %>%               # population for the SAME year
  mutate(gdp_pc = gdp_usd / population) %>%
  select(wb_name, iso3, wb_year = year, gdp_pc, population)

full <- egdi %>%
  left_join(name_map, by = "country") %>%
  left_join(indicators, by = "wb_name")

cat("Countries in EGDI:", nrow(full), "\n")
cat("Missing GDP per capita or population:", sum(is.na(full$gdp_pc) | is.na(full$population)), "\n")
print(full %>% filter(is.na(gdp_pc) | is.na(population)) %>% select(country, region, egdi) %>% arrange(desc(egdi)), n = 30)
cat("Data year used (World Bank):\n"); print(table(full$wb_year, useNA = "ifany"))

# Analysis table: complete cases, log-transformed predictors, binary outcome
df <- full %>%
  filter(!is.na(gdp_pc), !is.na(population)) %>%
  mutate(log_gdp = log10(gdp_pc), log_pop = log10(population),
         region = relevel(factor(region), ref = "Africa"),
         responded = factor(ifelse(msq_responded, "Yes", "No"), levels = c("No", "Yes")),
         very_high = as.integer(egdi >= 0.75))
n <- nrow(df); cat("\nCountries in the analysis:", n, "\n")
india <- df %>% filter(country == "India")

# ---------- 2. DESCRIPTIVE STATISTICS ----------
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
desc <- tibble(variable = c("EGDI", "GDP per capita (US$)", "log10 GDP per capita", "Population", "log10 Population"),
               mean = c(mean(df$egdi), mean(df$gdp_pc), mean(df$log_gdp), mean(df$population), mean(df$log_pop)),
               median = c(median(df$egdi), median(df$gdp_pc), median(df$log_gdp), median(df$population), median(df$log_pop)),
               sd = c(sd(df$egdi), sd(df$gdp_pc), sd(df$log_gdp), sd(df$population), sd(df$log_pop)),
               skewness = c(skew(df$egdi), skew(df$gdp_pc), skew(df$log_gdp), skew(df$population), skew(df$log_pop)))
print(desc %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), width = Inf)

# ---------- 3. CONFIDENCE INTERVALS ----------
ci_mean <- function(x) { m <- mean(x); se <- sd(x) / sqrt(length(x)); tq <- qt(0.975, length(x) - 1)
  c(mean = m, lower = m - tq * se, upper = m + tq * se) }
cat("\n95% CI for mean EGDI (all countries):\n"); print(round(ci_mean(df$egdi), 4))
region_ci <- df %>% group_by(region) %>%
  summarise(n = n(), mean = mean(egdi), lower = ci_mean(egdi)[2], upper = ci_mean(egdi)[3], .groups = "drop") %>%
  arrange(desc(mean))
print(region_ci %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

# Bootstrap CI for the median EGDI (no distribution assumption)
boot_med <- replicate(5000, median(sample(df$egdi, replace = TRUE)))
cat("\nBootstrap 95% CI for median EGDI:", round(quantile(boot_med, c(0.025, 0.975)), 4), "\n")

# ---------- 4. HYPOTHESIS TESTS ----------
# 4.1 Do regions differ? ANOVA (+ Tukey) and the non-parametric Kruskal-Wallis
aov1 <- aov(egdi ~ region, data = df); print(summary(aov1))
cat("Kruskal-Wallis:\n"); print(kruskal.test(egdi ~ region, data = df))
tuk <- TukeyHSD(aov1)$region; cat("\nTukey HSD (adjusted p < 0.05 shown):\n")
print(round(tuk[tuk[, "p adj"] < 0.05, ], 4))
eta2 <- summary(aov1)[[1]]$`Sum Sq`[1] / sum(summary(aov1)[[1]]$`Sum Sq`)
cat("Eta-squared (share of EGDI variance explained by region):", round(eta2, 3), "\n")

# 4.2 Do countries that responded to the questionnaire score differently? (Welch t-test + effect size)
tt <- t.test(egdi ~ responded, data = df); print(tt)
sp <- sqrt(((table(df$responded)[1] - 1) * var(df$egdi[df$responded == "No"]) +
            (table(df$responded)[2] - 1) * var(df$egdi[df$responded == "Yes"])) / (n - 2))
cat("Cohen's d:", round((mean(df$egdi[df$responded == "Yes"]) - mean(df$egdi[df$responded == "No"])) / sp, 2), "\n")

# 4.3 Is region associated with EGDI group? (chi-square with simulated p because of small counts)
print(chisq.test(table(df$region, df$egdi >= 0.75), simulate.p.value = TRUE, B = 5000))

# 4.4 Normality of EGDI within regions (checks the ANOVA assumption)
print(df %>% group_by(region) %>% summarise(shapiro_p = round(shapiro.test(egdi)$p.value, 3), .groups = "drop"))

# ---------- 5. CORRELATION ANALYSIS ----------
ct_p <- cor.test(df$egdi, df$log_gdp);                      print(ct_p)
ct_s <- cor.test(df$egdi, df$log_gdp, method = "spearman", exact = FALSE); cat("Spearman rho:", round(ct_s$estimate, 3), "\n")
cat("Correlation EGDI vs raw GDP per capita (not logged):", round(cor(df$egdi, df$gdp_pc), 3), "\n")
cat("Correlation EGDI vs log10 population:", round(cor(df$egdi, df$log_pop), 3), "\n")

# ---------- 6. PREDICTIVE MODELS ----------
# 6.1 Model comparison data: simple, multiple and non-linear linear models
m1 <- lm(egdi ~ log_gdp, data = df)
m2 <- lm(egdi ~ log_gdp + region + responded + log_pop, data = df)
m3 <- lm(egdi ~ poly(log_gdp, 2) + region + responded + log_pop, data = df)
forms <- list("M1 GDP only" = egdi ~ log_gdp,
              "M2 GDP + region + response + pop" = egdi ~ log_gdp + region + responded + log_pop,
              "M3 as M2 with quadratic GDP" = egdi ~ poly(log_gdp, 2) + region + responded + log_pop)

# 10-fold cross-validation (RMSE and MAE), repeated 5 times
cv_lm <- function(formula, data, k = 10, reps = 5) {
  res <- replicate(reps, { f <- sample(rep(1:k, length.out = nrow(data))); pr <- numeric(nrow(data))
    for (i in 1:k) { fit <- lm(formula, data[f != i, ]); pr[f == i] <- pmin(pmax(predict(fit, data[f == i, ]), 0), 1) }
    c(rmse = sqrt(mean((data$egdi - pr)^2)), mae = mean(abs(data$egdi - pr))) })
  rowMeans(res) }
cv <- t(sapply(forms, cv_lm, data = df))
fits <- list(m1, m2, m3)
comparison <- tibble(model = names(forms), adj_r2 = sapply(fits, function(m) summary(m)$adj.r.squared),
                     aic = sapply(fits, AIC), cv_rmse = cv[, "rmse"], cv_mae = cv[, "mae"])
cat("\nMODEL COMPARISON\n"); print(comparison %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
cat("Baseline RMSE (predict the mean for everyone):", round(sd(df$egdi), 4), "\n")

# 6.2 Chosen model: M2 (summary, CIs)
print(summary(m2)); cat("\n95% CI for coefficients:\n"); print(round(confint(m2), 4))
cat("\nEffect of a 10% rise in GDP per capita (model M2):",
    round(coef(m2)["log_gdp"] * log10(1.10), 4), "EGDI points\n")
cat("Effect of a 10-fold rise in GDP per capita (model M2):", round(coef(m2)["log_gdp"], 4), "\n")

# 6.3 Hold-out validation (80% train / 20% test)
idx   <- sample(seq_len(n), size = round(0.8 * n))
train <- df[idx, ]; test <- df[-idx, ]
m2_tr <- lm(formula(m2), data = train)
test$pred <- pmin(pmax(predict(m2_tr, test), 0), 1)
rmse <- sqrt(mean((test$egdi - test$pred)^2)); mae <- mean(abs(test$egdi - test$pred))
r2_test <- 1 - sum((test$egdi - test$pred)^2) / sum((test$egdi - mean(test$egdi))^2)
cat("\nHOLD-OUT TEST (n =", nrow(test), "): RMSE", round(rmse, 4), " MAE", round(mae, 4), " R2", round(r2_test, 3), "\n")
cat("Share of test countries predicted within 0.10 of actual:", round(mean(abs(test$egdi - test$pred) <= 0.10) * 100, 1), "%\n")

# 6.4 Assumption checks for M2
df$fitted <- fitted(m2); df$resid <- resid(m2); df$rstud <- rstudent(m2); df$cooks <- cooks.distance(m2)
cat("\nShapiro-Wilk test of residuals: p =", round(shapiro.test(df$resid)$p.value, 3), "\n")
aux <- summary(lm(I(resid^2) ~ fitted, data = df))$r.squared             # Koenker (studentised) Breusch-Pagan
bp_p <- pchisq(n * aux, df = 1, lower.tail = FALSE); cat("Breusch-Pagan (Koenker) test for constant variance: p =", round(bp_p, 3), "\n")
X <- model.matrix(m2)[, -1]
vif <- sapply(seq_len(ncol(X)), function(j) 1 / (1 - summary(lm(X[, j] ~ X[, -j]))$r.squared))
names(vif) <- colnames(X); cat("Variance inflation factors:\n"); print(round(vif, 2))
cat("Observations with Cook's distance > 4/n:", sum(df$cooks > 4 / n), "\n")
cat("Observations with |studentised residual| > 2:", sum(abs(df$rstud) > 2), "\n")

# 6.5 Prediction for India and the largest over/under performers
india <- df %>% filter(country == "India")      # refresh so India carries fitted values and residuals
india_pred <- predict(m2, india, interval = "prediction", level = 0.95)
cat("\nIndia actual EGDI:", india$egdi, "| predicted:", round(india_pred[1], 3),
    "| 95% prediction interval:", round(india_pred[2], 3), "to", round(india_pred[3], 3), "\n")
cat("India GDP per capita (US$):", round(india$gdp_pc), " (year", india$wb_year, ")\n")
over  <- df %>% arrange(desc(resid)) %>% slice_head(n = 8) %>% select(country, region, gdp_pc, egdi, fitted, resid)
under <- df %>% arrange(resid)       %>% slice_head(n = 8) %>% select(country, region, gdp_pc, egdi, fitted, resid)
cat("\nOver-performers (EGDI higher than model predicts):\n"); print(over %>% mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("\nUnder-performers (EGDI lower than model predicts):\n"); print(under %>% mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("Rank of India's residual among", n, "countries (1 = largest positive):", which(df$country[order(-df$resid)] == "India"), "\n")

# ---------- 7. LOGISTIC REGRESSION: PREDICT "VERY HIGH" STATUS ----------
lg <- glm(very_high ~ log_gdp + region + responded, data = df, family = binomial)
print(summary(lg))
or <- exp(cbind(OR = coef(lg), confint.default(lg))); cat("\nOdds ratios with 95% CI (Wald):\n"); print(round(or, 3))
lg_tr <- glm(formula(lg), data = train, family = binomial)
test$p_vh <- predict(lg_tr, test, type = "response"); test$pred_vh <- as.integer(test$p_vh >= 0.5)
cm <- table(Actual = test$very_high, Predicted = test$pred_vh); print(cm)
cat("Hold-out accuracy:", round(mean(test$very_high == test$pred_vh), 3),
    " | Baseline (always predict majority):", round(max(table(test$very_high)) / nrow(test), 3), "\n")
auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cat("Hold-out AUC:", round(auc(test$very_high, test$p_vh), 3), "\n")
df$p_vh_all <- predict(lg, df, type = "response")
cat("India probability of Very High (model):", round(df$p_vh_all[df$country == "India"], 3), "\n")

# 10-fold cross-validated accuracy and AUC for the logistic model
f <- sample(rep(1:10, length.out = n)); cvp <- numeric(n)
for (i in 1:10) cvp[f == i] <- predict(glm(formula(lg), data = df[f != i, ], family = binomial), df[f == i, ], type = "response")
cat("10-fold CV accuracy:", round(mean((cvp >= 0.5) == df$very_high), 3), " AUC:", round(auc(df$very_high, cvp), 3), "\n")

# ---------- 8. VISUALISATIONS ----------
theme_set(theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 13),
          plot.title.position = "plot", plot.caption = element_text(colour = "grey45", size = 8), panel.grid.minor = element_blank()))
src <- "Sources: UN E-Government Survey 2024; World Bank GDP and population (latest year 2019-2023)."
pal_reg <- c(Africa = "#D55E00", Americas = "#E69F00", Asia = "#009E73", Europe = "#0072B2", Oceania = "#CC79A7")
sv <- function(p, f, w = 8, h = 4.8) ggsave(paste0("figures/", f), p, width = w, height = h, dpi = 200, bg = "white")

# 8.1 Why log-transform GDP per capita?
p1 <- df %>% select(`GDP per capita (US$)` = gdp_pc, `log10 GDP per capita` = log_gdp) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(value)) + geom_histogram(bins = 25, fill = "#0072B2", colour = "white") +
  facet_wrap(~name, scales = "free") +
  labs(title = "GDP per capita is heavily skewed; the log scale is far more symmetric",
       x = NULL, y = "Countries", caption = src)
sv(p1, "w4_log_transform.png", 8, 3.8)

# 8.2 Relationship EGDI ~ log GDP per capita
p2 <- ggplot(df, aes(10^log_gdp, egdi)) +
  geom_point(aes(colour = region), size = 2.2, alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey20", se = TRUE, linewidth = 0.8) +
  geom_point(data = india, shape = 21, size = 5, stroke = 1.3, colour = "#CC0066", fill = NA) +
  annotate("text", x = india$gdp_pc * 0.82, y = india$egdi + 0.05, label = "India", colour = "#CC0066", fontface = "bold", hjust = 1) +
  scale_x_log10(labels = label_dollar()) + scale_colour_manual(values = pal_reg, name = "Region") +
  labs(title = "Richer countries score higher, but income is not the whole story",
       subtitle = sprintf("EGDI 2024 vs GDP per capita (log scale); r = %.2f, 95%% CI %.2f to %.2f",
                          ct_p$estimate, ct_p$conf.int[1], ct_p$conf.int[2]),
       x = "GDP per capita (US$, log scale)", y = "EGDI 2024", caption = src) +
  theme(legend.position = "bottom")
sv(p2, "w4_scatter_gdp.png")

# 8.3 Confidence intervals by region
p3 <- ggplot(region_ci, aes(reorder(region, mean), mean, colour = region)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), linewidth = 1, size = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = mean(df$egdi), linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = pal_reg) + coord_flip() +
  labs(title = "Regional mean EGDI with 95% confidence intervals",
       subtitle = "Dashed line = overall mean; intervals that do not overlap indicate clear differences",
       x = NULL, y = "Mean EGDI", caption = src)
sv(p3, "w4_region_ci.png", 8, 4)

# 8.4 Diagnostic plots for M2 (4 panels)
diag_df <- bind_rows(
  tibble(panel = "Residuals vs fitted", x = df$fitted, y = df$resid, country = df$country),
  tibble(panel = "Normal Q-Q of residuals", x = qqnorm(df$rstud, plot.it = FALSE)$x, y = sort(df$rstud)[rank(df$rstud)], country = df$country),
  tibble(panel = "Scale-location", x = df$fitted, y = sqrt(abs(df$rstud)), country = df$country),
  tibble(panel = "Cook's distance", x = seq_len(n), y = df$cooks, country = df$country))
p4 <- ggplot(diag_df, aes(x, y)) + geom_point(alpha = 0.7, colour = "#0072B2", size = 1.6) +
  geom_smooth(data = ~ filter(.x, panel %in% c("Residuals vs fitted", "Scale-location")), method = "loess", formula = y ~ x,
              se = FALSE, colour = "#D55E00", linewidth = 0.7) +
  geom_abline(data = tibble(panel = "Normal Q-Q of residuals"), aes(intercept = 0, slope = 1), linetype = "dashed") +
  geom_hline(data = tibble(panel = "Residuals vs fitted"), aes(yintercept = 0), linetype = "dashed") +
  geom_hline(data = tibble(panel = "Cook's distance"), aes(yintercept = 4 / n), linetype = "dashed", colour = "#D55E00") +
  facet_wrap(~ panel, scales = "free") +
  labs(title = "Model diagnostics for the multiple regression (M2)", x = NULL, y = NULL, caption = src)
sv(p4, "w4_diagnostics.png", 9, 6)

# 8.5 Predicted vs actual on the hold-out test set
p5 <- ggplot(test, aes(pred, egdi, colour = region)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 3) + scale_colour_manual(values = pal_reg, name = "Region") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Predicted vs actual EGDI on unseen\ncountries (hold-out test set)",
       subtitle = sprintf("n = %d, RMSE = %.3f, MAE = %.3f, R-squared = %.2f", nrow(test), rmse, mae, r2_test),
       x = "Predicted EGDI", y = "Actual EGDI", caption = src) + theme(legend.position = "bottom")
sv(p5, "w4_pred_vs_actual.png", 6.5, 6)

# 8.6 Model comparison by cross-validated error
p6 <- comparison %>% mutate(model = factor(model, levels = rev(model))) %>%
  ggplot(aes(model, cv_rmse)) + geom_col(fill = "#0072B2", width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", cv_rmse)), hjust = -0.15) +
  geom_hline(yintercept = sd(df$egdi), linetype = "dashed", colour = "#D55E00") +
  annotate("text", x = 0.55, y = sd(df$egdi), label = "baseline (predict the mean) ", hjust = 1, vjust = -0.5, colour = "#D55E00", size = 3.5) +
  coord_flip(ylim = c(0, 0.27)) +
  labs(title = "Each added layer lowers prediction error;\nall models beat the baseline",
       subtitle = "10-fold cross-validated RMSE (lower is better)", x = NULL, y = "RMSE", caption = src)
sv(p6, "w4_model_comparison.png", 8, 3.6)

# 8.7 Over- and under-performers relative to income
extr <- bind_rows(over, under, india %>% select(country, region, gdp_pc, egdi, fitted, resid)) %>%
  distinct(country, .keep_all = TRUE) %>%
  mutate(kind = ifelse(resid > 0, "Higher than predicted", "Lower than predicted"),
         country = stringr::str_wrap(country, 26))
p7 <- ggplot(extr, aes(reorder(country, resid), resid, fill = kind)) + geom_col(width = 0.7) +
  geom_col(data = filter(extr, country == "India"), fill = "#CC0066", width = 0.7) +
  scale_fill_manual(values = c("Higher than predicted" = "#009E73", "Lower than predicted" = "#D55E00"), name = NULL) +
  coord_flip() + labs(title = "Who beats or lags what their income\nand region predict?",
       subtitle = "Residual of model M2 (actual minus predicted EGDI); India in magenta", x = NULL, y = "Residual (EGDI points)", caption = src) +
  theme(legend.position = "bottom")
sv(p7, "w4_residual_ranking.png", 8, 6)

# 8.8 Coefficient plot (M2) with 95% CIs
cf <- as.data.frame(cbind(est = coef(m2), confint(m2))); names(cf) <- c("est", "lo", "hi"); cf$term <- rownames(cf)
cf <- cf %>% filter(term != "(Intercept)") %>%
  mutate(term = recode(term, log_gdp = "10x GDP per capita", regionAmericas = "Americas (vs Africa)", regionAsia = "Asia (vs Africa)",
                       regionEurope = "Europe (vs Africa)", regionOceania = "Oceania (vs Africa)", respondedYes = "Responded to questionnaire",
                       log_pop = "10x population"))
p8 <- ggplot(cf, aes(reorder(term, est), est)) + geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = lo, ymax = hi), colour = "#0072B2", linewidth = 0.9) + coord_flip() +
  labs(title = "What moves EGDI? Coefficients of the multiple regression",
       subtitle = "Change in EGDI (0-1 scale) with 95% CI;\nintervals that cross zero are not significant",
       x = NULL, y = "Estimated effect on EGDI", caption = src)
sv(p8, "w4_coefficients.png", 8, 4.2)

# 8.9 ROC curve for the logistic model (hold-out)
roc_pts <- function(y, p) { th <- sort(unique(c(0, p, 1)), decreasing = TRUE)
  tibble(th = th, tpr = sapply(th, function(t) mean(p[y == 1] >= t)), fpr = sapply(th, function(t) mean(p[y == 0] >= t))) }
roc <- roc_pts(test$very_high, test$p_vh)
p9 <- ggplot(roc, aes(fpr, tpr)) + geom_abline(linetype = "dashed", colour = "grey50") +
  geom_line(colour = "#0072B2", linewidth = 1.1) + coord_equal() +
  annotate("text", x = 0.62, y = 0.2, label = sprintf("Hold-out AUC = %.2f", auc(test$very_high, test$p_vh)), size = 4.2) +
  labs(title = "Logistic model: separating Very High\nfrom other countries",
       subtitle = "ROC curve on unseen test countries", x = "False positive rate", y = "True positive rate", caption = src)
sv(p9, "w4_roc.png", 6, 5.4)

# ---------- 9. SAVE ----------
write_csv(df %>% select(country, region, egdi, group, gdp_pc, population, wb_year, fitted, resid, p_vh_all), "data/egdi_model_data.csv")
cat("\nDone.\n")
