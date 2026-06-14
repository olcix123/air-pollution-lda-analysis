rm(list = ls())

library(MASS)
library(mvnormtest)
library(biotools)

data <- read.table("Amostra_Grupo 4.txt", header = TRUE, sep = "\t")

Year <- as.factor(data$X6)
X <- data[, c("X1", "X2", "X3", "X4", "X5")]
lda_data <- data.frame(Year, X)

cat("\nFigure 1 - A priori probabilities\n")
print(round(prop.table(table(Year)), 4))

cat("\nMultivariate normality tests\n")
for (g in levels(Year)) {
  cat("\n", g, "\n", sep = "")
  print(mshapiro.test(t(X[Year == g, ])))
}

cat("\nTable 2 - Box's M test\n")
print(boxM(X, Year))

cat("\nFigure 2 - Group means\n")
group_means <- aggregate(X, by = list(Year = Year), mean)
group_means[, -1] <- round(group_means[, -1], 2)
print(group_means)

k <- length(levels(Year))
p <- ncol(X)
cat("\nNumber of discriminant functions\n")
cat("m =", min(k - 1, p), "\n")

lda_model <- lda(Year ~ ., data = lda_data, na.action = "na.omit")
lda_pred <- predict(lda_model)

cat("\nFigure 3 - Coefficients of the discriminant functions\n")
print(round(lda_model$scaling, 5))

eigenvalues <- lda_model$svd^2
explained <- eigenvalues / sum(eigenvalues)
canonical_correlations <- sqrt(eigenvalues / (1 + eigenvalues))

cat("\nTable 5 - Importance of the discriminant functions\n")
importance <- data.frame(
  Function = c("LD1", "LD2"),
  Eigenvalue = round(eigenvalues, 4),
  Between_group_variability_explained = paste0(round(explained * 100, 2), "%"),
  Canonical_correlation = round(canonical_correlations, 4)
)
print(importance)

n <- nrow(X)
wilks_lambda <- prod(1 / (1 + eigenvalues))
chi_square <- -((n - 1) - ((p + k) / 2)) * log(wilks_lambda)
df <- p * (k - 1)
p_value <- pchisq(chi_square, df = df, lower.tail = FALSE)

cat("\nFigure 4 - Wilks' Lambda test\n")
cat("Wilks' Lambda =", round(wilks_lambda, 4), "\n")
cat("Chi-square =", round(chi_square, 4), "\n")
cat("df =", df, "\n")
cat("p-value =", format.pval(p_value, digits = 4), "\n")

cat("\nFigure 5 - First ten discriminant scores\n")
scores <- data.frame(Year = Year, lda_pred$x)
print(data.frame(Year = scores$Year[1:10], round(scores[1:10, -1], 4)))

plot(scores$LD1, scores$LD2,
     col = as.numeric(Year),
     pch = 19,
     xlab = "LD1",
     ylab = "LD2",
     main = "Linear Discriminant Analysis")
legend("topright", legend = levels(Year), col = seq_along(levels(Year)), pch = 19)

cat("\nFigure 7 / Table 6 - Classification matrix, resubstitution method\n")
conf_resub <- table(Original = Year, Predicted = lda_pred$class)
print(conf_resub)

cat("\nFigure 8 / Table 7 - Correct and incorrect classification percentages by group\n")
correct_by_group <- diag(conf_resub)
total_by_group <- rowSums(conf_resub)
classification_results <- data.frame(
  Group = names(correct_by_group),
  Correctly_classified = paste0(correct_by_group, "/", total_by_group),
  Correct_classification = paste0(sprintf("%.2f", correct_by_group / total_by_group * 100), "%"),
  Incorrect_classification = paste0(sprintf("%.2f", (1 - correct_by_group / total_by_group) * 100), "%")
)
print(classification_results)
cat("Total correct classification =", sprintf("%.2f%%", sum(diag(conf_resub)) / sum(conf_resub) * 100), "\n")

lda_cv <- lda(Year ~ ., data = lda_data, na.action = "na.omit", CV = TRUE)

cat("\nFigure 9 / Table 8 - Classification matrix, leave-one-out cross-validation\n")
conf_cv <- table(Original = Year, Predicted = lda_cv$class)
print(conf_cv)
cat("Total correct classification =", sprintf("%.2f%%", sum(diag(conf_cv)) / sum(conf_cv) * 100), "\n")

cat("\nFigure 10 - Summary\n")
cat("LD1 explains", sprintf("%.2f%%", explained[1] * 100), "of the between-group variability.\n")
cat("Total correct classification rate:", sprintf("%.2f%%", sum(diag(conf_resub)) / sum(conf_resub) * 100), "\n")
cat("Leave-one-out cross-validation rate:", sprintf("%.2f%%", sum(diag(conf_cv)) / sum(conf_cv) * 100), "\n")
