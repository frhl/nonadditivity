library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(41)

# --- 1. Define the Orthogonal Encoding Function ---
get_orthogonal_codes <- function(genotypes) {
  n <- length(genotypes)
  r <- sum(genotypes == 0) / n 
  h <- sum(genotypes == 1) / n 
  a <- sum(genotypes == 2) / n 
  
  dom_codes <- numeric(n)
  dom_codes[genotypes == 0] <- -h * a
  dom_codes[genotypes == 1] <- 2 * a * r
  dom_codes[genotypes == 2] <- -h * r
  
  dom_codes <- scale(dom_codes)
  add_codes <- scale(genotypes)
  
  return(list(add=add_codes, dom=dom_codes))
}

# ==============================================================================
# SIMULATION 1: Robustness / False Positive Rates
# ==============================================================================

N <- 5000         
MAF <- 0.3        
beta_add_values <- seq(0, 0.5, length.out = 20) 

# Create vectors to store results
sim1_list <- list()

for(b_add in beta_add_values) {
  fp_recessive <- 0
  fp_orthogonal <- 0
  
  for(i in 1:100) { 
    g <- rbinom(N, 2, MAF)
    # True Non-Additive effect is ZERO
    y <- b_add * g + rnorm(N)
    
    # Model A: Standard Recessive
    g_rec <- ifelse(g == 2, 1, 0)
    pval_rec <- summary(lm(y ~ g_rec))$coefficients[2,4]
    
    # Model B: Orthogonal
    codes <- get_orthogonal_codes(g)
    pval_orth <- summary(lm(y ~ codes$add + codes$dom))$coefficients[3,4]
    
    if(pval_rec < 0.05) fp_recessive <- fp_recessive + 1
    if(pval_orth < 0.05) fp_orthogonal <- fp_orthogonal + 1
  }
  
  sim1_list[[length(sim1_list) + 1]] <- data.frame(
    Beta_Add = b_add,
    Recessive = fp_recessive / 100,
    Orthogonal = fp_orthogonal / 100
  )
}

# Bind and Pivot to Long Format for ggplot
results <- do.call(rbind, sim1_list) %>%
  pivot_longer(cols = c("Recessive", "Orthogonal"), 
               names_to = "Model", 
               values_to = "FPR")

p1 <- ggplot(results, aes(x = Beta_Add, y = FPR, color = Model, shape = Model)) +
  # 1. The Threshold Line (Subtle)
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  # 2. The Data (Refined thickness)
  geom_line(linewidth = 0.7) +
  geom_point(size = 2, stroke = 0.5, fill = "white") + # Filled points look cleaner
  
  # 3. Professional Color Palette (Okabe-Ito compliant) & Shapes
  scale_color_manual(values = c("Recessive" = "#D55E00", "Orthogonal" = "#0072B2")) +
  scale_shape_manual(values = c("Recessive" = 19, "Orthogonal" = 21)) + # Solid vs Hollow circle
  
  # 4. Tighter Axes
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  
  # 5. Labels
  labs(
    title = "False Positive Rate vs. Additive Effect Size",
    x = "Magnitude of True Additive Effect",
    y = "False Positive Rate (Type I Error)"
  ) +
  
  # 6. The "Classic" Theme adjustments
  theme_bw(base_size = 12) + 
  theme(
    # Text hierarchy
    plot.title = element_text(face = "bold", size = 12, color = "black"),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9),
    
    # Legend placement (Inset is often cleaner for publication)
    legend.position = c(0.8, 0.5), # Adjust x,y coordinates to put legend INSIDE plot
    legend.background = element_rect(color = "black", size = 0.2), # Box around legend
    legend.title = element_blank()
  )

print(p1)

# ==============================================================================
# SIMULATION 2: Parameter Recovery
# ==============================================================================

set.seed(123)
N_genes <- 1000
N_samples <- 2000

true_beta_D <- numeric(N_genes)
est_beta_D  <- numeric(N_genes)

for(i in 1:N_genes) {
  b_A <- rnorm(1, mean=0, sd=0.2)
  b_D <- 0.8 * b_A + rnorm(1, 0, 0.05) 
  
  true_beta_D[i] <- b_D
  
  g <- rbinom(N_samples, 2, 0.2)
  codes <- get_orthogonal_codes(g) 
  
  y <- b_A * codes$add + b_D * codes$dom + rnorm(N_samples)
  fit <- lm(y ~ codes$add + codes$dom)
  est_beta_D[i] <- coef(fit)[3]
}

sim2_df <- data.frame(True = true_beta_D, Est = est_beta_D)

# Calculate correlation for the subtitle
corr_val <- round(cor(sim2_df$True, sim2_df$Est - sim2_df$True), 3)

p2 <- ggplot(sim2_df, aes(x = True, y = Est)) +
  geom_abline(slope = 1, intercept = 0, color = "gray30", linetype = "dotted") +
  
  # Use a hollow shape with transparency for density
  geom_point(alpha = 0.5, color = "#0072B2", shape = 1, size = 1.5) +
  
  labs(
    title = "Recovery of Non-Additive Effects",
    subtitle = paste0("Pearson correlation: ", round(cor(sim2_df$True, sim2_df$Est), 3)),
    x = "True Effect Size (Simulated)",
    y = "Estimated Effect Size (Orthogonal Model)"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9)
  )

print(p2)

# Verify correlation output
cat("Correlation between True Beta_D and Error:", corr_val, "\n")