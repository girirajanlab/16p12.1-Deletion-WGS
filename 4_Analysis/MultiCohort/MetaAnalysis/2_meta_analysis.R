library(meta)

# Perform a meta analyses

# Input and output files
FILENAMES <- c('/path/to/DD/UKB/stats.csv',
               '/path/to/UKB/MyCode/AoU/stats.csv',
               '/path/to/DD/SPARK/stats.csv') # Use the outputs of script 1_compile_data.py
OUTPUT_FILES <- c('/path/to/DD/UKB/meta.csv',
                  '/path/to/UKB/MyCode/AoU/meta.csv',
                  '/path/to/DD/SPARK/meta.csv') # These output files are presented in Table S5H

for (i in 1:3) {
  filename <- FILENAMES[i]
  df <- read.csv(filename)
  
  vars <- unique(df$Variant)
  phenotypes <- unique(df$Phenotype)
  
  stats <- data.frame(matrix(ncol=12, nrow=0))
  if (i==2) {
    stats <- data.frame(matrix(ncol=15, nrow=0))
  }
  for (v in vars) {
    for (p in phenotypes) {
      subdf <- df[(df$Phenotype==p) & (df$Variant==v), ]
      if ((length(rownames(subdf)))==0) {
        next
      }
      m.cont <- metacont(n.e = Case.N,
                         mean.e = Case.mean,
                         sd.e = Case.SD,
                         n.c = Control.N,
                         mean.c = Control.mean,
                         sd.c = Control.SD,
                         studlab = Cohort,
                         data = subdf,
                         sm = "SMD",
                         method.smd = "Cohen",
                         common = FALSE,
                         random = TRUE,
                         method.tau = "REML",
                         method.random.ci = "HK",
                         title = paste(v, p))
      
      # Parse output
      out <- data.frame(matrix(ncol=12, nrow=1))
      g1 <- m.cont$studlab[1]
      g2 <- m.cont$studlab[2]
      colnames(out) <- c('Variant', 'Phenotype',
                         paste(g1, '_effect', sep=''), paste(g1, '_CI_low', sep=''), paste(g1, '_CI_upp', sep=''),
                         paste(g2, '_effect', sep=''), paste(g2, '_CI_low', sep=''), paste(g2, '_CI_upp', sep=''),
                         'Random_effect', 'Random_CI_low', 'Random_CI_upp', 'p_value')
      if (i==2) {
        out <- data.frame(matrix(ncol=15, nrow=1))
        g3 <- m.cont$studlab[3]
        colnames(out) <- c('Variant', 'Phenotype',
                         paste(g1, '_effect', sep=''), paste(g1, '_CI_low', sep=''), paste(g1, '_CI_upp', sep=''),
                         paste(g2, '_effect', sep=''), paste(g2, '_CI_low', sep=''), paste(g2, '_CI_upp', sep=''),
                         paste(g3, '_effect', sep=''), paste(g3, '_CI_low', sep=''), paste(g3, '_CI_upp', sep=''),
                         'Random_effect', 'Random_CI_low', 'Random_CI_upp', 'p_value')
      }
      
      out$Variant <- v
      out$Phenotype <- p
      out[, paste(g1, '_effect', sep='')] <- m.cont$TE[1]
      out[, paste(g1, '_CI_low', sep='')] <- m.cont$lower[1]
      out[, paste(g1, '_CI_upp', sep='')] <- m.cont$upper[1]
      out[, paste(g2, '_effect', sep='')] <- m.cont$TE[2]
      out[, paste(g2, '_CI_low', sep='')] <- m.cont$lower[2]
      out[, paste(g2, '_CI_upp', sep='')] <- m.cont$upper[2]
      out$Random_effect <- m.cont$TE.random
      out$Random_CI_low <- m.cont$lower.random
      out$Random_CI_upp <- m.cont$upper.random
      out$p_value <- m.cont$pval.random
      
      if (i==2) {
        out[, paste(g3, '_effect', sep='')] <- m.cont$TE[3]
        out[, paste(g3, '_CI_low', sep='')] <- m.cont$lower[3]
        out[, paste(g3, '_CI_upp', sep='')] <- m.cont$upper[3]
      }
      
      stats <- rbind(stats, out)
    }
  }
  # FDR correction
  stats['FDR']=p.adjust(stats$p_value, method='BH')
  
  # Save
  outfile <- OUTPUT_FILES[i]
  write.csv(stats, outfile, row.names = F)
}
