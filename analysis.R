# Computational Biomedicine | Tutorial 3 | 16 November 2026
# Somatic mutations and survival analysis (TCGA AML / LAML)
# Run entire script with Source from RStudio. No separate data download needed.
# Output files will be created in Week3_R_outputs under your working directory.

options(repos = c(CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
if (!requireNamespace('maftools', quietly = TRUE)) BiocManager::install('maftools', ask = FALSE, update = FALSE)
for (pkg in c('survival', 'ggfortify', 'ggplot2')) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages(library(maftools))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(ggplot2))

out <- 'Week3_R_outputs'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
export_csv <- function(x, name) write.csv(as.data.frame(x), file.path(out, name), row.names = FALSE)
plot_pdf <- function(name, expr, width = 11, height = 7) {
  grDevices::pdf(file.path(out, name), width = width, height = height, onefile = TRUE)
  tryCatch(force(expr), finally = grDevices::dev.off())
}

# 1.1 Import MAF and clinical annotations bundled with maftools.
maf_path <- system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
clin_path <- system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
stopifnot(nzchar(maf_path), nzchar(clin_path), file.exists(maf_path), file.exists(clin_path))
laml <- maftools::read.maf(maf = maf_path, clinicalData = clin_path)
cat('\n===== MAF OBJECT =====\n'); print(laml)
cat('\n===== MAF STRUCTURE =====\n'); str(laml, max.level = 1)
ss <- as.data.frame(maftools::getSampleSummary(laml))
gs <- as.data.frame(maftools::getGeneSummary(laml))
clinical <- as.data.frame(maftools::getClinicalData(laml))
mut <- as.data.frame(laml@data)
export_csv(ss, '01_sample_summary.csv')
export_csv(gs, '02_gene_summary.csv')
export_csv(clinical, '03_clinical_data.csv')
writeLines(c('MAF object contains variant-level data and mutation/clinical summaries.',
             paste('Number of samples with variant calls:', length(unique(mut$Tumor_Sample_Barcode))),
             paste('Number of variant records:', nrow(mut)),
             paste('Clinical records:', nrow(clinical)),
             paste('Clinical columns:', paste(names(clinical), collapse = ', ')),
             paste('Variant columns:', paste(names(mut), collapse = ', '))),
           con = file.path(out, '00_data_inspection.txt'))

# 1.2 Mutation distribution and top-five oncoplot.
plot_pdf('04_mutation_summary.pdf', {
  maftools::plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
})
plot_pdf('05_oncoplot_top5.pdf', {
  maftools::oncoplot(maf = laml, top = 5)
}, width = 12, height = 7.5)

# 1.2 Transitions and transversions (SNP substitutions; not the GWAS SNP testing step).
titv_result <- maftools::titv(maf = laml, plot = FALSE, useSyn = TRUE)
plot_pdf('06_titv.pdf', {
  maftools::plotTiTv(res = titv_result)
}, width = 11, height = 8)

# 1.2 FLT3 lollipop with automatic choice of protein annotation column.
aa_candidates <- c('Protein_Change', 'AAChange', 'HGVSp_Short', 'HGVSp')
aa_found <- aa_candidates[aa_candidates %in% names(mut)]
if (length(aa_found) > 0) {
  aa_col <- aa_found[1]
  cat('Using amino-acid change column:', aa_col, '\n')
  tryCatch(plot_pdf('07_FLT3_lollipop.pdf', {
    maftools::lollipopPlot(maf = laml, gene = 'FLT3', AACol = aa_col,
                           showMutationRate = FALSE)
  }, width = 12, height = 6), error = function(e) {
    warning('Lollipop plot could not be produced: ', conditionMessage(e))
  })
} else {
  warning('No amino-acid annotation field found for lollipopPlot; check getFields(laml).')
}

# 2.1 Prepare survival time and event status; ensure right censoring.
required_cols <- c('Tumor_Sample_Barcode', 'days_to_last_followup',
                   'Overall_Survival_Status', 'FAB_classification')
missing_cols <- setdiff(required_cols, names(clinical))
if (length(missing_cols)) stop('Missing clinical columns: ', paste(missing_cols, collapse = ', '))
dat <- clinical[, required_cols]
dat$days_to_last_followup <- suppressWarnings(as.numeric(as.character(dat$days_to_last_followup)))
raw_event <- toupper(trimws(as.character(dat$Overall_Survival_Status)))
# The example file uses 1 = death, 0 = alive / censored.
dat$event <- ifelse(raw_event %in% c('1', 'DEAD', 'DECEASED', 'TRUE'), 1L,
                    ifelse(raw_event %in% c('0', 'ALIVE', 'LIVING', 'FALSE'), 0L, NA_integer_))
dat$FAB_classification <- trimws(as.character(dat$FAB_classification))
dat$FAB_classification[dat$FAB_classification %in% c('', 'NA', 'Unknown')] <- NA_character_
# Must not silently interpret missing follow-up/status or negative time.
dat <- dat[!is.na(dat$days_to_last_followup) & dat$days_to_last_followup >= 0 &
           !is.na(dat$event), , drop = FALSE]
stopifnot(nrow(dat) > 0)
cat('Survival records:', nrow(dat), ' events:', sum(dat$event), '\n')
export_csv(dat, '08_survival_data_clean.csv')
surv_obj <- survival::Surv(time = dat$days_to_last_followup, event = dat$event)
fit_all <- survival::survfit(surv_obj ~ 1, data = dat)
p_all <- ggfortify::autoplot(fit_all, conf.int = TRUE) +
  labs(title = 'Overall survival in TCGA AML cohort', x = 'Follow-up time (days)',
       y = 'Estimated survival probability') + theme_minimal(base_size = 13)
ggsave(file.path(out, '09_KM_overall.pdf'), p_all, width = 9, height = 6)

# FAB groups: exclude only missing FAB values, make reference explicit and reproducible.
fab <- dat[!is.na(dat$FAB_classification), , drop = FALSE]
fab$FAB_classification <- factor(fab$FAB_classification)
# Use M0 if present, otherwise the alphabetically first available group.
ref_fab <- if ('M0' %in% levels(fab$FAB_classification)) 'M0' else levels(fab$FAB_classification)[1]
fab$FAB_classification <- relevel(fab$FAB_classification, ref = ref_fab)
cat('Cox reference FAB group:', ref_fab, '\n')
fit_fab <- survival::survfit(Surv(days_to_last_followup, event) ~ FAB_classification, data = fab)
p_fab <- ggfortify::autoplot(fit_fab, conf.int = FALSE) +
  labs(title = 'Kaplan-Meier by FAB subtype', x = 'Follow-up time (days)',
       y = 'Estimated survival probability', color = 'FAB subtype') +
  theme_minimal(base_size = 12)
ggsave(file.path(out, '10_KM_FAB.pdf'), p_fab, width = 11, height = 7)
cox_fab <- survival::coxph(Surv(days_to_last_followup, event) ~ FAB_classification, data = fab)
cs <- summary(cox_fab)
cox_fab_tab <- data.frame(term = rownames(cs$coefficients),
                          hazard_ratio = cs$conf.int[, 'exp(coef)'],
                          CI95_low = cs$conf.int[, 'lower .95'],
                          CI95_high = cs$conf.int[, 'upper .95'],
                          p_value = cs$coefficients[, 'Pr(>|z|)'], row.names = NULL)
cox_fab_tab <- cox_fab_tab[order(cox_fab_tab$p_value), , drop = FALSE]
export_csv(cox_fab_tab, '11_Cox_FAB_results.csv')
ph_fab <- survival::cox.zph(cox_fab)
writeLines(c(paste('FAB reference:', ref_fab),
             paste('Lowest p-value versus reference:', cox_fab_tab$term[1]),
             paste('p-value:', signif(cox_fab_tab$p_value[1], 5)),
             '', 'Cox PH diagnostic:', capture.output(print(ph_fab))),
           file.path(out, '12_FAB_model_interpretation.txt'))

# 2.2 FLT3 mutation status, based on presence of FLT3 in the mutation table.
flt3_ids <- unique(as.character(mut$Tumor_Sample_Barcode[mut$Hugo_Symbol == 'FLT3']))
dat$FLT3_status <- factor(ifelse(as.character(dat$Tumor_Sample_Barcode) %in% flt3_ids,
                                 'Mutated', 'Not_mutated'),
                          levels = c('Not_mutated', 'Mutated'))
export_csv(dat, '13_survival_with_FLT3.csv')
fit_flt3 <- survival::survfit(Surv(days_to_last_followup, event) ~ FLT3_status, data = dat)
p_flt3 <- ggfortify::autoplot(fit_flt3, conf.int = FALSE) +
  labs(title = 'Kaplan-Meier by FLT3 mutation status', x = 'Follow-up time (days)',
       y = 'Estimated survival probability', color = 'FLT3') + theme_minimal(base_size = 13)
ggsave(file.path(out, '14_KM_FLT3.pdf'), p_flt3, width = 9, height = 6)
cox_flt3 <- survival::coxph(Surv(days_to_last_followup, event) ~ FLT3_status, data = dat)
fs <- summary(cox_flt3)
flt3_tab <- data.frame(comparison = rownames(fs$coefficients),
                       hazard_ratio = fs$conf.int[, 'exp(coef)'],
                       CI95_low = fs$conf.int[, 'lower .95'],
                       CI95_high = fs$conf.int[, 'upper .95'],
                       p_value = fs$coefficients[, 'Pr(>|z|)'], row.names = NULL)
export_csv(flt3_tab, '15_Cox_FLT3_results.csv')
logrank <- survival::survdiff(Surv(days_to_last_followup, event) ~ FLT3_status, data = dat)
logrank_p <- stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
writeLines(c('FLT3 mutation status, defined by any FLT3 variant call in the MAF',
             capture.output(print(table(dat$FLT3_status))),
             paste('Log-rank p:', signif(logrank_p, 5)),
             paste('Cox HR (Mutated vs Not_mutated):', signif(flt3_tab$hazard_ratio[1], 4)),
             paste('Cox p:', signif(flt3_tab$p_value[1], 5)),
             '', 'Cox proportional-hazards diagnostic:',
             capture.output(print(survival::cox.zph(cox_flt3)))),
           file.path(out, '16_FLT3_interpretation.txt'))

cat('\nDONE! Files saved in: ', normalizePath(out), '\n', sep = '')
cat('Interpretation caution: associations are not proof of causation; Cox PH is an assumption.\n')
