shell.executable("bash")

PANDOC_DIR = "/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"

rule all:
    input:
        # Main text figures
        "figures/fig2_growth_by_temp_hab.png",       # Fig 2: growth ~ max temp x habitat
        "figures/fig4_abscort_by_temp_hab.png",      # Fig 3: cort ~ max temp x habitat
        "figures/fig5_provis_by_temp_hab.png",       # Fig 4: provisioning ~ temp x habitat
        # Supplement figures
        "figures/dag3.png",                           # Fig S1: DAG
        "figures/max-weightedmean_outside.png",      # Fig S2: temp anomaly by land cover
        "figures/fig3_survival_by_temp_hab.png",     # Fig S3: survival ~ max temp x habitat
        "figures/figs2_s1cort_by_priordayt_hab.png", # Fig S4: baseline cort ~ prior-day temp
        "figures/fig6_growth_by_temp_hab.png",       # Fig S5: growth ~ provisioning + cort
        "figures/s1bypriordaymaxhhixhab_TRES.png",   # Fig S6: TRES baseline cort ~ heat index
        # Supplement tables (all from publication_tables.Rmd)
        "figures/full_samp_size_tbl.html",
        "figures/full_landcover_contrast_tbl.html",
        "figures/trend_tbl.html",
        "figures/int_tbl.html",
        "figures/growth_coef.html",
        "figures/survival_coef.html",
        "figures/s1_coef.html",
        "figures/abs_coef.html",
        "figures/s2_coef.html",
        "figures/provis_coef.html",
        "figures/g_s1_provis_coef.html",
        "figures/g_abs_provis_coef.html",
        "figures/val_tbl.html",
        "figures/tukeytable_tempbylanduse.html",
        "figures/tresprovistrend.html",
        "figures/val_tmaxxjday_tbl.html",


rule model_growth:
    input:
        data   = "data/growth_cort_provis_manytempmeasures.rds",
        provis = "data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds",
        script = "code/model_growth.R",
    output:
        "data/models_growth.RData",
    resources:
        runtime = "4h",
    shell:
        "Rscript code/model_growth.R"


rule model_cort:
    input:
        data   = "data/growth_cort_provis_manytempmeasures.rds",
        script = "code/model_cort.R",
    output:
        "data/models_cort.RData",
    resources:
        runtime = "4h",
    shell:
        "Rscript code/model_cort.R"


rule model_survival:
    input:
        data   = "data/growth_cort_provis_manytempmeasures.rds",
        surv   = "data/survival_attempt.rds",
        script = "code/model_survival.R",
    output:
        "data/models_survival.RData",
    resources:
        runtime = "2h",
    shell:
        "Rscript code/model_survival.R"


rule model_provis:
    input:
        data   = "data/growth_cort_provis_manytempmeasures.rds",
        provis = "data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds",
        script = "code/model_provis.R",
    output:
        "data/models_provis.RData",
    resources:
        runtime = "2h",
    shell:
        "Rscript code/model_provis.R"


rule model_seasonal:
    input:
        data   = "data/growth_cort_provis_manytempmeasures.rds",
        surv   = "data/survival_attempt.rds",
        provis = "data/provis_manytempmeasures.rds",
        script = "code/seasonal_sensitivity.R",
    output:
        "data/models_seasonal.RData",
    resources:
        runtime = "4h",
    shell:
        "Rscript code/seasonal_sensitivity.R"


rule figures:
    input:
        models_growth   = "data/models_growth.RData",
        models_cort     = "data/models_cort.RData",
        models_survival = "data/models_survival.RData",
        models_provis   = "data/models_provis.RData",
        temp            = "data/temp.rds",
        script          = "code/figures.R",
    output:
        "figures/fig2_growth_by_temp_hab.png",
        "figures/fig6_growth_by_temp_hab.png",
        "figures/fig4_abscort_by_temp_hab.png",
        "figures/figs2_s1cort_by_priordayt_hab.png",
        "figures/s1bypriordaymaxhhixhab_TRES.png",
        "figures/fig3_survival_by_temp_hab.png",
        "figures/fig5_provis_by_temp_hab.png",
        "figures/dag3.png",
        "figures/max-weightedmean_outside.png",
    resources:
        runtime = "30m",
    shell:
        "Rscript code/figures.R"


rule tables:
    input:
        models_growth   = "data/models_growth.RData",
        models_cort     = "data/models_cort.RData",
        models_survival = "data/models_survival.RData",
        models_provis   = "data/models_provis.RData",
        models_seasonal = "data/models_seasonal.RData",
        temp_data       = "data/temp.rds",
        script          = "code/publication_tables.Rmd",
    output:
        "figures/full_samp_size_tbl.html",
        "figures/full_landcover_contrast_tbl.html",
        "figures/trend_tbl.html",
        "figures/int_tbl.html",
        "figures/growth_coef.html",
        "figures/survival_coef.html",
        "figures/s1_coef.html",
        "figures/abs_coef.html",
        "figures/s2_coef.html",
        "figures/provis_coef.html",
        "figures/g_s1_provis_coef.html",
        "figures/g_abs_provis_coef.html",
        "figures/val_tbl.html",
        "figures/tukeytable_tempbylanduse.html",
        "figures/tresprovistrend.html",
        "figures/val_tmaxxjday_tbl.html",
    resources:
        runtime = "4h",
    shell:
        "export RSTUDIO_PANDOC='{PANDOC_DIR}' && Rscript -e \"rmarkdown::render('code/publication_tables.Rmd')\"".format(PANDOC_DIR=PANDOC_DIR)
