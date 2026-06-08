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


rule growth_figures:
    input:
        data = "data/growth_cort_provis_manytempmeasures.rds",
        script = "code/growth_graphs.R",
        model_eval = "code/growth_model_eval.R",
    output:
        fig2  = "figures/fig2_growth_by_temp_hab.png",
        figs5 = "figures/fig6_growth_by_temp_hab.png",
    resources:
        runtime = "2h",
    shell:
        "Rscript code/growth_graphs.R"


rule cort_figures:
    input:
        data = "data/growth_cort_provis_manytempmeasures.rds",
        script = "code/cort_graphs.R",
        model_eval = "code/cort_model_eval.R",
    output:
        fig3  = "figures/fig4_abscort_by_temp_hab.png",
        figs4 = "figures/figs2_s1cort_by_priordayt_hab.png",
        figs6 = "figures/s1bypriordaymaxhhixhab_TRES.png",
    resources:
        runtime = "2h",
    shell:
        "Rscript code/cort_graphs.R"


rule survival_figure:
    input:
        data = "data/survival_attempt.rds",
        script = "code/survival_model_eval.R",
    output:
        "figures/fig3_survival_by_temp_hab.png",
    resources:
        runtime = "1h",
    shell:
        "Rscript code/survival_model_eval.R"


rule provis_figure:
    input:
        data = "data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds",
        script = "code/provis_model_eval.R",
    output:
        "figures/fig5_provis_by_temp_hab.png",
    resources:
        runtime = "1h",
    shell:
        "Rscript code/provis_model_eval.R"


rule dag_figure:
    input:
        data = "data/growth.rds",
        script = "code/dag.Rmd",
    output:
        "figures/dag3.png",
    resources:
        runtime = "10m",
    shell:
        "export RSTUDIO_PANDOC='{PANDOC_DIR}' && Rscript -e \"rmarkdown::render('code/dag.Rmd')\"".format(PANDOC_DIR=PANDOC_DIR)


rule temp_anomaly_figure:
    input:
        temp   = "data/temp.rds",
        morph  = "data/morph.rds",
        canopy = "data/canopy_cover_hemispheR.rds",
        growth = "data/growth.rds",
        script  = "code/analysis.Rmd",
        helpers = "code/helper_functions.R",
    output:
        "figures/max-weightedmean_outside.png",
    resources:
        runtime = "30m",
    shell:
        "export RSTUDIO_PANDOC='{PANDOC_DIR}' && Rscript -e \"rmarkdown::render('code/analysis.Rmd')\"".format(PANDOC_DIR=PANDOC_DIR)


rule tables:
    input:
        growth_data       = "data/growth_cort_provis_manytempmeasures.rds",
        survival_data     = "data/survival_attempt.rds",
        provis_data       = "data/provis_with_attempt_1h_combined_mobilenetv3-original_dataset.h5.rds",
        provis_other_data = "data/provis_manytempmeasures.rds",
        growth_script     = "code/growth_model_eval.R",
        survival_script   = "code/survival_model_eval.R",
        cort_script       = "code/cort_model_eval.R",
        provis_script     = "code/provis_model_eval.R",
        growth_other      = "code/growth_model_eval_othertempmeasures.R",
        survival_other    = "code/survival_model_eval_othertempmeasures.R",
        provis_other      = "code/provis_model_eval_othertempmeasures.R",
        script            = "code/publication_tables.Rmd",
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
    resources:
        runtime = "4h",
    shell:
        "export RSTUDIO_PANDOC='{PANDOC_DIR}' && Rscript -e \"rmarkdown::render('code/publication_tables.Rmd')\"".format(PANDOC_DIR=PANDOC_DIR)
