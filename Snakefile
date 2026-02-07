rule runall:
    input:
        "results/tables/by-field.tex",
        "results/tables/by-outflow.tex",
        "results/tables/by-outflow-full.tex",
        "results/tables/offsets.tex",
        "results/figs_for_paper/appendix-1.pdf",
        "results/figs_for_paper/appendix-2.pdf",
        "results/figs_for_paper/appendix-3.pdf",
        "results/figs_for_paper/appendix-4.pdf",
        "results/figs_for_paper/appendix-5.pdf",
        "results/figs_for_paper/appendix-6.pdf",
        "results/figs_for_paper/appendix-7.pdf",
        "results/figs_for_paper/appendix-8.pdf",
        "results/figs_for_paper/appendix-9.pdf",
        "results/figs_for_paper/appendix-10.pdf",
        "results/figs_for_paper/appendix-11.pdf",
        "results/figs_for_paper/appendix-12.pdf",
        "results/figs_for_paper/appendix-13.pdf",
        "results/figs_for_paper/fig_1.pdf",
        "results/figs_for_paper/fig_2.pdf",
        "results/m0_outflow_maps.pdf",
        "results/all_m8_maps.pdf",
        "results/stat_test/histogram.pdf",
        "results/stat_test/DeltaPA_cumulat_deg.pdf",
        "results/master_reference.pdf"

rule prepare_tables:
    input:
        "data/input/reynolds2024_perseus_1.txt",
        "data/input/tobin2018_perseus_2.txt",
        "data/input/tobin2022_orion.txt",
        "data/input/notes.csv",
        "data/input/offsets.csv",
        "data/input/tobin2022_orion_pairings.txt",
        "data/input/tobin2022_perseus_pairings.txt"
    output:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv",
        "data/output/data_by_field.csv",
        "data/output/data_by_outflow.csv",
        "results/tables/by-field.tex",
        "results/tables/by-outflow.tex",
        "results/tables/by-outflow-full.tex",
        "results/tables/offsets.tex"
    shell:
        "python3 scripts/prepare_tables.py"

rule perform_stat_test:
    input:
        "data/output/outflow_data.csv"
    output:
        "results/stat_test/DeltaPA_cumulat_deg.pdf",
        "results/stat_test/histogram.pdf"
    shell:
        "python3 scripts/stat_test-combined.py"

rule make_figs_for_paper:
    input:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv"
    output:
        "results/figs_for_paper/appendix-1.pdf",
        "results/figs_for_paper/appendix-2.pdf",
        "results/figs_for_paper/appendix-3.pdf",
        "results/figs_for_paper/appendix-4.pdf",
        "results/figs_for_paper/appendix-5.pdf",
        "results/figs_for_paper/appendix-6.pdf",
        "results/figs_for_paper/appendix-7.pdf",
        "results/figs_for_paper/appendix-8.pdf",
        "results/figs_for_paper/appendix-9.pdf",
        "results/figs_for_paper/appendix-10.pdf",
        "results/figs_for_paper/appendix-11.pdf",
        "results/figs_for_paper/appendix-12.pdf",
        "results/figs_for_paper/appendix-13.pdf",
        "results/figs_for_paper/fig_1.pdf",
        "results/figs_for_paper/fig_2.pdf"
    shell:
        "python3 scripts/make_figs_for_paper.py"

rule make_all_m0_outflow_maps:
    input:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv"
    output:
        "results/m0_outflow_maps.pdf"
    shell:
        "python3 scripts/make_all_m0_outflow_maps.py"

rule make_all_m8_maps:
    input:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv"
    output:
        "results/all_m8_maps.pdf"
    shell:
        "python3 scripts/make_all_m8_maps.py"

rule master_reference_doc:
    input:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv"
    output:
        "results/master_reference.pdf"
    shell:
        "python3 scripts/master_reference_doc.py"