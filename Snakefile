rule runall:
    input:
        # directory("results/stat_test"),
        # directory("results/histogram"),
        # "results/m0_outflow_maps.pdf",
        # "results/all_m8_maps.pdf"
        # directory("results/figs_for_paper")

rule prepare_tables:
    input:
        "data/input/reynolds2024_perseus_1.txt"
        "data/input/tobin2018_perseus_2.txt",
        "data/input/tobin2022_orion.txt",
        "data/input/notes.csv",
        "data/input/tobin2022_orion_pairings.txt",
        "data/input/tobin2022_perseus_pairings.txt"
    output:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv",
        "data/output/data_by_field.csv",
        "data/output/data_by_outflow.csv",
        "results/tables/by-field.tex",
        "results/tables/by-outflow.tex"
    shell:
        "python3 scripts/prepare_tables.py"

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

rule make_figs_for_paper:
    input:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv"
    output:
        directory("results/figs_for_paper")
    shell:
        "python3 scripts/make_figs_for_paper.py"

rule perform_stat_test:
    input:
        "data/output/outflow_data.csv"
    output:
        directory("results/stat_test")
    shell:
        "python3 scripts/stat_test.py"

rule make_histogram:
    input:
        "data/output/outflow_data.csv"
    output:
        directory("results/histogram")
    shell:
        "python3 scripts/make_histogram.py"