rule runall:
    input:
        # directory("results/m8_maps"),
        # directory("results/m0_outflow_contours"),
        directory("results/m8_outflow_plots"),
        directory("results/m0_outflow_plots"),
        directory("results/stat_test"),
        directory("results/histogram"),
        directory("results/m0_multi_outflow_plots")

rule prepare_tables:
    input:
        "data/input/distances.txt",
        "data/input/notes.csv"
    output:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv"
    shell:
        "python3 scripts/prepare_tables.py"

rule make_all_m8_maps:
    input:
        "data/output/source_info.csv"
    output:
        directory("results/m8_maps")
    shell:
        "python3 scripts/make_all_m8_maps.py"

rule make_all_m8_outflow_plots:
    input:
        "data/output/outflow_data.csv"
    output:
        directory("results/m8_outflow_plots"),
    shell:
        "python3 scripts/make_all_m8_outflow_plots.py"

rule make_all_m0_outflow_plots:
    input:
        "data/output/outflow_data.csv"
    output:
        directory("results/m0_outflow_plots"),
    shell:
        "python3 scripts/make_all_m0_outflow_plots.py"

rule make_all_m0_outflow_contours:
    input:
        "data/output/outflow_data.csv"
    output:
        directory("results/m0_outflow_contours"),
    shell:
        "python3 scripts/make_all_m0_outflow_contours.py"

rule make_all_m0_multi_outflow_plots:
    input:
        "data/output/outflow_data.csv",
        "data/output/source_info.csv"
    output:
        directory("results/m0_multi_outflow_plots"),
    shell:
        "python3 scripts/make_all_m0_multi_outflow_plots.py"

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