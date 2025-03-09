rule runall:
    input:
        directory("results/m8_maps"),
        directory("results/m8_outflow_plots")

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
        "results/outflow_vs_separation_angles.csv"
    shell:
        "python3 scripts/make_all_m8_outflow_plots.py"