import os

plot_config = config["plots"]
input_timesteps = ""
for timestep in HMM["timesteps"]:
    input_timesteps += timestep + ","
input_references = ""
for reference in HMM["references"]:
    input_references += reference + ","
for replicate in HMM["replicates"]:
    os.makedirs(f"{out_fld}/genomewide_recombination/{replicate}/", exist_ok=True)
    os.makedirs(f"{out_fld}/coverage_arrays/{replicate}/", exist_ok=True)


rule plot_coverage_dynamics:
    input:
        hybrid_ref=rules.hybrid_ref.output.hybrid_ref,
        coverage_folder=directory(out_fld + "/coverage_arrays/{replicate}/"),
        wait=rules.HMM_all.output.finish,
    output:
        plots=out_fld + "/plots/coverage_dynamics/coverage_{replicate}.pdf",
    params:
        timesteps=input_timesteps,
        references=input_references,
        coverage_threshold=plot_config["coverage_threshold"],
    conda:
        "../conda_envs/sci_py.yml"
    shell:
        """
        python scripts/time_dynamics_coverage.py \
            --hybrid_ref {input.hybrid_ref} \
            --coverage {input.coverage_folder} \
            --timesteps {params.timesteps} \
            --references {params.references} \
            --coverage_threshold {params.coverage_threshold} \
            --out {output.plots}
        """


rule plot_recombination_dynamics:
    input:
        recombination_folder=directory(out_fld + "/genomewide_recombination/{replicate}/"),
        coverage_folder=directory(out_fld + "/coverage_arrays/{replicate}/"),
        wait=rules.HMM_all.output.finish,
    output:
        plots=out_fld + "/plots/recombination_dynamics/recombination_{replicate}.pdf",
    params:
        timesteps=input_timesteps,
        references=input_references,
        coverage_threshold=plot_config["coverage_threshold"],
    conda:
        "../conda_envs/sci_py.yml"
    shell:
        """
        python scripts/time_dynamics_recombination.py \
            --recombination {input.recombination_folder} \
            --coverage {input.coverage_folder} \
            --timesteps {params.timesteps} \
            --references {params.references} \
            --coverage_threshold {params.coverage_threshold} \
            --out {output.plots}
        """


rule unique_plot:
    input:
        hybrid_ref=rules.hybrid_ref.output.hybrid_ref,
        recombination_folder=directory(out_fld + "/genomewide_recombination/{replicate}/"),
        coverage_folder=directory(out_fld + "/coverage_arrays/{replicate}/"),
        wait=rules.HMM_all.output.finish,
    output:
        plots=out_fld + "/plots/unique_plots/unique_{replicate}.pdf",
    params:
        timesteps=input_timesteps,
        references=input_references,
        coverage_threshold=plot_config["coverage_threshold"],
    conda:
        "../conda_envs/sci_py.yml"
    shell:
        """
        python scripts/unique_plot.py \
            --hybrid_ref {input.hybrid_ref} \
            --recombination {input.recombination_folder} \
            --coverage {input.coverage_folder} \
            --timesteps {params.timesteps} \
            --references {params.references} \
            --coverage_threshold {params.coverage_threshold} \
            --out {output.plots}
        """

rule plot_all:
    input:
        coverage_dynamics=expand(rules.plot_coverage_dynamics.output.plots, replicate=HMM["replicates"]),
        recombination_dynamics=expand(rules.plot_recombination_dynamics.output.plots, replicate=HMM["replicates"]),
        unique_plots=expand(rules.unique_plot.output.plots, replicate=HMM["replicates"]),
        finish=rules.HMM_all.output.finish,
    shell:
        """
        rm {input.finish}
        """
