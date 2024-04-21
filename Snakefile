configfile: "config.yml"

#reads = 'data/reads_for_test/test_{population}_{timestep}.fastq.bgz'
reads = 'data/reads/{population}_{timestep}.fastq.gz'
references = 'data/references.fasta'

rule msa:
    input:
        reference = references
    output:
        msa = 'results/msa/msa_refs.fasta'
    conda:
        'conda_envs/read_mapping.yml'
    shell:
        """
        mafft --auto \
            {input.reference} \
            > {output.msa}
        """

rule hybrid_ref:
    input:
        msa = rules.msa.output.msa
    output:
        hybrid_ref = 'results/msa/hybrid_ref.fasta'
    conda:
        'conda_envs/scientific_python.yml'
    shell:
        """
        python scripts/hybrid_reference.py \
            --msa {input.msa} \
            --out {output.hybrid_ref}
        """

rule read_mapping:
    input:
        fastq=reads,
        ref=rules.hybrid_ref.output.hybrid_ref
    output:
        sam = 'results/alignments/{population}/{timestep}/{population}_{timestep}.sam',
        bam = 'results/alignments/{population}/{timestep}/{population}_{timestep}.bam',
        bai = 'results/alignments/{population}/{timestep}/{population}_{timestep}.bam.bai'
    conda:
        'conda_envs/read_mapping.yml'
    shell:
        """
        minimap2 -a -x map-ont {input.ref} {input.fastq} > {output.sam}
        samtools sort {output.sam} > {output.bam}
        samtools index {output.bam}
        """

rule evidence_arrays:
    input:
        bam=rules.read_mapping.output.bam,
        msa = rules.msa.output.msa
    output:
        evidences='results/evidence_arrays/{population}/{population}_{timestep}.tsv',
        coverage='results/coverage_arrays/{population}/{population}_{timestep}.npz'
    conda:
        'conda_envs/scientific_python.yml'
    params:
        length_threshold = config["length_threshold"]
    shell:
        """
        python scripts/extract_evidence_arrays.py \
            --bam {input.bam} \
            --msa_refs {input.msa} \
            --evidences_out {output.evidences} \
            --coverage_out {output.coverage} \
            --length_threshold {params.length_threshold}
        """

rule recombination_array:
    input:
        evidences=rules.evidence_arrays.output.evidences,
        msa = rules.msa.output.msa
    output:
        recombination='results/genomewide_recombination_arrays/{population}/{population}_{timestep}.npz'
    conda:
        'conda_envs/scientific_python.yml'
    params:
        cores = config["cores"],
        initial_probability = config["initial_probability"]["A"]+","+config["initial_probability"]["B"],
        transition_probability = config["transition_probability"]["A"]["A"]+","+config["transition_probability"]["A"]["B"]+"/"+config["transition_probability"]["B"]["A"]+","+config["transition_probability"]["B"]["B"],
        emission_probability = config["emission_probability"]["A"][0]+","+config["emission_probability"]["A"][1]+","+config["emission_probability"]["A"][2]+"/"+config["emission_probability"]["B"][0]+","+config["emission_probability"]["B"][1]+","+config["emission_probability"]["B"][2]
    shell:
        """
        python scripts/hmm_recombination_array.py \
            --evidences {input.evidences} \
            --msa_refs {input.msa} \
            --out {output.recombination} \
            --cores {params.cores} \
            --initial_p {params.initial_probability}\
            --transition_p {params.transition_probability}\
            --emission_p {params.emission_probability}
        """

rule plot_recombination_array:
    input:
        recombination=rules.recombination_array.output.recombination,
        coverage=rules.evidence_arrays.output.coverage
    output:
        plots='results/plots/{population}/{population}_{timestep}.png'
    conda:
        'conda_envs/scientific_python.yml'
    shell:
        """
        python scripts/plot_recombination_array.py \
            --recombination {input.recombination} \
            --coverage {input.coverage} \
            --out {output.plots}
        """

rule all:
    input:
        plots=expand(rules.plot_recombination_array.output.plots, population="P2", timestep="7")