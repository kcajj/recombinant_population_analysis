

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
        evidences='results/evidence_arrays/{population}/{timestep}/{population}_{timestep}.tsv',
        coverage='results/coverage_arrays/{population}/{timestep}/{population}_{timestep}.npz'
    conda:
        'conda_envs/scientific_python.yml'
    params:
        length_threshold=5000
    shell:
        """
        python scripts/extract_evidence_arrays.py \
            --bam {input.bam} \
            --msa_refs {input.msa} \
            --evidences_out {output.evidences} \
            --coverage_out {output.coverage} \
            --length_threshold {params.length_threshold}
        """

rule all:
    input:
        evidence_arrays=expand(rules.evidence_arrays.output.evidences, population="P2", timestep="7")