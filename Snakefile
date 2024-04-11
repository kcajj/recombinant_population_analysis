

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

rule all:
    input:
        mapped_reads=expand(rules.read_mapping.output.bai, population="P2", timestep="7")