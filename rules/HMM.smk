alignments_config = config["alignments"]
HMM_config = config["HMM_parameters"]

references = in_fld + "/references.fasta",
reads = in_fld + "/reads/{replicate}_{timestep}.fastq.gz"
#reads = 'data/reads/{population}_{timestep}.fastq.gz'

rule msa:
    input:
        reference = references
    output:
        msa = out_fld + '/msa/msa_refs.fasta'
    conda:
        '../conda_envs/read_mapping.yml'
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
        hybrid_ref = out_fld + '/msa/hybrid_ref.fasta'
    conda:
        '../conda_envs/sci_py.yml'
    shell:
        """
        python scripts/hybrid_reference.py \
            --msa {input.msa} \
            --out {output.hybrid_ref}
        """

rule read_mapping:
    input:
        fastq = reads,
        ref = rules.hybrid_ref.output.hybrid_ref
    output:
        sam = out_fld + '/alignments/{replicate}/{timestep}.sam',
        bam = out_fld + '/alignments/{replicate}/{timestep}.bam',
        bai = out_fld + '/alignments/{replicate}/{timestep}.bam.bai'
    conda:
        '../conda_envs/read_mapping.yml'
    shell:
        """
        minimap2 -a -x map-ont {input.ref} {input.fastq} > {output.sam}
        samtools sort {output.sam} > {output.bam}
        samtools index {output.bam}
        """

rule evidence_arrays:
    input:
        bam = rules.read_mapping.output.bam,
        msa = rules.msa.output.msa
    output:
        evidences = out_fld + '/evidence_arrays/{replicate}/{timestep}.tsv',
    conda:
        '../conda_envs/sci_py.yml'
    params:
        length_threshold = alignments_config["length_threshold"]
    shell:
        """
        python scripts/extract_evidence_arrays.py \
            --bam {input.bam} \
            --msa_refs {input.msa} \
            --evidences_out {output.evidences} \
            --length_threshold {params.length_threshold}
        """

rule prediction_arrays:
    input:
        evidences = rules.evidence_arrays.output.evidences,
        hybrid_ref = rules.hybrid_ref.output.hybrid_ref
    output:
        predictions = out_fld + '/prediction_arrays/{replicate}/{timestep}.tsv',
        coverage = out_fld + '/coverage_arrays/{replicate}/{timestep}.npz'
    conda:
        '../conda_envs/sci_py.yml'
    params:
        cores = HMM_config["cores"],
        initial_probability = HMM_config["initial_probability"]["A"]+","+HMM_config["initial_probability"]["B"],
        transition_probability = HMM_config["transition_probability"]["A"]["A"]+","+HMM_config["transition_probability"]["A"]["B"]+"/"+HMM_config["transition_probability"]["B"]["A"]+","+HMM_config["transition_probability"]["B"]["B"],
        emission_probability = HMM_config["emission_probability"]["A"][0]+","+HMM_config["emission_probability"]["A"][1]+","+HMM_config["emission_probability"]["A"][2]+"/"+HMM_config["emission_probability"]["B"][0]+","+HMM_config["emission_probability"]["B"][1]+","+HMM_config["emission_probability"]["B"][2]
    shell:
        """
        python scripts/hmm_prediction_arrays.py \
            --evidences {input.evidences} \
            --hybrid_ref {input.hybrid_ref} \
            --predictions_out {output.predictions} \
            --coverage_out {output.coverage} \
            --cores {params.cores} \
            --initial_p {params.initial_probability}\
            --transition_p {params.transition_probability}\
            --emission_p {params.emission_probability}
        """

rule genomewide_recombination_array:
    input:
        predictions = rules.prediction_arrays.output.predictions,
        hybrid_ref = rules.hybrid_ref.output.hybrid_ref
    output:
        genomewide_recombination = out_fld + '/genomewide_recombination/{replicate}/{timestep}.npz',
        genomewide_recombination_01 = out_fld + '/genomewide_recombination/{replicate}/{timestep}_01.npz',
        genomewide_recombination_10 = out_fld + '/genomewide_recombination/{replicate}/{timestep}_10.npz'
    conda:
        '../conda_envs/sci_py.yml'
    shell:
        """
        python scripts/genomewide_recombination.py \
            --predictions {input.predictions} \
            --hybrid_ref {input.hybrid_ref} \
            --out {output.genomewide_recombination} \
            --out_01 {output.genomewide_recombination_01} \
            --out_10 {output.genomewide_recombination_10}
        """

rule HMM_all:
    input:
        predictions = expand(rules.prediction_arrays.output.predictions, replicate=HMM["replicates"], timestep=HMM["timesteps"]),
        genomewide_recombination = expand(rules.genomewide_recombination_array.output.genomewide_recombination, replicate=HMM["replicates"], timestep=HMM["timesteps"])
    output:
        finish = out_fld + '/HMM_done.txt'
    shell:
        """
        touch {output.finish}
        """