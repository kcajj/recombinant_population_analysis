import numpy as np
import pysam

if __name__ == "__main__":
    population='P2'
    timestep='7'
    #bam_file = f'results/alignments/{population}_{timestep}.bam'
    bam_file = f'data/test/hybrid_test_{population}_{timestep}.bam'

    imported = pysam.AlignmentFile(bam_file, mode = 'rb')  # 'rb' ~ read bam
    coverage_arrays = imported.count_coverage("hybrid_ref")
    coverage=np.zeros(len(coverage_arrays[0]))
    for i in range(len(coverage_arrays[0])):
        for j in range(4):
            coverage[i]+=coverage_arrays[j][i]

    #np.savez("results/coverage_arrays/{population}_{timestep}.npz",coverage)
    np.savez(f"data/test/coverage_test_{population}_{timestep}.npz",coverage)