import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pysam

def get_hyb_ref_map(bam_file):
    map_hyb_ref={}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary):
                alignment=read.get_aligned_pairs()
                for hyb,ref in alignment:
                    if hyb==None:
                        continue
                    if hyb in map_hyb_ref.keys():
                        if map_hyb_ref[hyb]==None:
                            map_hyb_ref[hyb]=ref
                        else:
                            continue
                    else:
                        map_hyb_ref[hyb]=ref
    return map_hyb_ref

if __name__ == "__main__":
    populations=["P2","P3"]
    timesteps=["1","3","5","7"]

    recombination_folder="results/genomewide_recombination"
    coverage_folder="results/coverage_arrays"
    #refs_msa_path=
    #output_path="thesis/plots"

    threshold_frequency=0.01

    bas51_map=get_hyb_ref_map("thesis/alignments/hybrid_on_bas51.bam")
    bas54_map=get_hyb_ref_map("thesis/alignments/hybrid_on_bas54.bam")

    print(" ")
    print("bas51")
    print(bas51_map[98100])

    print(" ")
    print("hybrid ref")
    print(list(bas51_map.keys())[list(bas51_map.values()).index(95262)])
    print(list(bas51_map.keys())[list(bas51_map.values()).index(123100)])
    print(list(bas51_map.keys())[list(bas51_map.values()).index(125434)])



    print(" ")
    print("bas54")
    print(bas54_map[96213])
    print(bas54_map[124146])
    print(bas54_map[126481])