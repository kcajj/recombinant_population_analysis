- [ ] finding recombinant reads

    - [x] viterbi algorithm
        - [x] generate data on the basis of given probabilities
        - [x] understand viterbi script deeply
        - [x] make predictions on that data

    - [x] test on Aionostat data
        - [x] clones
            - [x] extract evidences from a clone sequence
            - [x] make a prediction with the viterbi algorithm
        - [x] reads
            - [x] extract evidences from a read
            - [x] make a prediction with the viterbi algorithm

    - [ ] parameter estimation
        - [ ] emission probabilities
            - [x] reads (analyse the reads of pure phages sequenced, infer the probabilities corresponding to each hidden state)
                - [x] build the msa of each read with the references and keep track of the frequency of evidences
            - [ ] clones
        - [ ] transition probabilities
            - [ ] maximum likelihood

    - [x] alignment between references
        - [x] snps and indels frequency

    - [x] gaps
        - [x] what happens if we don't consider gaps
        - [x] what happens if we consider gaps
        - [x] what is the best option? do not consider gaps

    - [ ] running on the whole dataset
        - [x] script to build msa of each read with references (using MAFFT) and analyse it
            - [x] map the reference index to the msa index
            - [x] cut the msa
            - [x] add the read to the msa
            - [x] we are using the fastest msa method possible. are we making a lot of error? probably yes
            - [x] we are letting MAFFT trim the reads to keep the length of the cut alignment the same. are we trimming a lot? we don't care
        - [ ] make the algorithm faster
            - [x] hybrid reference approach
                - [x] create the hybrid reference
                - [x] align the reads to the hybrid reference
                - [x] build the msa starting from minimap alignment (use method .get_aligned_pairs())
                - [x] extract the evidences with the msa and make the prediction with the viterbi algorithm
            - [x] set up a comparison framework between mafft method and fast method
                - [x] create a subset of 100 reads on which to test the methods
                    - [x] reads
                    - [x] parameter estimation
                - [x] run the two methods on the test set and see the differences
                    - [x] reads
                    - [x] parameter estimation
            - [ ] run on the whole dataset
                - [x] reads
                    - [x] set up a read_length threshold to make the algorithm faster
                - [ ] parameter estimation
        - [x] extract and store the recombination data
            - [x] array with position of recombination events in the reference genome
    - [x] representing whole-population recombination data

- [x] pipeline
    - [x] plots of directional recombination
        - [ ] put a threshold on the coverage
    - [x] plots of coverage

- [ ] extract the ranking of sites
- [ ] time dynamics plots
- [ ] correct coordinates