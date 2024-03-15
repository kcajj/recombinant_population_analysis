- [ ] finding recombinant reads

    - [ ] viterbi algorithm
        - [x] generate data on the basis of given probabilities
        - [ ] understand viterbi script deeply
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
            - [ ] reads (analyse the reads of pure phages sequenced, infer the probabilities corresponding to each hidden state)
                - [x] build the msa of each read with the references and keep track of the frequency of evidences
            - [ ] clones
        - [ ] transition probabilities
            - [ ] maximum likelihood

    - [x] alignment between references
        - [x] snps and indels frequency

    - [x] gaps
        - [x] what happens if we don't consider gaps
        - [x] what happens if we consider gaps
        - [x] what is the best option?

    - [ ] running on the whole dataset
        - [ ] script to build msa of each read with references and analyse it
            - [x] map the reference index to the msa index
            - [x] cut the msa
            - [x] add the read to the msa
            - [ ] we are using the fastest msa method possible. are we making a lot of error?
            - [ ] we are letting MAFFT trim the reads to keep the length of the cut alignment the same. are we trimming a lot?
        - [ ] make the algorithm faster
        - [ ] extract and store the recombination data
            - [ ] position of recombination on the reference genome

- [ ] representing whole-population recombination data