#Original function from project03.ipynb that we need to write code for
def GibbsMotifFinder (seqs, k, seed=None):
    '''
    Function to find a pfm from a list of strings using a Gibbs sampler
    
    Args: 
        seqs (str list): a list of sequences, not necessarily in same lengths
        k (int): the length of motif to find
        seed (int, default=None): seed for np.random

    Returns:
        pfm (numpy array): dimensions are 4xlength
    '''
    # Use rng to make random samples/selections/numbers
    # Example: randint = rng.integer(1, 10)
    random.seed(seed)
    rng = np.random.default_rng(seed)

    pass

#Pseudocode
#FUNCTION GibbsMotifFinder(seqs, k, seed)

    #Setup (before the loop)
        #Initialize random seed:  Already done in the code above
        #Motifs = a list that will hold the current guessed chunk (string) for each sequence
        #FOR each sequence in seqs:
            #pick a random starting position between 0 and sequence length - k
            #grab the k-length chunk starting at that position
            #save this chunk as the current guess for this sequence

    #Main Loop (5 basics steps of Gibbs Sampling)
        #REPEAT up to 10,000 times:
            #Randomly pick one sequence to leave out
            #Build a consensus (PFM/PWM) from all guesses in Motifs except sequence i
                #pfm = build_pfm(all Motifs except Motifs[i])
                #pwm = build_pwm(pfm)
            #Score every possible position in the "left out" sequence i, checking both strands:
                #FOR each valid starting position in sequence i:
                    #forward_chunk = the k-length chunk at this position
                    #reverse_chunk = reverse_complement(forward_chunk)
                    #forward_score = score_kmer(forward_chunk, pwm)
                    #reverse_score = score_kmer(reverse_chunk, pwm)
                    #keep the higher of forward_score/reverse_score and remember which it is
                        #this gives a list of best scores (and matching chunks), one per position
            
            #Convert scores into weights, since scores are log2-based and can be negative
                #weights = 2^(score) for each score in the list (undoes log2)
                #turns scores back into positive usable probability
            
            #Pick a new position for sequence i, randomly weighted by score
                #use np.random.choice, weighted by the converted weights and not picking single highest scoring position
            
            #Update Motifs[i] to the chunck (forward or reverse complement) that corresponds to newly picked position

    #After the Loop
        #Build one final PFM table using all Motifs 
            #final_pfm = build_pfm(Motifs)
        #Return final_PFM



#Original Program from The NRF1 Driver Program, need to write pseudocode for "TODO" portion
# Here we test your Gibbs sampler.
# You do not need to edit this or the section below. This is the Driver program

#read promoters, store in a list of strings
nrf1_file="<path/to/file.fa>" 

nrf1_peaks = []

#TODO: Data ingest of nrf1_peaks

# Run the gibbs sampler:
nrf1_pfm = GibbsMotifFinder(nrf1_peaks, 10)

# Plot the final pfm that is generated: 
seqlogo.seqlogo(seqlogo.CompletePm(pfm = nrf1_pfm.T))


#PSeudocode 
#missing portions of NRF1 Driver Program
    #nrf1_file needs the placeholder path replaced with real one
    #"TODO" 
        #Data ingest of nrf1_peaks
            #FOR each (name, seq) entry returned by get_fasta(nrf1_file):
            #add seq directly to nrf1_peaks (no filtering needed)