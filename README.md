# BINF6250 - Project 03

# Introduction
This project implements ```GibbsMotifFinder()```, a Gibbs Sampling algorithm that discovers a short, shared DNA motif across a list of sequences without knowing in advance where that motif falls within each one. The function works by repeatedly holding out one sequence at a time, building a consensus (a Position Frequency Matrix and Position Weight Matrix) from every other sequence's current guess, scoring every possible position in the held-out sequence on both the forward and reverse-complement strands, and then updating that sequence's guess via weighted random selection, favoring higher-scoring positions without ever being locked into always picking the single best score. This loop repeats for a fixed number of rounds, after which a final consensus is built from every sequence's last guess and returned.

The function was tested on two datasets of different scale: a B. subtilis genome, pre-filtered down to ~837 promoter sequences known to contain the Shine-Dalgarno motif (AGGAGG), and a larger, unfiltered set of 90,061 real NRF1 ChIP-seq peak sequences.

The program ran successfully and all required functionality described in the assignment instructions was completed. ```GibbsMotifFinder()``` runs correctly on both datasets, and the NRF1 data-ingest step (reading nrf1_gibbs.fa into a list) was implemented as well. In the course of testing at the assignment specified scale (10,000 rounds), we identified that the algorithm's convergence behavior scales with dataset size, detailed below in Evidence and Reasoning and Reflection. This was not a flaw in the implementation. It reflects a real property of the algorithm as specified.

# Pseudocode
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

#NRF1 driver program psuedocode
    #Data ingest of nrf1_peaks
        #FOR each (name, seq) entry returned by get_fasta(nrf1_file):
            #add seq directly to nrf1_peaks (no filtering needed)

# Results Output
Driver Program Results:
 [100  97   0   1   0   1   0   0 102  86]
 [106 128   0 832 836   5 833 836 157 223]
 [240 160   0   0   0   0   1   0 206 205]]
12.54257153610685

NRF1 Driver Program Results:
 [24705 26885 25540 25618 26974 25988 24796 25378 25447 24288]
 [27611 25652 25988 26627 25359 24721 26184 26144 24764 27040]
 [18415 17989 18939 19074 18659 19434 19652 18780 19154 18835]]
0.16053914705479966

NRF1 dataset (90,061 sequences, k=10, 10,000 rounds): IC remained low (under 0.2) and was climbing only very slowly even after several thousand rounds, in contrast to the promoter dataset's strong, confident result. See Evidence and Reasoning below for why.

# Evidence and Reasoning
The promoter dataset's PFM shows near-unanimous agreement at several positions (834-836 out of ~837 sequences), and an IC of 12.57 out of a theoretical max of 20 for k=10 is strong evidence the algorithm correctly converges on a real, non-random motif when given enough rounds relative to its dataset size, consistent with the known Shine-Dalgarno pre-filter applied to that data.

The NRF1 dataset did not show this same convergence within the specified 10,000 rounds. We worked out why mathematically rather than assuming a bug. The algorithm updates exactly one randomly chosen sequence per round. With 10,000 rounds spread across 837 promoter sequences, each sequence gets updated about 12 times on average, which is enough to converge. Spread across 90,061 NRF1 sequences, each sequence gets updated only about 0.11 times on average, meaning the large majority of sequences never leave their random initial guess by the time the loop ends. This fully explains the low IC without needing to assume anything is wrong with the implementation itself, and we verified that the code runs correctly and without error on both datasets at full scale.

We did not perform additional automated testing beyond this comparison and the math check above, so we can't rule out that NRF1's real motif signal is also inherently weaker or noisier than the artificially clean Shine-Dalgarno practice case.  The too few round-count limitation is the explanation. Based on the math, we believe the IC could be brought up to a comparable level either by running substantially more rounds for the NRF1 dataset specifically (over 1,000,000), or by restructuring the loop to update a batch of several sequences per round instead of just one, so that the total number of individual sequence updates scales with dataset size rather than staying fixed at one per round.

# Successes
Our main success was producing a correct, fully working implementation that behaves as expected on a clean, well-understood dataset, and then using that same implementation to discover a non-obvious insight about how the algorithm's design assumptions interact with dataset size. Rather than treating the NRF1 dataset's low IC as a bug to chase, working through the update-count math let us explain it and turn it into a real discussion point about the algorithm's scalability.

Yulia's contribution to the NRF1 driver program was a big win. She added periodic sanity-check print statements every 500 rounds, which let us watch a run that took over two hours, tracking the IC score's progress along the way rather than waiting blind for a single final result and questioning if the result was a due to a bug. Given how long the NRF1 run took, this made a real, practical difference in being able to observe and reason about what was happening in real time instead of guessing afterward.
  
# Struggles
A significant early struggle was etting a working environment set up. Sorting out which packages (bamnostic, seqlogo) were required versus optional, and discovering that seqlogo's visual plotting depends on an external system tool (Ghostscript) not included with the Python package itself. Per the assignment's own note that this library is optional, we worked around this by verifying results directly through additional print statements in the ```PFM``` and ```pfm_ic()``` rather than the visual plot.

We also hit a FileNotFoundError after unzipping the provided data. The files had been automatically decompressed during download, so our code's paths (expecting .gz extensions) no longer matched what was on the disk. This turned out to be a simple path mismatch rather than file corruption, which we confirmed by inspecting the decompressed file's contents directly before assuming anything was broken and adjusted filenames.
  
# Personal Reflections
## Group Leader: Graziano Peregrino 


## Other member: Ildiko Polyak
The most interesting moment of the project for me was working on the NRF1 "Challenge Yourself" section. Initially, we assumed 10,000 rounds would be more than enough, maybe even overkill, for a Gibbs sampler to converge. Instead, running it on the full 90,061-sequence NRF1 dataset, the IC score barely moved even after thousands of rounds, nowhere near the strong convergence we'd seen on the smaller promoter dataset. That contradiction between our assumption and the result is what made me want to understand why, rather than just assume something was broken.

Working through the math, I realized the algorithm updates exactly one randomly chosen sequence per round, so the real question isn't "how many rounds," but "how many times does each individual sequence  get touched." With promoters, that ratio was generous; with NRF1, it meant the majority of sequences never left their random starting guess. That reframing, from "rounds" to "updates per sequence," was the key insight.

That curiosity led me to look into how real ChIP-seq motif discovery is handled, given that real datasets are often even larger than our 90,061-sequence practice set. What I found was that this isn't really a solved problem so much as a managed one. People do run these algorithms for very large numbers of iterations, frequently with periodic print statements and live progress plots similar to what Yulia added to our NRF1 driver, and the run is manually stopped once convergence is visually apparent, rather than relying on a fixed round count decided in advance (perhaps this is different in experienced labs and beter known sequences). In other words, what we had started doing, watching the IC progress and deciding by eye to understand what was hapenning, is a simplified version of the real practice, not just a workaround for a homework assignment. I found that an interesting discovery.

## Other member: Julianne Murthy


# Generative AI Appendix
Anthropic Claude opus 5.5 was used for this project.
Each team member used Claude differently throughout the project:
One member used Claude to explain Gibbs sampling theory, the project in its entirety, and help in understanding both the pseudocode and code.
One member used Claude to help write code only for converting the log2 score into a weighted random score to pick a new position in the ```GibbsMotifFinder()``` function.
One member used Claude to help insert more print statements throughtout the code in appropriate places that would help monitor the function of the program while it worked.