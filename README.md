# BINF6250 - Project 03

# Introduction
This project main focus is to identify conserved DNA motifs using Gibbs sampling and also apply a probabilistic motif finding algorithm to biological sequence data and also identify a recurring pattern that may represent a potential transcription factor binding motif.  

# Pseudocode
## Gibbs Motif Finder Pseudocode

```text
FUNCTION GibbsMotifFinder(seqs, k, seed):

    # Setup
    SET random seed

    Motifs ← empty list

    FOR each sequence IN seqs:
        start ← random integer from 0 to LENGTH(sequence) - k
        motif ← sequence[start : start + k]
        APPEND motif TO Motifs

    # Main Gibbs Sampling Loop
    REPEAT up to 10,000 times:

        # Step 1: Randomly choose one sequence to leave out
        i ← random integer from 0 to LENGTH(seqs) - 1

        # Step 2: Build PFM and PWM from all motifs except Motifs[i]
        RemainingMotifs ← all motifs in Motifs except motif at index i

        pfm ← build_pfm(RemainingMotifs, k)
        pwm ← build_pwm(pfm)

        # Step 3: Score every possible k-mer in the left-out sequence
        candidate_motifs ← empty list
        candidate_scores ← empty list

        sequence ← seqs[i]

        FOR position FROM 0 TO LENGTH(sequence) - k:

            forward_chunk ← sequence[position : position + k]
            reverse_chunk ← reverse_complement(forward_chunk)

            forward_score ← score_kmer(forward_chunk, pwm)
            reverse_score ← score_kmer(reverse_chunk, pwm)

            IF forward_score >= reverse_score:
                APPEND forward_chunk TO candidate_motifs
                APPEND forward_score TO candidate_scores

            ELSE:
                APPEND reverse_chunk TO candidate_motifs
                APPEND reverse_score TO candidate_scores

        # Step 4: Convert log2 scores into positive weights
        weights ← empty list

        FOR each score IN candidate_scores:
            weight ← 2 ^ score
            APPEND weight TO weights

        # Normalize weights into probabilities
        total_weight ← SUM(weights)

        probabilities ← each weight / total_weight

        # Step 5: Sample a new motif using the probabilities
        chosen_index ← weighted random choice from
                       0 to LENGTH(candidate_motifs) - 1
                       using probabilities

        Motifs[i] ← candidate_motifs[chosen_index]

    # Final Result
    final_pfm ← build_pfm(Motifs, k)

    RETURN final_pfm
```
## NRF1 Driver Program Pseudocode

```text

nrf1_file ← path to NRF1 FASTA file

# Data Ingest: NRF1 Peaks
nrf1_peaks ← empty list

FOR each (name, seq) returned by get_fasta(nrf1_file):
    # No sequence filtering is required
    APPEND seq TO nrf1_peaks
```
# Successes
- One if the main successes while working on this project was understanding how Gibbs sampling can be applied on the motif discovery, Being able to organize the algorithm and get a good output.
- Another success I would bring was how the team put together a great project and work together to get the challenge part of the homework done even though to run the dataset took more than 2 hours but we were able to pull a correct output.
  
# Struggles
- One of the main struggles was organizing the Gibbs sampling loop and understanding which values to update during each iteration.
- Another point that we had was working with reverse complements and analyzing the output to make sure it was following the correct path.
- Working on the additional challenge brought us a bit of work since the package suggested was not providing an output that would bring some insights.
  
# Personal Reflections
## Group Leader: Graziano Peregrino 

## Other member: Ildiko Polyak

## Other member: Julianne Murthy

# Generative AI Appendix

The appendix entry must contain:
Description of which generative AI was used and its version.

The entire prompt that was used to generate the content.

An explanation of how it was used .


A justification for why generative AI was used.


