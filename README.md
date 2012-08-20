# Protein-protein Interaction #

This repository will contain some of the code that I wrote during my
research term under Dr. Luis Rueda at the University of Windsor.

Areas include:
* PPI related data structures
* FASTA file format parsing
* PDB file parsing
* Finding protein complex interation interfaces
* Gibbs Sampling implementation for "one occurence per sequnece" (OOPS)

# Overall Layout #
* `include/`
  * Contains data structurs, algorithms and utility classes 
* `lib/`
  * Contains pre-compiled libraries such as the FASTA parser
* `gibbs_opps/`
  * Contains code implementing an OOPS model Gibbs sampler 
* `fasta_parser/`
  * Contains the source for the fasta parser
* `pdb_processing/`
  * Contains sample code for loading pdb files, data sets
  * Contains benchmark code for finding interface residues
