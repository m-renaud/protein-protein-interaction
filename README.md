# Protein-protein Interaction #

This repository will contain some of the code that I wrote during my
research term under Dr. Luis Rueda at the University of Windsor.


### Areas include:###
* PPI related data structures
* FASTA file format parsing
* PDB file parsing
* Finding protein complex interation interfaces
* Gibbs Sampling implementation for "one occurence per sequnece" (OOPS)


## Overall Layout ##
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

## Requirements ##
* C++11 compiler
  * GCC version 4.7 or higher
  * Clange version 3.1 or higher

## Building Instructions ##
Before you start, run the `configure` file in the `motif/` directory
to build all necessary libraries

In each directory, there is a Makefile that can be used to build the
code in that directory. If your C++11 compiler is not the default, you
can set it with `CXX=your_c++11_compiler make`


## `include/` Directory ##
* Contains all header files used through the project
  * `bio_utils.hxx`contains routines for:
	* loading binary data set files
    * Removing unneeded chains
  * `coordinate.hxx` contains the coordinate class for representing positions
  * `fasta_parser.hxx`
	* Contains the grammar file used by Boost.Spirit for parsing FASTA files
  * `fasta_types.hxx` and `pdb_types.hxx`
	* Data structures for FASTA and PDB files
  * `fasta_utils.hxx` contains the function prototype for the FASTA
  parser
  * `gibbs_algorithm.hxx` contains the start of an abstracted Gibbs
    sampler, although not even close to completing
  * `gibbs_sampler.hxx`contains the class template for implementing a
    Gibbs sampling model
  * `pdb_utils.hxx` contains routines to:
	* Load PDB files
	* Computing residue centres
	* Finding interface residues
	* Loading a PDB dataset
	* Printing interacting residues
  * `utilities.hxx`
	* Contains several utility function used internally

## `lib/` Directory ##
* Contains pre-compiled fasta parser
  * Boost.Spirit takes a long time to compile, so only build it once
* When you run `configure`, this library file is build and placed in
  the correct directory
* If your default compiler does not support C++11, set it temporarily
  with `CXX=your_c++11_compiler ./configure`

## `gibbs_oops/` Directory ##
This directory contains code that implements an OOPS model for Gibbs
sampling.
* Build the code with `CXX=your_c++11_compiler make`
* Run it with `./main <dataset_file> <dataset_director>`

## `fasta_parser/` Directory ##
This directory contains the library file for the FASTA parser.
* The Gibbs OOPS implementation requires this file to be built
* This file is build automatically when `configure` is run

## `pdb_processing/` Dirctory ##
This directory contains:
* Code for loading PDB files
* Code for loading data sets
* Code for finding interface residues - `time_executions.cxx`
  * A bench mark is done comparing the Naive algorithm and the Residue
    Centre algorithm 
