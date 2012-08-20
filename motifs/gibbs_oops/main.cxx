#include <iostream>
#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <sstream>


#include <bio_utils.hxx>
#include <fasta_utils.hxx>
#include <pdb_utils.hxx>
#include <gibbs_sampler.hxx>


// argv[1] - dataset file
// argv[2] - dataset directory
int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    std::cerr << "Usage: " << argv[0]
	      << " (default | <dataset_file> <dataset_directory>)\n";
    exit(1);
  }

  // String representing the path to the dataset file and the directory
  std::string dataset_filename;
  std::string dataset_directory;
  std::string fasta_file_directory;
  std::string pdb_file_directory;

  int cycles = 100;

  if(std::string{argv[1]} == "default")
  {
    std::stringstream ss;
    ss << argv[2];
    ss >> cycles;

    dataset_filename = "../../data_sets/zhuListNonObligate.txt";
    dataset_directory = "../../data_sets/zhuNonObligate";
  }
  else
  {
    dataset_filename = argv[1];
    dataset_directory = argv[2];
  }

  fasta_file_directory = dataset_directory + "FASTAFiles/";
  pdb_file_directory = dataset_directory + "PDBFiles/";

  // String to hold the path to the protein complex currently being processed
  std::string protein_complex_path;

  // Holds a representation of the binary dataset file
  mrr::bio::binary_dataset dataset;

  // Holds the list of complexes based on the dataset file
  std::vector<mrr::bio::fasta_complex> fasta_complexes;
  std::vector<mrr::bio::pdb_complex> pdb_complexes;

  std::vector<std::string> all_sequences;
  std::set<char> amino_acid_set;
  std::vector<char> amino_acids;

  // Fill up the binary_dataset datastructure from the dataset file
  dataset = mrr::bio::process_dataset_file(dataset_filename);

  // For each complex in the dataset, parse the fasta file and put it in a vector
  for(mrr::bio::binary_complex const& bin_complex : dataset.complexes)
  {
    mrr::bio::fasta_complex fasta_protein_complex;

    mrr::bio::parse_fasta_file(
      fasta_file_directory + bin_complex.name + ".fasta",
      fasta_protein_complex
    );

    // Put the complex on the back of the vector
    fasta_complexes.push_back(fasta_protein_complex);

    std::cout << "Parse of " << bin_complex.name << " complete\n";

  } // for(binary_complex : dataset)


  // Display the parsed complexes and their sequences
  for(mrr::bio::fasta_complex const& complex: fasta_complexes)
  {
    std::cout << "Complex: " << complex.pdb_id << '\n';
    for(mrr::bio::fasta_complex::chain_type const& chain : complex.chains)
    {
      std::cout << "    Adding sequence of length: " << chain.sequence.size()
		<< std::endl;

      if((chain.sequence.size() >= 100) && (all_sequences.size() < 30))
      {
	all_sequences.push_back(chain.sequence.substr(0,100));
	for(auto aa : *std::prev(std::end(all_sequences)))
	  amino_acid_set.insert(std::toupper(aa));
      }
    }
  }

  std::cout << "\n\nNumber of sequences = " << all_sequences.size() << '\n';


  // Put the amino acids into a vector.
  std::move(
    begin(amino_acid_set), end(amino_acid_set),
    std::back_inserter(amino_acids)
  );

  // Create a Gibbs sampler and run it.
  mrr::bio::gibbs_sampler gibbs(7, amino_acids, all_sequences);

  gibbs.run_algorithm(cycles);


#if 0
  // Load the information from the pdb_files;
  mrr::bio::pdb_load_dataset(
    dataset_filename,
    pdb_file_directory,
    pdb_complexes
  );

  for(auto& complex : pdb_complexes)
  {
    std::cout << "Finding interface for complex " << complex.name << '\n';
    mrr::bio::find_interacting_residues(complex, 7);
  }
#endif

  return 0;
}
