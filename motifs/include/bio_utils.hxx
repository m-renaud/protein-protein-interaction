#ifndef BIO_UTILS_HXX_
#define BIO_UTILS_HXX_

#include <algorithm>
#include <array>
#include <vector>
#include <fstream>
#include <iostream>

namespace mrr {
namespace bio {

// Structure to repersent a binary interaction complex 
struct binary_complex
{
  std::string name;
  std::array<char,2> chains;
};

// Structure to represent a binary dataset file
struct binary_dataset
{
  std::vector<binary_complex> complexes;
};


// This function takes a filename, namely a dataset filename, and converts
// it to a binary_dataset structure to be used in other algorithms.
binary_dataset process_dataset_file(std::string const& filename)
{
  std::ifstream is(filename.c_str(), std::ios::binary);

  std::string line;
  std::string complex_name;
  std::array<char,2> chains;
  
  binary_dataset dataset;
  

  if(is.is_open())
  {
    while(std::getline(is, line))
    {
      complex_name = line.substr(0,4);
      chains[0] = line[5];
      chains[1] = line[7];
      dataset.complexes.push_back(binary_complex{complex_name, chains});
    }
  }

  return dataset;
}


// This function takes a fasta_complex and a list of chain names and removes
// all chains from the fasta complex that aren't in the chain_list.
template <
  typename FastaType,
  typename SequenceType,
  typename ChainListType
>
void remove_unneeded_chains(FastaType& complex, ChainListType const& chain_list)
{
  complex.chains.erase(
    std::remove_if(
      begin(complex.chains), end(complex.chains),
      [&](typename FastaType::chain_type const& chain)
      {
	return std::find(
	  begin(chain_list), end(chain_list),
	  chain.chain_name
	) != end(chain_list);
      }
    ), end(complex.chains)
  );
}


} // namespace bio
} // namespace mrr


#endif // #ifndef BIO_UTILS_HXX_
