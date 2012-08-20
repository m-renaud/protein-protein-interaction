#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <pdb_types.hxx>
#include <pdb_utils.hxx>

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//
int main(int argc, char* argv[])
{

  // Variables.
  mrr::bio::pdb_traits::complex_type complex;
  bool successful_parse;

  // Make sure a filename was entered.
  if(argc < 2)
  {
    std::cerr << "Enter a filename" << std::endl;
    return 1;
  }

  // Parse the pdb file
  successful_parse = mrr::bio::load_pdb_file(argv[1], complex);
  mrr::bio::compute_residue_centres(complex);
  mrr::bio::find_interacting_residues(complex, 7);

  // Optionally display parsed info.
  if(argc > 2 && std::string(argv[2]) == "-v")
    std::cout << std::endl << complex << std::endl;

  // Display the number of interface residues per chain
  std::cout << "Interface residues per chain: " << std::endl;
  for(char c : complex.chain_names)
    std::cout << "Chain " << c << ": " << complex.number_of_interface_residues(c)
              << std::endl
    ;

  endl(std::cout);

  std::cout << "Number of residues: "
            << complex.number_of_residues()
            << std::endl;

  std::cout << "Number of interface residues: "
            << complex.number_of_interface_residues()
            << std::endl;


  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  mrr::bio::print_interacting_residues(complex);

}
