//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Copyright 2012 Matthew R. Renaud
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This file is part of Protein-protein-interaction.
//
// Protein-protein-interaction is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or(at your option) any later version.
//
// Protein-protein-interaction is distributed in the hope that it will
// be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Protein-protein-interaction. If not, see
// <http://www.gnu.org/licenses/>.
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
