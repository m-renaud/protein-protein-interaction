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
#include <iomanip>
#include <chrono>
#include <fstream>
#include <vector>
#include <sstream>

#include <pdb_types.hxx>
#include <pdb_utils.hxx>
#include <utilities.hxx>

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
//
int main(int argc, char* argv[])
{
  using watch_type = mrr::util::stopwatch<std::chrono::high_resolution_clock>;

  // Variables.
  watch_type watch;
  mrr::bio::pdb_traits::complex_type complex, complex2;

  double exe_time;
  double naive_exe_time;

  // Make sure a filename was entered.
  if(argc < 2)
  {
    std::cerr << "Enter a filename" << std::endl;
    return 1;
  }

  mrr::bio::load_pdb_file(argv[1], complex);
  complex2 = complex;

  // Time the algorithm
  watch.start();
  mrr::bio::compute_residue_centres(complex);
  mrr::bio::find_interacting_residues(complex, 7);
  exe_time = watch.lap();

  // Time the naive algorithm
  mrr::bio::find_interacting_residues_naive(complex2, 7);
  naive_exe_time = watch.lap();

#if 1
  // Optionally display parsed info.
  if(argc > 2 && std::string(argv[2]) == "-v")
    std::cout << std::endl << complex << std::endl;

  // Display execution time
  std::cout << "Naive Execution Time: \t" << naive_exe_time << " seconds\n"
            << "Central Point Execution Time: \t" << exe_time << " seconds\n"
            << "Time Ratio: " << (naive_exe_time / exe_time) << std::endl
  ;

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
            << '\t'
            << complex2.number_of_interface_residues()
            << std::endl;
#endif

  std::cout << std::setw(5) << complex.name << " & "
            << std::setw(10) << naive_exe_time << " & "
            << std::setw(10) << exe_time << " & "
            << std::setw(10) << (naive_exe_time / exe_time) << " \\\\\n";
}
