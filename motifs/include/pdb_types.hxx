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

#ifndef MRR_PDB_TYPES_HXX_
#define MRR_PDB_TYPES_HXX_

#include <set>
#include <string>
#include <map>

#include <coordinate.hxx>

namespace mrr {
namespace bio {

// Prototypes
struct pdb_atom;
struct pdb_residue;
struct pdb_chain;
struct pdb_complex;


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
struct pdb_atom
{
  int number;
  int residue_number;
  std::string element;
  std::string compound;
  char chain_name;
  mrr::coordinate<double> position;
  double occupancy;
  double temperature_factor;
  std::string element_name;
};

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
struct pdb_residue
{
  using atom_list_type = std::vector<pdb_atom>;

  int number;
  std::string name;
  atom_list_type atoms;
  mrr::coordinate<double> centre;
  double max_distance;
  bool on_interface;
  std::vector<pdb_residue*> interface_residues;
};

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
struct pdb_chain
{
  using residue_list_type = std::vector<pdb_residue>;

  std::string name;
  residue_list_type residues;
};



//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
struct pdb_complex
{
  using chain_name_type = char;
  using chain_list_type = std::map<chain_name_type, pdb_chain>;

  std::string name;
  chain_list_type chains;
  std::set<chain_name_type> chain_names;

  //m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  std::size_t number_of_residues() const
  {
    std::size_t num_residues = 0;

    for(std::pair<const chain_name_type, mrr::bio::pdb_chain> const& p : chains)
      num_residues += p.second.residues.size();

    return num_residues;
  }

  //m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  std::size_t number_of_residues(chain_name_type const& name) const
  {
    return chains.at(name).residues.size();
  }

  //m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  int number_of_interface_residues() const
  {
    std::size_t interface_residues = 0;

    for(std::pair<const chain_name_type, mrr::bio::pdb_chain> const& p : chains)
      for(mrr::bio::pdb_residue const& res : p.second.residues)
        if(res.on_interface)
          ++interface_residues;

    return interface_residues;
  }

  //m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  int number_of_interface_residues(chain_name_type const& name) const
  {
    std::size_t interface_residues = 0;

    for(mrr::bio::pdb_residue const& res : chains.at(name).residues)
      if(res.on_interface)
        ++interface_residues;

    return interface_residues;
  }


};

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
struct pdb_traits
{
  using complex_type = pdb_complex;

  using chain_type = pdb_chain;
  using chain_list_type = typename pdb_complex::chain_list_type;

  using residue_type = pdb_residue;
  using residue_list_type = typename pdb_chain::residue_list_type;

  using atom_type = pdb_atom;
  using atom_list_type = typename pdb_residue::atom_list_type;
};


} // namespace bio
} // namespace mrr



//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Ostream overloads
std::ostream& operator <<(std::ostream& os, mrr::bio::pdb_atom const&);
std::ostream& operator <<(std::ostream& os, mrr::bio::pdb_residue const&);
std::ostream& operator <<(std::ostream& os, mrr::bio::pdb_chain const&);
std::ostream& operator <<(std::ostream& os, mrr::bio::pdb_complex const&);




//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
std::ostream& operator <<(std::ostream& os, mrr::bio::pdb_atom const& atom)
{
  os << "\t\t\tAtom #" << atom.number << std::endl
     << "\t\t\tElement: " << atom.element << std::endl
     << "\t\t\tCompound: " << atom.compound << std::endl
     << "\t\t\tChain: " << atom.chain_name << std::endl
     << "\t\t\tResidue #" << atom.residue_number<< std::endl
     << "\t\t\tCoondinates: " << atom.position;

  return os;
}

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
std::ostream& operator <<(std::ostream& os, mrr::bio::pdb_residue const& res)
{
  os << "\t\tResidue Name: " << res.name << std::endl
     << "\t\tResidue Num:  " << res.number << std::endl
     << "\t\tCentre: " << res.centre << std::endl
     << "\t\t.................................."
     << std::endl;

  for(mrr::bio::pdb_atom const& a : res.atoms)
    os << a << std::endl << std::endl;

  return os;
}

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
std::ostream& operator <<(std::ostream& os, mrr::bio::pdb_chain const& chain)
{
  os << "\tChain Name: " << chain.name << std::endl
     << "\t------------------------------------------"
     << std::endl;

  for(mrr::bio::pdb_residue const& r : chain.residues)
    os << r << std::endl;

  return os;
}

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
std::ostream& operator <<(std::ostream& os, mrr::bio::pdb_complex const& complex)
{
  os << "Complex Name: " << complex.name << std::endl
     << "=================================================="
     << std::endl;

  for(std::pair<char, mrr::bio::pdb_chain> const& c : complex.chains)
    os << c.second << std::endl;

  return os;
}

#endif // #ifndef MRR_PDB_TYPES_HXX_
