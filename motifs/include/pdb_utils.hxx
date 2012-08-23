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

#ifndef MRR_PDB_UTILS_HXX_
#define MRR_PDB_UTILS_HXX_

#include <thread>
#include <future>

#include <pdb_types.hxx>
#include <utilities.hxx>

namespace mrr {
namespace bio {

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Prototypes

template <typename Traits = pdb_traits>
bool load_pdb_file(
  std::string const& filename,
  typename Traits::complex_type& complex
);

template <typename Traits = pdb_traits>
void compute_residue_centres(typename Traits::complex_type& complex);

template <typename Traits = pdb_traits>
void find_interacting_residues(
  typename Traits::complex_type& complex,
  double threshold
);

template <typename Traits = pdb_traits>
void find_interacting_residues_naive(
  typename Traits::complex_type& complex,
  double threshold
);

template <typename Traits = pdb_traits>
void pdb_load_dataset(
  std::string const& dataset_filename,
  std::string const& dir,
  std::vector<typename Traits::complex_type>& complex_list
);


template <typename Traits = pdb_traits>
void print_interacting_residues(
  typename Traits::complex_type const& complex,
  std::ostream& os = std::cout
);


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a pdb file and a pdb_complex data structure and
// fills the data structure with the info from the file.
template <typename Traits = pdb_traits>
bool load_pdb_file(
  std::string const& filename,
  typename Traits::complex_type& complex
)
{
  using atom_type = typename Traits::atom_type;
  using residue_type = typename Traits::residue_type;
  using chain_type = typename Traits::chain_type;

  // For file and line operations...
  std::ifstream pdb_file(filename.c_str(), std::ios::binary);
  std::string line;
  std::stringstream ss;

  // Keep track of the current residue and chain being processed
  int current_residue_number;
  std::string current_compound; // ie residue
  char current_chain;

  // Flag so that we don't push back an empty vector of atoms the first time
  bool first_atom = true;

  // Temporary atom, residue and chain for processing
  atom_type atom;
  residue_type residue;
  chain_type chain;


  //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // If the file is successfully opened
  if(pdb_file.is_open())
  {
    // Extract the header
    std::vector<std::string> fields;
    std::getline(pdb_file, line);
    ss.str(line);

    fields.assign(
      std::istream_iterator<std::string>(ss),
      std::istream_iterator<std::string>()
    );

    complex.name = *--std::end(fields);

    // Clear the flags
    ss.clear();

    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // Loop through the file line by line
    while(std::getline(pdb_file, line))
    {
      // We only care about the lines that start with ``ATOM"
      if(mrr::util::has_prefix(line, "ATOM"))
      {
        // Put the line into a stringstream and extract the info into an atom
        ss.str(line);

        ss >> mrr::util::ignore<std::string>
           >> atom.number >> atom.element >> atom.compound // compound ie. residue name
           >> atom.chain_name >> atom.residue_number
           >> atom.position.x >> atom.position.y >> atom.position.z
           >> atom.occupancy >> atom.temperature_factor >> atom.element_name
        ;


        // If we are not on the first atom in a chain and we are
        // starting a new residue, add the old one to the chain and clear
        if(!first_atom && (current_residue_number != atom.residue_number))
        {
          residue.number = current_residue_number;
          residue.name = current_compound;
          residue.on_interface = false;
          chain.residues.push_back(std::move(residue));
          residue.atoms.clear();
        }

        // Set the current chain and residue values
        current_chain = atom.chain_name;
        current_residue_number = atom.residue_number;
        current_compound = atom.compound;

        // Put the atom at the back of the residue vector
        residue.atoms.push_back(std::move(atom));

        // We are no longer on the first atom if we've gotten to here
        first_atom = false;

      } // if(has_prefix "ATOM")

      // If it has the prefix ``TER", then we are at the end of the chain
      if(mrr::util::has_prefix(line, "TER"))
      {
        // std::cout << "End of chain " << current_chain << std::endl;
        // std::cout << "Number of residues: " << chain.residues.size() << std::endl;

        residue.number = current_residue_number;
        residue.name = current_compound;
        residue.on_interface = false;
        chain.residues.push_back(std::move(residue));
        residue.atoms.clear();

        chain.name = current_chain;
        complex.chains.insert(std::make_pair(current_chain, std::move(chain)));
        chain.residues.clear();

        first_atom = true;

      } // if(has_prefix "TER")

    } // while(std::getline(pdb_file, line))

    // Make a set of the chain names
    for(std::pair<const char, mrr::bio::pdb_chain> const& p : complex.chains)
      complex.chain_names.insert(p.first);

    return true;

  } // if(pdb_file is open)
  else
  {
    std::clog << "Cannot open pdb file: " << filename << std::endl;
    return false;
  }

} // bool load_pdb_file()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function computes the centre of the residues in the complex
// as well as setting the maximum distance from any atom in the
// residue to the centre.
template <typename Traits = pdb_traits>
void compute_residue_centres(typename Traits::complex_type& complex)
{
  mrr::bio::pdb_traits::atom_list_type::size_type num_atoms;
  double dist;
  double max_dist = 0;

  // Loop through the chains in the complex
  for(std::pair<const char,mrr::bio::pdb_chain>& chain: complex.chains)
  {
    // Loop through the residues
    for(mrr::bio::pdb_residue& residue: chain.second.residues)
    {
      // Keep adding the position of the current atom to the centre
      for(mrr::bio::pdb_atom& atom: residue.atoms)
      {
        residue.centre.x += atom.position.x;
        residue.centre.y += atom.position.y;
        residue.centre.z += atom.position.z;
      }

      // Find the number of atoms
      num_atoms = residue.atoms.size();

      // Divide the coordinates by the number of atoms
      residue.centre.x /= num_atoms;
      residue.centre.y /= num_atoms;
      residue.centre.z /= num_atoms;

      // Loop through the atoms and find the maximum distance from the centre
      // of the residue
      for(mrr::bio::pdb_atom const& atom: residue.atoms)
      {
        dist = mrr::distance(atom.position, residue.centre);
        max_dist = std::max(dist, max_dist);
      }

      // Set the max distance and reset it to zero for the next residue
      residue.max_distance = max_dist;
      max_dist = 0;

    } // for(residue : chain.second.residues)

  } // for(chain: complex.chains)

} // void compute_residue_centres()

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
template <typename Traits = pdb_traits>
void find_interacting_residues(
  typename Traits::complex_type& complex,
  double threshold
)
{
  // Loop through the chains
  for(std::pair<const char, mrr::bio::pdb_chain>& chain: complex.chains)
  {
    // Loop through the residues of the current chain
    for(mrr::bio::pdb_residue& residue: chain.second.residues)
    {
      // Loop through the residues of all the other chains
      for(std::pair<const char, mrr::bio::pdb_chain>& other_chain: complex.chains)
      {
        // If on a different chain
        if(other_chain.second.name != chain.second.name)
        {
          // Loop through the residues, compare their distance to the current residue
          for(mrr::bio::pdb_residue& other_res: other_chain.second.residues)
          {
            // If the distance is less than 7A + the distances from centre,
            // check the distances of the atoms themselves
            auto dist = mrr::distance(residue.centre, other_res.centre);

            if(dist < threshold + residue.max_distance + other_res.max_distance)
            {
              // Loop through the atoms of the current residue
              for(mrr::bio::pdb_atom const& a1 : residue.atoms)
              {
                // Loop through the atoms of the other residues
                for(mrr::bio::pdb_atom const& a2 : other_res.atoms)
                {
                  // Compare their distances
                  if(mrr::distance(a1.position, a2.position) < threshold)
                  {
                    residue.on_interface = true;
                    residue.interface_residues.push_back(&other_res);

                    // Get out of the loops, we have already found that the
                    // residue is on the interface.  This is the only valid
                    // use of a goto statement...and even this bothers me :P
                    goto end_of_interface_atom_check;  //===\/
                                                       //   \/
                  } // if(distance < threshold)        //   \/
                } // for(a2 : other_res.atoms)         //   \/
              } // for(a1 : residue.atoms)             //   \/
                                                       //   \/
                                                       //   \/
            end_of_interface_atom_check:  // <==============<<
              ;


            } // if(dist < threshold + max_distances)

          } // for(residue: other_chain.second.residues)

        } // if(different chain)

      } // for(pair other_chain : complex.chains)

    } // for(residue : chain.second)

  } // for(pair chain : complex.chains)

} // void find_interacting_residues()


#if 0
void helper(
  mrr::bio::pdb_complex& complex,
  double threshold,
  std::pair<const char, mrr::bio::pdb_chain>& chain,
  mrr::bio::pdb_residue& residue
)
{
  // Loop through the residues of all the other chains
  for(std::pair<const char, mrr::bio::pdb_chain>& other_chain: complex.chains)
  {
    // If on a different chain
    if(other_chain.second.name != chain.second.name)
    {
      // Loop through the residues and compare their distance to the current residue
      for(mrr::bio::pdb_residue& other_res: other_chain.second.residues)
      {
        // If the distance is less than 7A, display the interacting residues
        auto dist = mrr::distance(residue.centre, other_res.centre);

        if(dist < threshold + residue.max_distance + other_res.max_distance)
        {

          for(mrr::bio::pdb_atom const& a1 : residue.atoms)
          {
            for(mrr::bio::pdb_atom const& a2 : other_res.atoms)
            {
              if(mrr::distance(a1.position, a2.position) < threshold)
              {
                residue.on_interface = true;
                residue.interface_residues.push_back(&other_res);
                goto end_of_interface_atom_check;  //===\/
              }                                    //   \/
            }                                      //   \/
          }                                        //   \/
                                                   //   \/
                                                   //   \/
        end_of_interface_atom_check:  // <==============<<
          ;

        } // if(dist < threshold + max_distances)

      } // for(residue: other_chain.second.residues)

    } // if(different chain)

  } // for(pair other_chain : complex.chains)

} // helper()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void find_interacting_residues_async(
  mrr::bio::pdb_complex& complex,
  double threshold
)
{
  std::vector<std::future<bool> > fin;
  // Loop through the chains
  for(std::pair<const char, mrr::bio::pdb_chain>& chain: complex.chains)
  {
    // Loop through the residues of the current chain
    for(mrr::bio::pdb_residue& residue: chain.second.residues)
    {
      fin.push_back(
        std::async([&]{ helper(complex, threshold, chain, residue); return true;})
      );

    } // for(residue : chain.second)

  } // for(pair chain : complex.chains)

  for(auto& x: fin)
    x.get();

} // void find_interacting_residues_async()
#endif

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This is a naive implementation of an algorithm to find the residues that are
// on the interface.
template <typename Traits = pdb_traits>
void find_interacting_residues_naive(
  typename Traits::complex_type& complex,
  double threshold
)
{
  // Loop through the chains
  for(std::pair<const char, mrr::bio::pdb_chain>& chain: complex.chains)
  {
    // Loop through the residues of the current chain
    for(mrr::bio::pdb_residue& residue: chain.second.residues)
    {
      for(mrr::bio::pdb_atom& atom: residue.atoms)
      {
        // Loop through the residues of all the other chains
        for(std::pair<const char, mrr::bio::pdb_chain>& other_chain: complex.chains)
        {
          // If on a different chain
          if(other_chain.second.name != chain.second.name)
          {
            // Loop through the residues in that chain
            for(mrr::bio::pdb_residue& other_res: other_chain.second.residues)
            {
              // Loop through the atoms from the other residues
              for(mrr::bio::pdb_atom const& atom_2 : other_res.atoms)
              {
                // If the distance is within the threshold, they are interacting
                if(mrr::distance(atom.position, atom_2.position) < threshold)
                {
                  residue.on_interface = true;
                  residue.interface_residues.push_back(&other_res);
                }

              } // for(atom_2 : other_res.atoms)

            } // for(residue: other_chain.second.residues)

          } // if(different chain)

        } // for(pair other_chain : complex.chains)

      } // for(atom: residue.atoms)

    } // for(residue : chain.second)

  } // for(pair chain : complex.chains)

} // void find_interacting_residues_naive()



//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function loads a dataset into a list of protein complexes.
template <typename Traits>
void pdb_load_dataset(
  std::string const& dataset_filename,
  std::string const& dir,
  std::vector<typename Traits::complex_type>& complex_list
)
{
  using complex_type = typename Traits::complex_type;

  // Open the dataset file for reading
  std::ifstream dataset_file(dataset_filename, std::ios::binary);

  std::string line;
  std::string pdb_id;
  std::string path;

  // If the file is open
  if(dataset_file.is_open())
  {
    // Process each line in the dataset file
    while(std::getline(dataset_file, line))
    {
      // Temporary complex to parse the file into
      complex_type complex;

      // Get the path to the file based on the pdb id and the directory provided
      pdb_id = line.substr(0,4);
      path = dir + pdb_id + ".pdb";

      std::cout << "Processing " << path << std::endl;

      // Load the pdb file, compute the centres and find the interacting residues
      mrr::bio::load_pdb_file(path, complex);
      mrr::bio::compute_residue_centres(complex);
      mrr::bio::find_interacting_residues(complex, 7);

      // Push the complex onto the vector
      complex_list.push_back(complex);
    }
    std::cout << "Done processing dataset" << std::endl;

  }

} // void load_dataset()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function takes a protein complex and an ostream and print all of the
// residues on the interface, consecutive ones being grouped together.
template <typename Traits = pdb_traits>
void print_interacting_residues(
  typename Traits::complex_type const& complex,
  std::ostream& os = std::cout
)
{
  bool prev_on_interface;

  // Loop through the chains in the complex
  for(std::pair<const char, mrr::bio::pdb_chain> const& chain : complex.chains)
  {
    prev_on_interface = true;

    // Print the chain name
    os << "==================================================\n"
       << "Chain " << chain.first << ":\n";

    // Loop through the residues
    for(mrr::bio::pdb_residue const& res : chain.second.residues)
    {
      // If the residue is on the interface, print it's name
      if(res.on_interface)
      {
        os << res.name << std::endl;
        prev_on_interface = true;
      }
      else if(!res.on_interface && prev_on_interface)
      {
        prev_on_interface = false;
        os << "\n";
      }

    } // for(res : chain.second.residues)

  } // for(pair : complex.chains)

} // void print_interacting_residues()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

} // namespace bio
} // namespace mrr

#endif // #ifndef MRR_PDB_UTILS_HXX_
