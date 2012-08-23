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

#ifndef MRR_FASTA_TYPES_HXX_
#define MRR_FASTA_TYPES_HXX_

#include <string>
#include <vector>
#include <iostream>
#include <boost/fusion/include/adapt_struct.hpp>

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

namespace mrr {
namespace bio {

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

struct fasta_chain
{
  using sequence_type = std::string;

  std::string pdb_id;
  std::string chain_name;
  sequence_type sequence;
};

struct fasta_complex
{
  using chain_type = mrr::bio::fasta_chain;
  using chain_list_type = std::vector<chain_type>;
  std::string pdb_id;
  chain_list_type chains;
};

struct fasta_traits
{
  using sequence_type = mrr::bio::fasta_chain::sequence_type;
  using chain_type = mrr::bio::fasta_chain;
  using chain_list_type = mrr::bio::fasta_complex::chain_list_type;
  using start_state_return_type = mrr::bio::fasta_complex;
};


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
template <typename Ch, typename Tr>
std::basic_ostream<Ch,Tr>& operator <<(
  std::basic_ostream<Ch,Tr>& os,
  mrr::bio::fasta_chain const& fc
)
{
  os << '(' << fc.pdb_id << ',' << fc.chain_name << ',' << fc.sequence << ')';
  return os;
}

template <typename Ch, typename Tr>
std::basic_ostream<Ch,Tr>& operator <<(
  std::basic_ostream<Ch,Tr>& os,
  std::vector<mrr::bio::fasta_chain> const& v
)
{
  os << '[';
  for (auto i = std::begin(v), iEnd = std::end(v); i != iEnd; ++i)
  {
    os << *i;

    // Output a command between elements...
    auto j = i;
    ++j;
    if (j != iEnd)
      os << ',';

  }
  os << ']';
  return os;
}


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

} // namespace bio
} // namespace mrr

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

BOOST_FUSION_ADAPT_STRUCT(
  mrr::bio::fasta_chain,
  (std::string, pdb_id)
  (std::string, chain_name)
  (std::string, sequence)
)

BOOST_FUSION_ADAPT_STRUCT(
  mrr::bio::fasta_complex,
  (mrr::bio::fasta_complex::chain_list_type, chains)
)

#endif // #ifndef MRR_FASTA_TYPES_HXX_
