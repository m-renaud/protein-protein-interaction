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

#ifndef MRR_FASTA_PARSER_UTILS_HXX_
#define MRR_FASTA_PARSER_UTILS_HXX_

#include <fstream>
#include <string>
#include <vector>

#include <boost/spirit/include/support_istream_iterator.hpp>

#include <fasta_utils.hxx>
#include <fasta_parser.hxx>

namespace mrr {
namespace bio {


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool parse_fasta_file(
  std::string const& file_name,
  // std::vector<mrr::fasta_chain>& data
  mrr::bio::fasta_traits::start_state_return_type& data
)
{
  namespace sp = boost::spirit;

  std::ifstream is(file_name.c_str(), std::ios::binary);

  if(is.is_open())
  {
    is.unsetf(std::ios::skipws);
    sp::basic_istream_iterator<char> begin(is);
    sp::basic_istream_iterator<char> end;

    bool success = parse_fasta_file(begin, end, data);
    data.pdb_id = data.chains.begin()->pdb_id;
    return success;
  }
  else
  {
    return false;
  }

}

} // namespace bio
} // namespace mrr


#endif // #ifndef MRR_FASTA_PARSER_UTILS_HXX_
