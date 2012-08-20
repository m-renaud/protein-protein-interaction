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
