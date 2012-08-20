#ifndef MRR_FASTA_UTILS_HXX_
#define MRR_FASTA_UTILS_HXX_

#include <string>
#include <vector>

#include <fasta_types.hxx>

namespace mrr {
namespace bio {


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
bool parse_fasta_file(
  std::string const& file_name,
  // std::vector<mrr::fasta_chain>& data
  mrr::bio::fasta_traits::start_state_return_type& data
);


} // namespace bio
} // namespace mrr


#endif // #ifndef MRR_FASTA_UTILS_HXX_
