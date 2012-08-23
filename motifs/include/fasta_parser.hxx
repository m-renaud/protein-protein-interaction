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

#ifndef MRR_FASTA_PARSER_HXX_
#define MRR_FASTA_PARSER_HXX_

#define START_STATE_RULE_NAME   protein_complex;

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/repository/include/qi_flush_multi_pass.hpp>

#include <fasta_types.hxx>

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
namespace mrr {
namespace bio {


namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace spirit = boost::spirit;

//using return_type = std::vector<mrr::fasta_chain>;

template <typename Iterator>
struct white_space
  : qi::grammar<Iterator>
{
  white_space()
    : white_space::base_type(start_state)
  {
    using qi::lit;
    using qi::ascii::char_;
    using qi::ascii::space;

    start_state =
      space
      //| lit("/*") >> *(char_ - lit("*/")) >> lit("*/")
    ;
  }

  qi::rule<Iterator> start_state;
};

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
template <typename Iterator, typename Traits = fasta_traits>
struct fasta_grammar
  : qi::grammar<Iterator, typename Traits::start_state_return_type(), white_space<Iterator> >
{
  using start_state_return_type = typename Traits::start_state_return_type;
  using chain_list_type = typename Traits::chain_list_type;
  using chain_type = typename Traits::chain_type;
  using ws_ = white_space<Iterator>;

  fasta_grammar() : fasta_grammar::base_type(protein_complex)
  {
    using qi::lit;
    using qi::_val;
    using qi::double_;
    using qi::_1;
    using qi::space;
    using ascii::char_;

    pdb_id %=
      +(char_ - (space | char_(':')))
    ;

    chain_name %=
      +(char_ - (space | char_('|')))
    ;

    sequence %=
      +(char_ - char_('>'))
    ;

    chain %=
      lit('>') >> pdb_id
      >> lit(':') >> chain_name
      >> lit("|PDBID|CHAIN|SEQUENCE")
      >> sequence
    ;

    chains %= +chain;

    protein_complex = chains;

#if BOOST_SPIRIT_DEBUG
  #define SRN(s)  s.name(#s);
  #define SDR(n)  debug(n);
  #define S(n)    SRN(n); SDR(n);
    S(sequence);
    S(pdb_id);
    S(chain_name);
    S(chain);
    S(protein_complex);
#endif
  }

  qi::rule<Iterator, std::string()> pdb_id;
  qi::rule<Iterator, std::string()> chain_name;
  qi::rule<Iterator, std::string(), ws_> sequence;
  qi::rule<Iterator, chain_type(), ws_> chain;
  qi::rule<Iterator, chain_list_type(), ws_> chains;
  qi::rule<Iterator, start_state_return_type(), ws_> protein_complex;
};

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
template <typename Iter, typename Traits = fasta_traits>
inline bool parse_fasta_file(
  Iter& first, Iter& last,
  typename Traits::start_state_return_type& data
)
{
  return boost::spirit::qi::phrase_parse(
    first, last, fasta_grammar<Iter>(), white_space<Iter>(), data
  );
}

} // namespace bio
} // namespace mrr

#endif // #ifndef MRR_FASTA_PARSER_HXX_
