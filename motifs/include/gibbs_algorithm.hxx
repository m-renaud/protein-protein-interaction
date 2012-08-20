#ifndef MRR_GIBBS_ALGORITHM_HXX_
#define MRR_GIBBS_ALGORITHM_HXX_


template <typename Model, typename ModelTraits = gibbs_traits<Model> >
struct gibbs_algorithm
{
  using amino_acid_type = typename ModelTraits::amino_acid_type;
  using amino_acid_list_type = typename ModelTraits::amino_acid_list_type;

  using sequence_type = typename ModelTraits::sequence_type;
  using sequence_list_type = typename ModelTraits::sequence_list_type;

  using profile_type = typename ModelTraits::profile_type;
  using sorted_profile_type = typename ModelTraits::sorted_profile_type; \


  /////////////////////////////////////////////////////////////////////////////
  // Variables ------------------------------------------------------------
  std::size_t k;  // Length of the k-mer
  std::size_t t;  // The number of sequences
  std::size_t l;  // The length of the sequences
  std::size_t h;  // The kmer to leave out when computing profile X


  sequence_list_type sequences;

  // The list of k-mers
  sequence_list_type sh_kmer_list;

  // Container for the k-mers used for computing X
  sequence_list_type sites;

  profile_type X; // Profile generated from selected k-mers
  profile_type Q; // Background frequencies

  amino_acid_list_type amino_acids;


  /////////////////////////////////////////////////////////////////////////////
  // Random number stuff --------------------------------------------------
  std::mt19937 seed_engine;
  std::uniform_int_distribution<std::size_t> seed_distribution;
  std::function<std::size_t()> seed_generator;

  std::mt19937 site_pos_engine;
  std::uniform_int_distribution<std::size_t> site_pos_distribution;
  std::function<std::size_t()> site_generator;

  std::mt19937 sequence_selector_engine;
  std::uniform_int_distribution<std::size_t> sequence_selector_distribution;
  std::function<std::size_t()> sequence_selector;


  //m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  gibbs_algorithm()
  {

  }


  //m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  void run_algorithm(long int cycles)
  {
    using std::begin;

  } // void run_algorithm()

}; // struct gibbs_algorithm

#endif // #ifndef MRR_GIBBS_ALGORITHM_HXX_
