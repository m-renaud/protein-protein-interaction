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

#ifndef MRR_GIBBS_SAMPLER_HXX_
#define MRR_GIBBS_SAMPLER_HXX_

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <random>
#include <cmath>


namespace mrr {
namespace bio {

// Class to simulate a gibbs sampling algorithm.
struct gibbs_sampler
{
  /////////////////////////////////////////////////////////////////////////////
  // Type Aliases ---------------------------------------------------------
  using amino_acid_type = char;
  using amino_acid_list_type = std::vector<amino_acid_type>;

  using sequence_type = std::string;
  using sequence_list_type = std::vector<sequence_type>;

  using profile_type = std::pair<
      std::vector<std::map<amino_acid_type, std::size_t> >,
      std::vector<std::map<amino_acid_type, float> >
    >
  ;

  using sorted_profile_type =
    std::vector<std::vector<std::pair<amino_acid_type, double> > >
  ;


  /////////////////////////////////////////////////////////////////////////////
  // Variables ------------------------------------------------------------
  sequence_list_type sequences;

  std::size_t k;  // Length of the k-mer
  std::size_t t;  // The number of sequences
  std::size_t l;  // The length of the sequences
  std::size_t h;  // The kmer to leave out when computing profile X

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



  /////////////////////////////////////////////////////////////////////////////
  // Constructors --------------------------------------------------------
  //
  // Takes:
  //   - The size of the k-mer (std::size_t)
  //   - A list of the amino acids (vector<char>)
  //   - The list of sequences (vector<string>)
  gibbs_sampler(
    std::size_t const& k_,
    amino_acid_list_type const& amino_acids_,
    sequence_list_type const& sequences_
  )
    : sequences(sequences_), k(k_),
      t(sequences_.size()), l(begin(sequences_)->size()),

      amino_acids(amino_acids_),

      seed_engine(std::time(0)), seed_distribution(),
      seed_generator(std::bind(seed_distribution, seed_engine)),

      site_pos_distribution(0, l - 1 - k),
      sequence_selector_distribution(0, t - 1)
  {
    site_pos_engine.seed(seed_generator());
    sequence_selector_engine.seed(seed_generator());

    site_generator = std::bind(site_pos_distribution, site_pos_engine);
    sequence_selector = std::bind(
      sequence_selector_distribution, sequence_selector_engine
    );
  }



  ////////////////////////////////////////////////////////////////////////////
  // Function Prototypes -------------------------------------------------
  profile_type build_profile(sequence_list_type& kmers, int h);

  profile_type compute_background_frequencies(
    sequence_list_type const& sequences,
    std::size_t h
  );

  profile_type compute_background_frequencies_2(
    sequence_list_type const& sequences,
    std::size_t h
  );

  double probability_given_profile(
    sequence_type const& sequence,
    profile_type const& profile
  );

  void display_profile(profile_type const& profile);

  void display_sorted_profile(sorted_profile_type const& profile);

  double profile_score(profile_type const& profile);

  std::pair<bool,double> compute_site_match(profile_type const& profile, std::string const& site);

  auto find_sites(profile_type const& profile, double threshold)
    -> std::vector<std::pair<std::size_t, std::string::size_type> >;

  auto sort_profile(profile_type const& profile)
    -> std::vector<std::vector<std::pair<amino_acid_type, double> > >;

  void run_algorithm(long int cycles);

  double weight(sequence_type const& kmer);

  void sort_sh();

  int find_position_in_sh_kmer_list();

}; // struct gibbs_sampler


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Input:  A list of kmers, the position to exclude
// Output: A profile representing build from the list of kmers
auto gibbs_sampler::build_profile(sequence_list_type& kmers, int h = -1)
  -> profile_type
{
  std::size_t n = kmers.size(); // The number of sequences

  // If h is >= zero, a sequence has been removed
  if(h >= 0)
    --n;

  // Table to hold the frequency and probability of each AA.
  std::vector<std::map<amino_acid_type, std::size_t> > freq_table(k);
  std::vector<std::map<amino_acid_type, float> > prob_table(k);

  // Initialize the frequencies to zero.
  for(auto& x : freq_table)
    for(auto res: amino_acids)
      x[res] = 0;


  // Find the frequencies for each of the amino acids at each position...
  for(std::size_t pos = 0; pos < k; ++pos)
  {
    for(std::size_t seq = 0; seq < n; ++seq)
    {
      if(static_cast<int>(seq) == h)
        continue;
      ++freq_table[pos][kmers[seq][pos]];
    }
  }

  // Add pseudo counts to profile
  for(std::size_t pos = 0; pos < k; ++pos)
    for(auto& aa_freq_pair : freq_table[pos])
      ++aa_freq_pair.second;


  // Find the probability of a given residue at a position...
  for(std::size_t i = 0; i < k; ++i)
    for(auto& x: freq_table[i])
      prob_table[i][x.first] =
        static_cast<float>(x.second)/static_cast<float>(n + amino_acids.size());

  return std::make_pair(freq_table, prob_table);

} // build_profile()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function will compute the background frequencies from the list of
// sequences excluding the sequence h. It calculates this based on all kmers
// in the input sequences.
//
auto gibbs_sampler::compute_background_frequencies(
  sequence_list_type const& sequences,
  std::size_t h
) -> profile_type
{
  std::size_t n = 0;               // Number of k-mers

  // Table to hold the frequency and probability of each AA.
  std::vector<std::map<amino_acid_type, std::size_t> > freq_table(k);
  std::vector<std::map<amino_acid_type, float> > prob_table(k);


  // Initialive freqencies to zero.
  for(auto& x : freq_table)
    for(auto const& res: amino_acids)
      x[res] = 0;

  // Fill in the freq table
  for(std::size_t i = 0; i < t; ++i)
  {
    if(i == h)
      continue;

    for(std::size_t j = 0; j < l - k; ++j)
    {
      // Find every substring of length k
      std::string str = sequences[i].substr(j,k);
      ++n;

      // For each residue, update the frequency table
      for(std::size_t pos = 0; pos < k; ++pos)
        ++freq_table[pos][str[pos]];

    }
  }

  // Find the probability of a given residue at a position...
  for(std::size_t i = 0; i < k; ++i)
    for(auto& x: freq_table[i])
      prob_table[i][x.first] =
        static_cast<float>(x.second)/static_cast<float>(n);

  return std::make_pair(freq_table, prob_table);

} // compute_background_frequencies()



//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Find the probability that a kmer fits a given profile
//
auto gibbs_sampler::probability_given_profile(
  sequence_type const& sequence,
  profile_type const& profile
) -> double
{
  using std::begin;

  double result = 1;
  auto position_iter = begin(profile.second);

  for(amino_acid_type const& residue : sequence)
     result *= position_iter++->at(residue);

  return result;

} // probability_given_profile()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Self explanatory, displays a given profile
//
auto gibbs_sampler::display_profile(profile_type const& profile) -> void
{
  // Display the heading to the tables.
  std::cout << '\n' << std::left << "  Position  |  ";
  for(auto& x: amino_acids)
  {
    std::cout << std::setw(8) << x;
  }

  std::cout << '\n'
            << "  --------------------------------------------"
            << std::endl;

  std::size_t count = 1;

  // Display the freq_table...
  for(auto& m: profile.first)
  {
    std::cout << "    " << std::setw(8) << count++ << "|  ";
    for(auto& x: m)
      std::cout << std::setw(8) << x.second;
    endl(std::cout);
  }

  std::cout << std::fixed;

  // Display the probabilities...
  std::cout << "\nProbabilities:\n";

  // Reset count
  count = 1;

  std::cout << '\n' << std::left << "  Position  |  ";

  for(auto& x: amino_acids)
    std::cout << std::setw(8) << x;

  std::cout << "\n"
            << "  --------------------------------------------"
            << std::endl;

  for(auto& m: profile.second)
  {
    std::cout << "    " << std::setw(8) << count++ << "|  ";
    for(auto& x: m)
      std::cout << std::setw(8) << std::setprecision(4) << x.second;
    endl(std::cout);
  }

} // void display_profile()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Display the sorted profile nicely formatted
//
void gibbs_sampler::display_sorted_profile(sorted_profile_type const& sorted_profile)
{
  for(auto const& position : sorted_profile)
  {
    for(auto const& p : position)
    {
      std::cout << p.first << ": " << std::setw(7) << std::setprecision(4)
                << p.second << " ";
    }
    endl(std::cout);
  }

} // display_sorted_profile()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Compute the score of a profile
//
double gibbs_sampler::profile_score(profile_type const& profile)
{
  double score ;

  score = k * std::log(amino_acids.size());

  // NOTE:
  //  - position_list is a vector<map<residue,prob> >
  //  - probability is a pair<amino_acid, probability>

  for(auto const& position_list : profile.second)
    for(auto const& probability : position_list)
    {
      if(probability.second == 0)
        continue;
      else
        score += probability.second * std::log(probability.second);
    }

  score /= k;

  return score;
}


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Given a possible site, compute the value of the match based on a profile
//
std::pair<bool,double> gibbs_sampler::compute_site_match(
  profile_type const& profile, std::string const& site
)
{
  double match = 1;
  std::string::size_type pos;

  // For each position in the site, multiply the probability of finding
  // that residue by the current value

  // NOTE:  profile.second[pos].find(site[pos])->second;
  //        profile.second[pos] - returns a map from residue to prob
  //        find(site[pos])     - returns a pair of residue to prob
  //        ->second            - returns the probability
  for(pos = 0; pos < site.size(); ++pos)
  {
    if(profile.second[pos].find(site[pos])->second < 0.05)
      return std::make_pair(false, 0);
    match *= profile.second[pos].find(site[pos])->second;
  }

  return std::make_pair(true,match);

} // compute_site_match()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Returns a list of positions into the sequence list where sites are
// located based on the threshold provided
//
auto gibbs_sampler::find_sites(profile_type const& profile, double threshold)
  -> std::vector<std::pair<std::size_t,std::string::size_type> >
{
  using position_pair_type = std::pair<std::size_t, std::string::size_type>;
  std::vector<position_pair_type> position_list;
  std::pair<bool,double> match;

  std::size_t seq_list_size = sequences.size();

  for(std::size_t seq_num = 0; seq_num < seq_list_size; ++seq_num)
  {
    std::size_t seq_size = sequences[seq_num].size();

    for(std::string::size_type pos = 0; pos < seq_size - k; ++pos)
    {
      match = compute_site_match(profile, sequences[seq_num].substr(pos, k));
      if(match.first && match.second >= threshold)
      {
        position_list.push_back(std::make_pair(seq_num, pos));

        // Optionally skip to the end of the kmer to avoid overlap.
        // Add k-1 and the for loop will take care of the last one
        pos += k-1;
      }

    }
  }

  return position_list;

} // find_sites()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Sort the residues in a profile from highest to lowest probability
//
auto gibbs_sampler::sort_profile(profile_type const& profile)
  -> std::vector<std::vector<std::pair<amino_acid_type, double> > >
{
  std::vector<std::pair<amino_acid_type, double> > sorted_position;
  std::vector<std::vector<std::pair<amino_acid_type, double> > > sorted_profile;

  // For each position in the probability table, put the pairs into a list
  // then sort them
  for(auto const& position : profile.second)
  {
    sorted_position.clear();
    sorted_position.assign(begin(position), end(position));

    std::sort(
      begin(sorted_position), end(sorted_position),
      [](std::pair<amino_acid_type, double> a, std::pair<amino_acid_type, double> b)
      {
        return a.second > b.second;
      }
    );

    sorted_profile.push_back(sorted_position);

  } // for()

  return sorted_profile;

} // sort_profile();



//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Driver function for the Gibbs Sampling algorithm
void gibbs_sampler::run_algorithm(long int cycles)
{
  using std::begin;

  // Empty the list of sites
  sites.clear();

  // 1) Randomly select an k-mer a_i in each input sequence
  for(std::size_t i = 0; i < t; ++i)
    sites.push_back(sequences[i].substr(site_generator(), k));

  for(int loop_count = 1; loop_count <= cycles; ++loop_count)
  {

    // 2) Randomly select one input sequence s_h
    h = sequence_selector();

    // 3) Build a 20xl profile X from a_1,...,a_h-1,a_h+1,...,a_t

    // Build a profile from the sites, excluding the one at position h
    X = build_profile(sites, (loop_count == 1) ? h : -1 );

    // display_profile(X);

    // 4) Compute background frequencies Q from input
    //    sequences s1,...,s_h-1,s_h+1,...,s_t
    Q = compute_background_frequencies(sequences, h);

    // display_profile(Q);

    // 5/6) - Compute the weight for each kmer in s_h
    sort_sh(); // Sort the list of kmers in s_h

    int pos = find_position_in_sh_kmer_list();

    sites[h] = sh_kmer_list[pos];

    // Every 100 iterations, display the current profile and sites
    if(loop_count % 100 == 0)
    {
      display_profile(X);
      auto sorted_profile = sort_profile(X);
      display_sorted_profile(sorted_profile);
      endl(std::cout);
      std::cout << "Score = " << profile_score(X) << std::endl;
      auto position_list = find_sites(X, 0.00000000078125);

      for(std::pair<std::size_t, std::string::size_type> const& p: position_list)
      {
        std::cout << std::setw(5) << p.first  << "   "
                  << std::setw(5) << p.second << "   "
                  << sequences[p.first].substr(p.second, k) << "   "
                  << std::setprecision(100)
                  << compute_site_match(
                    X, sequences[p.first].substr(p.second, k)
                  ).second
                  << std::endl;
      }

    } // if(loop_count % 100 == 0)

  } // for(loop_count)

} // run_algorithm()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This functions computes the weight of a given kmer based on the current
// profile and the background frequencies.
//
double gibbs_sampler::weight(sequence_type const& kmer)
{
  double prob_given_X, prob_given_Q;

  prob_given_X = probability_given_profile(kmer, X);
  prob_given_Q = probability_given_profile(kmer, Q);

  return prob_given_X/prob_given_Q;

} // weight()


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function sorts the kmers in sequence 'h' which can then be used by
// find_position_in_sh_kmer_list()
//
void gibbs_sampler::sort_sh()
{
  double sum_of_weight = 0;

  // Empty the previous list of kmers
  sh_kmer_list.clear();

  for(std::size_t i = 0; i < l - k; ++i)
  {
    sh_kmer_list.push_back(sequences[h].substr(i,k));
    sum_of_weight += weight(sequences[h].substr(i,k));
  }

  // Sort the list of kmers in sequence h based on their weight
  std::sort(
    begin(sh_kmer_list), end(sh_kmer_list),
    [&](std::string const& s1, std::string const& s2)
    {
      return this->weight(s1) < this->weight(s2);
    }
  );
}


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This function selects a position in sequence 'h' to be used in the kmer list
// for computing the profile
//
int gibbs_sampler::find_position_in_sh_kmer_list()
{
  double continuous_prob_sum;

  continuous_prob_sum = std::accumulate(
    begin(sh_kmer_list), end(sh_kmer_list),
    0,
    [&](double d,  sequence_type const& seq)
    {
      double weight = this->weight(seq);
      return d + weight;
    }
  );

  std::mt19937 prob_engine(seed_generator());
  std::uniform_real_distribution<double>  prob_distribution(0, continuous_prob_sum);
  std::function<double()> prob_generator {
    std::bind(
      prob_distribution, prob_engine
    )};

  double prob = prob_generator();

  double min = 0;

  for(std::size_t i = 0; i < sh_kmer_list.size(); ++i)
  {
    if(prob >= min && prob < min + weight(sh_kmer_list[i]))
      return i;
    min += weight(sh_kmer_list[i]);
  }

  return sh_kmer_list.size() - 1;

} // find_position_in_sh_kmer_list()


} // namespace bio
} // namespace mrr


#endif // #ifndef MRR_GIBBS_SAMPLER_HXX_
