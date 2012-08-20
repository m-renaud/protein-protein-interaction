#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <pdb_types.hxx>
#include <pdb_utils.hxx>


int main(int argc, char* argv[])
{
  std::vector<mrr::bio::pdb_traits::complex_type> complex_list;

  if(argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <dataset_file> <directory>\n";
    quick_exit(1);
  }

  mrr::bio::pdb_load_dataset(argv[1], argv[2], complex_list);

  return 0;
}
