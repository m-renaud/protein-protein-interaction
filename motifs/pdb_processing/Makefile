CXX_OPTS=-std=c++11 -I../include

all: main time_executions process_dataset

main: main.cxx #pdb_types.hxx pdb_algs.hxx utility.hxx
	$(CXX) $(CXX_OPTS) $< -o $@

time_executions: time_executions.cxx #pdb_types.hxx pdb_algs.hxx utility.hxx
	$(CXX) $(CXX_OPTS) $< -o $@

process_dataset: process_dataset.cxx #pdb_types.hxx pdb_algs.hxx utility.hxx
	$(CXX) $(CXX_OPTS) $< -o $@

clean:
	rm -f main time_executions process_dataset
