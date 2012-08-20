INCLUDE_DIR=../include
LIB_DIR=../lib

CXX_OPTS=-std=c++11 -I$(INCLUDE_DIR)

DEPS=fasta_parser.hxx fasta_utils.hxx
DEPS=$(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

############################################################
all: fasta_parser_utils.o archive

fasta_parser_utils.o: fasta_parser_utils.cxx $(DEPS)
	$(CXX) $(CXX_OPTS) -c $<

archive:
	ar cr $(LIB_DIR)/libfasta_parser.a fasta_parser_utils.o

.PHONY: clean

clean:
	rm -vf *.o
