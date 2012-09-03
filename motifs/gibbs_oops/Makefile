# Compilation command:
#   g++ -std=c++11 -I../include -L../lib main.cxx -lfasta_parser

# colors
C1=\e[1;39m
C2=\e[1;33m
C3=\e[1;32m
C4=\e[1;36m
C5=\e[1;34m
C6=\e[1;35m
C7=\e[1;31m
GR=\e[1;30m
BL=\e[1;8m
HI=\e[1;7m
UN=\e[1;0m
WA=\e[1;9m
UW=\e[1;11m
RST=\e[0m


LIB_DIR=../lib
INCLUDE_DIR=../include

LIBS=-lfasta_parser

CXX_OPTS=-std=c++11 -I$(INCLUDE_DIR) -L$(LIB_DIR) -Wall -Wextra

_DEPS=bio_utils.hxx fasta_utils.hxx pdb_utils.hxx gibbs_sampler.hxx
DEPS=$(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

########################################
all: options main

options:
	@echo -e "\n${C1}Build Options:\n------------------------------\n"
	@echo -e "${C3}CXX           ${C1}= ${C3}$(CXX)"
	@echo -e "${C3}CXX_OPTS      ${C1}= ${C5}$(CXX_OPTS)"
	@echo -e "${C3}LIB_DIR       ${C1}= ${C5}$(LIB_DIR)"
	@echo -e "${C3}INCLUDE_DIR   ${C1}= ${C5}$(INCLUDE_DIR)"
	@echo -e "${C3}LIBS          ${C1}= ${C5}$(LIBS)"
	@echo -e "${C3}DEPENDENCIES  ${C1}= ${C2}$(DEPS)"
	@echo -e "${RST}"

main: main.cxx $(DEPS)
	@echo -en "${RST}"
	$(CXX) $(CXX_OPTS) -o $@ $< $(LIBS)
#	@echo -e "${C3}$(CXX)" "${C5}-o" "${C6}$(@)${C1}" "${C2}$(<) $(DEPS)"
	@echo -e "${RST}"

colours:
	@echo -e "${C1}C1\n"
	@echo -e "${C2}C2\n"
	@echo -e "${C3}C3\n"
	@echo -e "${C4}C4\n"
	@echo -e "${C5}C5\n"
	@echo -e "${C6}C6\n"
	@echo -e "${C7}C7\n"
	@echo -e "${GR}GR\n"
	@echo -e "${BL}BL\n"
	@echo -e "${HI}HI\n"
	@echo -e "${UN}UN\n"
	@echo -e "${WA}WA\n"
	@echo -e "${UW}UW\n"
	@echo -e "${RST}"

.PHONY: clean

clean:
	rm -f *.o
	rm -f main
