OBJ = \
./source/needleman_wunsch.o \
./source/needleman_wunsch_affine_gaps_viterbi.o \
./source/reader.o \
./source/score_profile_similarity.o \
./source/score_profile_similarity_linear_normalized.o \
./source/score_sequence_similarity.o \
./source/score_sequence_similarity_profile_dependent.o \
./source/sequence.o \
./source/smith_waterman.o \
./source/string_functions.o \
./main/align_pairs.o

OBJ_PARA = \
./HMMER/src/generic_decoding.c \
./source/needleman_wunsch.cc \
./source/needleman_wunsch_affine_gaps_para.cc \
./source/reader.cc \
./source/score_profile_similarity.cc \
./source/score_profile_similarity_linear_normalized.cc \
./source/score_sequence_similarity.cc \
./source/score_sequence_similarity_profile_dependent.cc \
./source/sequence.cc \
./source/smith_waterman.cc \
./source/string_functions.cc \
./main/align_pairs.cpp

CXX = g++
CPP = g++-6

NAME = alignme2.0_mac.exe
NAME_PARA = alignme2.0_mac_parall.exe

CPPFLAGS = -v -Wall -O0 -g3 -fmessage-length=0 -Wno-deprecated -I./HMMER/include -L./HMMER/include -I./HMMER/src -L./HMMER/src -I./HMMER/easel -L./HMMER/easel -lhmmer -leasel -lm
PARA_FLAGS = -fopenmp -I/Users/edoardosarti/programs/openmp/openmp/runtime/build/runtime/src/ -L/Users/edoardosarti/programs/openmp/openmp/runtime/build/runtime/src/


single: $(OBJ)
	$(CXX) ${CPPFLAGS} $(OBJ) -o $(NAME)

parallel: $(OBJ)
	$(CPP) $(PARA_FLAGS) $(OBJ_PARA) -o $(NAME_PARA)

clean:
	rm -f $(NAME)
	rm -f $(OBJ)

