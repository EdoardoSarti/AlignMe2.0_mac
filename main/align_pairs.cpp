/////////////////////////////////////////////////////////////////////////
//!
//  Copyright (C) 2010 by Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//  AlignMe@rzg.mpg.de
//
//  AlignMe is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  AlignMe is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//!
//!
//!  AlignMe - a program to identify and align similar proteins with or
//! without sequence similarity. It allows to use a flexible number of
//! scales or profiles for identifying similarity.
//!
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @version: 1.1
//! @date: 16.12.2011
/////////////////////////////////////////////////////////////////////////


// CHANGELOG
// Version 1.1
// - added option to generate alignments based on PSSM
// - weights can now be adjusted automatically (<< INCLUDED???????????????????)
// - added option to set the format of the alignment output to clustalW or fasta-format
// - gaps in profiles were assigned 0 as gap-value - now have the gap-value ?0
// - sequence alignment of MSA vs MSA can now be printed out (<<< CHECK!!!)
// - fixed a bug concerning incompatibilities with Ubuntu
// - added some error messages that were missing
// - improved documentation
//
// Version 1.0
//	- first release


#include <cmath>
#include <set>


#include "../include/definitions.h"
#include "../include/macro_functions_read_write.h"


#ifdef POSITION_SPECIFIC_SIMILARITY
#include "../include/score_position_specific_similarity.h"
#endif

#include "../include/amino_acid.h"
#include "../include/sum_function.t.h"
#include "../include/needleman_wunsch.h"
#include "../include/needleman_wunsch_affine_gaps.h"
#include "../include/smith_waterman.h"
#include "../include/command_line_manager.h"
#include "../include/string_functions.h"
#include "../include/score_sequence_similarity.h"
#include "../include/score_profile_similarity.h"
#include "../include/score_sequence_similarity_profile_dependent.h"
#include "../include/std_functions.h"
#include "../include/profile_step_function.h"
#include "../include/constant_function.h"



#include "../include/alignment_variables.h"
#include "../include/allowed_flags.h"
#include "../include/alignment_functions.h"
#include "../include/write_help.h"
#include "../include/msa_averager.h"
#include "../include/alignment_write.h"
#include "../include/reader.h"
#include "../include/commands_defined.h"


// NOTES:
// - fasta header for MSAs should contain only single strings




// reading in the commandline. argc = number of parameters of the commandline, argv = string-array with the commandlineparameters
int main( const int argc, const char * argv[])
{
  ///////////////////////////////////
  ///  DECLARATIONS
  ///////////////////////////////////
  DebugWrite( argv[0]);
  DebugWrite( "Welcome to the Debugger of AlignMe1.1");



  std::map<std::string,std::string>
	  allowed_flags;

  DefineAllowedFlagsForAlignPairs( allowed_flags);



  CommandLineManager
	  cmd( argc, argv, allowed_flags);

  AlignmentVariables
	  vars;


  SetDefaults( vars);


  if( cmd.IsFlagSet( "help") || cmd.IsFlagSet( "-help"))
  {
	  WriteHelpForAlignPairs();
  }


  ReadCommandsForAlignPairs( cmd, vars);



  ///////////////////////////////////////////
  /// READ VALUES FROM COMMANDLINE OPTIONS
  ///////////////////////////////////////////


  ReadSequencesProfilesAndScalesForAlignPairs( vars);



  /////////////////////////////////////////////
  ////  BUILD MATRIX
  ///////////////////////////////////////////

  vars.dynamic_programing_matrix.Set( vars.first_sequence.size() + 1, vars.second_sequence.size() + 1);


  /////////////////////////////////////////////
  ////  PERFORM ALIGNMENT
  ///////////////////////////////////////////

  
  // 5.4.17 Ok until here

  std::ofstream
	  out;


  if( vars.algorithm == "global_linear")
  {

#ifdef SECURE

	  if( cmd.IsFlagSet( "thresholds_for_penalties"))
	  {
		  std::cerr << "ERROR: the 'thresholds_for_penalties' flag has no meaning when the algorithm is set to 'global_linear', it only makes sense with 'global_affine'\n" << "\n";
		  exit(-1);
	  }
#endif
	  AlignPairsGlobal (vars);
  }
  else if( vars.algorithm == "local")
  {

#ifdef SECURE
	  if( cmd.IsFlagSet( "thresholds_for_penalties"))
	  {
		  std::cerr << "ERROR: the 'thresholds_for_penalties' flag has no meaning when the algorithm is set to 'local', it only makes sense with 'global_affine'\n" << "\n";
		  exit(-1);
	  }
#endif
	  AlignPairsLocal (vars);
  }
  else if( vars.algorithm == "global_affine")
  {
	  AlignPairsGlobalAffine (vars);
  }

  //out.close();
  //out.clear();
  // 5.4.17 Problems before this point

  /////////////////////////////////////
  /// WRITE ALIGNMENT
  ////////////////////////////////////

  WriteAlignment( cmd, vars);


#ifdef SECURE
  if( cmd.IsFlagSet("completed"))
  {
	  out.open( cmd.GetArguments("completed")[0].c_str());
	  out << "COMPLETED CORRECTLY" << std::endl;
	  out.close(); out.clear();
  }
#endif
/*
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *fwd     = NULL;
  P7_GMX         *bck     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fsc, bsc;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  sq     = esl_sq_CreateDigital(abc);
  if (esl_sqfile_Open(seqfile, format, NULL, &sqfp) != eslOK) p7_Fail("Failed to open sequence file %s", seqfile);
  if (esl_sqio_Read(sqfp, sq)                       != eslOK) p7_Fail("Failed to read sequence");
  esl_sqfile_Close(sqfp);

  bg = p7_bg_Create(abc);
*/


  return 0;
}

