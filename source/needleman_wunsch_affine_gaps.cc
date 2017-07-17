/////////////////////////////////////////////////////////////////////////
//!
//  This file is part of the AlignMe program
//
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
//!  Needleman Wunsch algorithm with affine gap penalties
//! (different penalties for opening and continuing gaps).
//! As described in "Biological sequence analysis"
//! from Durbin, Eddy, Krogh, Mitchison
//!
//!
//! @author: Marcus Stamm, Kamil Khafizov,
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////


#include "../include/needleman_wunsch_affine_gaps.h"
#include <limits>
#include <climits>

// default constructor
NeedlemanWunschAffineGaps::NeedlemanWunschAffineGaps()
  : m_Matrix( 0, 0)
{}


// construct from data
// needleman_wunsch( gap_opening_penalty, gap_extension_penalty, termini_gap_opening_penalty, termini_gap_extension_penalty,
// first_sequence, second_sequence, scores);
// assign these values directly to class members with an initialization-list

NeedlemanWunschAffineGaps::NeedlemanWunschAffineGaps
( 
		 const ShPtr< Function< std::vector< double>, double> > &GAP_OPENING_PENALTY_FCT,
		 const ShPtr< Function< std::vector< double>, double> > &GAP_EXTENSION_PENALTY_FCT,
		 const double &TERMIN_GAP_OPENING_PENALTY,
		 const double &TERMIN_GAP_EXTENSION_PENALTY,
		 const Sequence &FIRST_SEQUENCE,
		 const Sequence &SECOND_SEQUENCE,
		 const ShPtr< Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double> > &SCORE,
                 const std::pair< ShPtr< Vec3>, ShPtr< Vec3> > &HMMProfilePair,
		 Matrix< DynamicProgrammingMatrixElement> &MATRIX
):
	  m_GapOpeningPenaltyFunction( GAP_OPENING_PENALTY_FCT),   // build member: call constructor of double and initialize with value GAP_OPENING_PENALTY
	  m_GapExtensionPenaltyFunction( GAP_EXTENSION_PENALTY_FCT),
	  m_TerminiGapOpeningPenalty( TERMIN_GAP_OPENING_PENALTY),
	  m_TerminiGapExtensionPenalty( TERMIN_GAP_EXTENSION_PENALTY),
	  m_FirstSequence( FIRST_SEQUENCE),
	  m_SecondSequence( SECOND_SEQUENCE),
	  m_Score( SCORE),
          m_HMMProfilePair( HMMProfilePair),
	  m_Matrix( MATRIX)
{
	m_HMMProfile = m_HMMProfilePair.first;
	m_HMMForward = m_HMMProfilePair.second;
}


std::vector< double> CalculateUncertainty( ShPtr< Vec3> &FOWBAW, ShPtr< Vec3> &FOW, std::vector< double> &UNCERT)
{
	int N = FOWBAW->PMX.size();   //Length of Sequence
	int M = FOWBAW->PMX[0].size();  //# of nodes in HMM
	UNCERT.resize(N);
	for (int i = 0; i<N; i++)
	{
		UNCERT[i] = 0.0;
		for (int j = 0; j<M; j++)
		{
//			std::cout << __FUNCTION__ << "   Building UNCERT" << std::endl;
//			std::cout << __FUNCTION__ << FOW->PMX.size() << std::endl;
			UNCERT[i] += std::abs(FOWBAW->PMX[i][j] - FOW->PMX[i][j])/M;    //This is a trivial way to estimate the error... Have a better idea!
		}
	}
	return UNCERT;
}


void NeedlemanWunschAffineGaps::CalculateMatrix()
{
	const std::vector< std::vector < double> >*	pTMX = &(m_HMMProfile->TMX);
	const std::vector < std::pair< int, int> >*	pMMX[3] = {&(m_HMMProfile->MMX[0]), &(m_HMMProfile->MMX[1]), &(m_HMMProfile->MMX[2])};
	const std::map< char, std::vector< double> >*   pEMI = &(m_HMMProfile->EMI);
	size_t v_size = pTMX->size();
	
	DynamicProgrammingMatrixElement
		element( v_size),
		* element_ptr( new DynamicProgrammingMatrixElement(v_size)),
		element_1( v_size),
		element_2( v_size),
		element_3( v_size);

	double
		local_score;
	std::vector< double>
		profiles;

        std::pair< std::vector< double>, double> max_value;

//	for( size_t i = 0; i < 9; ++i)
//	{
//		m_Matrix( 0, 0).AddValue(i, 0);
//	}

	// each dynamic-programing-matrix-element has a 9 element affine-path-way member
	// first three represent the previous diagonal (times its three possible previous elements)
	// second three go horizontal previous element, last three vertical previous element
	// each triplet follows the same structure

        for( size_t i = 0; i < m_Matrix.GetNumberOfRows(); ++i)  // first sequence
                for( size_t j = 0; j < m_Matrix.GetNumberOfColumns(); ++j)  // second sequence
		{
//			std::cout << __FUNCTION__ << "   entry #  " << i << "  "<< j << std::endl;
			m_Matrix(i, j).SetPriors( pTMX->size());
//			std::cout << __FUNCTION__ << "   Setting Matrix   " << i << "   " << j << std::endl;
			m_Matrix(i, j).SetMatrices( pTMX, pMMX, pEMI, &(m_HMMProfile->PMX[0]));
		}

	::CalculateUncertainty( m_HMMProfile, m_HMMForward, m_uncertainty);

	std::cout << __FUNCTION__ << "   Setting Matrix: done!" << std::endl;
//	std::vector< double> nullV;
//	nullV.resize(pMMX[0]->size());
//	for (size_t i = 0; i < pMMX[0]->size(); i++)
//		nullV[i] = 0.0;

	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() -1; ++i)  // first sequence
	{
/*		const std::vector< double> * p_pmx3m, p_pmx2m;
		const char aatype3m, aatype2m;
		const double indet3m, indet2m;
		if ( i > 1)
		{
			aatype2m = m_SecondSequence[i -2].GetType();
			p_pmx2m = &m_HMMProfile->PMX[i -2];
			indet2m = m_uncertainty[i -2];
			if (i > 2)
			{
				aatype3m = m_SecondSequence[i -3].GetType();
				p_pmx3m = &m_HMMProfile->PMX[i -3];
				indet3m = m_uncertainty[i -3];
			}
			else
			{
				aatype3m = 'A';
				p_pmx3m = &(nullV);
				indet3m = 0.0;
			}
		}
		else
		{
			aatype2m = 'A';
			p_pmx2m = &(nullV);
			indet2m = 0.0;
		}
*/		
		const std::vector< double> * p_pmxm = &m_HMMProfile->PMX[i -1];
		const std::vector< double> * p_pmx = &m_HMMProfile->PMX[i];
		const char aatypem = m_SecondSequence[i -1].GetType();
		const char aatype = m_SecondSequence[i].GetType();
		const double indetm = m_uncertainty[i -1];
		const double indet = m_uncertainty[i];

		for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() -1; ++j)  // second sequence
		{	
//			std::cout << __FUNCTION__ << "   " << i << "   " << j << std::endl;
			//std::cout << __FUNCTION__ << "   entry ##  " << m_Matrix.GetNumberOfRows() << "  "<< m_Matrix.GetNumberOfColumns() << std::endl;
			// DIAGONAL
			element = m_Matrix( i -1, j -1); // access previous diagonal (aligned) element
			local_score =  m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[ j - 1])); // similarity score, note that i is matrix index, matrix has zero in first element, sequence does not, that's why i-1 and j-1 are used here
			//std::cout << __FUNCTION__ << "   here  " << p_pmxm[0] << std::endl;
                        max_value = element.BestSubPathWayScore( 0, p_pmxm, aatypem, indetm);
                        max_value.second += local_score;
			m_Matrix( i, j).AddValue( 0, max_value); // fill array containing affine path ways, first three values from diagonal step
//			m_EndElement.AddValue(i-1, j-1, 0, max_value);
//			std::cout << __FUNCTION__ << "   here  " << std::endl;
                        max_value = element.BestSubPathWayScore( 1, p_pmxm, aatypem, indetm);
                        max_value.second += local_score;
			m_Matrix( i, j).AddValue( 1, max_value);
//			m_EndElement.AddValue(i-1, j-1, 1, max_value);
                        max_value = element.BestSubPathWayScore( 2, p_pmxm, aatypem, indetm);
                        max_value.second += local_score;
			m_Matrix( i, j).AddValue( 2, max_value);
//			m_EndElement.AddValue(i-1, j-1, 2, max_value);
//			std::cout << __FUNCTION__ << "   DIAG  " << std::endl;

			// HORIZONTAL
			element = m_Matrix( i -1, j);
			max_value = element.BestSubPathWayScore( 0, p_pmxm, aatypem, indetm);
                        profiles = MinPerElement( m_SecondSequence[ j].GetProfiles(),m_SecondSequence[ j - 1].GetProfiles());
                        max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 3, max_value);
//			m_EndElement.AddValue(i-1, j, 0, max_value);
			max_value = element.BestSubPathWayScore( 1, p_pmxm, aatypem, indetm);
                        max_value.second -= m_GapExtensionPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 4, max_value);
//			m_EndElement.AddValue(i-1, j, 1, max_value);
			max_value = element.BestSubPathWayScore( 2, p_pmxm, aatypem, indetm);
                        max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 5, max_value);
//			m_EndElement.AddValue(i-1, j, 2, max_value);
//			std::cout << __FUNCTION__ << "   HOR  " << std::endl;

			// VERTICAL
//			std::cout << __FUNCTION__ << "   here  " << std::endl;
			element = m_Matrix( i, j -1);
			max_value = element.BestSubPathWayScore( 0, p_pmx, aatype, indet);
                        profiles = MinPerElement( m_FirstSequence[ i].GetProfiles(), m_FirstSequence[ i - 1].GetProfiles());
                        max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 6, max_value);
//			m_EndElement.AddValue(i, j-1, 0, max_value);
			max_value = element.BestSubPathWayScore( 1, p_pmx, aatype, indet);
                        max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 7, max_value);
//			m_EndElement.AddValue(i, j-1, 1, max_value);
			max_value = element.BestSubPathWayScore( 2, p_pmx, aatype, indet);
                        max_value.second -= m_GapExtensionPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 8, max_value);
//			m_EndElement.AddValue(i, j-1, 2, max_value);
//			std::cout << __FUNCTION__ << "   VER  " << std::endl;
		}
	}
	std::cout << __FUNCTION__ << "   Matrix core: done!" << std::endl;

	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() - 1; ++i)
	{
		const std::vector< double> * p_pmxm = &m_HMMProfile->PMX[i -1];
		const std::vector< double> * p_pmx = &m_HMMProfile->PMX[i];
		const char aatypem = m_SecondSequence[i -1].GetType();
		const char aatype = m_SecondSequence[i].GetType();
		const double indetm = m_uncertainty[i -1];
		const double indet = m_uncertainty[i];

		element_ptr   = &m_Matrix( i, m_Matrix.GetNumberOfColumns() - 1);
		element_1 = m_Matrix( i - 1, m_Matrix.GetNumberOfColumns() - 1 - 1);
		element_2 = m_Matrix( i - 1, m_Matrix.GetNumberOfColumns() - 1);
		element_3 = m_Matrix( i, m_Matrix.GetNumberOfColumns() - 1 - 1);
		local_score = m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[  m_Matrix.GetNumberOfColumns()-1 - 1]));

                max_value = element_1.BestSubPathWayScore(0, p_pmxm, aatypem, indetm);
                max_value.second += local_score;
		element_ptr->AddValue( 0, max_value);
//		m_EndElement.AddValue(i-1, m_Matrix.GetNumberOfColumns()-1-1, 0, max_value);
                max_value = element_1.BestSubPathWayScore(1, p_pmxm, aatypem, indetm);
                max_value.second += local_score;
                element_ptr->AddValue( 1, max_value);
//		m_EndElement.AddValue(i-1, m_Matrix.GetNumberOfColumns()-1-1, 1, max_value);
                max_value = element_1.BestSubPathWayScore(2, p_pmxm, aatypem, indetm);
                max_value.second += local_score;
                element_ptr->AddValue( 2, max_value);
//		m_EndElement.AddValue(i-1, m_Matrix.GetNumberOfColumns()-1-1, 2, max_value);

                max_value = element_2.BestSubPathWayScore(0, p_pmxm, aatypem, indetm);
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 3, max_value);
//		m_EndElement.AddValue(i-1, m_Matrix.GetNumberOfColumns()-1, 0, max_value);
                max_value = element_2.BestSubPathWayScore(1, p_pmxm, aatypem, indetm);
                max_value.second -= m_TerminiGapExtensionPenalty;
                element_ptr->AddValue( 4, max_value);
//		m_EndElement.AddValue(i-1, m_Matrix.GetNumberOfColumns()-1, 1, max_value);
                max_value = element_2.BestSubPathWayScore(2, p_pmxm, aatypem, indetm);
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 5, max_value);
//		m_EndElement.AddValue(i-1, m_Matrix.GetNumberOfColumns()-1, 2, max_value);

		profiles = MinPerElement(m_FirstSequence[i].GetProfiles(),m_FirstSequence[i-1].GetProfiles() );
                max_value = element_3.BestSubPathWayScore(0, p_pmx, aatype, indet);
                max_value.second -= ( *m_GapOpeningPenaltyFunction)( profiles);
                element_ptr->AddValue( 6, max_value);
//		m_EndElement.AddValue(i, m_Matrix.GetNumberOfColumns()-1-1, 0, max_value);
                max_value = element_3.BestSubPathWayScore(1, p_pmx, aatype, indet);
                max_value.second -= ( *m_GapOpeningPenaltyFunction)( profiles);
                element_ptr->AddValue( 7, max_value);
//		m_EndElement.AddValue(i, m_Matrix.GetNumberOfColumns()-1-1, 1, max_value);
                max_value = element_3.BestSubPathWayScore(2, p_pmx, aatype, indet);
                max_value.second -= ( *m_GapExtensionPenaltyFunction)( profiles);
                element_ptr->AddValue( 8, max_value);
//		m_EndElement.AddValue(i, m_Matrix.GetNumberOfColumns()-1-1, 2, max_value);
	}

	std::cout << __FUNCTION__ << "   Matrix last row: done!" << std::endl;

	int nrows = m_Matrix.GetNumberOfRows() -1;
	const std::vector< double> * p_pmxm = &m_HMMProfile->PMX[nrows -1];
	const std::vector< double> * p_pmx = &m_HMMProfile->PMX[nrows];
	const char aatypem = m_SecondSequence[nrows -1].GetType();
	const char aatype = m_SecondSequence[nrows].GetType();
	const double indetm = m_uncertainty[nrows -1];
	const double indet = m_uncertainty[nrows];

	for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() - 1; ++j)
	{
		element_ptr   = &m_Matrix( m_Matrix.GetNumberOfRows() - 1, j);
		element_1 = m_Matrix( m_Matrix.GetNumberOfRows() - 1 - 1, j - 1);
		element_2 = m_Matrix( m_Matrix.GetNumberOfRows() - 1 - 1, j);
		element_3 = m_Matrix( m_Matrix.GetNumberOfRows() - 1, j - 1);
		local_score = m_Score->operator()( std::make_pair( m_FirstSequence[  m_Matrix.GetNumberOfRows()-1 - 1], m_SecondSequence[ j - 1]));

                max_value = element_1.BestSubPathWayScore(0, p_pmxm, aatypem, indetm);
                max_value.second += local_score; 
                element_ptr->AddValue( 0, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1-1, j-1, 0, max_value);
                max_value = element_1.BestSubPathWayScore(1, p_pmxm, aatypem, indetm);
                max_value.second += local_score; 
                element_ptr->AddValue( 1, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1-1, j-1, 1, max_value);
                max_value = element_1.BestSubPathWayScore(2, p_pmxm, aatypem, indetm);
                max_value.second += local_score; 
                element_ptr->AddValue( 2, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1-1, j-1, 2, max_value);

                profiles = MinPerElement(m_SecondSequence[ j].GetProfiles(),m_SecondSequence[ j - 1].GetProfiles());
                max_value = element_2.BestSubPathWayScore(0, p_pmxm, aatypem, indetm);
                max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 3, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1-1, j, 0, max_value);
                max_value = element_2.BestSubPathWayScore(1, p_pmxm, aatypem, indetm);
                max_value.second -= m_GapExtensionPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 4, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1-1, j, 1, max_value);
                max_value = element_2.BestSubPathWayScore(2, p_pmxm, aatypem, indetm);
                max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 5, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1-1, j, 2, max_value);

                max_value = element_3.BestSubPathWayScore(0, p_pmx, aatype, indet); 
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 6, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1, j-1, 0, max_value);
                max_value = element_3.BestSubPathWayScore(1, p_pmx, aatype, indet); 
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 7, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1, j-1, 1, max_value);
                max_value = element_3.BestSubPathWayScore(2, p_pmx, aatype, indet);
                max_value.second -= m_TerminiGapExtensionPenalty;
                element_ptr->AddValue( 8, max_value);
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows()-1, j-1, 2, max_value);
	}

	std::cout << __FUNCTION__ << "   Matrix last column: done!" << std::endl;

	element_ptr   = &m_Matrix( m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns() -1);
	element_1     = m_Matrix( m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1-1);
	element_2 	  = m_Matrix( m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1);
	element_3     = m_Matrix( m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns()-1-1);
	local_score   = m_Score->operator()( std::make_pair( m_FirstSequence[  m_Matrix.GetNumberOfRows()-1 - 1], m_SecondSequence[ m_Matrix.GetNumberOfColumns()-1 - 1]));

        max_value = element_1.BestSubPathWayScore(0, p_pmxm, aatypem, indetm);
        max_value.second += local_score; 
        element_ptr->AddValue( 0, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1-1, 0, max_value);
        max_value = element_1.BestSubPathWayScore(1, p_pmxm, aatypem, indetm);
        max_value.second += local_score; 
        element_ptr->AddValue( 1, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1-1, 1, max_value);
        max_value = element_1.BestSubPathWayScore(2, p_pmxm, aatypem, indetm);
        max_value.second += local_score; 
        element_ptr->AddValue( 2, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1-1, 2, max_value);

        max_value = element_2.BestSubPathWayScore(0, p_pmxm, aatypem, indetm);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 3, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1, 0, max_value);
        max_value = element_2.BestSubPathWayScore(1, p_pmxm, aatypem, indetm);
        max_value.second -= m_TerminiGapExtensionPenalty;
        element_ptr->AddValue( 4, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1, 1, max_value);
        max_value = element_2.BestSubPathWayScore(2, p_pmxm, aatypem, indetm);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 5, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1, 2, max_value);

        max_value = element_3.BestSubPathWayScore(0, p_pmx, aatype, indet);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 6, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns()-1-1, 0, max_value);
        max_value = element_3.BestSubPathWayScore(1, p_pmx, aatype, indet);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 7, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns()-1-1, 1, max_value);
        max_value = element_3.BestSubPathWayScore(2, p_pmx, aatype, indet);
        max_value.second -= m_TerminiGapExtensionPenalty;
        element_ptr->AddValue( 8, max_value);
//	m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns()-1-1, 2, max_value);

//	for (int k = 0; k < 3; ++k)
//		m_EndElement.AddValue(m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns()-1, k, element_ptr->BestSubPathWayScore(k, p_pmx, aatype, indet));

	/*
	// Records the end node: searches for the highest "exit score" in the matrix.
	// The exit score must also be greater than the "evolution score" of that entry.
	// In case the exit score of a non-terminal entry is the greatest and the entry is chosen as terminal one,
	//  does this invalidate the construction of the rest of the matrix? No, since it is as if we carried on the
	//  evolution of the second-best pathway in that position, the best one being the exit. This is consistent
	//  with the dynamic programming phylosophy (and the matrix can be used if global alignment is forced).
	double 
		best = -std::numeric_limits<int>::max(),
		evolution_score,
		exit_score;
	for(unsigned int i=0; i<m_Matrix.GetNumberOfRows(); i++)
	{
		for(unsigned int j=0; j<m_Matrix.GetNumberOfRows(); j++)
		{
			exit_score = m_Matrix(i,j).EndNode();
			evolution_score = m_Matrix(i,j).GetScore();
			if(best < exit_score && exit_score > evolution_score)
			{
				best = m_Matrix(i,j).EndNode();
				m_EndCoord = std::make_pair(i, j);
			}
		}
	}
	*/

	std::cout << __FUNCTION__ << "   Matrix: done!" << std::endl;
}


std::pair< double, std::vector< std::pair< int, int> > >
NeedlemanWunschAffineGaps:: TraceBack() const
{
	//Declarations
	std::vector< std::vector < double> >*
		// Pointer to the transition matrix. Size: (3M+2, 3M+2) where M = number of nodes of the HMM
		pTMX = &(m_HMMProfile->TMX);

        std::vector<std::pair<int, int> > *
		// Array of 4 pointers to the mask matrices. Since they are diagonal, they are stored as vectors. Size: 3M+2 where M = number of nodes of the HMM  
		pMMX[3] = {&(m_HMMProfile->MMX[0]), &(m_HMMProfile->MMX[1]), &(m_HMMProfile->MMX[2])};

	size_t 
		v_size = pTMX->size();    // 3M+2 where M = number of nodes of the HMM

        DynamicProgrammingMatrixElement
                element( v_size);    // Constructor initializes 12 Prior vectors. Size: 3M+2 where M = number of nodes of the HMM

	double 
		best = -std::numeric_limits< double>::max();

	int
		maxint = std::numeric_limits<int>::max(),    // Short name for calling the maximum integer
		position = 0, 
		group;

	std::list<std::pair<int, int> > 
		alignment;    // List of pairs containing the alignments. Indices start from . The length of the list is not known a priori.
	std::pair< int, int> 
		endpoint;    // Pair of indices defining the endpoint. The endpoint is not known a priori when HMM is present.

//----This part needs structure------------
	// if HMM...
	//endpoint = m_EndCoord;
	// otherwise
	endpoint = std::make_pair(m_Matrix.GetNumberOfRows() - 1, m_Matrix.GetNumberOfColumns() - 1);
//-----------------------------------------

	alignment.push_front(endpoint);

	{
		// Use copy constructor?
		DynamicProgrammingMatrixElement element = m_Matrix(endpoint.first, endpoint.second);

		for (size_t i = 0; i < 9; i++)
		{
			if ( element.GetAffinePathWay( i) > best)
			{
				best = element.GetAffinePathWay( i);
				position = i;
			}
		}
	}

	group =  int (double (position) / 3.0);
	position = position - 3*group;

	std::cout << __FUNCTION__ << " here " << std::endl;

	int next_ind, indchk;
        size_t i = endpoint.first;
        size_t j = endpoint.second;

	std::cout << __FUNCTION__ << " here " << std::endl;

        do
	{
		if (group == 0)
		{
			std::cout << __FUNCTION__ << "  M  " << i-1 << "   " << j-1 << std::endl;
			alignment.push_front(std::make_pair(i - 1, j - 1));
			element = m_Matrix(i - 1, j - 1);
			--i;
			--j;
		}
		else if (group == 1)
		{
			std::cout << __FUNCTION__ << "  I  " << i-1 << "   " << j << std::endl;
			alignment.push_front(std::make_pair(i - 1, maxint));
			element = m_Matrix(i - 1, j);
			--i;
		}
		else if (group == 2)
		{
			std::cout << __FUNCTION__ << "  D  " << i << "   " << j-1 << std::endl;
			alignment.push_front(std::make_pair(maxint, j - 1));
			element = m_Matrix(i, j - 1);
			--j;
		}

		std::cout << __FUNCTION__ << " here " << std::endl;
		group = position;
		position = element.GetBestPath( group);
		std::cout << __FUNCTION__ << position << std::endl;
	} while ( i > 0 || j > 0);

	double 
		diffscore = m_Matrix(i, j).GetValue();

	std::vector< std::pair< int, int> > 
		converted( alignment.begin(), alignment.end());

	return std::make_pair( best-diffscore, converted);
}


const Matrix< DynamicProgrammingMatrixElement> &
NeedlemanWunschAffineGaps::GetMatrix() const
{
	return m_Matrix;
}
