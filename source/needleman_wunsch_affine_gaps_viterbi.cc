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
	std::vector< double> old_prior;
	std::vector< double> null_prior;
	double                                          d_inf = std::numeric_limits< double>::infinity();
	
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


	old_prior.resize(pTMX->size());
	null_prior.resize(pTMX->size());
	for( size_t i = 0; i < null_prior.size(); ++i)
	{
		null_prior[i] = -d_inf;
	}
        for( size_t i = 0; i < m_Matrix.GetNumberOfRows(); ++i)  // first sequence
                for( size_t j = 0; j < m_Matrix.GetNumberOfColumns(); ++j)  // second sequence
		{
//			std::cout << __FUNCTION__ << "   entry #  " << i << "  "<< j << std::endl;
			m_Matrix(i, j).SetPriors( pTMX->size());
//			std::cout << __FUNCTION__ << "   Setting Matrix   " << i << "   " << j << std::endl;
			m_Matrix(i, j).SetMatrices( pTMX, pMMX, pEMI, &(m_HMMProfile->PMX[0]));
		}
	std::cout << __FUNCTION__ << "   Setting Matrix: done!" << std::endl;

	// First element
	max_value.first = m_Matrix( 0, 0).GetPrior( 0);
	max_value.second = 0.0;
	for (int k=0; k<9; ++k)
	{
		m_Matrix( 0, 0).AddValue( k, max_value);
	}

	for (int k=0; k<9; ++k)
	{
		if (k > 2 && k < 6)
		{
			max_value.second = 0.0;
		}
		else
		{
			max_value.second = -d_inf;
		}
		m_Matrix( 0, 1).AddValue( k, max_value);	
	}

	for (int k=0; k<9; ++k)
	{
		if (k > 5)
		{
			max_value.second = 0.0;
		}
		else
		{
			max_value.second = -d_inf;
		}
		m_Matrix( 1, 0).AddValue( k, max_value);	
	}

	// First row \ {(0,0), (0,1)}
	for( size_t j = 2; j < m_Matrix.GetNumberOfColumns(); ++j)
	{
		for (int k=0; k<9; ++k)
		{
			if (k == 4)
			{
				max_value.second = -(m_TerminiGapOpeningPenalty + (j-1)*m_TerminiGapExtensionPenalty);
			}
			else
			{
				max_value.second = -d_inf;
			}
			m_Matrix( 0, j).AddValue( k, max_value);
		}
	}

	// First column \ {(0,0), (1,0)}
	for( size_t i = 2; i < m_Matrix.GetNumberOfRows(); ++i)
	{
		for (int k=0; k<9; ++k)
		{
			if (k == 8)
			{
				max_value.second = -(m_TerminiGapOpeningPenalty + (i-1)*m_TerminiGapExtensionPenalty);
			}
			else
			{
				max_value.second = -d_inf;
			}
			m_Matrix( i, 0).AddValue( k, max_value);
		}
	}
	

	// Core
	char aatypem, aatype;
	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() -1; ++i)  // first sequence
	{
		aatypem = 'X';
		for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() -1; ++j)  // second sequence
		{
			aatype = m_SecondSequence[j -1].GetType();

			std::cout << __FUNCTION__ << "   " << i << "   " << j << " type(i-1) " << aatypem << " type(i) " << aatype << std::endl;

			// DIAGONAL
			element = m_Matrix( i -1, j -1); // access previous diagonal (aligned) element

			local_score =  m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[ j - 1])); // similarity score, note that i is matrix index, matrix has zero in first element, sequence does not, that's why i-1 and j-1 are used here
			max_value = element.BestSubPathWayScore_PureV( 0, aatype, 0);
                        max_value.second += local_score;
			m_Matrix( i, j).AddValue( 0, max_value);
			max_value = element.BestSubPathWayScore_PureV( 1, aatype, 0);
                        max_value.second += local_score;
			m_Matrix( i, j).AddValue( 1, max_value);
			max_value = element.BestSubPathWayScore_PureV( 2, aatype, 0);
                        max_value.second += local_score;
			m_Matrix( i, j).AddValue( 2, max_value);

			// HORIZONTAL
			element = m_Matrix( i, j-1);

			max_value = element.BestSubPathWayScore_PureV( 0, aatype, 1);
                        profiles = MinPerElement( m_SecondSequence[ j].GetProfiles(),m_SecondSequence[ j - 1].GetProfiles());
                        max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 3, max_value);
			max_value = element.BestSubPathWayScore_PureV( 1, aatype, 1);
                        max_value.second -= m_GapExtensionPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 4, max_value);
			max_value = element.BestSubPathWayScore_PureV( 2, aatype, 1);
                        max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 5, max_value);

			// VERTICAL
			element = m_Matrix( i-1, j);

			max_value = element.BestSubPathWayScore_PureV( 0, aatype, 2);
                        profiles = MinPerElement( m_FirstSequence[ i].GetProfiles(), m_FirstSequence[ i - 1].GetProfiles());
                        max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 6, max_value);
			max_value = element.BestSubPathWayScore_PureV( 1, aatype, 2);
                        max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 7, max_value);
			max_value = element.BestSubPathWayScore_PureV( 2, aatype, 2);
                        max_value.second -= m_GapExtensionPenaltyFunction->operator ()( profiles);
			m_Matrix( i, j).AddValue( 8, max_value);

			for (int k=0; k<9; ++k)
			{
				std::cout << " Affine Pathway " << k << " " << m_Matrix( i, j).GetAffinePathWay(k) << std::endl;			
//				std::cout << " Prior " << k << " " << m_Matrix( i, j).GetPrior(k) << std::endl;
			}

			aatypem = aatype;
		}
	}
	std::cout << __FUNCTION__ << "   Matrix core: done!" << std::endl;

	aatypem = 'X';
	size_t i_last = m_Matrix.GetNumberOfRows() - 1; // Last row index
	for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() - 1; ++j)
	{
		std::cout << __FUNCTION__ << " " << j << std::endl;
		aatype = m_SecondSequence[j -1].GetType();

		element_ptr   = &m_Matrix( i_last, j);
		element_1 = m_Matrix( i_last - 1, j - 1);
		element_2 = m_Matrix( i_last, j - 1);
		element_3 = m_Matrix( i_last - 1, j);
		local_score = m_Score->operator()( std::make_pair( m_FirstSequence[ i_last - 1], m_SecondSequence[ j - 1]));

		max_value = element_1.BestSubPathWayScore_PureV( 0, aatype, 0);
                max_value.second += local_score;
		element_ptr->AddValue( 0, max_value);
		max_value = element_1.BestSubPathWayScore_PureV( 1, aatype, 0);
                max_value.second += local_score;
                element_ptr->AddValue( 1, max_value);
		max_value = element_1.BestSubPathWayScore_PureV( 2, aatype, 0);
                max_value.second += local_score;
                element_ptr->AddValue( 2, max_value);

		std::cout << __FUNCTION__ << " " << j << std::endl;
		max_value = element_2.BestSubPathWayScore_PureV( 0, aatype, -1);
		std::cout << __FUNCTION__ << " " << j << std::endl;
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 3, max_value);
		max_value = element_2.BestSubPathWayScore_PureV( 1, aatype, -1);
                max_value.second -= m_TerminiGapExtensionPenalty;
                element_ptr->AddValue( 4, max_value);
		max_value = element_2.BestSubPathWayScore_PureV( 2, aatype, -1);
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 5, max_value);

		std::cout << __FUNCTION__ << " " << j << std::endl;
		profiles = MinPerElement(m_SecondSequence[j].GetProfiles(),m_SecondSequence[j-1].GetProfiles() );
		max_value = element_3.BestSubPathWayScore_PureV( 0, aatype, 2);
                max_value.second -= ( *m_GapOpeningPenaltyFunction)( profiles);
                element_ptr->AddValue( 6, max_value);
		max_value = element_3.BestSubPathWayScore_PureV( 1, aatype, 2);
                max_value.second -= ( *m_GapOpeningPenaltyFunction)( profiles);
                element_ptr->AddValue( 7, max_value);
		max_value = element_3.BestSubPathWayScore_PureV( 2, aatype, 2);
                max_value.second -= ( *m_GapExtensionPenaltyFunction)( profiles);
                element_ptr->AddValue( 8, max_value);

		aatype = aatypem;
	}

	std::cout << __FUNCTION__ << "   Matrix last row: done!" << std::endl;

	size_t j_last = m_Matrix.GetNumberOfColumns() -1;
	aatypem = m_SecondSequence[j_last -2].GetType();
	aatype = m_SecondSequence[j_last -1].GetType();
	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() - 1; ++i)
	{
		element_ptr   = &m_Matrix( i, j_last);
		element_1 = m_Matrix( i - 1, j_last - 1);
		element_2 = m_Matrix( i, j_last - 1);
		element_3 = m_Matrix( i - 1, j_last);
		local_score = m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[ j_last - 1]));

		max_value = element_1.BestSubPathWayScore_PureV( 0, aatype, 0);
                max_value.second += local_score; 
                element_ptr->AddValue( 0, max_value);
		max_value = element_1.BestSubPathWayScore_PureV( 1, aatype, 0);
                max_value.second += local_score; 
                element_ptr->AddValue( 1, max_value);
		max_value = element_1.BestSubPathWayScore_PureV( 2, aatype, 0);
                max_value.second += local_score; 
                element_ptr->AddValue( 2, max_value);

                profiles = MinPerElement(m_FirstSequence[ i].GetProfiles(),m_FirstSequence[ i - 1].GetProfiles());
		max_value = element_2.BestSubPathWayScore_PureV( 0, aatype, 1);
                max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 3, max_value);
		max_value = element_2.BestSubPathWayScore_PureV( 1, aatype, 1);
                max_value.second -= m_GapExtensionPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 4, max_value);
		max_value = element_2.BestSubPathWayScore_PureV( 2, aatype, 1);
                max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 5, max_value);

		
		max_value = element_3.BestSubPathWayScore_PureV( 0, aatype, -2);
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 6, max_value);
		max_value = element_3.BestSubPathWayScore_PureV( 1, aatype, -2);
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 7, max_value);
		max_value = element_3.BestSubPathWayScore_PureV( 2, aatype, -2);
                max_value.second -= m_TerminiGapExtensionPenalty;
                element_ptr->AddValue( 8, max_value);
	}

	std::cout << __FUNCTION__ << "   Matrix last column: done!" << std::endl;

	element_ptr   = &m_Matrix( m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns() -1);
	element_1     = m_Matrix( m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns() -1 -1);
	element_2 	  = m_Matrix( m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns() -1 -1);
	element_3     = m_Matrix( m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns() -1);
	local_score   = m_Score->operator()( std::make_pair( m_FirstSequence[  m_Matrix.GetNumberOfRows() -1 -1], m_SecondSequence[ m_Matrix.GetNumberOfColumns()-1 - 1]));

	max_value = element_1.BestSubPathWayScore_PureV( 0, aatype, 0);
        max_value.second += local_score; 
        element_ptr->AddValue( 0, max_value);

	max_value = element_1.BestSubPathWayScore_PureV( 1, aatype, 0);
        max_value.second += local_score; 
        element_ptr->AddValue( 1, max_value);

	max_value = element_1.BestSubPathWayScore_PureV( 2, aatype, 0);
        max_value.second += local_score; 
        element_ptr->AddValue( 2, max_value);


	max_value = element_2.BestSubPathWayScore_PureV( 0, aatype, -1);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 3, max_value);

	max_value = element_2.BestSubPathWayScore_PureV( 1, aatype, -1);
        max_value.second -= m_TerminiGapExtensionPenalty;
        element_ptr->AddValue( 4, max_value);

	max_value = element_2.BestSubPathWayScore_PureV( 2, aatype, -1);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 5, max_value);


	max_value = element_3.BestSubPathWayScore_PureV( 0, aatype, -2);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 6, max_value);

	max_value = element_3.BestSubPathWayScore_PureV( 1, aatype, -2);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 7, max_value);

	max_value = element_3.BestSubPathWayScore_PureV( 2, aatype, -2);
        max_value.second -= m_TerminiGapExtensionPenalty;
        element_ptr->AddValue( 8, max_value);

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
		best = -std::numeric_limits< double>::max(),
		best2 = -std::numeric_limits< double>::max();

	int
		maxint = std::numeric_limits<int>::max(),    // Short name for calling the maximum integer
		position = 0,
		start,
		end, 
		group;

	std::list<std::pair<int, int> > 
		alignment;    // List of pairs containing the alignments. Indices start from . The length of the list is not known a priori.
	std::pair< int, int> 
		endpoint;    // Pair of indices defining the endpoint. The endpoint is not known a priori when HMM is present.



	endpoint = std::make_pair(m_Matrix.GetNumberOfRows() - 1, m_Matrix.GetNumberOfColumns() - 1);
	std::cout << " ENDPOINT " << endpoint.first << " " << endpoint.second << std::endl;

//	alignment.push_front(endpoint);   //This is an error: it elongates "alignment" by one, making the function WriteAlignedSequencesInClustalwFormat in include/alignment_write.h fail.

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

	int next_ind, indchk;
        int i = endpoint.first;
        int j = endpoint.second;

	double diffscore;

        while (i > 0 || j > 0)
	{
		if (group == 0)
		{
			std::cout << __FUNCTION__ << "  M  " << i-1 << "   " << j-1 << std::endl;
			alignment.push_front(std::make_pair(i - 1, j - 1));
			--i;
			--j;
		}
		else if (group == 1)
		{
			std::cout << __FUNCTION__ << "  I  " << i << "   " << j-1 << std::endl;
			alignment.push_front(std::make_pair(maxint, j -1));
			--j;
		}
		else if (group == 2)
		{
			std::cout << __FUNCTION__ << "  D  " << i-1 << "   " << j << std::endl;
			alignment.push_front(std::make_pair(i -1, maxint));
			--i;
		}

		group = position;
                start = position * 3;
                end = (position * 3) + 3;
//              std::cout << "group: " << group << " position: " << position << " start: " << start << " end: " << end << " i: " << i << " j: " << j << std::endl;

                best2 = -std::numeric_limits< double>::max();

                for (int count = start; count < end; ++count)
                {
			std::cout << " candidate " << count << " " << m_Matrix(i, j).GetAffinePathWay( count) << std::endl;
                        if (m_Matrix(i, j).GetAffinePathWay( count) > best2)
                        {
                                best2 = m_Matrix(i, j).GetAffinePathWay( count);
                                position = count;
                        }
                }
		std::cout << " best " << position << " " << best2 << std::endl;
//              std::cout << "pos: " << position << std::endl;
                position -= start;
	}

	std::vector< std::pair< int, int> > 
		converted( alignment.begin(), alignment.end());

	for(std::list<std::pair<int, int> >::iterator it = alignment.begin(); it != alignment.end(); ++it)
//		std::cout << __FUNCTION__  << " " << it->first << " " << it->second << std::endl;

	return std::make_pair( best-diffscore, converted);
}


const Matrix< DynamicProgrammingMatrixElement> &
NeedlemanWunschAffineGaps::GetMatrix() const
{
	return m_Matrix;
}
