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
#include <omp.h>

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
                 const ShPtr< Vec3> &HMMProfile,
		 Matrix< DynamicProgrammingMatrixElement> &MATRIX
 )
:
	  m_GapOpeningPenaltyFunction( GAP_OPENING_PENALTY_FCT),   // build member: call constructor of double and initialize with value GAP_OPENING_PENALTY
	  m_GapExtensionPenaltyFunction( GAP_EXTENSION_PENALTY_FCT),
	  m_TerminiGapOpeningPenalty( TERMIN_GAP_OPENING_PENALTY),
	  m_TerminiGapExtensionPenalty( TERMIN_GAP_EXTENSION_PENALTY),
	  m_FirstSequence( FIRST_SEQUENCE),
	  m_SecondSequence( SECOND_SEQUENCE),
	  m_Score( SCORE),
          m_HMMProfile( HMMProfile),
	  m_Matrix( MATRIX)
{}


void NeedlemanWunschAffineGaps::CalculateMatrix()
{
	std::vector< std::vector < double> >*	pTMX = &(m_HMMProfile->TMX);
	std::vector < double>*	pMMX[4] = {&(m_HMMProfile->MMX[0]), &(m_HMMProfile->MMX[1]), &(m_HMMProfile->MMX[2]), &(m_HMMProfile->MMX[3])};

	DynamicProgrammingMatrixElement
		element( pTMX, pMMX),
		* element_ptr( new DynamicProgrammingMatrixElement(pTMX, pMMX)),
		element_1( pTMX, pMMX),
		element_2( pTMX, pMMX),
		element_3( pTMX, pMMX);

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
			//std::cout << __FUNCTION__ << "   entry #  " << i << "  "<< j << std::endl;
			m_Matrix(i, j).SetPriors( pMMX[0]->size());
			m_Matrix(i, j).SetMatrices( pTMX, pMMX);
		}

	std::cout << __FUNCTION__ << "   Setting Matrix: done!" << std::endl;

	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() -1; ++i)  // first sequence
	{
		for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() -1; ++j)  // second sequence
		{
			std::cout << __FUNCTION__ << "   " << i << "   " << j << std::endl;
			omp_set_num_threads(9);
#pragma omp parallel private (element, max_value, local_score)
			{	
				int ID = omp_get_thread_num();
				//std::cout << __FUNCTION__ << "   entry ##  " << m_Matrix.GetNumberOfRows() << "  "<< m_Matrix.GetNumberOfColumns() << std::endl;
				// DIAGONAL
				if (ID < 3)
				{
					element = m_Matrix( i -1, j -1); // access previous diagonal (aligned) element
					local_score =  m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[ j - 1])); // similarity score, note that i is matrix index, matrix has zero in first element, sequence does not, that's why i-1 and j-1 are used here
					max_value = element.BestSubPathWayScore( ID, m_HMMProfile->PMX[i -1]);
					max_value.second += local_score;
					m_Matrix( i, j).AddValue( ID, max_value); // fill array containing affine path ways, first three values from diagonal step
				}
				else if (ID > 2 && ID < 6)
				{
					// HORIZONTAL
					element = m_Matrix( i -1, j);
					max_value = element.BestSubPathWayScore( ID-3, m_HMMProfile->PMX[i -1]);
					profiles = MinPerElement( m_SecondSequence[ j].GetProfiles(),m_SecondSequence[ j - 1].GetProfiles());
					max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
					m_Matrix( i, j).AddValue( ID, max_value);
				}
				else if (ID > 5)
				{
					// VERTICAL
					element = m_Matrix( i, j -1);
					max_value = element.BestSubPathWayScore( ID-6, m_HMMProfile->PMX[i]);
					profiles = MinPerElement( m_FirstSequence[ i].GetProfiles(), m_FirstSequence[ i - 1].GetProfiles());
					max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
					m_Matrix( i, j).AddValue( ID, max_value);
				}
			}
		}
	}

	std::cout << __FUNCTION__ << "   Matrix core: done!" << std::endl;

	for( size_t i = 1; i < m_Matrix.GetNumberOfRows() - 1; ++i)
	{
		element_ptr   = &m_Matrix( i, m_Matrix.GetNumberOfColumns() - 1);
		element_1 = m_Matrix( i - 1, m_Matrix.GetNumberOfColumns() - 1 - 1);
		element_2 = m_Matrix( i - 1, m_Matrix.GetNumberOfColumns() - 1);
		element_3 = m_Matrix( i, m_Matrix.GetNumberOfColumns() - 1 - 1);
		local_score = m_Score->operator()( std::make_pair( m_FirstSequence[ i - 1], m_SecondSequence[  m_Matrix.GetNumberOfColumns()-1 - 1]));

                max_value = element_1.BestSubPathWayScore(0, m_HMMProfile->PMX[i -1]);
                max_value.second += local_score;
		element_ptr->AddValue( 0, max_value);
                max_value = element_1.BestSubPathWayScore(1, m_HMMProfile->PMX[i -1]);
                max_value.second += local_score;
                element_ptr->AddValue( 1, max_value);
                max_value = element_1.BestSubPathWayScore(2, m_HMMProfile->PMX[i -1]);
                max_value.second += local_score;
                element_ptr->AddValue( 2, max_value);

                max_value = element_2.BestSubPathWayScore(0, m_HMMProfile->PMX[i -1]);
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 3, max_value);
                max_value = element_2.BestSubPathWayScore(1, m_HMMProfile->PMX[i -1]);
                max_value.second -= m_TerminiGapExtensionPenalty;
                element_ptr->AddValue( 4, max_value);
                max_value = element_2.BestSubPathWayScore(2, m_HMMProfile->PMX[i -1]);
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 5, max_value);

		profiles = MinPerElement(m_FirstSequence[i].GetProfiles(),m_FirstSequence[i-1].GetProfiles() );
                max_value = element_3.BestSubPathWayScore(0, m_HMMProfile->PMX[i]);
                max_value.second -= ( *m_GapOpeningPenaltyFunction)( profiles);
                element_ptr->AddValue( 6, max_value);
                max_value = element_3.BestSubPathWayScore(1, m_HMMProfile->PMX[i]);
                max_value.second -= ( *m_GapOpeningPenaltyFunction)( profiles);
                element_ptr->AddValue( 7, max_value);
                max_value = element_3.BestSubPathWayScore(2, m_HMMProfile->PMX[i]);
                max_value.second -= ( *m_GapExtensionPenaltyFunction)( profiles);
                element_ptr->AddValue( 8, max_value);
	}

	std::cout << __FUNCTION__ << "   Matrix last row: done!" << std::endl;

        int nrows = m_Matrix.GetNumberOfRows() -1;
	for( size_t j = 1; j < m_Matrix.GetNumberOfColumns() - 1; ++j)
	{
		element_ptr   = &m_Matrix( m_Matrix.GetNumberOfRows() - 1, j);
		element_1 = m_Matrix( m_Matrix.GetNumberOfRows() - 1 - 1, j - 1);
		element_2 = m_Matrix( m_Matrix.GetNumberOfRows() - 1 - 1, j);
		element_3 = m_Matrix( m_Matrix.GetNumberOfRows() - 1, j - 1);
		local_score = m_Score->operator()( std::make_pair( m_FirstSequence[  m_Matrix.GetNumberOfRows()-1 - 1], m_SecondSequence[ j - 1]));

                max_value = element_1.BestSubPathWayScore(0, m_HMMProfile->PMX[nrows -1]);
                max_value.second += local_score; 
                element_ptr->AddValue( 0, max_value);
                max_value = element_1.BestSubPathWayScore(1, m_HMMProfile->PMX[nrows -1]);
                max_value.second += local_score; 
                element_ptr->AddValue( 1, max_value);
                max_value = element_1.BestSubPathWayScore(2, m_HMMProfile->PMX[nrows -1]);
                max_value.second += local_score; 
                element_ptr->AddValue( 2, max_value);

                profiles = MinPerElement(m_SecondSequence[ j].GetProfiles(),m_SecondSequence[ j - 1].GetProfiles());
                max_value = element_2.BestSubPathWayScore(0, m_HMMProfile->PMX[nrows -1]);
                max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 3, max_value);
                max_value = element_2.BestSubPathWayScore(1, m_HMMProfile->PMX[nrows -1]);
                max_value.second -= m_GapExtensionPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 4, max_value);
                max_value = element_2.BestSubPathWayScore(2, m_HMMProfile->PMX[nrows -1]);
                max_value.second -= m_GapOpeningPenaltyFunction->operator ()( profiles);
                element_ptr->AddValue( 5, max_value);

                max_value = element_3.BestSubPathWayScore(0, m_HMMProfile->PMX[nrows]); 
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 6, max_value);
                max_value = element_3.BestSubPathWayScore(1, m_HMMProfile->PMX[nrows]); 
                max_value.second -= m_TerminiGapOpeningPenalty;
                element_ptr->AddValue( 7, max_value);
                max_value = element_3.BestSubPathWayScore(2, m_HMMProfile->PMX[nrows]);
                max_value.second -= m_TerminiGapExtensionPenalty;
                element_ptr->AddValue( 8, max_value);
	}

	std::cout << __FUNCTION__ << "   Matrix last column: done!" << std::endl;

	element_ptr   = &m_Matrix( m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns() -1);
	element_1     = m_Matrix( m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1-1);
	element_2 	  = m_Matrix( m_Matrix.GetNumberOfRows() -1 -1, m_Matrix.GetNumberOfColumns()-1);
	element_3     = m_Matrix( m_Matrix.GetNumberOfRows() -1, m_Matrix.GetNumberOfColumns()-1-1);
	local_score   = m_Score->operator()( std::make_pair( m_FirstSequence[  m_Matrix.GetNumberOfRows()-1 - 1], m_SecondSequence[ m_Matrix.GetNumberOfColumns()-1 - 1]));

        max_value = element_1.BestSubPathWayScore(0, m_HMMProfile->PMX[nrows -1]);
        max_value.second += local_score; 
        element_ptr->AddValue( 0, max_value);
        max_value = element_1.BestSubPathWayScore(1, m_HMMProfile->PMX[nrows -1]);
        max_value.second += local_score; 
        element_ptr->AddValue( 1, max_value);
        max_value = element_1.BestSubPathWayScore(2, m_HMMProfile->PMX[nrows -1]);
        max_value.second += local_score; 
        element_ptr->AddValue( 2, max_value);

        max_value = element_2.BestSubPathWayScore(0, m_HMMProfile->PMX[nrows -1]);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 3, max_value);
        max_value = element_2.BestSubPathWayScore(1, m_HMMProfile->PMX[nrows -1]);
        max_value.second -= m_TerminiGapExtensionPenalty;
        element_ptr->AddValue( 4, max_value);
        max_value = element_2.BestSubPathWayScore(2, m_HMMProfile->PMX[nrows -1]);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 5, max_value);

        max_value = element_3.BestSubPathWayScore(0, m_HMMProfile->PMX[nrows]);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 6, max_value);
        max_value = element_3.BestSubPathWayScore(1, m_HMMProfile->PMX[nrows]);
        max_value.second -= m_TerminiGapOpeningPenalty;
        element_ptr->AddValue( 7, max_value);
        max_value = element_3.BestSubPathWayScore(2, m_HMMProfile->PMX[nrows]);
        max_value.second -= m_TerminiGapExtensionPenalty;
        element_ptr->AddValue( 8, max_value);


	double best = -std::numeric_limits<int>::max();
	for(unsigned int i=0; i<m_Matrix.GetNumberOfRows() -1; i++)
	{
		for(unsigned int j=0; j<m_Matrix.GetNumberOfRows() -1; j++)
		{
			if(best < m_Matrix(i,j).EndNode())
			{
				best = m_Matrix(i,j).EndNode();
				m_EndCoord = std::make_pair(i, j);
			}
		}
	}

	std::cout << __FUNCTION__ << "   Matrix: done!" << std::endl;
}


std::pair< double, std::vector< std::pair< int, int> > >
NeedlemanWunschAffineGaps:: TraceBack() const
{
	std::vector< std::vector < double> >*   pTMX = &(m_HMMProfile->TMX);
        std::vector < double>*  pMMX[4] = {&(m_HMMProfile->MMX[0]), &(m_HMMProfile->MMX[1]), &(m_HMMProfile->MMX[2]), &(m_HMMProfile->MMX[3])};

        DynamicProgrammingMatrixElement
                element( pTMX, pMMX);

	double 
		best = -std::numeric_limits< double>::max();

	int
		max = std::numeric_limits<int>::max(),
		position = 0,
		group;

	std::list<std::pair<int, int> > alignment;
	std::pair< int, int> endpoint;

	// if HMM...
	endpoint = m_EndCoord;
	// otherwise
	// endpoint = std::make_pair(m_Matrix.GetNumberOfRows() - 1, m_Matrix.GetNumberOfColumns() - 1)

	alignment.push_front(endpoint);

	{
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

	int next_ind, indchk;
        size_t i = endpoint.first;
        size_t j = endpoint.second;

        do
	{
		if (group == 0)
		{
			alignment.push_front(std::make_pair(i - 1, j - 1));
			element = m_Matrix(i - 1, j - 1);
			--i;
			--j;
		}
		else if (group == 1)
		{
			alignment.push_front(std::make_pair(i - 1, max));
			element = m_Matrix(i - 1, j);
			--i;
		}
		else if (group == 2)
		{
			alignment.push_front(std::make_pair(max, j - 1));
			element = m_Matrix(i, j - 1);
			--j;
		}

		next_ind = element.GetBestPath( group);
		group = int (double (next_ind) / 4.0);
		indchk = next_ind % 4;
	} while (indchk != 0 );

	double diffscore = m_Matrix(i, j).GetAffinePathWay( group);

	std::vector< std::pair< int, int> > converted( alignment.begin(), alignment.end());
	return std::make_pair( best-diffscore, converted);

/*
	size_t i = endpoint.first;
	size_t j = endpoint.second;
	group =  int (double (position) / 3.0);
	position -= 3 * group;

	while (i > 0 || j > 0)
	{
		if (group == 0)
		{
			alignment.push_front(std::make_pair(i - 1, j - 1));
			--i;
			--j;
		}

		else if (group == 1)
		{
			alignment.push_front(std::make_pair(i - 1, max));
			--i;
		}

		else if (group == 2)
		{
			alignment.push_front(std::make_pair(max, j - 1));
			--j;
		}

		group = position;
		start = position * 3;
		end = (position * 3) + 3;

		best2 = -std::numeric_limits< double>::max();

		for (int count = start; count < end; ++count)
		{
			if (m_Matrix(i, j).GetAffinePathWay( count) > best2)
			{
				best2 = m_Matrix(i, j).GetAffinePathWay( count);
				position = count;
			}
		}
		position -= start;
	}

	std::vector< std::pair< int, int> >
		converted( alignment.begin(), alignment.end());
	return std::make_pair( best, converted);  // return last matrix element as first of pair
*/
}

/*
std::pair< double, std::vector< std::pair< int, int> > >
NeedlemanWunschAffineGaps:: HMMTraceBack() const
{

        double
                tot_best = -std::numeric_limits< double>::max(),
                best = -std::numeric_limits< double>::max(),
                best2;

        int
                max = std::numeric_limits<int>::max(),
                position = 0,
                group,
                start,
                end;

        std::list<std::pair<int, int> >
                alignment;


        for (unsigned int a = m_Matrix.GetNumberOfRows() - 1; a>0; a++)
        {
            for (unsigned int b = m_Matrix.GetNumberOfColumns() - 1; b>0; b++)
            {

                std::list<std::pair<int, int> > alignment_try;
                DynamicProgrammingMatrixElement element = m_Matrix(a,b);

                for (size_t i = 0; i < 3; i++)
                {
                        if ( element.GetTotalscore(i) > best)
                        {
                                best = element.GetTotalscore(i);
                                group = i;
                        }
                }

                size_t i = a;
                size_t j = b;

                while (group != 4)
                {
                    if (group == 0)
                    {
                        alignment_try.push_front(std::make_pair(i - 1, j - 1));
                        --i;
                        --j;
                    }
                    else if (group == 1)
                    {
                        alignment_try.push_front(std::make_pair(i - 1, max));
                        --i;
                    }
                    else if (group == 2)
                    {
                        alignment_try.push_front(std::make_pair(max, j - 1));
                        --j;
                    }

                    best2 = -std::numeric_limits< double>::max();

                    for (int count = 0; count < 3; ++count)
                    {
                        if (m_Matrix(i, j).GetTotalscore(i) > best2)
                        {
                                best2 = m_Matrix(i, j).GetTotalscore( count);
                                group = i;
                        }
                    }
                }

                if (best > tot_best)
                {
                    tot_best = best;
                    alignment = alignment_try;
                }
            }
        }

        std::vector< std::pair< int, int> > converted( alignment.begin(), alignment.end());
        return std::make_pair( tot_best, converted);  // return last matrix element as first of pair
}
*/

const Matrix< DynamicProgrammingMatrixElement> &
NeedlemanWunschAffineGaps::GetMatrix() const
{
	return m_Matrix;
}


