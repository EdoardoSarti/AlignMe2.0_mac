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
//!  NEEDS CLEANUP
//! Beside the similarity value in a matrix element the DynamicProgrammingMatrixElement
//! contains links to previous elements for determining the gap penalty.
//!
//!
//! @author: Rene Staritzbichler, Kamil Khafizov, Marcus Stamm
//! @date: 18.3.2010
/////////////////////////////////////////////////////////////////////////

#ifndef DYNAMIC_PROGRAMMING_MATRIX_ELEMENT_H
#define DYNAMIC_PROGRAMMING_MATRIX_ELEMENT_H

#include <cstddef>

#include <boost/shared_ptr.hpp>

#include "pearson_correlation.h"


class DynamicProgrammingMatrixElement
{
private:

	double                                  		m_Value;			//!<  the total similarity score at this element of the dynamic programming matrix
	std::pair< size_t, size_t>              		m_PreviousElementIndices;	//!<  indices of the previous element, the one with the best score
	ShPtr< PearsonCorrelation>              		m_Correlation;
	std::vector<double>                     		m_AffinePathWays;		//!<  needed for Needleman Wunsch Affine Gap (2 shells)
        std::vector< double>                    		m_Prior[9];
	std::vector< double>                                    m_null_prior;
	std::vector< double>					m_init_prior;
        double                                  		m_Totalscore[3];	
        int                                     		m_Best_i[3] = {-1, -1, -1};
	const std::vector< std::vector < double> >* 		m_pTMX;
	const std::vector < std::pair< int, int> >* 		m_pMMX[3];
	const std::map< char, std::vector< double> >*		m_pEMI;
	std::vector < double>					m_new_B;
	int							m_nrows;
	int							m_ncols;
	double							d_inf = std::numeric_limits< double>::infinity();

public:

	// Default constructor
	DynamicProgrammingMatrixElement()
	: m_Value( 0.0),
	m_PreviousElementIndices(),
	m_Correlation(),
        m_AffinePathWays( 9, 0.0)
	{};

        // Matrices constructor
        DynamicProgrammingMatrixElement( size_t v_size)
        : m_Value( 0.0),
        m_PreviousElementIndices(),
        m_Correlation(),
        m_AffinePathWays( 9, 0.0)
        {
		for(int i=0; i<9; i++)
		{
			m_Prior[i].resize(v_size);
			m_Prior[i][0] = 0.0;
			for(unsigned int j=1; j<m_Prior[i].size(); j++)
			{
				m_Prior[i][j] = -d_inf;
			}
		}
	}

        DynamicProgrammingMatrixElement( size_t v_size, size_t nrows, size_t ncols)
        : m_Value( 0.0),
        m_PreviousElementIndices(),
        m_Correlation(),
	m_nrows( nrows),
	m_ncols( ncols),
        m_AffinePathWays( 3*nrows*ncols, 0.0)
        {}
	

	// Destructor
	~DynamicProgrammingMatrixElement(){/*std::cout << __FUNCTION__ << std::endl;*/}


	std::vector< double> GetPrior( const int &INDEX)
	{
		return m_Prior[INDEX];
	} 

	// Adds the total similarity VALUE to all path ways
	void AddValue(const double &VALUE)
	{
		for( std::vector<double>::iterator itr = m_AffinePathWays.begin(); itr != m_AffinePathWays.end(); ++itr)
		{
			*itr += VALUE;
		}

		m_Value += VALUE;
	}


	// Adds a prior to the prior list and updates one of the nine pathways
	// m_AffinePathWays accepts indeces from 0 to 8, m_Prior from 0 to 11.
	// Correspondance:
	// m_AffinePathWays[0] -- m_Prior[0]
	// m_AffinePathWays[1] -- m_Prior[1]
	// m_AffinePathWays[2] -- m_Prior[2]
	// m_AffinePathWays[3] -- m_Prior[4]
	// m_AffinePathWays[4] -- m_Prior[5]
	// m_AffinePathWays[5] -- m_Prior[6]
	// m_AffinePathWays[6] -- m_Prior[8]
	// m_AffinePathWays[7] -- m_Prior[9]
	// m_AffinePathWays[8] -- m_Prior[10]
	void AddValue(const int &INDEX, const std::pair< std::vector < double>, double> &PAIR)
	{
		m_Prior[INDEX] = PAIR.first;
		m_AffinePathWays[INDEX] += PAIR.second;
	}


	// Returns all affine pathways
	const std::vector<double> &GetAffinePathWays() const
	{
		return m_AffinePathWays;
	}


	// Returns the pathway relative to the corresponding index
	const double &GetAffinePathWay( const int &INDEX) const
	{
		return m_AffinePathWays[INDEX];
	}


	// Returns the logarithm of the inner product of two vectors
	double logmultv( const std::vector< double> FIRST, const std::vector< double> SECOND) const
	{
		if(FIRST.size() == SECOND.size())
		{
            		double m = 0.;
            		for( unsigned int i=0; i<FIRST.size(); i++)
			{
				m += FIRST[i]*SECOND[i];
			}

			return logg(m);
		}
		else
		{
			std::cout << "logmultv: vectors are not the same length. First: " << FIRST.size() << "   Second: " << SECOND.size() << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	double logsum(double a, double b)
	{
		return log(exp(a) + exp(b));
	}

	double logdiff(double a, double b)
	{
		return log(exp(a) - exp(b));
	}
	


	// Returns the Kullback-Leibler Divergence of the two vectors (it is not symmetric wrt the input variables)
	double KLD( const std::vector< double> FIRST, const std::vector< double> SECOND) const
	{
		if(FIRST.size() == SECOND.size())
		{
            		double m = 0.;
            		for( unsigned int i=0; i<FIRST.size(); i++)
			{
				if (FIRST[i] > 0 && SECOND[i] > 0) // The KLD is defined only for P/=0 and Q/=0
				{
					m += FIRST[i] * logg(FIRST[i]/SECOND[i]);
				}
			}
//			std::cout << "  new prior  " << FIRST << "  profile  " << SECOND << std::endl;

			return m;
		}
		else
		{
//			std::cout << "KLD: vectors are not the same length. First: " << FIRST.size() << "   Second: " << SECOND.size() << std::endl;
			exit(EXIT_FAILURE);
		}
	}


	// Returns the logarithm of a number. When the argument is zero, it returns the maximum negative value.
	double logg( const double num) const
	{
		if(num == 0)
		{
			return -d_inf;
		}
		else
		{
			return log(num);
		}
	}


	// (Original) From the current element, it searches one of the three neighbor elements determined by MAIN_PATH_WAY_INDEX (0..2) for the best total path way (0..8)
	// MAIN_PATH_WAY_INDEX = 0: Diagonal element subspace
	// MAIN_PATH_WAY_INDEX = 1: Horizontal element subspace
	// MAIN_PATH_WAY_INDEX = 2: Vertical element subspace
	const double BestSubPathWayScore( const size_t &MAIN_PATH_WAY_INDEX) const
	{
		double best = m_AffinePathWays[ 3 * MAIN_PATH_WAY_INDEX];
		for( size_t i = 3 * MAIN_PATH_WAY_INDEX + 1; i < 3 * MAIN_PATH_WAY_INDEX + 3; i++)
		{
			if( m_AffinePathWays[i] > best)
			{
				best = m_AffinePathWays[i];                   
			}
		}

		return best;
	}


	void SetPriors( const int DIM)
	{
//		std::cout << __FUNCTION__ << "   Setting Priors" << std::endl;
		for(int i=0; i<9; i++)
		{
			m_Prior[i].resize(DIM);
			m_Prior[i][0] = 0.0;
                        for(unsigned int j=1; j<DIM; j++)
                        {
                                m_Prior[i][j] = -d_inf;
                        }			
		}
		m_null_prior.resize(DIM);
		for(unsigned int j=0; j<DIM; j++)
		{
			m_null_prior[j] = -d_inf;
		}
		m_init_prior = m_null_prior;
		m_init_prior[0] = 0.0;
	}

	void SetMatrices( const std::vector< std::vector < double> >* TMX, const std::vector < std::pair< int, int> >* MMX[3], const std::map< char, std::vector< double> >* EMI, const std::vector< double>* pB)
	{
//		std::cout << __FUNCTION__ << "   Getting TMX size" << std::endl;
		int M = TMX->size();
//		std::cout << __FUNCTION__ << "   Wiring TMX" << std::endl;
		m_pTMX = TMX;
		for(int i=0; i<3; i++)
		{
//			std::cout << __FUNCTION__ << "   Filling pMMX" << std::endl;
			m_pMMX[i] = MMX[i];
		}
//		std::cout << __FUNCTION__ << "   Wiring EMI" << std::endl;
		m_pEMI = EMI;
//		std::cout << __FUNCTION__ << "   Resizing m_new_B   " << M << std::endl;
		m_new_B.resize(M);
		
		for (int i=0; i<M; i++)
		{
//			std::cout << __FUNCTION__ << "   Filling m_new_B" << std::endl;
			m_new_B[i] = 0.0;
			for (int j=0; j<1; j++)
			{
				//This is run on the whole matrix because there is no list for the begin state - maybe implement it here, hardcoded? This should be the only place
				m_new_B[i] += (*pB)[i] + m_pTMX->at(i)[j];
			}
		}
//		std::cout << __FUNCTION__ << "   Finished   " << M << std::endl;
	}	


	double HMMCriterion(std::vector< double> prior, std::vector< double> profile, double indet)
	{
		// This estimator should answer to the question: Is "prior" compatible with the evolution described by "profile"?
		//  Note that the evolution of "prior" is more similar to a forward process, while "profile" went through a
		//  forward-backward. But this is the point: independently of the process (which is not even a forward: it is a
		//  Needlema-Wunsch, if anything), we want the vector which is best mimicking the profile.
		return - KLD(prior, profile);    // The more the two vectors are orthogonal, the less it will count.
	}


	std::pair< std::vector < double>, double> BestSubPathWayScore_NewV( const int &MAIN_PATH_WAY_INDEX, char aa_type, const int &CURRENT_PATH_WAY_INDEX, std::vector < std::pair < int, double > > Vt_el)
	{
		// MAIN CYCLE: go through the elements of the vector (now, there should always be only ONE) and make the HMM states in the prior evolve according to them (they must be the terminal point) and their quality (M and I must have different maps, i.e. for a diagonal CURRENT move, if the CURRENT state is I13, the final state must be I13 (and so the TMX must end in I13), while if it is M15 must be M15. For a vertical/horizontal move, if the CURRENT state is I13, the final state must NOT evolve, but if it is M15, the final state must be M15).
		int M = m_pTMX->size();
		std::vector< double> new_prior, best_new_prior;
		new_prior.resize(M);
		best_new_prior.resize(M);
		double best = -d_inf;
		double emi_modifier = 1.0;

		int start = 3 * MAIN_PATH_WAY_INDEX;
		int end = 3 * MAIN_PATH_WAY_INDEX + 3;
		bool last_row = false;
		bool last_column = false;

		double temp, bestemp;

		if(CURRENT_PATH_WAY_INDEX != 0)   // Because Insertion scores are put to 0 (are equal to null model)
		{
			emi_modifier = 0.0;
		}

		int mev[3];
		if(CURRENT_PATH_WAY_INDEX == 0)
		{
			mev[0] = 3;
			mev[1] = 2;
			mev[2] = 1;
		}
		else if (abs(CURRENT_PATH_WAY_INDEX) == 1)
		{
			mev[0] = 1;
			mev[1] = 0;
			mev[2] = 2;
		}
		else if (abs(CURRENT_PATH_WAY_INDEX) == 2)
		{
			mev[0] = 5;
			mev[1] = 4;
			mev[2] = 3;
		}

		int iev[3];
		if(CURRENT_PATH_WAY_INDEX == 0)
		{
			iev[0] = 1;
			iev[1] = 0;
			iev[2] = 2;
		}
		else
		{
			iev[0] = -1;
			iev[1] = -1;
			iev[2] = -1;
		}

		int *ev;
		int bestq;
		best_new_prior = m_null_prior;
		std::cout << "VECTOR SIZE " << Vt_el.size() << std::endl;

		for( int i = start; i < end; ++i)
		{
			double contribution = 0.0;
			bool once_here = false;
			new_prior = m_null_prior;

			if ((CURRENT_PATH_WAY_INDEX > -1) || (CURRENT_PATH_WAY_INDEX == -1 && i != start+1) || (CURRENT_PATH_WAY_INDEX == -2 && i != end-1))
			{
//				for(std::vector < std::pair< int, int> >::const_iterator it = (*m_pMMX[abs(CURRENT_PATH_WAY_INDEX)]).begin(); it != (*m_pMMX[abs(CURRENT_PATH_WAY_INDEX)]).end(); ++it)
				for (std::vector < std::pair < int, double > >::iterator iv = Vt_el.begin(); iv != Vt_el.end(); ++iv)
				{
					std::cout << "endpoint " << iv->first << std::endl;
					if (iv->first < 4)
					{
						if ((iv->first == 1 && CURRENT_PATH_WAY_INDEX == 0) || (iv->first == 2 && CURRENT_PATH_WAY_INDEX == 1))
						{
							std::cout << "staring " << iv->first - 1 << " TMX " << m_pTMX->at(iv->first - 1)[iv->first] << std::endl;
							new_prior[iv->first] = m_Prior[i][iv->first - 1] + m_pTMX->at(iv->first - 1)[iv->first];
						}
						else if (iv->first == 3 && CURRENT_PATH_WAY_INDEX == 2)
						{
							std::cout << "no change" << std::endl;
							new_prior = m_Prior[i];
						}
						continue;
					}

					if (iv->first % 3 == 1)  // CURRENT HMM STATE == M
					{
						ev = &mev[0];
						new_prior[iv->first] = m_Prior[i][iv->first - ev[MAIN_PATH_WAY_INDEX]] + m_pTMX->at(iv->first - ev[MAIN_PATH_WAY_INDEX])[iv->first];
						std::cout << new_prior[iv->first] << " " << m_Prior[i][iv->first - ev[MAIN_PATH_WAY_INDEX - 3]] << std::endl;

						bestq = -1;
						if (CURRENT_PATH_WAY_INDEX != 1 && MAIN_PATH_WAY_INDEX != 1)
						{
							bestemp = new_prior[iv->first];
							for (int q = iv->first - ev[MAIN_PATH_WAY_INDEX] -3; q > 0; q -= 3)
							{
								temp = m_Prior[i][q] + m_pTMX->at(q)[iv->first];
								if (temp > bestemp)
								{
									bestemp = temp;
									bestq = q;
								}
							}
							temp = m_Prior[i][0] + m_pTMX->at(0)[iv->first];
							if (temp > bestemp)
							{
								bestemp = temp;
								bestq = 0;
							}
							new_prior[iv->first] = bestemp;
						}
						std::cout << new_prior[iv->first] << " BEST q: " << bestq << " VALUE " << bestemp << std::endl;
					}
					else if (iv->first % 3 == 2)   // CURRENT HMM STATE == I
					{
						ev = &iev[0];
						new_prior[iv->first] = m_Prior[i][iv->first - ev[MAIN_PATH_WAY_INDEX]] + m_pTMX->at(iv->first - ev[MAIN_PATH_WAY_INDEX])[iv->first];
					}

//					std::cout << m_pTMX->size() << std::endl;
//					std::cout << iv->first << " " << iv->first - ev[MAIN_PATH_WAY_INDEX] << " " << m_pTMX->at(iv->first - ev[MAIN_PATH_WAY_INDEX])[iv->first] << std::endl;
//					new_prior[k+c[MAIN_PATH_WAY_INDEX]] = logsum(new_prior[k+c[MAIN_PATH_WAY_INDEX]], m_Prior[i][k] + m_pTMX->at(k)[k+c[MAIN_PATH_WAY_INDEX]]);
					
//					new_prior[iv->first] = m_Prior[i][iv->first - ev[MAIN_PATH_WAY_INDEX]] + m_pTMX->at(iv->first - ev[MAIN_PATH_WAY_INDEX])[iv->first];

					// VERY QUESTIONABLE LINES: WHY A SUM??___________
//					if (iv->first % 3 == 1 && iv->first < M-1)
//						new_prior[iv->first] = logsum(new_prior[iv->first], m_Prior[i][0] + m_pTMX->at(0)[iv->first]);
					// ________________________________________________


					contribution += iv->second * exp(new_prior[iv->first])/(1.0 + exp(new_prior[iv->first]));
					std::cout << "Contribution " << contribution << std::endl;
				}
			}
			else
			{
				new_prior = m_Prior[i];
			}
			std::cout << "Contribution from PathWay" << m_AffinePathWays[i] << std::endl;
			contribution += m_AffinePathWays[i];
			if( contribution > best)
			{
				best = contribution;
				m_Best_i[MAIN_PATH_WAY_INDEX] = i;
	 			best_new_prior = new_prior;
			}
		}
		for (std::vector < std::pair < int, double > >::iterator iv = Vt_el.begin(); iv != Vt_el.end(); ++iv)
		{
			if (iv->first % 3 == 1 && iv->first < M-1)
				best_new_prior[iv->first] += emi_modifier*m_pEMI->at(aa_type)[iv->first]; //Emissions are put AFTER decision! (See Viterbi algorithm)
		}
		return std::make_pair(best_new_prior, best);
	}


	std::pair< std::vector < double>, double> BestSubPathWayScore_PureV( const int &MAIN_PATH_WAY_INDEX, char aa_type, const int &CURRENT_PATH_WAY_INDEX)
	{
		int M = m_pTMX->size();
		std::vector< double> new_prior, best_new_prior;
		new_prior.resize(M);
		best_new_prior.resize(M);
		double best = -d_inf;
		double emi_modifier = 1.0;

		int start = 3 * MAIN_PATH_WAY_INDEX;
		int end = 3 * MAIN_PATH_WAY_INDEX + 3;
		bool last_row = false;
		bool last_column = false;


		if(CURRENT_PATH_WAY_INDEX != 0)   // Because Insertion scores are put to 0 (are equal to null model)
		{
			emi_modifier = 0.0;
		}

		int c[3];
		if(CURRENT_PATH_WAY_INDEX == 0)
		{
			c[0] = 3;
			c[1] = 2;
			c[2] = 1;
		}
		else if (abs(CURRENT_PATH_WAY_INDEX) == 1)
		{
			c[0] = 1;
			c[1] = 0;
			c[2] = 2;
		}
		else if (abs(CURRENT_PATH_WAY_INDEX) == 2)
		{
			c[0] = 5;
			c[1] = 4;
			c[2] = 3;
		}

		best_new_prior = m_null_prior;
		for( int i = start; i < end; ++i)
		{
			double contribution = 0.0;
			bool once_here = false;
			new_prior = m_null_prior;

			if ((CURRENT_PATH_WAY_INDEX > -1) || (CURRENT_PATH_WAY_INDEX == -1 && i != start+1) || (CURRENT_PATH_WAY_INDEX == -2 && i != end-1))
			{
//				for(std::vector < std::pair< int, int> >::const_iterator it = (*m_pMMX[abs(CURRENT_PATH_WAY_INDEX)]).begin(); it != (*m_pMMX[abs(CURRENT_PATH_WAY_INDEX)]).end(); ++it)
				for(int k = 1+MAIN_PATH_WAY_INDEX; k < M-4; k+=3)
				{
//					new_prior[it->second] = logsum(new_prior[it->second], m_Prior[i][it->first] + m_pTMX->at(it->first)[it->second]);
//					std::cout << k << " " << k+c[MAIN_PATH_WAY_INDEX] << " " << m_pTMX->at(k)[k+c[MAIN_PATH_WAY_INDEX]] << std::endl;
//					new_prior[k+c[MAIN_PATH_WAY_INDEX]] = logsum(new_prior[k+c[MAIN_PATH_WAY_INDEX]], m_Prior[i][k] + m_pTMX->at(k)[k+c[MAIN_PATH_WAY_INDEX]]);
					new_prior[k+c[MAIN_PATH_WAY_INDEX]] = m_Prior[i][k] + m_pTMX->at(k)[k+c[MAIN_PATH_WAY_INDEX]];
				}
				if(CURRENT_PATH_WAY_INDEX == 0)
				{
					for(int kk = 1; kk < M-1; kk += 3)
					{
//						std::cout << 0 << " " << kk << " " << m_pTMX->at(0)[kk] << std::endl;
						new_prior[kk] = logsum(new_prior[kk], m_Prior[i][0] + m_pTMX->at(0)[kk]);
					}
				}
				for (int j = 1+abs(CURRENT_PATH_WAY_INDEX); j < M-1; j+=3)
				{
					//double primo = (new_prior[j] == -d_inf) ? 0.0 : new_prior[j];
					//double secondo = (m_Prior[i][j] == -d_inf) ? 0.0 : m_Prior[i][j];
					//std::cout << CURRENT_PATH_WAY_INDEX << " " << new_prior[j] << " " << m_Prior[i][j] << std::endl;
//					contribution += (new_prior[j] > m_Prior[i][j]) ? 1.0 : 0.0;
//					contribution += (new_prior[j] > m_Prior[i][j]) ? exp(new_prior[j])/(1.0 + exp(new_prior[j])) : 0.0;
					contribution += exp(new_prior[j])/(1.0 + exp(new_prior[j]));// - exp(m_Prior[i][j])/(1.0 + exp(m_Prior[i][j]));
				}
//				new_prior = m_Prior[i];
//				std::cout << CURRENT_PATH_WAY_INDEX << " " << i << " " << new_prior << " " << m_Prior[i] << std::endl;
			}
			else
			{
				new_prior = m_Prior[i];
			}
//			if (m_AffinePathWays[i] == -d_inf)
				contribution += m_AffinePathWays[i];
			if( contribution > best)
//			if( i == 0)
			{
				best = contribution;
				m_Best_i[MAIN_PATH_WAY_INDEX] = i;
	 			best_new_prior = new_prior;
			}
		}
		for (int j = 1+abs(CURRENT_PATH_WAY_INDEX); j < M; j+=3)
		{
			best_new_prior[j] += emi_modifier*m_pEMI->at(aa_type)[j]; //Emissions are put AFTER decision! (See Viterbi algorithm)
		}
		return std::make_pair(best_new_prior, best);
	}


	std::pair< std::vector < double>, double> BestSubPathWayScore_PV( const int &MAIN_PATH_WAY_INDEX, char aa_type, const int &CURRENT_PATH_WAY_INDEX)
	{
		int M = m_pTMX->size();
		std::vector< double> new_prior, best_new_prior;
		new_prior.resize(M);
		best_new_prior.resize(M);
		double best = -d_inf;
		double emi_modifier = 1.0;
		//m_Best_i[MAIN_PATH_WAY_INDEX] = -1; // No choice is given because no score greater than -inf was found. This implies we are on a starting border of the NW matrix.

		int start = 3 * MAIN_PATH_WAY_INDEX;
		int end = 3 * MAIN_PATH_WAY_INDEX + 3;
		bool last_row = false;
		bool last_column = false;


		if(CURRENT_PATH_WAY_INDEX != 0)   // Because Insertion scores are put to 0 (are equal to null model)
		{
			emi_modifier = 0.0;
		}

		best_new_prior = m_null_prior;
		for( int i = start; i < end; ++i)
		{
			double contribution = 0.0;
			bool once_here = false;
			new_prior = m_null_prior;

//			for(std::vector < std::pair< int, int> >::const_iterator it = (*m_pMMX[i - 3*MAIN_PATH_WAY_INDEX]).begin(); it != (*m_pMMX[i - 3*MAIN_PATH_WAY_INDEX]).end(); it++)    // This is wrong: the priors you're taking are already evolved, you have to make them do the last step. So, when you ask the question "which one of these three paths is better?", you THEN take the winner and evolve it by the rules. But it is smarter doing the comparison among the three priors "as if" they were already evolved (example: let's say we are in the H branch and the winner before evolution seems V. Then you evolve it and it goes to null. In hindsight, would you have chosen it?). So we evolve the 3 priors by the same matrix, then we compare them. 
			if ((CURRENT_PATH_WAY_INDEX > -1) || (CURRENT_PATH_WAY_INDEX == -1 && i != start+1) || (CURRENT_PATH_WAY_INDEX == -2 && i != end-1))
			{
				for(std::vector < std::pair< int, int> >::const_iterator it = (*m_pMMX[abs(CURRENT_PATH_WAY_INDEX)]).begin(); it != (*m_pMMX[abs(CURRENT_PATH_WAY_INDEX)]).end(); ++it)
				{
//					std::cout << __FUNCTION__ << " " << new_prior[it->first] << std::endl;
					new_prior[it->second] = logsum(new_prior[it->second], m_Prior[i][it->first] + m_pTMX->at(it->first)[it->second]);
//					std::cout << __FUNCTION__ << " " << i << "  " << it->first << " " << it->second << " " << m_Prior[i][it->first] << " " << m_pTMX->at(it->first)[it->second] << " " << new_prior[it->second] << std::endl;
				}
				contribution += exp(new_prior[0])/(1.0 + exp(new_prior[0]));
				for (int j = 1+abs(CURRENT_PATH_WAY_INDEX); j < M; j+=3)
				{
//					std::cout << __FUNCTION__ << " PRIOR " << j << " " << m_Prior[i][j] << std::endl;	
//					std::cout << __FUNCTION__ << " UPDATED PRIOR " << j << " " << new_prior[j] << std::endl;
					contribution += exp(new_prior[j])/(1.0 + exp(new_prior[j])); //Logistic function. WARNING: THERE IS NO LOG-ODD PRIOR YET!
				}
				contribution += exp(new_prior[M-1])/(1.0 + exp(new_prior[M-1]));
			}
			else
			{
				new_prior = m_Prior[i];
			}
//			std::cout << CURRENT_PATH_WAY_INDEX << " " << i << " " << new_prior << " " << m_Prior[i] << std::endl;
//			std::cout << __FUNCTION__ << i << "  " << m_Prior[i].size() << std::endl;
//			std::cout << __FUNCTION__ << " CONTRIBUTION BEFORE AFFINEPATHS " << contribution << std::endl;
			if (m_AffinePathWays[i] == -d_inf)
				contribution += m_AffinePathWays[i];
//			std::cout << __FUNCTION__ << " TOTAL SCORE (IF ONLY HMM MATCHES MUST BE NON-NEGATIVE) " << i << "  " << m_AffinePathWays[i] << " " << contribution << std::endl;

			if( contribution > best)
			{
//				std::cout << __FUNCTION__ << i << " contribution " << contribution << " best " << best << " " << m_AffinePathWays[i] << " " << new_prior << std::endl;
				best = contribution;
				m_Best_i[MAIN_PATH_WAY_INDEX] = i;
	 			best_new_prior = new_prior;
			}
			//std::cout << __FUNCTION__ << i << "  complete" << std::endl;
		}
		for (int j = 1; j < M; j+=3)
		{
			best_new_prior[j] += emi_modifier*m_pEMI->at(aa_type)[j]; //Emissions are put AFTER decision! (See Viterbi algorithm)
		}

		return std::make_pair(best_new_prior, best);
	}


	// Returns the total similarity value of this element
	const double &GetValue() const
	{
		return m_Value;
	}


	const double GetTotalScore( const int &GROUP) const
	{
		return m_Totalscore[GROUP];
	}


	const int GetBestPath( const int &GROUP) const
	{
		return m_Best_i[GROUP];
	}


	// Sets the total similarity value of this element
	void SetValue( const double &VALUE)
	{
		m_Value = VALUE;
	}


	// Returns the indices of the previous matrix element (where it comes from)
	const std::pair< size_t, size_t> &GetIndicesOfPreviousElement() const
	{
		return m_PreviousElementIndices;
	}


	std::ostream &Write( std::ostream &STREAM) const
	{
		STREAM  << m_Value;
		return STREAM;
	}


	const ShPtr< PearsonCorrelation> &GetCorrelation() const
	{
		return m_Correlation;
	}


	void SetCorrelation( const ShPtr< PearsonCorrelation> &CORRELATION)
	{
		m_Correlation = CORRELATION;
	}


	double operator *= ( const double &VALUE)
	{
		return (m_Value *= VALUE);
	}
};


inline
std::ostream &operator << ( std::ostream &STREAM, const DynamicProgrammingMatrixElement &ELEMENT)
{
	return ELEMENT.Write( STREAM);
}


inline
double operator *( const double &VALUE, const DynamicProgrammingMatrixElement &ELEMENT)
{
	return VALUE * ELEMENT.GetValue();
}


inline
double operator *( const DynamicProgrammingMatrixElement &ELEMENT, const double &VALUE)
{
	return VALUE * ELEMENT.GetValue();
}


#endif
