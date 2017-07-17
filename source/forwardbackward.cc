#include <math.h>
#include <limits>
#include <boost/shared_ptr.hpp>
#include <memory>

double LogSum( double a, double b)
{
    return log( exp( a) + exp( b));
};

double TMX(char X, char Y, int i, const HiddenMarkovModel &H)
{
        return H.GetState(i, std::make_pair(X,Y));
};

std::pair<ShPtr< Vec3>, ShPtr< Vec3> > ForwardBackward( const HiddenMarkovModel &HMM, ShPtr< Vec3> &PP, const Sequence &SEQ)
{
    double infin = std::numeric_limits< double>::max();
    unsigned int L = SEQ.size();
    int M = HMM.GetLastNode();

    // Multiple hit: 3 / L+3   One hit: 2 / L+2 (but we should not allow multiple hits!)
    double lmove = log(2.0/(L + 2.0)); 	// N, C, J move states
    double lloop = log(1.0 - 2.0/(L + 2.0));	// N, C, J loop states

    // The evolving vector is constituted by the Begin state, all HMM states and the End state, so
    //   it is equivalent to ev[i] = ( BVX[i], SMX[i], IMX[i], DMX[i], EVX[i]).
    // The prior (with respect to the L observations) is set as ev[0] = ( 1, 0, 0, ..., 0),
    //   because it is certain that we begin from the Begin node (so here there is no "least informative
    //   prior").
    // The matrices SMX, IMX and DMX have as zeroth node the prior, which is set to -inf (P=0). 
//     ShPtr< Vec3> FW( new Vec3(L+1,M+1)), BW( new Vec3(L+1,M+1)); 
    ShPtr< Vec3> FW( PP), BW( new Vec3(L+1,M+1));

    //Vec3(L+1,M)
    //-----FORWARD-----

    //Setting the prior (i=0 -> you sure are at the Begin node)
    FW->BVX[0] = lmove; //FW->BVX[0] = 0.0;
    for (int k = 0; k <= M; k++)
        FW->SMX[0][k] = FW->IMX[0][k] = FW->DMX[0][k] = -infin;
    FW->EVX[0] = -infin;
    //End
    std::cout << __FUNCTION__ << " BVX " << 0 << "   " << PP->BVX[0] << std::endl;
    std::cout << __FUNCTION__ << " SMX " << 0 << "   " << PP->SMX[0] << std::endl;	
    std::cout << __FUNCTION__ << " IMX " << 0 << "   " << PP->IMX[0] << std::endl;
    std::cout << __FUNCTION__ << " DMX " << 0 << "   " << PP->DMX[0] << std::endl;

    for (unsigned int i = 1; i <= L; i++)  //The observations (the aas) are L, the counter starts at 1
    {
        double sc;
        char t = SEQ[i-1].GetType(); //Real sequence here starts from 1, while it's stored from 0 in SEQ!!
        
        FW->SMX[i][0] = FW->IMX[i][0] = FW->DMX[i][0] = -infin; //Because the Begin node is mute, so these are definitions QUESTION: is Ins[0] wing-retracted with local alns?
        FW->EVX[i] = -infin; //This instead will be updated afterwards

        for (int k = 1; k < M; k++)
        {
            //Substitution
                                 //Subs-Subs
            sc = LogSum( LogSum( FW->SMX[i-1][k-1] + TMX('S','S',k-1,HMM),
                                 //Ins-Subs
                                 FW->IMX[i-1][k-1] + TMX('I','S',k-1,HMM)),
                                 //Del-Subs
                         LogSum( FW->DMX[i-1][k-1] + TMX('D','S',k-1,HMM), 
                                 //Begin-Subs
                                 FW->BVX[i-1] + TMX('B','S',k-1,HMM)));

            FW->SMX[i][k] = sc + HMM.GetSubs(k, t);

            //Insertion
                         //Subs-Ins
            sc = LogSum( FW->SMX[i-1][k] + TMX('S','I',k,HMM),
                         //Ins-Ins
                         FW->IMX[i-1][k] + TMX('I','I',k,HMM));

            FW->IMX[i][k] = sc + HMM.GetIns(k, t);

            //Deletion
                                   //Subs-Del 
            FW->DMX[i][k] = LogSum( FW->SMX[i][k-1] + TMX('S','D',k-1,HMM),
                                   //Del-Del
                                   FW->DMX[i][k-1] + TMX('D','D',k-1,HMM));

            //Update the End state //The profile is always considered as local (TMX(*,E) = 0.0)
            FW->EVX[i] = LogSum( LogSum( FW->SMX[i][k],
                                        FW->DMX[i][k]),
                                FW->EVX[i]);
        }

        sc = LogSum( LogSum( FW->SMX[i-1][M] + TMX('S','S',M-1,HMM),
                             FW->IMX[i-1][M] + TMX('I','S',M-1,HMM)),
                     LogSum( FW->DMX[i-1][M] + TMX('D','S',M-1,HMM),
                             FW->BVX[i-1] + TMX('B','S',M-1,HMM)));
        FW->SMX[i][M] = sc + HMM.GetSubs(M, t);

        FW->IMX[i][M] = -infin;

        FW->DMX[i][M] = LogSum( FW->SMX[i][M-1] + TMX('S','D',M-1,HMM),
                               FW->DMX[i][M-1] + TMX('D','D',M-1,HMM));

        FW->EVX[i] = LogSum( LogSum( FW->SMX[i][M],
                                    FW->DMX[i][M]),
                            FW->EVX[i]);

        FW->BVX[i] = 0.;

	std::cout << __FUNCTION__ << " BVX " << i << "   " << PP->BVX[i] << std::endl;
	std::cout << __FUNCTION__ << " SMX " << i << "   " << PP->SMX[i] << std::endl;	
	std::cout << __FUNCTION__ << " IMX " << i << "   " << PP->IMX[i] << std::endl;
	std::cout << __FUNCTION__ << " DMX " << i << "   " << PP->DMX[i] << std::endl;
    }

    exit(1);

    //-----BACKWARD-----

    BW->BVX[L] = -infin;
    BW->EVX[L] = 0.0;

    BW->SMX[L][M] = BW->DMX[L][M] = BW->EVX[L];
    BW->IMX[L][M] = -infin;

    for (int k = M-1; k >= 1; k--) {
        BW->SMX[L][k] = LogSum( BW->EVX[L] + 0.,
                               BW->DMX[L][k+1] + TMX('S','D',k,HMM));
        BW->DMX[L][k] = LogSum( BW->EVX[L] + 0.,
                               BW->DMX[L][k+1] + TMX('D','D',k,HMM));
        BW->IMX[L][k] = -infin;
    }

    for (int i = L-1; i >= 1; i--)
    {
        char t = SEQ[i-1].GetType();

        BW->BVX[i] = BW->SMX[i+1][1] + TMX('B','S',1,HMM) + HMM.GetSubs(1, t);
        for (int k = 2; k <= M; k++)
            BW->BVX[i] = LogSum( BW->BVX[i], 
                                BW->SMX[i+1][k] + TMX('B','S',k,HMM) + HMM.GetSubs(k, t));

        BW->EVX[i] = 0.0;

        BW->SMX[i][M] = BW->DMX[i][M] = BW->EVX[i];
        BW->IMX[i][M] = -infin; 

        for (int k = M-1; k >= 0; k--)
        { 
            BW->SMX[i][k] = LogSum( LogSum( BW->SMX[i+1][k+1] + TMX('S','S',k,HMM) + HMM.GetSubs(k+1, t),
                                          BW->IMX[i+1][k]   + TMX('S','I',k,HMM) + HMM.GetIns(k, t)),
                                  LogSum( BW->EVX[i],
                                          BW->DMX[i][k+1] + TMX('S','D',k,HMM)));

            BW->IMX[i][k] = LogSum( BW->SMX[i+1][k+1] + TMX('I','S',k,HMM) + HMM.GetSubs(k+1, t),
                                   BW->IMX[i+1][k]   + TMX('I','I',k,HMM) + HMM.GetIns(k, t));

            BW->DMX[i][k] = LogSum( BW->SMX[i+1][k+1] + TMX('D','S',k,HMM) + HMM.GetSubs(k+1, t),
                                   LogSum( BW->DMX[i][k+1] + TMX('D','D',k,HMM),
                                           BW->EVX[i]));
        }
    }

    char t1 = SEQ[1].GetType();
        
    BW->BVX[0] = BW->SMX[1][1] + HMM.GetState(0, std::make_pair('S','S')) + HMM.GetSubs(1, t1);
    for (int k = 2; k <= M; k++)
        BW->BVX[0] = LogSum( BW->BVX[0], 
                            BW->SMX[1][k] + BW->BVX[k-1] + HMM.GetSubs(k, t1));
    BW->EVX[0] = -infin;

    for (int k = M; k >= 1; k--)
        BW->SMX[0][k] = BW->IMX[0][k] = BW->DMX[0][k] = -infin;

    //-----COMBINE-----

    PP->EVX[0] = 0.0;
    PP->BVX[0] = 0.0;

    for (int k = 0; k <= M; k++)
        PP->SMX[0][k] = PP->IMX[0][k] = PP->DMX[0][k] = 0.0;

    for (int i = 1; i <= L; i++)
    {
        double denomS = 0.0, denomI = 0.0;
        PP->SMX[i][0] = PP->IMX[i][0] = PP->DMX[i][0] = 0.0;

        for (int k = 1; k < M; k++)
        {
            PP->SMX[i][k] = FW->SMX[i][k] + BW->SMX[i][k]; 
            denomS += PP->SMX[i][k];
            PP->IMX[i][k] = FW->IMX[i][k] + BW->IMX[i][k]; 
            denomI += PP->IMX[i][k];
            PP->DMX[i][k] = FW->DMX[i][k] + BW->DMX[i][k];
        }

        PP->SMX[i][M] = FW->SMX[i][M] + BW->SMX[i][M];
        denomS += PP->SMX[i][M];
        PP->IMX[i][M] = 0.;
        PP->DMX[i][M] = 0.;

        // order doesn't matter.  note that this whole function is trivially simd parallel
        PP->EVX[i] = 0.0;
        PP->EVX[i] = 0.0;

//	std::cout << __FUNCTION__ << "  " << i << " SMX " << PP->SMX[i] << "  IMX  " << PP->IMX[i] << std::endl;

        denomS = 1.0 / denomS;
	denomI = 1.0 / denomI;
        for (int k = 1; k < M; k++)
        {
            PP->SMX[i][k] *= denomS; 
            PP->IMX[i][k] *= denomI; 
        }
        PP->SMX[i][M]     *= denomS;

//	std::cout << __FUNCTION__ << "  " << i << " SMX " << PP->SMX[i] << "  IMX  " << PP->IMX[i] << std::endl;
    }

//    std::cout << __FUNCTION__ << FW->PMX.size() << std::endl;
	exit(1);
    return std::make_pair(PP, FW);
}
