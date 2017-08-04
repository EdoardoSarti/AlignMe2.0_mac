#include <math.h>
#include <limits>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/assign.hpp>

//#include "generic_decoding.c"
#include "p7_config.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "hmmer.h"
//#include "modelconfig.c"

#define STYLES     "--fs,--sw,--ls,--s"   /* Exclusive choice for alignment mode     */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "--fs",      eslARG_NONE,"default",NULL, NULL, STYLES,  NULL, NULL, "multihit local alignment",                         0 },
  { "--sw",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit local alignment",                           0 },
  { "--ls",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "multihit glocal alignment",                        0 },
  { "--s",       eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit glocal alignment",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of posterior decoding, generic implementation";

static void dump_matrix_csv(FILE *fp, P7_GMX *pp, int istart, int iend, int kstart, int kend);

//c++11
//std::map<char, int> ctoi = {{'a', 0}, {'c', 1}, {'d', 2}, {'e', 3}, {'f', 4}, {'g', 5}, {'h', 6}, {'i', 7}, {'k', 8}, {'l', 9}, {'m', 10}, {'n', 11}, {'p', 12}, {'q', 13}, {'r', 14}, {'s', 15}, {'t', 16}, {'v', 17}, {'w', 18}, {'y', 19}, {'A', 0}, {'C', 1}, {'D', 2}, {'E', 3}, {'F', 4}, {'G', 5}, {'H', 6}, {'I', 7}, {'K', 8}, {'L', 9}, {'M', 10}, {'N', 11}, {'P', 12}, {'Q', 13}, {'R', 14}, {'S', 15}, {'T', 16}, {'V', 17}, {'W', 18}, {'Y', 19}};
static std::map<char, int> ctoi = boost::assign::map_list_of ('a', 0) ('c', 1) ('d', 2) ('e', 3) ('f', 4) ('g', 5) ('h', 6) ('i', 7) ('k', 8) ('l', 9) ('m', 10) ('n', 11) ('p', 12) ('q', 13) ('r', 14) ('s', 15) ('t', 16) ('v', 17) ('w', 18) ('y', 19) ('A', 0) ('C', 1) ('D', 2) ('E', 3) ('F', 4) ('G', 5) ('H', 6) ('I', 7) ('K', 8) ('L', 9) ('M', 10) ('N', 11) ('P', 12) ('Q', 13) ('R', 14) ('S', 15) ('T', 16) ('V', 17) ('W', 18) ('Y', 19);
static std::map<int, char> itoc = boost::assign::map_list_of (0, 'A') (1, 'C') (2, 'D') (3, 'E') (4, 'F') (5, 'G') (6, 'H') (7, 'I') (8, 'K') (9, 'L') (10, 'M') (11, 'N') (12, 'P') (13, 'Q') (14, 'R') (15, 'S') (16, 'T') (17, 'V') (18, 'W') (19, 'Y');

//extern "C" P7_BG *p7_bg_Create(const ESL_ALPHABET *abc);

struct Vec3
{
    std::vector< std::vector< double> > IMX;
    std::vector< std::vector< double> > SMX;
    std::vector< std::vector< double> > DMX;
    std::vector< double> BVX;
    std::vector< double> EVX;
    int seq_length, model_length, active;
    std::vector< std::vector < double> > TMX; //profile matrix
    std::vector< std::vector < double> > PMX;
    std::map< char, std::vector< double > > EMI;
    std::vector < std::pair< int, int> > MMX[3]; //mask lists. 0:Subs, 1:Ins, 2:Del    
    std::vector < std::vector < std::pair < int, double > > > Vtrace;

    Vec3()
        : seq_length(0),
        model_length(0),
        active(0)
        {};

	// L := Number of columns in the NW = Sequence length + 1
	// M := Number of NODES = LastNode index + 1 (because there also is Node 0)
    Vec3(int L, int M)
    {
        SMX.resize(L);
        for (int i=0; i<L; i++)
            SMX[i].resize(M);
        IMX.resize(L);
        for (int i=0; i<L; i++)
            IMX[i].resize(M);
        DMX.resize(L);
        for (int i=0; i<L; i++)
            DMX[i].resize(M);
        BVX.resize(L);
        EVX.resize(L);

        seq_length = L;
        model_length = M;
        active = 1;
    };

    ~Vec3()
    {
        for (unsigned int i=0; i<SMX.size(); i++)
            SMX[i].clear();
        SMX.clear();

        for (int i=0; i<IMX.size(); i++)
            IMX[i].clear();
        IMX.clear();

        for (int i=0; i<DMX.size(); i++)
            DMX[i].clear();
        DMX.clear();

        BVX.clear();
        EVX.clear();
    }

    void copy( Vec3 C)
    {
        unsigned int L = C.seq_length;
        unsigned int M = C.model_length;

        SMX.resize(L);
        for (unsigned int i=0; i<L; i++)
        {
            SMX[i].resize(M);
            SMX[i] = C.SMX[i];
        }
        IMX.resize(L);
        for (unsigned int i=0; i<L; i++)
        {
            IMX[i].resize(M);
            IMX[i] = C.IMX[i];
        }
        DMX.resize(L);
        for (unsigned int i=0; i<L; i++)
        {
            DMX[i].resize(M);
            DMX[i] = C.DMX[i];
        }
        BVX.resize(L);
        BVX = C.BVX;
        EVX.resize(L);
        EVX = C.EVX;

        seq_length = L;
        model_length = M;
        active = 1;
    }

    void fillPMX()
    {
        unsigned int L = seq_length;
        unsigned int M = model_length;

        PMX.resize(L);
        for (unsigned int i=0; i<L; i++)
        {
            PMX[i].resize(3*M+2);
            PMX[i][0] = exp(-BVX[i]);
            for(unsigned int j=0; j<M; j++)
            {
                PMX[i][3*j+1] = exp(-SMX[i][j]);
                PMX[i][3*j+2] = exp(-IMX[i][j]);
                PMX[i][3*j+3] = exp(-DMX[i][j]);
            }
            PMX[i][3*M+1] = exp(-EVX[i]);
        }

	std::cout << "Number of HMM states: " << PMX[0].size() << std::endl;
    }
};



class HiddenMarkovModel
{
protected:
    std::vector< std::map< char, double> > Substitution;
    std::vector< std::map< char, double> > Insertion;
    std::map< char, std::vector< double > > Emissions;
    std::vector< std::map< std::pair< char, char>, double> > State; 
    std::vector< std::pair< char, char> > StateName; // SS, SI, SD, IS, II, DS, DD
    int LastNode;
    std::vector< char> aa_types;
    std::vector< double> occ;
    std::vector< std::vector < double> > TMX; //transition matrix
    std::vector < std::pair< int, int > > MMX[3]; //mask matrices. Only diagonal is saved. 0:Begin, 1:Ins, 2:Subs, 3:Del
    double d_inf = std::numeric_limits< double>::infinity();

public:
    HiddenMarkovModel()
    : Substitution(),
    Insertion(),
    Emissions(),
    State(),
    LastNode(0),
    aa_types()
    {
        StateName.resize(8);
        StateName[0] = std::make_pair('S','S');
        StateName[1] = std::make_pair('S','I');
        StateName[2] = std::make_pair('S','D');
        StateName[3] = std::make_pair('I','S');
        StateName[4] = std::make_pair('I','I');
        StateName[5] = std::make_pair('D','S');
        StateName[6] = std::make_pair('D','D');
        StateName[7] = std::make_pair('B','S');
    };

    ~HiddenMarkovModel()
	{
		for (unsigned int i=0; i<TMX.size(); i++)
		{
			TMX[i].clear();
		}
		TMX.clear();

		for (unsigned int i=0; i<3; i++)
		{
			MMX[i].clear();
		}

		StateName.clear();

		occ.clear();

		aa_types.clear();

		std::cout << State.size() << "\n";
                for (unsigned int i=0; i<State.size(); i++)
                {
			std::cout << State[i].at(std::make_pair('S','S')) << "\n";
			State[i].clear();
		}
		State.clear();

                for (unsigned int i=0; i<Substitution.size(); i++)
                {
                        Substitution[i].clear();
                }
                Substitution.clear();

                for (unsigned int i=0; i<Insertion.size(); i++)
                {
                        Insertion[i].clear();
                }
                Insertion.clear();

                for (unsigned int i=0; i<Emissions.size(); i++)
                {
                        Emissions[i].clear();
                }
                Emissions.clear();
	};

    void SetSubs( const int &NODE, const char &AA_TYPE, const double &VALUE)
    {
        Substitution.resize(LastNode+1);
        Substitution[NODE][AA_TYPE] = VALUE;
    }
 
    void SetIns( const int &NODE, const char &AA_TYPE, const double &VALUE)
    {
        Insertion.resize(LastNode+1);
        Insertion[NODE][AA_TYPE] = VALUE;
    }


    // For each aa type sets a vector (as long as all states in the HMM) with the emission probabilities of each state relative to the aa type (begin, end and delete states probabilities are indeed 0)

    void SetEmissions()
    {
        std::vector< char> aa_vec;
	for (size_t i = 0; i < LastNode+1; i++)
        {
            for(std::map<char, double>::iterator it = Substitution[i].begin(); it != Substitution[i].end(); ++it)
            {
                aa_vec.push_back(it->first);
            }
        }
	for (size_t i = 0; i < aa_vec.size(); i++)
        {
            char aa_type = aa_vec[i];
            Emissions[aa_type].resize(LastNode*3+2);
            Emissions[aa_type][0] = 0.0;
            for(unsigned int j=0; j<LastNode; j++)
            {
                Emissions[aa_type][3*j+1] = Substitution[j][aa_type];
                Emissions[aa_type][3*j+2] = Insertion[j][aa_type];
                Emissions[aa_type][3*j+3] = 0.0;  // Deletions
            }
            Emissions[aa_type][3*LastNode+1] = 0.0;
        }

    }

    void SetState( const int &NODE, const std::pair< char, char> &STATES, const double &VALUE)
    {
        State.resize(LastNode+1);
        State[NODE][STATES] = VALUE;
    }

    const int GetLastNode() const
    {
        return LastNode;
    }

    const double &GetSubs( const int &NODE, const char &AA_TYPE) const
    {
        return Substitution[NODE].at(AA_TYPE);
    }

    const double &GetIns( const int &NODE, const char &AA_TYPE) const
    {
        return Insertion[NODE].at(AA_TYPE);
    }

    const double &GetState( const int &NODE, const std::pair< char, char> &STATES) const
    {
        return State[NODE].at(STATES);
    }

    const std::pair< char, char> &GetStateName( const int &STATENUM) const
    {
        return StateName[STATENUM];
    }

    void LastNodeInc()
    {
        LastNode++;
    }

    double GetProb( const int &NODE, const std::pair< char, char> &STATES, const char &AA_TYPE) const
    {
        if(STATES.second == 'S') 
        {
            return State[NODE].at(STATES)*Substitution[NODE].at(AA_TYPE);
        }

        else if(STATES.second == 'I')
        {
            return State[NODE].at(STATES)*Insertion[NODE].at(AA_TYPE);
        }

        else if(STATES.second == 'D')
        {
            return State[NODE].at(STATES);
        }

        else
        {
            exit(1);
        }   
    }

    void SetAA( char AA)
    {
        aa_types.push_back(AA);
    }

    char GetAA( const int ind) const
    {
        return aa_types[ind];
    }

    void CalculateBeginStateTransitionProbabilities()
    {
        occ.resize(LastNode+1);
        double Z = 0.;

        //First, calculate the occupancy of each Substution state
        occ[0] = 0.;
        occ[1] = exp(-State[0].at(std::make_pair('S','I'))) + exp(-State[0].at(std::make_pair('S','S')));  //The zeroth node does not have a Del state
        //The occupancy of the present Subs state is the occupancy of the previous Subs state times the Trans Prob of it going to the present Subs state
        //  or to the previous Ins state (that will surely go sooner or later in the present Subs state), plus the occupancy of the previous Del state
        //  times the Trans Prob of it going to the present Subs state.
        for (int k = 2; k <= LastNode; k++)
            occ[k] = occ[k-1] * (exp(-State[k-1].at(std::make_pair('S','S'))) + exp(-State[k-1].at(std::make_pair('S','I')))) +
                     (1.0-occ[k-1]) * exp(-State[k-1].at(std::make_pair('D','S')));

        //Calculate a normalization for the occupancy, but I don't understand how it works...
        for (int k = 1; k <= LastNode; k++)
            Z += occ[k] * (double) (LastNode-k+1);

        //The Begin state Transition Probability to the k-th Subs state is just the normalized occupancy probability for node k.
        //  (i.e. how easy is to happen in that (Substitution) state)
        State[0][std::make_pair('B','S')] = 0.;  //I ADDED THIS LINE...
        for (int k = 1; k <= LastNode; k++)
            State[k][std::make_pair('B','S')] = -log(occ[k] / Z);
        return;
    }

    void CalculateEmissionScores()
    {
        //Substitution score
        for (int k = 1; k <= LastNode; k++) 
        {
            for (unsigned int i = 0; i < aa_types.size(); i++)
                Substitution[k][aa_types[i]] = Substitution[k][aa_types[i]] - Substitution[0][aa_types[i]]; //Subs state zero is the null model!
        }
        
        //Insertion score --- NOT USED IN REAL HMMer!!! THEY ACTUALLY SET ALL TO ZERO, AS IF INS[k] = SUB[0] !!!
        for (int k = 1; k < LastNode; k++)  //Ins state of last node missing
        {
            for (unsigned int i = 0; i < aa_types.size(); i++)
                Insertion[k][aa_types[i]] = Insertion[k][aa_types[i]] - Substitution[0][aa_types[i]]; //Subs state zero is the null model!
        }

        for (unsigned int i = 0; i < aa_types.size(); i++)
            Insertion[LastNode][aa_types[i]] = -d_inf;

        return;
    }

    void FillTMX( boost::shared_ptr< Vec3> &HMMP)
    {
        int M = LastNode;

        HMMP->TMX.resize(3*(M+1)+2); // 1 Begin state, M Subs states, M Ins states, M Del states, 1 End state
        for (int i=0; i<3*(M+1)+2; i++)
        {
            HMMP->TMX[i].resize(3*(M+1)+2);
            for (int j=0; j<3*(M+1)+2; j++)
                HMMP->TMX[i][j] = 0.0;
        }
            
	std::cout << "TMX Resized: " << HMMP->TMX.size() << std::endl;

        for (int j=0; j<M; j++)
            HMMP->TMX[0][j*3+1] = exp(-State[j][std::make_pair('B','S')]);
        for (int i=0; i<M-1; i++)
        {
            HMMP->TMX[i*3+1][i*3+4] = exp(-State[i][std::make_pair('S','S')]);
            HMMP->TMX[i*3+1][i*3+2] = exp(-State[i][std::make_pair('S','I')]);
            HMMP->TMX[i*3+1][i*3+6] = exp(-State[i][std::make_pair('S','D')]);
            HMMP->TMX[i*3+2][i*3+2] = exp(-State[i][std::make_pair('I','I')]);
            HMMP->TMX[i*3+2][i*3+4] = exp(-State[i][std::make_pair('I','S')]);
            HMMP->TMX[i*3+3][i*3+4] = exp(-State[i][std::make_pair('D','S')]);
            HMMP->TMX[i*3+3][i*3+6] = exp(-State[i][std::make_pair('D','D')]);
        }
        HMMP->TMX[3*M-2][3*M+1] = exp(-State[M][std::make_pair('S','S')]);
        HMMP->TMX[3*M-1][3*M+1] = exp(-State[M][std::make_pair('I','S')]);
        HMMP->TMX[3*M][3*M+1] = exp(-State[M][std::make_pair('D','S')]);
    }

    void FillMasks( boost::shared_ptr< Vec3> HMMP)
    {
        int M = LastNode;

	/* // OLD VERSION; WRONG: WE WANT PROBABILITIES TO GO INTO, NOT TO EXIT FROM
	for (int i=0; i<M-1; i++)
        {
            //Substitution
            HMMP->MMX[0].push_back(std::make_pair(0, i*3+1));    
            HMMP->MMX[0].push_back(std::make_pair(i*3+1, i*3+2));
            HMMP->MMX[0].push_back(std::make_pair(i*3+1, i*3+4));
            HMMP->MMX[0].push_back(std::make_pair(i*3+1, i*3+6));
            //Insertion (In HMMer model the (3M-1)th state is empty)
            HMMP->MMX[1].push_back(std::make_pair(i*3+2, i*3+2));
            HMMP->MMX[1].push_back(std::make_pair(i*3+2, i*3+4));
            //Deletion (In HMMer model the 3rd state is empty)
            HMMP->MMX[2].push_back(std::make_pair(i*3+3, i*3+4));
            HMMP->MMX[2].push_back(std::make_pair(i*3+3, i*3+6));
        }
        HMMP->MMX[0].push_back(std::make_pair(0, (M-1)*3+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M-2, 3*M+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M-1, 3*M+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M, 3*M+1));
	*/

	// NEWER VERSION
	for (int i=0; i<M-1; i++)
        {
            //Substitution
            HMMP->MMX[0].push_back(std::make_pair(0, i*3+1));
            HMMP->MMX[0].push_back(std::make_pair(i*3+1, i*3+4));
            HMMP->MMX[0].push_back(std::make_pair(i*3+2, i*3+4));
            HMMP->MMX[0].push_back(std::make_pair(i*3+3, i*3+4));
            //Insertion (In HMMer model the (3M-1)th state is empty)
            HMMP->MMX[1].push_back(std::make_pair(i*3+1, i*3+2));
            HMMP->MMX[1].push_back(std::make_pair(i*3+2, i*3+2));
            //Deletion (In HMMer model the 3rd state is empty)
            HMMP->MMX[2].push_back(std::make_pair(i*3+1, i*3+6));
            HMMP->MMX[2].push_back(std::make_pair(i*3+3, i*3+6));
        }
        HMMP->MMX[0].push_back(std::make_pair(0, (M-1)*3+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M-2, 3*M+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M-1, 3*M+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M, 3*M+1));

	/*
	// M9 VERSION TO BE IMPLEMENTED
	for (int i=0; i<M-1; i++)
        {
            //Substitution
            HMMP->MMX[0].push_back(std::make_pair(0, i*3+1));
            HMMP->MMX[0].push_back(std::make_pair(i*3+1, i*3+4));
            HMMP->MMX[0].push_back(std::make_pair(i*3+2, i*3+4));
            HMMP->MMX[0].push_back(std::make_pair(i*3+3, i*3+4));
            //Insertion (In HMMer model the (3M-1)th state is empty)
            HMMP->MMX[1].push_back(std::make_pair(i*3+1, i*3+2));
            HMMP->MMX[1].push_back(std::make_pair(i*3+2, i*3+2));
            //Deletion (In HMMer model the 3rd state is empty)
            HMMP->MMX[2].push_back(std::make_pair(i*3+1, i*3+6));
            HMMP->MMX[2].push_back(std::make_pair(i*3+3, i*3+6));
        }
        HMMP->MMX[0].push_back(std::make_pair(0, (M-1)*3+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M-2, 3*M+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M-1, 3*M+1));
        HMMP->MMX[0].push_back(std::make_pair(3*M, 3*M+1));
	*/
	

//	for(std::vector < std::pair< int, int> >::const_iterator it = HMMP->MMX[0].begin(); it != HMMP->MMX[0].end(); it++)
//		std::cout << __FUNCTION__ << " " << it->first << " " << it->second << std::endl;
//	exit(1);
    }

    void FillEMI( boost::shared_ptr< Vec3> &HMMP)
    {
        HMMP->EMI = Emissions;
    }


    void ConvertHMM( boost::shared_ptr< Vec3> &HMMProfile)
    {
        CalculateBeginStateTransitionProbabilities();
        CalculateEmissionScores();
        FillTMX( HMMProfile);
        FillMasks( HMMProfile);
	FillEMI( HMMProfile);
        HMMProfile->fillPMX();
        return;
    }

///*
    void fromHMMER(char *HMMfilename, char *fastafilename, boost::shared_ptr< Vec3> &HMMP)
    {
        int            argc     = 4;
        char           **argv   = new char*[4];
        // --- Check if this is true
        argv[0] = new char[1];
        argv[0][0] = 'X';
	argv[1] = "--sw";
        argv[2] = HMMfilename;
        argv[3] = fastafilename;

        // --- end
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
        P7_TRACE       *tr      = NULL;
        ESL_SQ         *sq      = NULL;
        ESL_SQFILE     *sqfp    = NULL;
        int             format  = eslSQFILE_UNKNOWN;
        float           sc, fsc, bsc;

        // Read in one query profile
        if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
        if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
        p7_hmmfile_Close(hfp);

        // Read in one target sequence
        sq     = esl_sq_CreateDigital(abc);
        if (esl_sqfile_Open(seqfile, format, NULL, &sqfp) != eslOK) p7_Fail("Failed to open sequence file %s", seqfile);
        if (esl_sqio_Read(sqfp, sq)                       != eslOK) p7_Fail("Failed to read sequence");
        esl_sqfile_Close(sqfp);

        // Configure a profile from the HMM
        bg = p7_bg_Create(abc);
        gm = p7_profile_Create(hmm->M, abc);

	p7_ProfileConfig(hmm, bg, gm, sq->n, p7_UNILOCAL);  // In the original algorithm it is p7_LOCAL. Be careful!

	/* Allocate matrix and a trace */
	fwd = p7_gmx_Create(gm->M, sq->n);
	float          *xmx     = fwd->xmx;   // For the XMX macro, but we don't use it
	tr  = p7_trace_Create();

	/* Run Viterbi; do traceback */
	p7_GViterbi (sq->dsq, sq->n, gm, fwd, &sc);
	p7_GTrace   (sq->dsq, sq->n, gm, fwd, tr);

	int M = gm->M;
	HMMP->Vtrace.resize(3*(M)+2);
	int state;
	for (int z = 0; z < tr->N; z++)
	{
		if (tr->i[z] > 0)
		{
			std::cout << __FUNCTION__ << " " << tr->i[z] << " " << p7_hmm_DecodeStatetype( tr->st[z]) << tr->k[z] << std::endl;
			if (strcmp(p7_hmm_DecodeStatetype( tr->st[z]),"M") == 0)
			{
				state = 1 + 3*tr->k[z];
				std::cout << "M" << state << std::endl;
			}
			else if (strcmp(p7_hmm_DecodeStatetype( tr->st[z]),"I") == 0 || strcmp(p7_hmm_DecodeStatetype( tr->st[z]),"N") == 0)
			{
				state = 2 + 3*tr->k[z];
			}
			else if (strcmp(p7_hmm_DecodeStatetype( tr->st[z]),"C") == 0)
			{
				// THIS and N are wrong, because I didn't integrate the N and C "LOOP" Tstates yet. They will be placed in (0,0) and (last,last) of the T matrix, and all the transition from B to Ms will be in fact transitions from N to B to Ms; whereas all transitions from Ms to E will be transitions from Ms to C to E.
				state = 3*(M)-1;
			}
			else
			{
				continue;
			}
			HMMP->Vtrace[tr->i[z]].push_back( std::make_pair( state, 1.0));  // Double is for generalization to FB (posterior probability of being in that state)
		}
	}

	std::cout << __FUNCTION__ << " " << sq->n << std::endl;

	std::cout << __FUNCTION__ << " HMMER part done! " << std::endl;

        // Now we convert the structures into my kind of structures
        HMMP->TMX.resize(3*(M)+2); // 1 Begin state, M Subs states, M Ins states, M Del states, 1 End state
        for (int i=0; i<3*(M)+2; i++)
        {
            HMMP->TMX[i].resize(3*(M)+2);
            for (int j=0; j<3*(M)+2; j++)
                HMMP->TMX[i][j] = -d_inf;
        }

	std::cout << __FUNCTION__ << " TMX Resized: " << HMMP->TMX.size() << std::endl;
        for (int j=0; j<M; j++)
	{
            std::cout << __FUNCTION__ << " " <<  p7P_TSC(gm, j, p7P_BM) << std::endl;
            HMMP->TMX[0][j*3+1] = p7P_TSC(gm, j, p7P_BM);
	}
	HMMP->TMX[0][M*3+1] = -d_inf;
        for (int i=0; i<M-1; i++)
        {
//	std::cout << __FUNCTION__ << " " << p7P_TSC(gm, 0, p7P_MM) << " " << p7P_TSC(gm, 0, p7P_MI) << " " <<  p7P_TSC(gm, 0, p7P_MD) << " " << p7P_TSC(gm, 0, p7P_II) << " " << p7P_TSC(gm, 0, p7P_IM) << " " << p7P_TSC(gm, 0, p7P_DM) << " " << p7P_TSC(gm, 0, p7P_DD) << " " << p7P_TSC(gm, 0, p7P_BM) << std::endl;
//            std::cout << __FUNCTION__ << " " << p7P_TSC(gm, i+1, p7P_MM) << " " << p7P_TSC(gm, i+1, p7P_MI) << " " << p7P_TSC(gm, i+1, p7P_MD) << " " << p7P_TSC(gm, i+1, p7P_II) << " " << p7P_TSC(gm, i+1, p7P_IM) << " " << p7P_TSC(gm, i+1, p7P_DM) << " " << p7P_TSC(gm, i+1, p7P_DD) << std::endl;
	    /* The Begin state only goes into the first Match state. This is the reason why the entire first row and column have to be initialized with the prior [1,0,...0]
	       This is also the reason why we have to stop the prior evolution in the last row/column */
//	std::cout << __FUNCTION__ << " " << p7P_TSC(gm, i+1, p7P_MM) << std::endl;
            HMMP->TMX[i*3+1][i*3+4] = p7P_TSC(gm, i+1, p7P_MM); 
            HMMP->TMX[i*3+1][i*3+2] = p7P_TSC(gm, i+1, p7P_MI);
            HMMP->TMX[i*3+1][i*3+6] = p7P_TSC(gm, i+1, p7P_MD);
            HMMP->TMX[i*3+2][i*3+2] = p7P_TSC(gm, i+1, p7P_II);
            HMMP->TMX[i*3+2][i*3+4] = p7P_TSC(gm, i+1, p7P_IM);
            HMMP->TMX[i*3+3][i*3+4] = p7P_TSC(gm, i+1, p7P_DM);
            HMMP->TMX[i*3+3][i*3+6] = p7P_TSC(gm, i+1, p7P_DD);
        }
//	std::cout << __FUNCTION__ << " " << 3*M-2 << " " << 3*M+1 << " " << p7P_TSC(gm, M, p7P_MM) << std::endl;
	//----- ARE THESE 3 LINES OK? IS THERE AN M+1th node for MM, IM and DM in tsc?
	// Answer: no. The transitions to the END state are all put to 0 if the model is local (as it should), -inf if it's global. This is also why the last row of the TMX matrix is 0.
        HMMP->TMX[3*M-2][3*M+1] = p7P_TSC(gm, M, p7P_MM);
        HMMP->TMX[3*M-1][3*M+1] = p7P_TSC(gm, M, p7P_IM);
        HMMP->TMX[3*M][3*M+1] = p7P_TSC(gm, M, p7P_DM);

	double accum, dm, md;
	for (int i=1; i<M-1; i++)
	{
		dm = HMMP->TMX[i*3+3][(i+1)*3+1];
		accum = 0.0;
		for (int j=i-1; j>0; j--)
		{
			md = HMMP->TMX[j*3+1][(j+1)*3+3];
			HMMP->TMX[j*3+1][(i+1)*3+1] = md + accum + dm;
			HMMP->TMX[j*3+1][i*3+3] = md + accum;

			accum += HMMP->TMX[j*3+3][(j+1)*3+3];
			HMMP->TMX[j*3+3][(i+1)*3+1] = accum + dm;
			HMMP->TMX[j*3+3][i*3+3] = accum;
		}
		HMMP->TMX[1][(i+1)*3+1] = md + accum + dm;
		HMMP->TMX[1][i*3+3] = md + accum;
	}
	//-----
        std::vector< char> aa_vec;

	std::cout << __FUNCTION__ << " Until here?  " << sq->n << std::endl;

        for(size_t j = 0; j < sq->n; ++j)
        {
//	    std::cout << __FUNCTION__ << " " << j+1 << std::endl;
//	    std::cout << __FUNCTION__ << " " << (int) sq->dsq[j+1] << std::endl;
//	    std::cout << __FUNCTION__ << " " << itoc[(int) sq->dsq[j+1]] << std::endl;
            aa_vec.push_back(itoc[(int) sq->dsq[j+1]]);
        }

	std::cout << __FUNCTION__ << " Until here?  " << std::endl;

	for (size_t i = 0; i < aa_vec.size(); i++)
        {
            char aa_type = aa_vec[i];
            HMMP->EMI[aa_type].resize(M*3+2);
            HMMP->EMI[aa_type][0] = 0.0;
            for(unsigned int j=0; j<M; j++)
            {
//		std::cout << __FUNCTION__ << " " << i << " " << aa_type << " " << 3*j+1 << " " << p7P_MSC(gm, j+1, sq->dsq[i+1]) << std::endl;
//		std::cout << __FUNCTION__ << " " << i << " " << 3*j+2 << " " << p7P_ISC(gm, j, sq->dsq[i+1]) << std::endl;
                HMMP->EMI[aa_type][3*j+1] = p7P_MSC(gm, j+1, sq->dsq[i+1]);
//                HMMP->EMI[aa_type][3*j+2] = p7P_ISC(gm, j, sq->dsq[i+1]);
		HMMP->EMI[aa_type][3*j+2] = 0; //p7P_MSC(gm, j+1, sq->dsq[i+1]);  //IN HMMER, ISCs are put to 0. Then, we use MSC. Also, I0 doesn't seem to exist, so boh.
                HMMP->EMI[aa_type][3*j+3] = 0;  // Deletions
            }
//            HMMP->EMI[aa_type][3*M+1] = p7P_MSC(gm, j+1, sq->dsq[i+1]);
	}


	/* //SENZA ALCUN SENSO
	int L = aa_vec.size();
	for (int i = 0; i < M; ++i)
	{
		double s = 0.0;
		for (int j = 0; j < M; ++j)
		{
			double h = -d_inf;
			for (int ii = 0; ii < M; ++ii)
			{
				h = log(exp(h) + exp(HMMP->TMX[j][ii]));
			}
			std::cout << "TMX " << h << std::endl;
			for (int k = 0; k < L; ++k)
			{
				char ck = aa_vec[k];
				s += exp(h + HMMP->EMI[ck][j]) / (1.0 + exp(h + HMMP->EMI[ck][j]));
			}
		}
		std::cout << "NORM " << s << std::endl;
	}
	exit(1);
	*/

	LastNode = M-1;
	return;
    }
//*/

    void HMMER_ConvertHMM( char *HMMfilename, char *fastafilename, boost::shared_ptr< Vec3> &HMMProfile)
    {
        fromHMMER( HMMfilename, fastafilename, HMMProfile);
	FillMasks( HMMProfile);
        HMMProfile->fillPMX();
        return;
    }
};
