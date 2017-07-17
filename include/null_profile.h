#include "function.t.h"
#include "sequence.h"

#include <list>


class ScoreNullProfile
: public Function< std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid>, double>
{
 protected:
  size_t    m_ProfileID;    //!< ID of the profile in the GeneralizedAminoAcid that an object of this class refers to, needed in operator()

 public:

  ScoreNullProfile(){};

  ~ScoreNullProfile(){};

double operator()( const std::pair< GeneralizedAminoAcid, GeneralizedAminoAcid> &AA)
{
#ifdef SECURE
        if (AA.first.GetProfiles().size() != 0)
        {
                std::cerr << "Null Profile is not null! \n";
                exit(-1);
        }
#endif
        return 0.0;  // ask substitution matrix/map for value of AA pair 'AB'
};

}; // end class
