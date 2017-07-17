#include <iostream>
#include <vector>

class Edo
{
private:
	std::vector< double>*	m_V;
public:
	std::vector< double>    VX[4];

	Edo(std::vector< double>* V)
	: m_V(V)
	{
		for (int vi=0; vi<4; vi++)
		{
			VX[vi].resize(100);
			for (int i=0; i<VX[vi].size(); i++)
				VX[vi][i] = i*10+vi;
		}
	}

	void read_V()
	{
		for (int i=0; i<m_V->size(); i++)
			std::cout << m_V->at(i) << std::endl;
	}
};

int main ()
{
	std::vector< double> vec;

	vec.resize(44);
	for (int i=0; i<vec.size(); i++)
		vec[i] = i*5;

	std::vector< double>* pvec = &vec;

	for (int i=0; i<pvec->size(); i++)
		std::cout << pvec->at(i) << std::endl;

	Edo* edoinst = new Edo(pvec);
	edoinst->read_V();


	//------------
	std::vector< double>* vvv[4] = {&(edoinst->VX[0]), &(edoinst->VX[1]), &(edoinst->VX[2]), &(edoinst->VX[3])};
	std::cout << "vvv[2][9]   " << (*vvv[2])[9] << std::endl;

	return(0);
}
