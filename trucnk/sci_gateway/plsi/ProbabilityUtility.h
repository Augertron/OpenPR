//author: Mingbo Wang
//revised by: Jia Wu

#pragma once

static unsigned int seed=1;

class ProbabilityUtility
{
public:

	ProbabilityUtility(void);
	~ProbabilityUtility(void);	
	double * getPSum(int size, int precision);
	int * getDiffValue(int maxNum,int count);
	int seed1;
};
