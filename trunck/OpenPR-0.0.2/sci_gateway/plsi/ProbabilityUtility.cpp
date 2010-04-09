//author: Mingbo Wang
//revised by: Jia Wu

#include <ctime>
#include <cstdlib>

#include "ArrayUtility.h"
#include "ProbabilityUtility.h"

ProbabilityUtility::ProbabilityUtility(void)
{
	seed=(int)time(NULL);
}

ProbabilityUtility::~ProbabilityUtility(void)
{
}

double * ProbabilityUtility::getPSum(int size, int precision)
{

	double * p=ArrayUtility<double>::initialArray(size,0);
	int mod=(int)pow(10.0,precision);
	srand(seed++);
	int tsize=0;
	double num=1;
	double sum=0;
	double count=0;
	while(tsize<size)
	{
		count=rand()%mod+1;
	
		sum+=count;
		p[tsize]=count;
		tsize++;
	}
	for(int i=0;i<size;i++)
	p[i]/=sum;
	return p;
}

int * ProbabilityUtility::getDiffValue( int size,int count)
{
	int * array=ArrayUtility<int>::initialArray(size,0);
	int i=0;
	for(;i<size;i++)
		array[i]=i;
	i=size;
	srand(seed++);

	int index;
	for(;i>1;i--)
	{
    		index=rand()%i;

		int temp=array[i-1];
		array[i-1]=array[index];
		array[index]=temp;
	}

	int * res=ArrayUtility<int>::initialArray(count,0);
	for( i=0;i<count;i++)
		res[i]=array[i];
	delete []array;
	return res;
}
