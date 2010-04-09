// author: Mingbo Wang
// revised by: Jia Wu

#include <iostream>
#include <fstream>

#include "ArrayUtility.h"
#include "PLSITrain.h"

using namespace std;

const double zeor = 1e-6;

PLSITrain::PLSITrain()
{
}

PLSITrain::~PLSITrain()
{
}

void PLSITrain::train(int docnum, int wordnum, int **data, int nz, int beta, int itenum, double frate, double eta, int iter,  double *pz, double **pwz, double **pdz, int *foldnums)
{
	this->nz = nz;
	this->beta = beta;
	this->itenum = itenum;
	this->eta = eta;
	this->iter = iter;
	foldnum = (int)docnum*frate;
	nd = docnum;
	nw = wordnum;	
	
	//EM method
	em(data, pz,pwz,pdz,foldnums);
}

//EM method
void PLSITrain::em(int **data, double *pz,double **pwz,double **pdz,int *foldnums)
{
	bool noItem=false;
	if(itenum==0)
		iter=false;
	if(itenum<0)
		noItem=true;
	int i=0;

	double * pzC=ArrayUtility<double>::initialArray(nz,0);
	double **pwzC=ArrayUtility<double>::initial2DArray(nw,nz,0);
	double **pdzC=ArrayUtility<double>::initial2DArray(nd,nz,0);

	double betapreplexity=0;	

	double lastPerplexity=0;
	double perplexity=0;
	double fper=-1;
	int items=0;

	while((items++)<itenum||noItem)
	{
		cout<<"itenum: "<<items<<endl;
		perplexity=emIteration(data,pz,pwz,pdz,pzC,pwzC,pdzC,nd,nw,nz,foldnums,foldnum,beta,iter);
		if(fabs(perplexity-lastPerplexity)<zeor)
		{
			cout<<lastPerplexity<<" outer equals "<<perplexity<<endl;

			if(fabs(perplexity-fper)<zeor)
			{
				cout<<fper<<" outer equals "<<perplexity<<endl;
				break;
			}
			else
			{
				fper=perplexity;
				lastPerplexity=-1;
				beta*=eta;
			}
		}
		else
		{
			cout<<lastPerplexity<<" "<<perplexity<<" not equlas"<<endl;
			lastPerplexity=perplexity;
			double * tempPZ=pz;
			double ** tempPWZ=pwz;
			double ** tempPDZ=pdz;
			pz=pzC;
			pwz=pwzC;
			pdz=pdzC;
			pzC=tempPZ;
			pwzC=tempPWZ;
			pdzC=tempPDZ;
			ArrayUtility<double>::setValue(pzC,0,nz);
			ArrayUtility<double>::setValue(pwzC,0,nw,nz);
			ArrayUtility<double>::setValue(pdzC,0,nd,nz);
		}
		
		cout<<endl;
	}

	ArrayUtility<double>::finalize2DArray(pwzC,nw);
	ArrayUtility<double>::finalize2DArray(pdzC,nd);
	delete []pzC;
}

// one em iteration
double PLSITrain::emIteration(int **data,double * pz,double **pwz,double **pdz,double *pzC,double **pwzC,double ** pdzC,int nd,int nw,int nz,int * foldnums,int foldnum,double beta,bool iter)
{	
	double preplexity=0;
	cout<<"beta: "<<beta<<endl;
	
	if(iter)
	{
		getchar();
		cout<<"infor:"<<endl;
		cout<<"pz"<<endl;
		ArrayUtility<double>::show(pz,nz);
		cout<<"pwz"<<endl;
		ArrayUtility<double>::show(pwz,nw,nz);
		cout<<"pdz"<<endl;
		ArrayUtility<double>::show(pdz,nd,nz);			
		cout<<"pzdw"<<endl;
	}

	Mstep_1(data,pzC,pwzC,pdzC,pz,pwz,pdz,nw,nd,nz,beta);
	cout<<"Mstep_1 completed"<<endl;
	
	if(iter)
	{
		cout<<"compare"<<endl;
		cout<<"pzC"<<endl;ArrayUtility<double>::show(pzC,nz);
		cout<<"pwzC"<<endl;ArrayUtility<double>::show(pwzC,nw,nz);
		cout<<"pdzC"<<endl;ArrayUtility<double>::show(pdzC,nd,nz);
	}

	double perplex=getPerplexity(data,foldnums,nw,nz,foldnum,pzC,pwzC,pdzC);
	cout<<"emIteration perplex"<<endl;

	return perplex;
}

//pzwd
void PLSITrain::buildNormal(int nz,int nd,int nw,double* pz,double **pwz,double ** pdz,double** normal,double beta)
{
	cout<<"build normal matrix"<<endl;
	for(int w=0;w<nw;w++)
	{
		for(int d=0;d<nd;d++)
		{
			for(int z=0;z<nz;z++)
			{
				normal[d][w]+=pz[z]*pow(pdz[d][z]*pwz[w][z],beta);
			}
		}	
	}
	cout<<"build normal completed."<<endl;
}

void PLSITrain::Mstep_1(int **data, double* pzC,double ** pwzC,double ** pdzC,double *pz,double ** pwz,double ** pdz,int nw,int nd,int nz,double beta)
{

	double ** normal=ArrayUtility<double>::initial2DArray(nd,nw,0);
	buildNormal(nz,nd,nw,pz,pwz,pdz,normal,beta);
	double pzdw;
	int  ndw;
	double value;
	double summaryA;

	int N;
	double p;
	for (int z=0;z<nz;z++)
	{		
		N=0;
		p=0;
		summaryA=0;
		
		for(int w=0;w<nw;w++)
		{
			for(int d=0;d<nd;d++)
			{
				pzdw=getPZDW(z,d,w,pz,pwz,pdz,normal,beta);
				ndw = getNDW(data, d, w);
				N+=ndw;
				value=ndw*pzdw;
				p+=value;
				pwzC[w][z]+=value;
				pdzC[d][z]+=value;
				summaryA+=value;			
			}			
		}
				
		if(summaryA!=0)
		{
			for(int w=0;w<nw;w++)
			{
				pwzC[w][z]/=summaryA;
			}
		
			for(int d=0;d<nd;d++)
			{
				pdzC[d][z]/=summaryA;
			}
		}
		else
			cout<<"error"<<z<<" "<<endl;
		pzC[z]=p;

		if(N!=0)
			pzC[z]/=N;
				
	}
	ArrayUtility<double>::finalize2DArray(normal,nd);
}

double  PLSITrain::getPZDW(int z,int d,int w,double *pz,double **pwz,double ** pdz,double** normal,double beta)
{
	double p=pz[z]*pow(pdz[d][z]*pwz[w][z],beta);
	
	if(normal[d][w]!=0)
		p/=normal[d][w];

	return p;	
}

double PLSITrain::getPerplexity(int **data, int *foldnums, int nw, int nz, int foldnum, double *pz, double **pwz, double **pdz)
{
	double fenzi=0;
	double fenmu=0;
	for(int w=0;w<nw;w++)
	{
		for(int i=0;i<foldnum;i++)
		{
			int count = getNDW(data, foldnums[i], w);
			double pdw=getPDW(nz,foldnums[i],w,pz,pwz,pdz);
			fenzi-=count*pdw;

			fenmu+=count;
		}
	}
	
	if(fenmu==0)
	{
		cout<<"zero!!"<<endl;
	}
	 
	return exp(fenzi/fenmu);
}

int PLSITrain::getNDW(int **data, int d, int w)
{
	if(d >= nd || w >= nw)
		return -1;
	
	return data[d][w];
}

double PLSITrain::getPDW(int nz,int d,int w,double* pz,double **pwz,double **pdz)
{
	double pd=0;
	double pwd=0;
	for(int z=0;z<nz;z++)
	{
		pd+=pz[z]*pdz[d][z];
		pwd+=pz[z]*pdz[d][z]*pwz[w][z];
	}

	if(pd==0)
		return 0;

	return pwd/pd;	
}
