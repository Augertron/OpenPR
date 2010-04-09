////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Usage:
// Author: Mingbo Wang (NLPR)
// Revised by: Jia Wu 
// Reference:
//
// PLSITrian::train参数说明
//
//  datafile: 样本的数据文件，文件每一行代表一个文档向量
//  issparse: 0-不采用稀疏方式；1-采用稀疏方式, 即wordid:count，其中wordid=0,1,2,...
//  outputpath: 最终结果(pz.txt, pzd.txt, pwz.txt, foldnums.txt)的输出路径，其中foldnums.txt保存了用于TEM中预留的文档id
//  nz: 最终topic个数
//  beta: 算法中beta的值
//  itenum: 最大迭代次数
//  frate: 预留文档比例
//  eta: beta的衰减系数
//  ite: 
//  modelpath: 前一次结果(pz.txt, pzd.txt, pwz.txt, foldnums.txt)的路径，当需要继续前一次结果进行迭代时指定，默认为空
////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

//#include <string>
//#include "DataSet.h"

//using namespace std;

class PLSITrain  
{
public:	

	PLSITrain();
	virtual ~PLSITrain();		

	void train(int docnum, int wordnum, int **data, int nz, int beta, int itenum, double frate, double eta, int iter,  double *pz, double **pwz, double **pdz, int *foldnums);

private:
	void em(int **data, double *pz,double **pwz, double **pdz, int *foldnums);
	double emIteration(int **data,double * pz,double **pwz,double **pdz,double *pzC,double **pwzC,double ** pdzC,int nd,int nw,int nz,int * foldnums,int foldnum,double beta,bool iter);	

	int getNDW(int **data, int, int);
	double getPDW(int nz,int d,int w,double* pz,double **pwz,double **pdz);
	double getPZDW(int z,int w,int d,double *pz,double **pwz,double ** pdz,double** normal,double beta);
    double getPerplexity(int **data, int *foldnums, int nw, int nz, int foldnum, double *pz, double **pwz, double **pdz);  
	void Mstep_1(int **data,double* pzC,double ** pwzC,double ** pdzC,double *pz,double ** pwz,double ** pdz,int nw,int nd,int nz,double beta);
	void buildNormal(int nz,int nd,int nw,double* pz,double **pwz,double ** pdz,double** normal,double beta);

	int foldnum;
	int itenum;
	bool iter;

	int nz;
	int nd;
	int nw;

	double beta;
	double eta;
	
//	int **data;
};
