#include "DataSet.h"

extern "C" 
{
	#include "stack-c.h"

	#include "sciprint.h"
	#include "Scierror.h"

	int int_readdata(char *fname) 
	{
		if(Rhs == 0)
		{
			sciprint("Usage:\n[docnum wordnum matrix] = plsiread(filename, issparse)\n");
			return -1;
		}
	  	int mR1 = 0, nR1 = 0, lR1 = 0;
	  	int mR2 = 0, nR2 = 0, lR2 = 0;
		int mL, nL;	
	  	
	    char *filename = NULL;
	    int sparse;
	    bool init = false;	    	  	
	
	  	CheckRhs(2,2); 	  	
	  	CheckLhs(3,3); 
	  	
	    DataSet dataset;
		dataset.clear();
	  	
	  	GetRhsVar(1, STRING_DATATYPE, &mR1, &nR1, &lR1); 
	  	filename = cstk(lR1);	  	

	  	GetRhsVar(2, MATRIX_OF_DOUBLE_DATATYPE, &mR2, &nR2, &lR2);
	  	sparse = (int)*stk(lR2); 
	  	
	  	//read data from datafile
		if(sparse)
			init = dataset.initData(filename);
		else
			init = dataset.initDataF(filename);
		
		if(!init)
		{
			Scierror(999, "%s: Data initialization from file %s failed.\r\n", fname, filename);
			return -1;
		}
	    
		//docnum
		double *docnum = new double[1];
		docnum[0] = (double)dataset.getDocNum();
		mL = 1;
		nL = 1;
		CreateVarFromPtr(Rhs+1, MATRIX_OF_DOUBLE_DATATYPE, &mL, &nL, &docnum);
		delete []docnum;		

		//wordnum
		double *wordnum = new double[1];
		wordnum[0] = (double)dataset.getWordNum();
		mL = 1;
		nL = 1;
		CreateVarFromPtr(Rhs+2, MATRIX_OF_DOUBLE_DATATYPE, &mL, &nL, &wordnum);
		delete []wordnum;
////		CreateVar(Rhs+2, MATRIX_OF_DOUBLE_DATATYPE, &mL2, &nL2, &lL2);
////		*stk(lL2) = dataset.getWordNum();

		//dwmatrix
		mL = dataset.getDocNum();
		nL = dataset.getWordNum();
		double *pData = new double[mL*nL];
//		mL = 3;
//		nL = 600;
//		double *pData = new double[mL*nL];
		
		//the data should be stored in column sequence
		int i, j;
/*		for(i = 0; i < mL; i++)
		{
			for(j = 0; j < nL; j++)
			{
				pData[i*nL+j] = (double)dataset.data[i][j];
			}
		}*/
		int k = 0;
		for(j = 0; j < nL; j++)
		{
			for(i = 0; i < mL; i++)
			{
				pData[k++] = (double)dataset.data[i][j];
			}
		} 
		
		
		//test
/*		for(i = 0; i < nL; i++)
		{
			cout<<pData[i]<<" ";
		}		
		cout<<endl;*/
		
		CreateVarFromPtr(Rhs+3, MATRIX_OF_DOUBLE_DATATYPE, &mL, &nL, &pData);
		delete []pData;		
		
		LhsVar(1) = Rhs+1;
		LhsVar(2) = Rhs+2;
		LhsVar(3) = Rhs+3;
  	
  	  	return 0;
	}
}
