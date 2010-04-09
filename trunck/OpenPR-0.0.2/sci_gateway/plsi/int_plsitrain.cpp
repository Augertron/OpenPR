//////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2009 OpenPR
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of OpenPR nor the names of its 
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL HOLDER AND CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
////////////////////////////////////////////////////////////////////////////

#include "PLSITrain.h"
#include "ArrayUtility.h"
#include "ProbabilityUtility.h"

#include <iostream>

// data stored in Scilab are in column sequence

extern "C"
{
	#include "stack-c.h"

	#include "sciprint.h"
	#include "Scierror.h"
	
	int int_plsitrain(char *fname)
	{
		if(Rhs == 0)
		{
			sciprint("Usage:\n[pz pwz pdz foldnums]=plsitrain(docnum,wordnum,matrix,nz,beta,itenum,frate,eta,iter[,pz,pwz,pdz,foldnums])\n");
			return -1;
		}

		CheckRhs(9, 13);
		CheckLhs(4, 4);
		
		if((Rhs == 9) || (Rhs == 13))
		{	
			int mR, nR, lR;
			int mL, nL, lL;
			
			int i, j, k;
			int docnum, wordnum;

			int nz, beta, itenum;
			double frate, eta;
			int iter;
			
			GetRhsVar(1, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
			docnum = (int)*stk(lR);
			
			GetRhsVar(2, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
			wordnum = (int)*stk(lR);
			
			double *pdata = ArrayUtility<double>::initialArray(docnum*wordnum, 0);
			int **p2Ddata = ArrayUtility<int>::initial2DArray(docnum, wordnum, 0);
			GetRhsVar(3, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);

			if((mR != docnum) || (nR != wordnum))
			{
				Scierror(999, "%s: size of matrix does not match doc number and word number.\r\n", fname);
				return -1;
			}
			pdata = stk(lR);  //pdata stored in column sequence			

			k = 0;
			for(j = 0; j < wordnum; j++)
			{
				for(i = 0; i < docnum; i++)
				{
					p2Ddata[i][j] = (int)pdata[k++];  //p2Ddata stored in row sequence
				}
			}			
			
			GetRhsVar(4, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
			nz = (int)*stk(lR);
			
			GetRhsVar(5, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
			beta = (int)*stk(lR);
			
			GetRhsVar(6, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
			itenum = (int)*stk(lR);
			
			GetRhsVar(7, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
			frate = *stk(lR);
			
			GetRhsVar(8, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
			eta = *stk(lR);
			
			GetRhsVar(9, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
			iter = (int)*stk(lR);
			
			double *pz = ArrayUtility<double>::initialArray(nz, 0);
			double *pwz = ArrayUtility<double>::initialArray(wordnum*nz, 0);
			double *pdz = ArrayUtility<double>::initialArray(docnum*nz, 0);
			int *foldnums = ArrayUtility<int>::initialArray((int)(docnum*frate), 0);
			
			double **pwz2D = ArrayUtility<double>::initial2DArray(wordnum, nz, 0);
			double **pdz2D = ArrayUtility<double>::initial2DArray(docnum, nz, 0);
			
			if(Rhs == 13)
			{			
				double *pz_tmp = ArrayUtility<double>::initialArray(nz, 0);
				double *pwz_tmp = ArrayUtility<double>::initialArray(wordnum*nz, 0);
				double *pdz_tmp = ArrayUtility<double>::initialArray(docnum*nz, 0);
				double *dfoldnums_tmp = ArrayUtility<double>::initialArray((double)(docnum*frate), 0);
			
				GetRhsVar(10, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
				pz_tmp = stk(lR);
				
				GetRhsVar(11, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
				pwz_tmp = stk(lR);
				
				GetRhsVar(12, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
				pdz_tmp = stk(lR);
				
				GetRhsVar(13, MATRIX_OF_DOUBLE_DATATYPE, &mR, &nR, &lR);
				dfoldnums_tmp = stk(lR);	

				//pz
				memcpy(pz, pz_tmp, nz*8);

				//pwz2D
				k = 0;
				for(j = 0; j < nz; j++)
				{
					for(i = 0; i < wordnum; i++)
					{
						pwz2D[i][j] = pwz_tmp[k++];
					}
				}
				//pdz2D
				k = 0;
				for(j = 0; j < nz; j++)
				{
					for(i = 0; i < docnum; i++)
					{
						pdz2D[i][j] = pdz_tmp[k++];
					}
				}
				
				//foldnums
				for(j = 0; j < ((int)(docnum*frate)); j++)
				{
					foldnums[j] = (int)dfoldnums_tmp[j];
				}
			}
			else
			{
				ProbabilityUtility util;	

				for(i = 0; i < nz; i++)
				{
					pz = util.getPSum(wordnum, 5);

					for(j = 0; j < wordnum; j++)
						pwz2D[j][i] = pz[j];	//pwz2D
				}
		
				for(i = 0; i < nz; i++)
				{
					pz = util.getPSum(docnum, 5);
					
					for(j = 0; j < docnum; j++)
						pdz2D[j][i] = pz[j];	//pdz2D
				}
		
				pz = util.getPSum(nz, 5);	//pz	
		
				foldnums = util.getDiffValue(docnum, (int)(docnum*frate));	//foldnums
			}

			PLSITrain train;
	    	train.train(docnum, wordnum, p2Ddata, nz, beta, itenum, frate, eta, iter, pz, pwz2D, pdz2D, foldnums);	

	    	//pz, pwz, pdz, foldnums should be stored in column sequence	    	
	    	k = 0;
	    	for(j = 0; j < nz; j++)
	    	{
	    		for(i = 0; i < wordnum; i++)
	    		{
	    			pwz[k++] = pwz2D[i][j];			//pwz
	    		}
	    	}
	    	k = 0;
	    	for(j = 0; j < nz; j++)
	    	{
	    		for(i = 0; i < docnum; i++)
	    		{
	    			pdz[k++] = pdz2D[i][j];			//pdz
	    		}
	    	}
	    	
	    	mL = 1;
	    	nL = nz;
			CreateVarFromPtr(Rhs+1, MATRIX_OF_DOUBLE_DATATYPE, &mL, &nL, &pz);
			
			mL = wordnum;
			nL = nz;
			CreateVarFromPtr(Rhs+2, MATRIX_OF_DOUBLE_DATATYPE, &mL, &nL, &pwz);

			mL = docnum;
			nL = nz;
			CreateVarFromPtr(Rhs+3, MATRIX_OF_DOUBLE_DATATYPE, &mL, &nL, &pdz);
			
			double *dfoldnums = new double[(int)(docnum*frate)];
			for(i = 0; i < ((int)(docnum*frate)); i++)
			{
				dfoldnums[i] = (double)foldnums[i];
			}
			mL = 1;
			nL = (int)(docnum*frate);
			CreateVarFromPtr(Rhs+4, MATRIX_OF_DOUBLE_DATATYPE, &mL, &nL, &dfoldnums);
					
			ArrayUtility<double>::finalize2DArray(pwz2D, wordnum);
			ArrayUtility<double>::finalize2DArray(pdz2D, docnum);
			
			delete []pz;
			delete []pwz;
			delete []pdz;
			delete []foldnums;
			delete []dfoldnums;
			
			ArrayUtility<int>::finalize2DArray(p2Ddata, docnum);
			
			LhsVar(1) = Rhs+1;
			LhsVar(2) = Rhs+2;
			LhsVar(3) = Rhs+3;
			LhsVar(4) = Rhs+4;
		}
		else
		{
			Scierror(999, "%s: 9 or 13 arguments should be passed.\r\n", fname);
			return -1;		
		}
		
		return 0;
	}
}
