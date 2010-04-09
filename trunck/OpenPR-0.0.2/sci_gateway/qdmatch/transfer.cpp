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

#include "transfer.h"


extern "C" {

	#include "stack-c.h"
	
}


#define MAT_WIDTH 5

/*vector<Match> SciMat2VecMatch(int nPos)
{
	int m, n, l;
	double *pmat;
	int i, j;
	
	GetRhsVar(nPos, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
	pmat = stk(l);

	vector<Match> matches;	
	for(i = 0; i < m; i++)
	{
		Point tmp_pt_t(pmat[m+i], pmat[i]);
		Point tmp_pt_f(pmat[(3*m)+i], pmat[(2*m)+i]);
		double tmp_dist = pmat[(4*m)+i];
		
		Match tmp_mh(tmp_pt_t, tmp_pt_f, tmp_dist);
		 
		matches.push_back(tmp_mh);	
	}
	
	return matches;
}*/

int SciMat2VecMatch(int nPos, vector<Match>& matches)
{
	int m, n, l;
	double *pmat;
	int i, j;
	
	GetRhsVar(nPos, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
	pmat = stk(l);

	for(i = 0; i < m; i++)
	{
		Point tmp_pt_t(pmat[m+i], pmat[i]);
		Point tmp_pt_f(pmat[(3*m)+i], pmat[(2*m)+i]);
		double tmp_dist = pmat[(4*m)+i];
		
		Match tmp_mh(tmp_pt_t, tmp_pt_f, tmp_dist);
		 
		matches.push_back(tmp_mh);	
	}	
	
	return 0;
}

int VecMatch2SciMat(const vector<Match> &matches, int nPos)
{                                                                                                            
	int matchNum = matches.size();
	double *pmat = new double[matchNum*MAT_WIDTH];
	
	int i = 0;
	vector<Match>::const_iterator j;
	
//	for(i = 0; i < matchNum; i++)
//	{
		for(j = matches.begin(); j != matches.end(); ++j)
		{
			pmat[i] = (*j).GetTarget().GetY();
			pmat[matchNum+i] = (*j).GetTarget().GetX();
			pmat[(2*matchNum)+i] = (*j).GetFound().GetY();
			pmat[(3*matchNum)+i] = (*j).GetFound().GetX();
			pmat[(4*matchNum)+i] = (*j).GetDist();	
			
			i++;				
		}
//	}
	
	int m, n;
	m = matchNum;
	n = MAT_WIDTH;
	CreateVarFromPtr(nPos, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &pmat);
	
	delete []pmat;
	
	return 0;
}
