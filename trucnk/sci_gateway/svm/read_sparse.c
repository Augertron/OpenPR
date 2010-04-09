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

/****************************************************************************************\
                                COPYRIGHT NOTICE
                                ----------------

  The code has been derived from libsvm library (version 2.88 for matlab interface)
  (http://www.csie.ntu.edu.tw/~cjlin/libsvm).

  Here is the orignal copyright:
------------------------------------------------------------------------------------------
    Copyright (c) 2000-2008 Chih-Chung Chang and Chih-Jen Lin
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.

    3. Neither name of copyright holders nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.


    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
\****************************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>

#include "mex.h"


#define max(x,y) (((x)>(y))?(x):(y))
#define min(x,y) (((x)<(y))?(x):(y))

void exit_with_help_read()
{
	mexPrintf(
	"Usage: [label_vector, instance_matrix] = readsparse(fname);\n"
	);
}

static void fake_answer(mxArray *plhs[])
{
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

static char *line;
static int max_line_len;

static char* readline(FILE *input)
{
	int len;
	
	if(fgets(line,max_line_len,input) == NULL)
		return NULL;

	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line, max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
			break;
	}
	return line;
}

// read in a problem (in svmlight format)
void read_problem(const char *filename, mxArray *plhs[])
{
	int max_index, min_index, inst_max_index, i;
	long elements, k;
	FILE *fp = fopen(filename,"r");
	int l = 0;
	char *endptr;
	int *ir, *jc;
	double *labels, *samples;
	
	if(fp == NULL)
	{
		mexPrintf("can't open input file %s\n",filename);
		fake_answer(plhs);
		return;
	}

	max_line_len = 1024;
	line = (char *) malloc(max_line_len*sizeof(char));

	max_index = 0;
	min_index = 1; // our index starts from 1
	elements = 0;
	while(readline(fp) != NULL)
	{
		char *idx, *val;
		// features
		int index = 0;

		inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
		strtok(line," \t"); // label
		while (1)
		{
			idx = strtok(NULL,":"); // index:value
			val = strtok(NULL," \t");
			if(val == NULL)
				break;

			errno = 0;
			index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || index <= inst_max_index)
			{
				mexPrintf("Wrong input format at line %d\n",l+1);
				fake_answer(plhs);
				return;
			}
			else
				inst_max_index = index;

			min_index = min(min_index, index);
			elements++;
		}
		max_index = max(max_index, inst_max_index);
		l++;
	}
	rewind(fp);

	// y
	plhs[0] = mxCreateDoubleMatrix(l, 1, mxREAL);
	// x^T
	if (min_index <= 0)
		plhs[1] = mxCreateSparse(max_index-min_index+1, l, elements, mxREAL);
	else
		plhs[1] = mxCreateSparse(max_index, l, elements, mxREAL);

	labels = mxGetPr(plhs[0]);
	samples = mxGetPr(plhs[1]);
	ir = mxGetIr(plhs[1]);
	jc = mxGetJc(plhs[1]);

	k=0;
	for(i=0;i<l;i++)
	{
		char *idx, *val, *label;
		jc[i] = k;

		readline(fp);

		label = strtok(line," \t");
		labels[i] = (int)strtol(label,&endptr,10);
		if(endptr == label)
		{
			mexPrintf("Wrong input format at line %d\n",i+1);
			fake_answer(plhs);
			return;
		}

		// features
		while(1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");
			if(val == NULL)
				break;

			ir[k] = (int) (strtol(idx,&endptr,10) - min_index); // precomputed kernel has <index> start from 0

			errno = 0;
			samples[k] = strtod(val,&endptr);
			if (endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
			{
				mexPrintf("Wrong input format at line %d\n",i+1);
				fake_answer(plhs);
				return;
			}
			++k;
		}
	}
	jc[l] = k;

	fclose(fp);
	free(line);

	{
		mxArray *rhs[1], *lhs[1];
		rhs[0] = plhs[1];
		if(mexCallSCILAB(1, lhs, 1, rhs, "transpose"))
		{
			mexPrintf("Error: cannot transpose problem\n");
			fake_answer(plhs);
			return;
		}
		plhs[1] = lhs[0];
	}
}

#ifdef __SCILAB__
void mex_read_sparse( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
#else
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] )
#endif		
{
	if(nrhs == 1)
	{
		char filename[256];

		mxGetString(prhs[0], filename, mxGetN(prhs[0]) + 1);

		if(filename == NULL)
		{
			mexPrintf("Error: filename is NULL\n");
			return;
		}

		read_problem(filename, plhs);
	}
	else
	{
		exit_with_help_read();
		fake_answer(plhs);
		return;
	}
}

