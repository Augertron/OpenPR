//////////////////////////////////////////////////////////////////////////////////////////////////
// Author:		Styx(Zhenhui Xu)
// Version:		0.1
// Date:		May 20, 2009 
// Description: SVD Decomposition algorithms.
// 
// Copyright(C) 2009 OpenPR
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright notice, this
//       list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright notice, 
//       this list of conditions and the following disclaimer in the documentation 
//       and/or other materials provided with the distribution.
//     * Neither the name of the OpenPR nor the names of its contributors may
//       be used to endorse or promote products derived from this software without 
//       specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
// SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////////////////////////////////

#include "NyxSvd.h"
#include <cmath>
#include <iostream>

/*********************************************************************************
*  Original SVD decomposition
*	a = u * w * v^T
*	a: m x n
*	u: m x m
*	v: n x n
*	q: n x 1  -->  w
*	eps: 1e-15 is a proper choice
*	tol: 0.1*eps is right
*	withu: whether u is needed to be computed 
*	withv: whether v is needed to be computed
*  PS: the singularities have not been ordered in this function.
*********************************************************************************/
int algol_svd(const NyxMat<double>& a, int withu, int withv, NyxMat<double>& u, 
			  NyxMat<double> &q, NyxMat<double> &v, double eps, double tol)
{
	int i,j,k,l,l1,iter,retval;
	double c,f,g,h,s,x,y,z;
	double *e;
	int m = a.GetRow();
	int n = a.GetCol();

	e = (double *)calloc(n,sizeof(double));
	retval = 0;

	/* Copy 'a' to 'u' */    
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++)
			u(i,j) = a(i,j);
	}
	/* Householder's reduction to bidiagonal form. */
	g = x = 0.0;    
	for (i=0;i<n;i++) {
		e[i] = g;
		s = 0.0;
		l = i+1;
		for (j=i;j<m;j++)
			s += (u(j,i)*u(j,i));
		if (s < tol)
			g = 0.0;
		else {
			f = u(i,i);
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
			u(i,i) = f - g;
			for (j=l;j<n;j++) {
				s = 0.0;
				for (k=i;k<m;k++)
					s += (u(k,i)*u(k,j));
				f = s / h;
				for (k=i;k<m;k++)
					u(k,j) += (f * u(k,i));
			} /* end j */
		} /* end s */
		q(i,0) = g;
		s = 0.0;
		for (j=l;j<n;j++)
			s += (u(i,j)*u(i,j));
		if (s < tol)
			g = 0.0;
		else {
			f = u(i,i+1);
			g = (f < 0) ? sqrt(s) : -sqrt(s);
			h = f * g - s;
			u(i,i+1) = f - g;
			for (j=l;j<n;j++) 
				e[j] = u(i,j)/h;
			for (j=l;j<m;j++) {
				s = 0.0;
				for (k=l;k<n;k++) 
					s += (u(j,k)*u(i,k));
				for (k=l;k<n;k++)
					u(j,k) += (s * e[k]);
			} /* end j */
		} /* end s */
		y = fabs(q(i,0)) + fabs(e[i]);                         
		if (y > x)
			x = y;
	} /* end i */

	/* accumulation of right-hand transformations */
	if (withv) {
		for (i=n-1;i>=0;i--) {
			if (g != 0.0) {
				h = u(i,i+1) * g;
				for (j=l;j<n;j++)
					v(j,i) = u(i,j)/h;
				for (j=l;j<n;j++) {
					s = 0.0;
					for (k=l;k<n;k++) 
						s += (u(i,k) * v(k,j));
					for (k=l;k<n;k++)
						v(k,j) += (s * v(k,i));

				} /* end j */
			} /* end g */
			for (j=l;j<n;j++)
				v(i,j) = v(j,i) = 0.0;
			v(i,i) = 1.0;
			g = e[i];
			l = i;
		} /* end i */

	} /* end withv, parens added for clarity */

	/* accumulation of left-hand transformations */
	if (withu) {
		for (i=n;i<m;i++) {
			for (j=n;j<m;j++)
				u(i,j) = 0.0;
			u(i,i) = 1.0;
		}
	}
	if (withu) {
		for (i=n-1;i>=0;i--) {
			l = i + 1;
			g = q(i,0);
			for (j=l;j<m;j++)  /* upper limit was 'n' */
				u(i,j) = 0.0;
			if (g != 0.0) {
				h = u(i,i) * g;
				for (j=l;j<m;j++) { /* upper limit was 'n' */
					s = 0.0;
					for (k=l;k<m;k++)
						s += (u(k,i)*u(k,j));
					f = s / h;
					for (k=i;k<m;k++) 
						u(k,j) += (f * u(k,i));
				} /* end j */
				for (j=i;j<m;j++) 
					u(j,i) /= g;
			} /* end g */
			else {
				for (j=i;j<m;j++)
					u(j,i) = 0.0;
			}
			u(i,i) += 1.0;
		} /* end i*/
	} /* end withu, parens added for clarity */

	/* diagonalization of the bidiagonal form */
	eps *= x;
	for (k=n-1;k>=0;k--) {
		iter = 0;
test_f_splitting:
		for (l=k;l>=0;l--) {
			if (fabs(e[l]) <= eps) goto test_f_convergence;
			if (fabs(q(l-1,0)) <= eps) goto cancellation;
		} /* end l */

		/* cancellation of e[l] if l > 0 */
cancellation:
		c = 0.0;
		s = 1.0;
		l1 = l - 1;
		for (i=l;i<=k;i++) {
			f = s * e[i];
			e[i] *= c;
			if (fabs(f) <= eps) goto test_f_convergence;
			g = q(i,0);
			h = q(i,0) = sqrt(f*f + g*g);
			c = g / h;
			s = -f / h;
			if (withu) {
				for (j=0;j<m;j++) {
					y = u(j,l1);
					z = u(j,i);
					u(j,l1) = y * c + z * s;
					u(j,i) = -y * s + z * c;
				} /* end j */
			} /* end withu, parens added for clarity */
		} /* end i */
test_f_convergence:
		z = q(k,0);
		if (l == k) goto convergence;

		/* shift from bottom 2x2 minor */
		iter++;
		if (iter > 30) {
			retval = k;
			break;
		}
		x = q(l,0);
		y = q(k-1,0);
		g = e[k-1];
		h = e[k];
		f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y);
		g = sqrt(f*f + 1.0);
		f = ((x-z)*(x+z) + h*(y/((f<0)?(f-g):(f+g))-h))/x;
		/* next QR transformation */
		c = s = 1.0;
		for (i=l+1;i<=k;i++) {
			g = e[i];
			y = q(i,0);
			h = s * g;
			g *= c;
			e[i-1] = z = sqrt(f*f+h*h);
			c = f / z;
			s = h / z;
			f = x * c + g * s;
			g = -x * s + g * c;
			h = y * s;
			y *= c;
			if (withv) {
				for (j=0;j<n;j++) {
					x = v(j,i-1);
					z = v(j,i);
					v(j,i-1) = x * c + z * s;
					v(j,i) = -x * s + z * c;
				} /* end j */
			} /* end withv, parens added for clarity */
			q(i-1,0) = z = sqrt(f*f + h*h);
			c = f/z;
			s = h/z;
			f = c * g + s * y;
			x = -s * g + c * y;
			if (withu) {
				for (j=0;j<m;j++) {
					y = u(j,i-1);
					z = u(j,i);
					u(j,i-1) = y * c + z * s;
					u(j,i) = -y * s + z * c;
				} /* end j */
			} /* end withu, parens added for clarity */
		} /* end i */
		e[l] = 0.0;
		e[k] = f;
		q(k,0) = x;
		goto test_f_splitting;
convergence:
		if (z < 0.0) {
			/* q[k] is made non-negative */
			q(k,0) = - z;
			if (withv) {
				for (j=0;j<n;j++)
					v(j,k) = -v(j,k);
			} /* end withv, parens added for clarity */
		} /* end z */
	} /* end k */

	free(e);
	return retval;
}


/***********************************************************************
 * Wrapper of the original SVD Decomposition.
 * Here the singularities have been arranged in descending order.
 ***********************************************************************/
int nyx_svd(const NyxMat<double>& a, int withu, int withv, 
			NyxMat<double>& u, NyxMat<double>& q, NyxMat<double> &v)
{
	int err;
	int m = a.GetRow();
	int n = a.GetCol();
	double EPSILON = 1e-15;

	if (withu)
		u = NyxMat<double>(m, m);
	else
		u = NyxMat<double>(m, n);  /* u不需要计算出来，因此只需要与a同维的u作为临时变量 */

	v = NyxMat<double>(n, n);
	q = NyxMat<double>(n, 1);
	err = algol_svd(a, withu, withv, u, q, v, EPSILON, 0.1*EPSILON);
	if (err != 0)
	{
		std::cout << "nyx_svd WARNING: svd在第 " << err << " 次迭代失败!" << std::endl;
	}

	/*  按奇异值大小重新排列 u, q, v  */
	arrange_uqv(withu, withv, u, q, v, m, n);
	
	return err;
}


/***********************************************************************
 * Arrange the vector q, in descending order. At the same time, u and v 
 * are arranged correspondently.
 ***********************************************************************/
void arrange_uqv(int withu, int withv, NyxMat<double>& u, NyxMat<double>& q, 
				 NyxMat<double> &v, int m, int n)
{
	/* 冒泡排序 */
	int i, j;
	int maxind;
	double tmp, maxval = 0;

	for (i = 0; i < n; ++i)
	{
		maxind = i;
		for (j = i+1; j < n; ++j)
		{
			if (q(j,0) > q(maxind,0))
				maxind = j;
		}
		
		/* 交换最大值到最前面 */
		if (maxind != i)
		{
			tmp = q(i,0);
			q(i,0) = q(maxind,0);
			q(maxind,0) = tmp;

			if (withu)
				u.ExchangeCols(i, maxind);

			if (withv)
				v.ExchangeCols(i, maxind);
		}
	}
}
