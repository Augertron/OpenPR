/***********************************************************************
 * SIVP - Scilab Image and Video Processing toolbox
 * Copyright (C) 2005  Shiqi Yu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 ***********************************************************************/


#include "common.h"

/**************************************************
 * nRow: the first dim
 * nCol: the second dim
 * pData: pData is the data for matrix, and it can be freed after this function is called.
 * nType: I_CHAR, I_INT16, I_INT32, I_UCHAR,  I_UINT16,  I_UINT32
 **************************************************/
int
Create2DIntMat(int nPos, int nRow, int nCol, void * pData, int nType)
{
  SciIntMat IntMat;

  //create SciIntMat for the image data
  IntMat.m = nRow;
  IntMat.n = nCol;
  IntMat.l = -1; /* dimensions matrix not in scilab stack */
  IntMat.it = nType;
  IntMat.D = pData;

  CreateVarFromPtr(nPos,"I",&(IntMat.m), &(IntMat.n), &IntMat ); 
  return TRUE;
}


/**************************************************
 * Create a 2D nRow * nCol float matrix(single precision, 4bytes)
 * nRow: the first dim
 * nCol: the second dim
 * pData: pData is the data for matrix, and it can be freed after this function is called.
 **************************************************/
int
Create2DFloatMat(int nPos, int nRow, int nCol, float * pData)
{
  CreateVarFromPtr(nPos,"r",&nRow, &nCol, &pData); 
  return TRUE;
}

/**************************************************
 * Create a 2D nRow * nCol double matrix(double precision, 8bytes)
 * nRow: the first dim
 * nCol: the second dim
 * pData: pData is the data for matrix, and it can be freed after this function is called.
 **************************************************/
int
Create2DDoubleMat(int nPos, int nRow, int nCol, double * pData)
{
  CreateVarFromPtr(nPos,"d",&nRow, &nCol, &pData); 
  return TRUE;
}

/**************************************************
 * Create a 3D nRow * nCol * nCh integer matrix
 * nRow: the first dim
 * nCol: the second dim
 * nCh:  the third dim
 * pData: pData is the data for matrix, and it can be freed after this function is called.
 * nType: I_CHAR, I_INT16, I_INT32, I_UCHAR,  I_UINT16,  I_UINT32
 **************************************************/
int 
Create3DIntMat(int nPos, int nRow, int nCol, int nCh, void* pData, int nType)
{
  static char *Str[]= { "hm","dims","entries"}; 
  int m1=1,n1=3;
  int mL=3,nL=1,lL, un=1;

  SciIntMat Dims;
  SciIntMat IntData;

  //create SciIntMat for dimentional information
  Dims.m = 1;
  Dims.n = 3;
  Dims.l = -1; /* dimensions matrix not in scilab stack */
  Dims.it = I_INT32;
  Dims.D = malloc(Dims.n * sizeof(int)); //dim data
  if (!Dims.D)
  {
    Scierror(999,"Unable to alloc memory for the image\n");
    return FALSE;
  }

  IC_INT32(Dims.D)[0] = nRow;
  IC_INT32(Dims.D)[1] = nCol;
  IC_INT32(Dims.D)[2] = nCh;

  //create SciIntMat for the image data
  IntData.m = nRow ;
  IntData.n = nCol * nCh;
  IntData.l = -1; /* dimensions matrix not in scilab stack */
  IntData.it = nType;
  IntData.D = pData;

  CreateVar(nPos,"m", &mL, &nL, &lL);
  CreateListVarFromPtr(nPos,1,"S", &m1, &n1, Str);
  CreateListVarFromPtr(nPos,2,"I",&(Dims.m), &(Dims.n), &Dims); 
  CreateListVarFromPtr(nPos,3,"I",&(IntData.m), &(IntData.n), &IntData ); 

/*   free memory */
  free(Dims.D);

  return TRUE;
}

/**************************************************
 * Create a 3D nRow * nCol * nCh float matrix (single precision, 4 bytes)
 * nRow: the first dim
 * nCol: the second dim
 * nCh:  the third dim
 * pData: pData is the data for matrix, and it can be freed after this function is called.
 **************************************************/
int 
Create3DFloatMat(int nPos, int nRow, int nCol, int nCh, float* pData)
{
  static char *Str[]= { "hm","dims","entries"}; 
  int m1=1,n1=3;
  int mL=3,nL=1,lL;
  int un=1;
  int nSize = nRow * nCol * nCh;

  SciIntMat Dims;

  //create SciIntMat for dimentional information
  Dims.m = 1;
  Dims.n = 3;
  Dims.l = -1; /* dimensions matrix not in scilab stack */
  Dims.it = I_INT32;
  Dims.D = malloc(Dims.n * sizeof(int)); //dim data
  if (!Dims.D)
  {
    Scierror(999,"Unable to alloc memory for the image\n");
    return FALSE;
  }

  IC_INT32(Dims.D)[0] = nRow;
  IC_INT32(Dims.D)[1] = nCol;
  IC_INT32(Dims.D)[2] = nCh;


  CreateVar(nPos,"m", &mL, &nL, &lL);
  CreateListVarFromPtr(nPos,1,"S", &m1, &n1, Str);
  CreateListVarFromPtr(nPos,2,"I",&(Dims.m), &(Dims.n), &Dims); 
  CreateListVarFromPtr(nPos,3,"r", &nSize, &un, &pData );

/*   free memory */
  free(Dims.D);

  return TRUE;
}

/**************************************************
 * Create a 3D nRow * nCol * nCh double matrix (double precision, 8 bytes)
 * nRow: the first dim
 * nCol: the second dim
 * nCh:  the third dim
 * pData: pData is the data for matrix, and it can be freed after this function is called.
 **************************************************/
int 
Create3DDoubleMat(int nPos, int nRow, int nCol, int nCh, double* pData)
{
  static char *Str[]= { "hm","dims","entries"}; 
  int m1=1,n1=3;
  int mL=3,nL=1,lL;
  int un=1;
  int nSize = nRow * nCol * nCh;

  SciIntMat Dims;

  //create SciIntMat for dimentional information
  Dims.m = 1;
  Dims.n = 3;
  Dims.l = -1; /* dimensions matrix not in scilab stack */
  Dims.it = I_INT32;
  Dims.D = malloc(Dims.n * sizeof(int)); //dim data
  if (!Dims.D)
  {
    Scierror(999,"Unable to alloc memory for the image\n");
    return FALSE;
  }

  IC_INT32(Dims.D)[0] = nRow;
  IC_INT32(Dims.D)[1] = nCol;
  IC_INT32(Dims.D)[2] = nCh;


  CreateVar(nPos,"m", &mL, &nL, &lL);
  CreateListVarFromPtr(nPos,1,"S", &m1, &n1, Str);
  CreateListVarFromPtr(nPos,2,"I",&(Dims.m), &(Dims.n), &Dims); 
  CreateListVarFromPtr(nPos,3,"d", &nSize, &un, &pData );

/*   free memory */
  free(Dims.D);

  return TRUE;
}


/************************************************************
 * convert IplImage to SCI matrix
************************************************************/
int IplImg2Mat(IplImage * pImage, int nPos)
{
  void * pMatData;
  int nBytes;
  int nType;
  
  if(pImage == NULL)
    return FALSE;

  //if bottom-left origin
  if(pImage->origin==1)
    {
      cvFlip(pImage, NULL, 0);
      pImage->origin=0;
    }

  /*how many bytes per pixel per channel*/
  nBytes = pImage->depth;
  if (nBytes > IPL_DEPTH_SIGN)
    nBytes -= IPL_DEPTH_SIGN;
  nBytes = nBytes >> 3;
  
  /*alloc memory for matrix data*/
  pMatData = malloc(pImage->width * pImage->height * pImage->nChannels * nBytes);
  if(pMatData == NULL)
    return FALSE;

  ImgData2MatData(pImage, pMatData);

  /*convert IplImage data type to scilab data type*/
  nType = IplType2SciType(pImage->depth);
  if (nType <= 0)
    {
      free(pMatData);
      return FALSE;
    }
  /*create matrix for scilab*/
  if(pImage->nChannels == 1)
    {
      switch(nType){
      case I_CHAR:
      case I_UCHAR:
      case I_INT16:
      case I_UINT16:
      case I_INT32:
	Create2DIntMat(nPos, pImage->height, pImage->width, pMatData, nType);
	break;
      case SIVP_FLOAT:
	Create2DFloatMat(nPos,pImage->height, pImage->width, pMatData);
	break;
      case SIVP_DOUBLE:
	Create2DDoubleMat(nPos,pImage->height, pImage->width, pMatData);
	break;
      }
    }
  else
    {
      switch(nType){
      case I_CHAR:
      case I_UCHAR:
      case I_INT16:
      case I_UINT16:
      case I_INT32:
	Create3DIntMat(nPos, pImage->height, pImage->width, pImage->nChannels, pMatData, nType);
	break;
      case SIVP_FLOAT:
	Create3DFloatMat(nPos,pImage->height, pImage->width, pImage->nChannels, pMatData);
	break;
      case SIVP_DOUBLE:
	Create3DDoubleMat(nPos,pImage->height, pImage->width, pImage->nChannels, pMatData);
	break;
      }
    }

  /*free matrix data*/
  free(pMatData);
  return TRUE;
}


/************************************************************
 * get image dimension information and the data address
 * input: nPos(the nPos'th  argument, the  argument shouble be a mlist)
 * 
************************************************************/
IplImage * CreateIplImgFromHm(int nPos)
{
  char ** pStr;
  int * pListHeader;
  int * pDataHeader;
  void * pData;
  int m, n, l;
  int m1, n1, m2, n2, m3, n3, l3;
  SciIntMat Dims;
  SciIntMat IntMat;
  
  int nWidth, nHeight, nCh=1;
  int iplType;
  IplImage * pImg = NULL;

  GetListRhsVar(nPos, 1 ,"S",&m1,&n1,&pStr);
  /*check whether the the argument is a hypermatrix*/
  if( m1 !=1 || n1 !=3)
    goto NOT_HM;
  if( strcmp(pStr[0], "hm") != 0 || strcmp(pStr[1], "dims") != 0 || strcmp(pStr[2], "entries") != 0)
    goto NOT_HM;
  
  /*get the dimension information, it's stored in the second element of the mlist*/

  GetListRhsVar(nPos,2,"I", &m2, &n2,&Dims);
  if( m2*n2 !=2 && m2*n2 !=3)
    goto NOT_HM;

  nHeight = IC_INT32(Dims.D)[0];
  nWidth  = IC_INT32(Dims.D)[1];
  if(m2*n2 == 3)
    nCh = IC_INT32(Dims.D)[2];
  else
    nCh = 1;


  /*get mlist data structure*/
  pListHeader = (int*)GetData(nPos);
  /*data is stored in the third element*/
  /*data struct is get first*/
  pDataHeader = (int*)( ((char*)pListHeader) + 24 + (pListHeader[4]-1)*8 );
  
  /*the next step is to get the data type and the data*/
  /*the data type is store in the first 4 bytes as a integer*/
  switch( pDataHeader[0]  ){
  case 1:  /*if the data is real*/
    iplType = IPL_DEPTH_64F;
    GetListRhsVar(nPos, 3, "d", &m3, &n3, &l3);
    pData = stk(l3);
    break;
  case 8:  /*integer*/
    GetListRhsVar(nPos, 3, "I", &m3, &n3, &IntMat);
    m3 = IntMat.m;
    n3 = IntMat.n;
    iplType = SciType2IplType(IntMat.it);
    if(iplType==0)
      {
	sciprint("This integer data type is not supported by SIVP. Integer type number: %d. \r\n", IntMat.it);
	goto EXIT_TAG;
      }
    pData = IntMat.D;
    break;
  default:
    sciprint("The data type of %d'th argument is %d. It can't be converted to an image.\r\n",nPos, pDataHeader[0] );
    goto EXIT_TAG;
  }

  /*check dimension */
  if(m3*n3 != nWidth * nHeight * nCh)
    {
      sciprint("Broken hypermatrix: The hypermatrix declares %d X %d X %d, but actually %d elements.\r\n", nHeight, nWidth, nCh, m3*n3);
      goto EXIT_TAG;
    }

  /*create a IplImage to receiving the data*/
  pImg = cvCreateImage(cvSize(nWidth, nHeight), iplType, nCh);
  if(pImg == NULL)
  {
    sciprint("Create IplImage for %d'th argument failed.\r\n", nPos);
    goto EXIT_TAG;
  }

  /*change data order and copy data*/
  MatData2ImgData(pImg, pData );

  
  FreeRhsSVar(pStr);
  return pImg;

 NOT_HM:
  sciprint("The %d'th argument is not a hypermatrix.\r\n", nPos);
  FreeRhsSVar(pStr);
  return NULL;
 EXIT_TAG:
  FreeRhsSVar(pStr);
  return NULL;
  
}

/************************************************************
 * convert SCI matrix to IplImage
************************************************************/
IplImage * Mat2IplImg(int nPos)
{
  IplImage * pImg;

  int mR1, nR1, lR1;
  SciIntMat IntMat;
  int iplType;

  switch(VarType(nPos)){
  case 1: /*real or complex constant matrix.*/
    GetRhsVar(nPos, "d", &mR1, &nR1, &lR1);
    pImg = cvCreateImage(cvSize(nR1, mR1), IPL_DEPTH_64F, 1);
    /*if can not create the IplImage*/
    if(pImg == NULL)
      {
	sciprint("Create IplImage for %d'th argument failed.\r\n", nPos);
	return NULL;
      }
    MatData2ImgData(pImg, stk(lR1) );
    return pImg;
    break;

    /*integer matrix*/
  case 8: 

    GetRhsVar(nPos, "I", &mR1, &nR1, &IntMat);
    iplType = SciType2IplType(IntMat.it);
    if(iplType==0)
      {
	sciprint("This integer data type is not supported by SIVP. Integer type number: %d. \r\n", IntMat.it);
	return NULL;
      }
    pImg = cvCreateImage(cvSize(nR1, mR1),iplType , 1);
    if(pImg == NULL)
      {
	sciprint("Create IplImage for %d'th argument failed.\r\n", nPos);
	return NULL;
      }
    MatData2ImgData(pImg, IntMat.D );
    return pImg;
    break;

  case 17:
    return CreateIplImgFromHm( nPos);
    break;
  default:
    sciprint("This data type can't be converted to an image.\r\n");
    return NULL;
  }
}

/************************************************************
 * change the data order from column-wise to row-wise
 ************************************************************/
int MatData2ImgData(IplImage * pImage, void * pMatData)
{
  //  IPL_DEPTH_8U, IPL_DEPTH_8S, IPL_DEPTH_16U,
  //IPL_DEPTH_16S, IPL_DEPTH_32S, IPL_DEPTH_32F and IPL_DEPTH_64F 
  int row, col, ch;
  long nCount = 0;
  int nBytes;

  char * pDst = NULL;
  char * pSrc = NULL;

  if (pImage == NULL || pMatData == NULL)
    return FALSE;

  pDst = (char*)(pImage->imageData);
  pSrc = (char*)pMatData;

  /*how many bytes per pixel per channel*/
  nBytes = pImage->depth;
  if (nBytes > IPL_DEPTH_SIGN)
    nBytes -= IPL_DEPTH_SIGN;
  nBytes = nBytes >> 3;

  for(ch = 0; ch < pImage->nChannels ; ch++) //the order of IplImage is BGR
    for(col =0; col < pImage->width; col++)
      for(row = 0; row < pImage->height; row++)
	{
	  memcpy(pDst + pImage->widthStep*row + (col*pImage->nChannels + (pImage->nChannels-ch-1))*nBytes, pSrc+nCount, nBytes );
	  nCount += nBytes;
	}


  return TRUE;
}

/************************************************************
 * change the data order from row-wise to column-wise
 ************************************************************/
int ImgData2MatData(IplImage * pImage, void * pMatData)
{
  //  IPL_DEPTH_8U, IPL_DEPTH_8S, IPL_DEPTH_16U,
  //IPL_DEPTH_16S, IPL_DEPTH_32S, IPL_DEPTH_32F and IPL_DEPTH_64F 
  int row, col, ch;
  long nCount = 0;
  int nBytes;

  char * pSrc = NULL;
  char * pDst = NULL;

  if (pImage == NULL || pMatData == NULL)
    return FALSE;

  pSrc = (char*)(pImage->imageData);
  pDst = (char*)pMatData;

  /*how many bytes per pixel per channel*/
  nBytes = pImage->depth;
  if (nBytes > IPL_DEPTH_SIGN)
    nBytes -= IPL_DEPTH_SIGN;
  nBytes = nBytes >> 3;


  for(ch = 0; ch < pImage->nChannels ; ch++) //the order of IplImage is BGR
    for(col =0; col < pImage->width; col++)
      for(row = 0; row < pImage->height; row++)
	{
	  memcpy(pDst+nCount, pSrc + pImage->widthStep*row + (col*pImage->nChannels + (pImage->nChannels-ch-1))*nBytes, nBytes );
	  nCount += nBytes;
	}


  return TRUE;
}

/************************************************************
 * unsigned int32 is not supported (no this type in OpenCV)
 ************************************************************/
int IplType2SciType(int IplType)
{
  switch(IplType)    {
  case IPL_DEPTH_8U:  return I_UCHAR;
  case IPL_DEPTH_8S:  return I_CHAR;
  case IPL_DEPTH_16U: return I_UINT16;
  case IPL_DEPTH_16S: return I_INT16;
  case IPL_DEPTH_32S: return I_INT32;
  case IPL_DEPTH_32F: return SIVP_FLOAT;
  case IPL_DEPTH_64F: return SIVP_DOUBLE;
  default: return 0;
  }
}

/************************************************************
 * unsigned int32 is not supported (no this type in OpenCV)
 ************************************************************/
int SciType2IplType(int SciType)
{
  switch(SciType)    {
  case I_UCHAR:  return IPL_DEPTH_8U;
  case I_CHAR:   return IPL_DEPTH_8S;
  case I_UINT16: return IPL_DEPTH_16U;
  case I_INT16:  return IPL_DEPTH_16S;
  case I_INT32:  return IPL_DEPTH_32S;
  case SIVP_FLOAT:  return IPL_DEPTH_32F;
  case SIVP_DOUBLE: return IPL_DEPTH_64F;
  default: return 0;
  }
}

/* convert data from row-wise to columnwise */
void img2mat(unsigned char* pSrc, unsigned char * pDst, int nWidth, int nHeight, int nCh)
{
  int row, col, ch;
  long nCount = 0;

  for(ch =0; ch < nCh; ch++)
    for(col =0; col < nWidth; col++)
      for(row = 0; row < nHeight; row++)
	{
	  *(pDst+nCount) = *(pSrc+(nWidth*nCh)*row+col*nCh+ch);
	  nCount++;
	}
}

/* convert data from columnwise to row-wise */
void mat2img(unsigned char * pMat, unsigned char *pImg, int nWidth, int nHeight, int nCh)
{
  int row, col, ch;
  long offset;
  long nCount = 0;

  for(row=0; row < nHeight; row++)
    for(col=0; col < nWidth; col++)
      for (ch=0; ch < nCh; ch++)
	{
	  offset = ch*(nWidth*nHeight) + col * nHeight + row;
	  pImg[nCount] = pMat[offset];
	  nCount++;
	}
}

/////add by Jia//////////////////////////////////////////////////////////////////////////////
/***********************************************************
* convert SCI matrix to OpenCV CvMat 
* change the data order from column-wise to row-wise
***********************************************************/ 
CvMat* SciMat2CvMat(int nPos)
{
	CvMat* pMat = NULL;
	int m, n, l;
    int i, j, k;
    
    GetRhsVar(nPos, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &l);
//    double* pSrc = stk(l);
    
    pMat = cvCreateMat(m, n, CV_32FC1);
    k = 0;
    for(i = 0; i < n; i++)
    {
    	for(j = 0; j < m; j++)
    	{
    		//pMat->data.fl[j*n+i] = (float)pSrc[k++];
			pMat->data.fl[j*n+i] = (float)stk(l)[k++];
    	}
    }
        
    return pMat;
}


/***********************************************************
* convert OpenCV CvMat (CV_32FC1 or CV_32SC1) to SCI matrix
* change the data order from row-wise to column-wise
***********************************************************/
void CvMat2SciMat(CvMat* pMat, int nPos)
{
	int rows = pMat->rows;
	int cols = pMat->cols;
	
	if(CV_MAT_TYPE(pMat->type) == CV_32FC1)
	{
		float *pData = pMat->data.fl;
		Create2DFloatMat(nPos, rows, cols, pData);		
	}
	
	if(CV_MAT_TYPE(pMat->type) == CV_32SC1)
	{
		int *pData = pMat->data.i;
		int nType = I_INT32;
		Create2DIntMat(nPos, rows, cols, pData, nType);
	}	
}
/////add by Jia//////////////////////////////////////////////////////////////////////////////



