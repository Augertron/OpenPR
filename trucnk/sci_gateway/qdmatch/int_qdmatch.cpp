#include "SparseMatch.h"
#include "QuasiMatch.h"

#include "transfer.h"

#include "common.h"

extern "C" {



	#include "sciprint.h"
	#include "Scierror.h"

	void chIpl2Img(Image& img, IplImage* iplImg);
	const int MAX_IMG_SIZE = 2000;

	int int_qdmatch(char *fname)
	{
		if(Rhs==0)
		{
			sciprint("Usage:\nmatchAB=qdmatch(imA,imB)\n[matchBC,matchABC]=qdmatch(imB,imC,matchAB)\n");
			return -1;	
		}

		CheckRhs(2, 3);
		CheckLhs(1, 2);
		
		if((Rhs==2)&&(Lhs==2))
		{
	     	Scierror(999, "%s: numbers of input and output arguments do not match.\r\n", fname);
			return -1;
		}
		if((Rhs==3)&&(Lhs==1))
		{
	     	Scierror(999, "%s: numbers of input and output arguments do not match.\r\n", fname);
			return -1;
		}
		
		/////////////////////////////////////////////////////////////////////////
		IplImage *pImgL = Mat2IplImg(1);
		IplImage *pImgR = Mat2IplImg(2);
		
		char *matchABFile;
		char *matchBCFile;
		char *matchABCFile;
		
		Image leftImg;
		Image rightImg;	
		//foundamental matrix
		double fundF[3][3];
		vector<Match> com3Match;
		vector<Match> *pmatchAB;
		vector<Match> matchAB;
		vector<Match> matchBC;	
		
		if(Rhs == 2)
		{
			matchABFile = "matchABFile";
			pmatchAB = 0;
		}
		else
		{
			matchBCFile = "matchBCFile";	
			matchABCFile = "matchABCFile";	
			
			SciMat2VecMatch(3, matchAB);		
			pmatchAB = &matchAB;				
		}	
		
		if((pImgL == NULL) || (pImgR == NULL))
		{
			Scierror(999, "%s: loading images failed.\r\n", fname);
			return -1;
		}		
		
		//the input images should be gray images
		if((pImgL->nChannels != 1) || (pImgR->nChannels != 1))
		{
			Scierror(999, "%s: input images should be gray images.\r\n", fname);
			return -1;	
		}
		
		//the input images should have the same size
		if((pImgL->width != pImgR->width) || (pImgL->height != pImgR->height))
		{
			Scierror(999, "%s: input images should have the same size.\r\n", fname);
			return -1;
		}
		
		//if input images are too large, scale them to a proper size
		int maxsize = pImgL->width>pImgL->height?pImgL->width:pImgL->height;
		if (maxsize > MAX_IMG_SIZE)
		{
			double scale = double(maxsize)/MAX_IMG_SIZE;
			IplImage *tmpImg;
			if (maxsize == pImgL->width)
			{
			//	cout << "Image resized to " << int(pImgL->height/scale)<< "x" << MAX_IMG_SIZE << endl;
				tmpImg = cvCreateImage(cvSize(MAX_IMG_SIZE, int(pImgL->height/scale)), IPL_DEPTH_8U, 1);
			}
			else
			{
			//	cout << "Image resized to " << MAX_IMG_SIZE<< "x" << int(pImgL->width/scale) << endl;
				tmpImg = cvCreateImage(cvSize(int(pImgL->width/scale), MAX_IMG_SIZE), IPL_DEPTH_8U, 1);	
			}
			
			cvResize(pImgL, tmpImg);
			//change IplImage to user defined format Image
			chIpl2Img(leftImg, tmpImg);

			cvResize(pImgR, tmpImg);
			chIpl2Img(rightImg, tmpImg);

			//release tmpImg
			cvReleaseImage(&tmpImg);
		}
		else
		{
			chIpl2Img(leftImg, pImgL);
			chIpl2Img(rightImg, pImgR);
		}
		
		cvReleaseImage(&pImgL);
		cvReleaseImage(&pImgR);
		
		int matchNum = SparseMatch(leftImg, rightImg, com3Match, fundF);
		if (matchNum < 10)
		{
			// matching error	
			Scierror(999, "%s: matching error.\r\n", fname);
			return -1;    
		}	
		sciprint("Sparse matching is over ...\n");

		//quasi-dense propagating
		//matrix fundM;
		//matchNum = QuasiMatch(leftImg, rightImg, com3Match, fundF, pmatch01, match12File);
		matchNum = QuasiMatch(leftImg, rightImg, com3Match, fundF, pmatchAB, matchBCFile);
		sciprint("Quasi matching is over ...\n");
		if (matchNum < 100)
		{
			// propagating error
			Scierror(999, "%s: propagating error.\r\n", fname);
			return -1;    
		}

		if (pmatchAB == 0)  
		{
			// write out the matches
			//WriteMatchFile(matchABFile, com3Match);
			//matchAB
			VecMatch2SciMat(com3Match, Rhs+1);
		}
		else
		{
			// match file contains matches between 3 images: image 0, 1 and 2
		//	WriteMatchFile(match012File, pmatch01, com3Match);
			//matchBC		
			ReadMatchFile(matchBC, matchBCFile);
			remove(matchBCFile);
			VecMatch2SciMat(matchBC, Rhs+1);		
		
			//WriteMatchFile(matchABCFile, pmatchAB, com3Match);
		
			//matchABC
			if(pmatchAB->size() != com3Match.size())
			{
				Scierror(999, "%s: The matches in (*pmatchAB) and com3Match should be of the same size!\r\n", fname);
				return -1;	
			}
			matchNum = com3Match.size();
			int width = 6;
			double *pmatch = new double[matchNum*width];
			
			int i;
			int m, n;
			
			//column-wise
			for(i = 0; i < matchNum; i++)
			{
				Match& refMatchAB = (*pmatchAB)[i];
				const Match& refMatchBC = com3Match[i];
						
				pmatch[i] = refMatchAB.GetTarget().GetY();
				pmatch[i+matchNum] = refMatchAB.GetTarget().GetX();
				pmatch[i+(2*matchNum)] = refMatchAB.GetFound().GetY();
				pmatch[i+(3*matchNum)] = refMatchAB.GetFound().GetX();
				pmatch[i+(4*matchNum)] = refMatchBC.GetFound().GetY();
				pmatch[i+(5*matchNum)] = refMatchBC.GetFound().GetX();	
			}
			
			m = matchNum;
			n = width;
			CreateVarFromPtr(Rhs+2, MATRIX_OF_DOUBLE_DATATYPE, &m, &n, &pmatch);
			
			delete []pmatch;			
		}	
		
		LhsVar(1) = Rhs+1;
		if(Rhs == 3)
		{
			LhsVar(2) = Rhs+2;
		}
				
		return 0;
	}

	void chIpl2Img(Image& img, IplImage* iplImg)
	{
		img = Image(iplImg->width, iplImg->height);
		for (int i = 0; i < iplImg->width; i++)
		{
			for (int j = 0; j < iplImg->height; j++)
			{
				img(i, j) = ((unsigned char)(iplImg->imageData[j*iplImg->widthStep+i]))/255.0;
			}
		}
	}	

}
