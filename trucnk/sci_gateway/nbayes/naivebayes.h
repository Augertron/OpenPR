/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

#ifndef __ML_H__
#define __ML_H__

// disable deprecation warning which appears in VisualStudio 8.0
#if _MSC_VER >= 1400
#pragma warning( disable : 4996 )
#endif

#ifndef SKIP_INCLUDES

  #include "cxcore.h"
  #include <limits.h>

  #if defined WIN32 || defined WIN64
    #include <windows.h>
  #endif

#else // SKIP_INCLUDES

  #if defined WIN32 || defined WIN64
    #define CV_CDECL __cdecl
    #define CV_STDCALL __stdcall
  #else
    #define CV_CDECL
    #define CV_STDCALL
  #endif

  #ifndef CV_EXTERN_C
    #ifdef __cplusplus
      #define CV_EXTERN_C extern "C"
      #define CV_DEFAULT(val) = val
    #else
      #define CV_EXTERN_C
      #define CV_DEFAULT(val)
    #endif
  #endif

  #ifndef CV_EXTERN_C_FUNCPTR
    #ifdef __cplusplus
      #define CV_EXTERN_C_FUNCPTR(x) extern "C" { typedef x; }
    #else
      #define CV_EXTERN_C_FUNCPTR(x) typedef x
    #endif
  #endif

  #ifndef CV_INLINE
    #if defined __cplusplus
      #define CV_INLINE inline
    #elif (defined WIN32 || defined WIN64) && !defined __GNUC__
      #define CV_INLINE __inline
    #else
      #define CV_INLINE static
    #endif
  #endif /* CV_INLINE */

  #if (defined WIN32 || defined WIN64) && defined CVAPI_EXPORTS
    #define CV_EXPORTS __declspec(dllexport)
  #else
    #define CV_EXPORTS
  #endif

  #ifndef CVAPI
    #define CVAPI(rettype) CV_EXTERN_C CV_EXPORTS rettype CV_CDECL
  #endif

#endif // SKIP_INCLUDES


#ifdef __cplusplus

// Apple defines a check() macro somewhere in the debug headers
// that interferes with a method definiton in this header
#undef check

/****************************************************************************************\
*                               Main struct definitions                                  *
\****************************************************************************************/

/* log(2*PI) */
#define CV_LOG2PI (1.8378770664093454835606594728112)

/* columns of <trainData> matrix are training samples */
#define CV_COL_SAMPLE 0

/* rows of <trainData> matrix are training samples */
#define CV_ROW_SAMPLE 1

#define CV_IS_ROW_SAMPLE(flags) ((flags) & CV_ROW_SAMPLE)

struct CvVectors
{
    int type;
    int dims, count;
    CvVectors* next;
    union
    {
        uchar** ptr;
        float** fl;
        double** db;
    } data;
};

#if 0
/* A structure, representing the lattice range of statmodel parameters.
   It is used for optimizing statmodel parameters by cross-validation method.
   The lattice is logarithmic, so <step> must be greater then 1. */
typedef struct CvParamLattice
{
    double min_val;
    double max_val;
    double step;
}
CvParamLattice;

CV_INLINE CvParamLattice cvParamLattice( double min_val, double max_val,
                                         double log_step )
{
    CvParamLattice pl;
    pl.min_val = MIN( min_val, max_val );
    pl.max_val = MAX( min_val, max_val );
    pl.step = MAX( log_step, 1. );
    return pl;
}

CV_INLINE CvParamLattice cvDefaultParamLattice( void )
{
    CvParamLattice pl = {0,0,0};
    return pl;
}
#endif

/* Variable type */
#define CV_VAR_NUMERICAL    0
#define CV_VAR_ORDERED      0
#define CV_VAR_CATEGORICAL  1

#define CV_TYPE_NAME_ML_SVM         "opencv-ml-svm"
#define CV_TYPE_NAME_ML_KNN         "opencv-ml-knn"
#define CV_TYPE_NAME_ML_NBAYES      "opencv-ml-bayesian"
#define CV_TYPE_NAME_ML_EM          "opencv-ml-em"
#define CV_TYPE_NAME_ML_BOOSTING    "opencv-ml-boost-tree"
#define CV_TYPE_NAME_ML_TREE        "opencv-ml-tree"
#define CV_TYPE_NAME_ML_ANN_MLP     "opencv-ml-ann-mlp"
#define CV_TYPE_NAME_ML_CNN         "opencv-ml-cnn"
#define CV_TYPE_NAME_ML_RTREES      "opencv-ml-random-trees"

class CV_EXPORTS CvStatModel
{
public:
    CvStatModel();
    virtual ~CvStatModel();

    virtual void clear();

    virtual void save( const char* filename, const char* name=0 );
    virtual void load( const char* filename, const char* name=0 );

    virtual void write( CvFileStorage* storage, const char* name );
    virtual void read( CvFileStorage* storage, CvFileNode* node );

protected:
    const char* default_model_name;
};


/****************************************************************************************\
*                                 Normal Bayes Classifier                                *
\****************************************************************************************/

/* The structure, representing the grid range of statmodel parameters.
   It is used for optimizing statmodel accuracy by varying model parameters,
   the accuracy estimate being computed by cross-validation.
   The grid is logarithmic, so <step> must be greater then 1. */
struct CV_EXPORTS CvParamGrid
{
    // SVM params type
    enum { SVM_C=0, SVM_GAMMA=1, SVM_P=2, SVM_NU=3, SVM_COEF=4, SVM_DEGREE=5 };

    CvParamGrid()
    {
        min_val = max_val = step = 0;
    }

    CvParamGrid( double _min_val, double _max_val, double log_step )
    {
        min_val = _min_val;
        max_val = _max_val;
        step = log_step;
    }
    //CvParamGrid( int param_id );
    bool check() const;

    double min_val;
    double max_val;
    double step;
};

class CV_EXPORTS CvNormalBayesClassifier : public CvStatModel
{
public:
    CvNormalBayesClassifier();
    virtual ~CvNormalBayesClassifier();

    CvNormalBayesClassifier( const CvMat* _train_data, const CvMat* _responses,
        const CvMat* _var_idx=0, const CvMat* _sample_idx=0 );

    virtual bool train( const CvMat* _train_data, const CvMat* _responses,
        const CvMat* _var_idx = 0, const CvMat* _sample_idx=0, bool update=false );

    virtual float predict( const CvMat* _samples, CvMat* results=0 ) const;
    virtual void clear();

    virtual void write( CvFileStorage* storage, const char* name );
    virtual void read( CvFileStorage* storage, CvFileNode* node );

//protected:
public:
    int     var_count, var_all;
    CvMat*  var_idx;
    CvMat*  cls_labels;
    CvMat** count;
    CvMat** sum;
    CvMat** productsum;
    CvMat** avg;
    CvMat** inv_eigen_values;
    CvMat** cov_rotate_mats;
    CvMat*  c;
};

#endif

#endif /*__ML_H__*/
/* End of file. */
