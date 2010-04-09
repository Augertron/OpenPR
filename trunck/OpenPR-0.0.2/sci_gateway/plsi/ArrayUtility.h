//author: Mingbo Wang

#pragma once

#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;

template <class T> class ArrayUtility
{
public:
	ArrayUtility(void);
	~ArrayUtility(void);
	static double * initialArray(int );
	static T * initialArray(int ,T );
	static T ** initial2DArray(int , int );
	static T ** initial2DArray(int , int ,T);
	static void show(T * p, int size);
	static void show(T ** p, int r,int c);
	static T*** initial3DArray(int,int,int,T);
	static bool findInArray(T *,T value,int );
	static void setValue(T *,T value,int size);
	static void setValue(T **,T value,int r,int c);
	static void finalize2DArray(T**,int r);
	static void show(ostream &out,T*,int);
	static void show2D(ostream &out,T**,int,int);
};

template <class T>
void ArrayUtility<T>::show( ostream &out,T * arr,int size )
{
	for(int i=0;i<size-1;i++)
	{
		out<<arr[i]<<" ";
	}
	if(size>0)
		out<<arr[size-1];
	out<<endl;
}

template <class T>
void ArrayUtility<T>::show2D( ostream &out,T ** arr,int r,int c )
{
	for(int i=0;i<r;i++)
		show(out,arr[i],c);
}

template <class T> ArrayUtility< T>::ArrayUtility(void)
{
}

template <class T> ArrayUtility< T>::~ArrayUtility(void)
{
}

template <class T>  T ** ArrayUtility< T>::initial2DArray(int r, int w)
{
	T ** p=new T*[r];
	for(int i=0;i<r;i++)
		p[i]=ArrayUtility<T>::initialArray(w);
	return p;
}

template <class T>  T ** ArrayUtility< T>::initial2DArray(int r, int w,T value)
{
	T ** p=new T*[r];
	for(int i=0;i<r;i++)
		p[i]=ArrayUtility<T>::initialArray(w,value);
	return p;
}

template <class T>  double * ArrayUtility< T>::initialArray(int r)
{
	double *p=new double[r];
	return p;
}

template <class T> T * ArrayUtility< T>::initialArray(int r,T value)
{
	T *p=new T[r];
	for(int i=0;i<r;i++)
		p[i]=value;
	return p;
}

template <class T> T *** ArrayUtility< T>::initial3DArray(int r,int c,int z,T value)
{
	T ***p=new T**[r];
	for(int i=0;i<r;i++)
		p[i]=initial2DArray(c,z,value);
	return p;
}

template<typename T> void ArrayUtility<T>::show(T * p, int size)
{
	for(int i=0;i<size;i++)
		cout<<p[i]<<" ";
	cout<<endl;
}

template<typename T>  void ArrayUtility<T>::show(T ** p, int r,int c)
{
	for(int i=0;i<r;i++)
		show(p[i],c);
}

template<typename T> bool ArrayUtility<T>::findInArray(T * data,T value,int size){
	for(int i=0;i<size;i++)
	{
		if(data[i]==value)
			return true;
	}
	return false;
}

template<typename T> void ArrayUtility<T>::setValue(T * p,T value,int size)
{
	for(int i=0;i<size;i++)
		p[i]=value;
}

template<typename T> void ArrayUtility<T>::setValue(T ** p,T value,int r,int c)
{
	for(int i=0;i<r;i++)
		setValue(p[i],value,c);
}

template<typename T> void ArrayUtility<T>::finalize2DArray(T ** p,int r)
{
	for(int i=0;i<r;i++)
		delete []p[i];
	delete []p;
}
