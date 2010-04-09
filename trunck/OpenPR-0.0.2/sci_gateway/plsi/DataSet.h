//author: Mingbo Wang
//revised by: Jia Wu

#pragma once

#include <string>

#include "ArrayUtility.h"

class DataSet
{
public:
	DataSet(void);
	~DataSet(void);
	int getCount(int d, int w);
	bool initData(string);
	int** data;
        int getDocNum();
	int getWordNum();
	void fold(int num,int **);
private:
	int wordnum;
	int docnum;
public:
	bool initDataF(string);
	DataSet(int,int);
	void setValue(int,int,int);
	void showData(void);
	bool isEmpty();
	bool clear();
};
