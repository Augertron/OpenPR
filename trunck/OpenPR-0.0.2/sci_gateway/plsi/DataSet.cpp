//author: Mingbo Wang
//revised by: Jia Wu

#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdlib>
#include <cmath>

#include "DataSet.h"

using namespace std;

DataSet::DataSet(void): docnum(0),wordnum(0),data(NULL)
{
}

DataSet::DataSet(int dnum, int wnum)
{
	docnum=dnum;
	wordnum=wnum;
        data=ArrayUtility<int>::initial2DArray(docnum,wordnum,0);
}
DataSet::~DataSet(void)
{
}

int DataSet::getCount(int d, int w)
{
	if(d>=docnum||w>=wordnum)
		return -1;
	return data[d][w];
}

int DataSet::getDocNum()
{
	return docnum;
}

int DataSet::getWordNum()
{
	return wordnum;
}

/*
docnum
wordnum
wordid:count
*/
bool DataSet::initData(string s)
{
        ifstream instream(s.c_str());
        if(instream==NULL)
        {
        	cout<<"file not find"<<endl;
		return false;
        }

        string content;// a line of the file
        string w;

        int wordid;
        int docid;
        int index;
        int count;

	getline(instream,content);
	wordnum=atoi(content.c_str());
	cout<<"wordnum is "<<wordnum<<endl;
	getline(instream,content);
	docnum=atoi(content.c_str());
	cout<<"docnum is "<<docnum<<endl;

	data=ArrayUtility<int>::initial2DArray(docnum,wordnum,0);
	while(getline(instream,content))
	{
		//cout<<1<<endl;
		istringstream ins(content);
		ins>>w;
		docid=atoi(w.c_str());

		while(ins>>w)
		{
			index=(int)w.find(':');
			string s1=w.substr(0,index);
			string s2=w.substr(index+1,w.size()-index);

			wordid=atoi(s1.c_str());
			count=atoi(s2.c_str());
			//cout<<wordid<<" "<<count<<endl;
			if(docid<docnum&&wordid<wordnum)
				data[docid][wordid]=count;
			else
			{
				cout<<"data error!"<<endl;
				return false;
			}
		}

	}
	
	cout<<"data initialized"<<endl;
	
	return true;
}

void DataSet::fold(int num,int ** data)
{
}

void DataSet::showData(void)
{
	ArrayUtility<int>::show(data,docnum,wordnum);
}

void DataSet::setValue(int r, int w, int value)
{
	if(r<docnum&&r>=0&&w<docnum&&w>=0)
		data[r][w]=value;
}

/*
docnum
wordnum
wordcount
pure vector
*/
bool DataSet::initDataF(string path)
{
	ifstream instream(path.c_str());
	if(instream==NULL)
	{		
		cout<<"file not find"<<endl;
		return false;
	}
	
	string content;// a line of the file
	string w;
	
	int wordid;
	int docid=0;

	int count;
	getline(instream,content);
	wordnum=atoi(content.c_str());
	cout<<"wordnum is: "<<wordnum<<endl;
	getline(instream,content);
	docnum=atoi(content.c_str());
	cout<<"docnum is: "<<docnum<<endl;	
	
	data=ArrayUtility<int>::initial2DArray(docnum,wordnum,0);
	while(getline(instream,content))
	{
		wordid=0;
		//cout<<1<<endl;
		istringstream ins(content);	
		
		while(ins>>w)
		{		
			count=atoi(w.c_str());
			//cout<<wordid<<" "<<count<<endl;
			if(docid<docnum&&wordid<wordnum)
			{
				//cout<<docid<<wordid<<endl;
				data[docid][wordid++]=count;
			}
			else
			{
				cout<<"data error!"<<endl;
				return false;
			}
		}
		docid++;
	}
	cout<<"data initialized"<<endl;
	return true;
}

bool DataSet::isEmpty()
{
	if(wordnum==0||docnum==0)
		return true;
	return false;
}

bool DataSet::clear()
{
	if(!isEmpty())
	{
		ArrayUtility<int>::finalize2DArray(data,docnum);
		wordnum=0;
		docnum=0;
		return true;
	}
	
	return false;
}
