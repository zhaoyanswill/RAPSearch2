#ifndef __MERGEUNIT_H__
#define __MERGEUNIT_H__

// this class is used to read each temp file and combine the result

#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include "sortUnit.h"
#include "cindex.h"

using namespace std;



class CMergeUnit
{
public:
	// open, alloc, read
	CMergeUnit(const char* szFile, int nSize, bool bAln);

	// close
	~CMergeUnit()
	{
		m_ifFile.close();
		remove(m_sFile.c_str());
	}

	// insert all hits <= sComplete into v, read next data block
	string Update(vector<CSortUnit>& v);
	void Update(int nID, vector<CSortUnit>& v);
	
	// get the ID of last element in this file
	int GetLast();


private:
	ifstream m_ifFile;
	vector<char> m_vPool;
	int m_nSize;
	int m_nLeft;
	string m_sMax;
	// point to the next query string
	vector<char>::iterator m_itEnd;
	string m_sFile;
	bool m_bAln;
	bool m_bDone;
	char m_aSpace[1];
	char m_aTab[1];
	char m_aIdicator[1];

	vector<CIndex> m_vIdx;
	int m_nCur;
};

#endif
