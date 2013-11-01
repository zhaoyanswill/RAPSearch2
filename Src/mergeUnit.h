#ifndef __MERGEUNIT_H__
#define __MERGEUNIT_H__

// this class is used to read each temp file and combine the result

#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include "sortUnit.h"
#include "cindex.h"



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
	void Update(int nID, std::vector<CSortUnit>& v);
	
	// get the ID of last element in this file
	int GetLast();


private:
	std::ifstream m_ifFile;
	std::vector<char> m_vPool;
	int m_nSize;
	int m_nLeft;
	// point to the next query string
	std::vector<char>::iterator m_itEnd;
	std::string m_sFile;
	bool m_bAln;

	std::vector<CIndex> m_vIdx;

	char m_aSpace[1];
	char m_aTab[1];
	std::string m_sLogE;
};

#endif
