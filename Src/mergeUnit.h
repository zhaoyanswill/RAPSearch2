#ifndef __MERGEUNIT_H__
#define __MERGEUNIT_H__

// this class is used to read each temp file and combine the result

#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "hitUnit.h"
#include "cindex.h"



class CMergeUnit
{
public:
	// open, alloc, read
	CMergeUnit(const char* szFile);

	// close
	~CMergeUnit()
	{
		m_ifFile.close();
		remove(m_sFile.c_str());
		remove((m_sFile+".idx").c_str());
	}

	// insert all hits <= sComplete into v, read next data block
	void Update(int nID, std::vector<CHitUnit>& v);
	
	// get the ID of last element in this file
	int GetLast();


private:
	std::ifstream m_ifFile;
	std::string m_sFile;
	std::vector<CIndex> m_vIdx;
};

#endif
