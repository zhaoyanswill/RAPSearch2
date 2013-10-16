// this class is responsible for merge the temp results files
// by yongzhao
//
//
//


#include "mergeUnit.h"
#include <algorithm>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
using namespace boost;


CMergeUnit::CMergeUnit(const char* szFile, int nSize, bool bAln)
{
	m_aSpace[0] = ' ';
	m_aTab[0] = '\t';
	m_aIdicator[0] = '>';
	// indicate that the file is processed
	m_sMax.assign("~");

	m_bDone = false;
	m_bAln = bAln;
	m_nSize = nSize;
	m_nSize = 100000000;
	m_nLeft = 0;

	m_ifFile.open(szFile);
	if (!m_ifFile.good())
	{
		cout << "can not open temp file..." << endl;
	}
	
	m_sFile.assign(szFile);
	m_vPool.assign(m_nSize, '\0');

	ifstream ifIdx((string(szFile)+".idx").c_str());
	archive::binary_iarchive ia(ifIdx);
	ia >> m_vIdx;
	ifIdx.close();
	remove((string(szFile)+".idx").c_str());
	m_nCur = 0;
}


string CMergeUnit::Update(vector<CSortUnit>& v)
{
	if (!m_ifFile.good())
	{
		return m_sMax;
	}

	char* pSep = NULL;
	if (true == m_bAln)
	{
		pSep = m_aSpace;
	}
	else
	{
		pSep = m_aTab;
	}

	m_ifFile.read(&m_vPool[m_nLeft], m_nSize-m_nLeft);
	int nRead = m_ifFile.gcount();
	m_itEnd = m_vPool.begin();
	if (m_nLeft+nRead < m_nSize)
	{
		m_vPool[m_nLeft+nRead] = '>';
		advance(m_itEnd, m_nLeft+nRead+1);
	}
	else
	{
		advance(m_itEnd, m_nLeft+nRead);
	}

	vector<char>::iterator itSt = find(m_vPool.begin(), m_itEnd, '>');
	vector<char>::iterator itEd = find(itSt+1, m_itEnd, '>');
	while (itEd != m_itEnd)
	{
		string sCur(itSt, find(itSt+1, itEd, *pSep));

		double d = 0.0;
		if (true == m_bAln)
		{
			vector<char>::iterator itEValue = find(itSt+1,  itEd, '(');
			advance(itEValue, 10);
			d = lexical_cast<double>(string(itEValue, find(itEValue, itEd, *pSep)));
			v.push_back(CSortUnit(sCur, d, string(itSt, itEd)));
		}
		else
		{
			vector<char>::iterator itEValue = find_end(itSt+1,  itEd, pSep, pSep+1);
			d = lexical_cast<double>(string(find_end(itSt+1, itEValue-1, pSep, pSep+1)+1, itEValue));
			// get rid of the first '>'
			v.push_back(CSortUnit(sCur, d, string(itSt+1, itEd)));
		}

		itSt = itEd;
		itEd = find(itSt+1, m_itEnd, '>');
	}

	vector<char>::iterator it = find_end(m_vPool.begin(), itSt-1, m_aIdicator, m_aIdicator+1);
	string sLast(it, find(it+1, itSt, ' '));

	copy(itSt, itEd, m_vPool.begin());
	m_nLeft = distance(itSt, itEd);

	if (m_nLeft+nRead < m_nSize)
	{
		return m_sMax;
	}
	else
	{
		return sLast;
	}
}



void CMergeUnit::Update(int nID, vector<CSortUnit>& v)
{
	/*
	if (m_nCur >= m_vIdx.size())
	{
		return;
	}
	*/

	if (m_vIdx[nID].m_nSize == 0)
	{
		return;
	}

	char* pSep = NULL;
	if (true == m_bAln)
	{
		pSep = m_aSpace;
	}
	else
	{
		pSep = m_aTab;
	}

	m_ifFile.seekg(m_vIdx[nID].m_llBeg, ios::beg);
	m_ifFile.read(&m_vPool[0], m_vIdx[nID].m_nSize);

	int nRead = m_ifFile.gcount();
	m_vPool[m_nLeft+nRead] = '>';
	m_itEnd = m_vPool.begin();
	advance(m_itEnd, nRead+1);

	vector<char>::iterator itSt = find(m_vPool.begin(), m_itEnd, '>');
	vector<char>::iterator itEd = find(itSt+1, m_itEnd, '>');
	while (itEd != m_itEnd)
	{
		string sCur(itSt, find(itSt+1, itEd, *pSep));

		double d = 0.0;
		if (true == m_bAln)
		{
			vector<char>::iterator itEValue = find(itSt+1,  itEd, '(');
			advance(itEValue, 10);
			d = lexical_cast<double>(string(itEValue, find(itEValue, itEd, *pSep)));
			v.push_back(CSortUnit(sCur, d, string(itSt, itEd)));
		}
		else
		{
			vector<char>::iterator itEValue = find_end(itSt+1,  itEd, pSep, pSep+1);
			d = lexical_cast<double>(string(find_end(itSt+1, itEValue-1, pSep, pSep+1)+1, itEValue));
			// get rid of the first '>'
			v.push_back(CSortUnit(sCur, d, string(itSt+1, itEd)));
		}

		itSt = itEd;
		itEd = find(itSt+1, m_itEnd, '>');
	}

	//++m_nCur;
}

int CMergeUnit::GetLast()
{
	return m_vIdx.size();
}
