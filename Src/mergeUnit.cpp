// this class is responsible for merge the temp results files
// by yongzhao
//
//
//


#include <algorithm>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include "mergeUnit.h"
using namespace std;
using namespace boost;



CMergeUnit::CMergeUnit(const char* szFile, int nSize, bool bAln)
{
	m_aSpace[0] = ' ';
	m_aTab[0] = '\t';
	m_sLogE = "(E-value";
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
}


void CMergeUnit::Update(int nID, vector<CSortUnit>& v)
{
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
		double d = 0.0;
		if (true == m_bAln)
		{
			vector<char>::iterator itEValue = search(itSt+1,  itEd, m_sLogE.begin(), m_sLogE.end());
			//vector<char>::iterator itEValue = find(itSt+1,  itEd, '(');
			advance(itEValue, 10);
			d = lexical_cast<double>(string(itEValue, find(itEValue, itEd, *pSep)));
			v.push_back(CSortUnit(d, itSt, itEd));
		}
		else
		{
			vector<char>::iterator itEValue = find_end(itSt+1,  itEd, pSep, pSep+1);
			d = lexical_cast<double>(string(find_end(itSt+1, itEValue-1, pSep, pSep+1)+1, itEValue));
			// get rid of the first '>'
			v.push_back(CSortUnit(d, itSt+1, itEd));
		}

		itSt = itEd;
		itEd = find(itSt+1, m_itEnd, '>');
	}

}

int CMergeUnit::GetLast()
{
	return m_vIdx.size();
}
