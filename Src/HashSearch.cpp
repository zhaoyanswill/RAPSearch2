
#include "HashSearch.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <climits>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/chrono/thread_clock.hpp>
#include "threadpool.hpp"
#include "weight.h"
#include "aa.h"
#include "n2a.h"
#include "mergeUnit.h"
using namespace std;
using namespace boost;
using namespace boost::threadpool;


boost::mutex muMonitor;

const ushort ONEBYTE = 15;
const ushort TWOBYTE = 255;
const ushort THRBYTE = 4095;
const ushort FOUBYTE = 65535;



CHashSearch::CHashSearch(int nThreadNum)
{
	// for any letter which is not in the 20 aa
	m_uMask = 10;
	m_uSeg = 8;
	fill_n(m_aChar2Code, 256, (m_uMask<<4));
	fill_n(m_aCode2Char, 256, m_uMask);
	fill_n(m_aCode2Ten, 256, m_uMask);
	// read group info in aa.h and build mapping array
	// defaultly use murphy10s
	for (int i = 0; i < 500; ++i)
	{
		if ('\0' == murphy10s[i][0])
		{
			break;
		}

		char* p = murphy10s[i];
		for (uint j = 0; j < strlen(p); ++j) 
		{
			// use new format to fix alignment break in SEGed region
			// the first four bits: group id
			// then one bit: SEGed? 1 : 0
			// then three bits: offset in a group
			uint unIdx = (i << 4) + j+1;
			m_aChar2Code[p[j]] = unIdx;
			m_aChar2Code[p[j]+32] = unIdx;	// add lower-case character
			m_aCode2Char[unIdx] = p[j];
			m_aCode2Ten[unIdx] = i;

			unIdx |= m_uSeg;	// set SEGed bit
			m_aChar2Code[p[j]+128] = unIdx;	// SEGed letter
		}
	}

	// build substitute matrix for compressed char
	fill_n((int*)m_aSubMatrix, 256*256, -5);
	for (uint i = 0; i < strlen(aAAAlph); ++i)
	{
		for (uint j = 0; j < strlen(aAAAlph); ++j)
		{
			if (i < 20 && j < 20)
			{
				m_aSubMatrix[m_aChar2Code[aAAAlph[i]]][m_aChar2Code[aAAAlph[j]]] = blosum62[i][j];
			}
		}
	}
	m_aChar2Code['.'] = (10 << 4);
	m_aCode2Char[(10<<4)] = '.';
	m_aCode2Ten[(10<<4)] = m_uMask;
	m_aCode2Char['-'] = '-';

	m_pBlastSig = NULL;
	m_pComptor = NULL;

	int LONGQUERY = 4096;
	if (0 == nThreadNum)
	{
		m_nThreadNum = 1;
	}
	else
	{
		m_nThreadNum = nThreadNum;
	}
	m_vTrace.assign(m_nThreadNum, vector<vector<char> >(LONGQUERY, vector<char>(LONGQUERY)));
	m_vETrace.assign(m_nThreadNum, vector<vector<char> >(LONGQUERY, vector<char>(LONGQUERY)));
	m_vDTrace.assign(m_nThreadNum, vector<vector<char> >(LONGQUERY, vector<char>(LONGQUERY)));
	m_vBlastPt.assign(m_nThreadNum, -1);

	m_unTotalSeeds = 0;
	m_unTotalQuery = 0;
	m_unTotalSubj = 0;

	m_bSeqType = false;

	// used for convert index from with-fram to non-frame
	m_nIdxScl = 1;
	m_nQueryType = 0;

	m_sOutBase = "";
	m_sOutput = "";
	m_sOutput.reserve(100000000);
	m_sM8 = "";
	//m_sM8.reserve(50000000);
	m_llOutCum = 0;
	m_llM8Cum = 0;
	m_nSeqBase = 0;
	
	// for test on gap extension
	m_unGapExt = 0;

	// hssp
	m_vCriteria.assign(100, 0);
	for (int i = 1; i < 100; ++i)
	{
		float f = 290.15 * pow(i, -0.562);
		f = f * i / 100;
		m_vCriteria[i] = (int)ceil(f);

	}

	m_unXmlSp = 0;
	m_unXmlCnt = 1;
	m_lnSeqNum = 0;
	m_lnTotalAa = 0;
	m_nStdout = 0;
	m_sQFile = "";
	m_sDFile = "";
	m_sStartTime = "";
	m_sLeft = "";
}


struct CompDbObj
{
	CompDbObj(VUCHAR& vSeqs, VUINT& vLens, uint& unMer) : m_vSeqs(vSeqs), m_vLens(vLens), m_unMer(unMer) {}
	bool operator() (const uint pos1, const uint pos2) const
	{
		// init paras
		int nIdx1 = pos1>>11;
		int nLen1 = m_vLens[nIdx1+1] - m_vLens[nIdx1];
		int nOff1 = (pos1&0x7ff) + m_unMer;
		uchar* p1 = &m_vSeqs[m_vLens[nIdx1]] + nOff1;

		int nIdx2 = pos2>>11;
		int nLen2 = m_vLens[nIdx2+1] - m_vLens[nIdx2];
		int nOff2 = (pos2&0x7ff) + m_unMer;
		uchar* p2 = &m_vSeqs[m_vLens[nIdx2]] + nOff2;

		// comp
		int nDiff = 0;
		if (nLen1-nOff1 >= 4 && nLen2-nOff2 >= 4)
		{
			nDiff = 4;
		}
		else
		{
			nDiff = (nLen1-nOff1) < (nLen2-nOff2) ? (nLen1-nOff1) : (nLen2-nOff2);
		}
		for (int i = 0; i < nDiff; ++i)
		{
			if ((*(p1+i)>>4) != (*(p2+i)>>4))
			{
				return (*(p1+i)>>4) < (*(p2+i)>>4);
			}
		}
		return ((nLen1-nOff1) < (nLen2-nOff2));
	}

	VUCHAR& m_vSeqs;
	VUINT& m_vLens;
	uint& m_unMer;
};


int CHashSearch::BuildDHash(const char* szFile, string& sOutFile, int nSplitNum, bool bFullId)
{
	/***************************************************************/
	// revise the size of database according to rapsearch
    long int lnSeqNum = 0;
	long int lnAaNum = 0;
    GuessTotSeq(szFile, lnSeqNum, lnAaNum);

	ifstream is(szFile);
	is.seekg(0, ios::end);
	long int lnFileSize = is.tellg();
	is.close();
	long int lnBlockSize = (long int)((1<<30) * 0.618);
	uint unBlockSize = m_unDSize = lnBlockSize;

	if (0 != nSplitNum)
	{
		unBlockSize = m_unDSize = lnFileSize/nSplitNum + 1;
	}

	/***************************************************************/
		
	// the para for seed variants in both fast and slow mode
	m_bFast = true;
    ifstream ifFile(szFile);
    if (!ifFile.good())
    {
        ifFile.close();
        cout << "can not open the file: " << szFile << endl;
        return -1;
    }

	ofstream of(sOutFile.c_str());
	ofstream ofInfo((sOutFile+".info").c_str());
	if (!of.good() || !ofInfo.good())
	{
		cout << "can not write files..." << endl;
		return -1;
	}
	archive::binary_oarchive oa(of);
	archive::binary_oarchive oaInfo(ofInfo);
	
	// container for db para
	MINDEX vHash(m_unTotalIdx, VUINT()); // all k-mer of database
	VUINT vLens;
	VUCHAR vSeqs;
	VNAMES vNames;
	vector<double> vFreq(strlen(murphy10r), 0);
	VUINT vWordCnts(m_unTotalIdx, 0);
	uint unMedian = 0;
	long int lnTotalAa = 0;

	// data pool for processing data file
	POOL vPool(m_unDSize, 0);
	int nBlock = 0;
	uint unLeft = 0;
    while (ifFile.good())
    {
		ifFile.read(&vPool[unLeft], m_unDSize-unLeft);
		int nRead = ifFile.gcount();
		ITER itStop = vPool.begin();
		if (unLeft+nRead < m_unDSize)
		{
			// for last block of data, give it a '>' to split the last sequence
			vPool[unLeft+nRead] = '>';
			advance(itStop, unLeft+nRead+1);
		}
		else
		{
			advance(itStop, unLeft+nRead);
		}


		// find a completed sequence
		ITER itSt = find(vPool.begin(), itStop, '>');
		ITER itBeg = find(itSt, itStop, '\n');
		ITER itEd = find(itBeg, itStop, '>');
		
		while (itEd != itStop)
		{
			vLens.push_back(0);
			while (itEd != itStop && vSeqs.size() < unBlockSize && vNames.size() < 2096152)	// 2^21=2097152, assume that the longest sequence is less than 1848000
			{
				++itBeg;

				// the lengths of some sequences are more than 2048
				ITER itEnd = remove(itBeg, itEd, '\r');
				itEnd = remove(itBeg, itEnd, '\n');
				int nLen = distance(itBeg, itEnd);
				if (nLen > 2048)
				{
					// store longer seq into several fragments with overlaps of size 200
					int nNum = (nLen-200) / 1848;
					if ((nLen-200) % 1848 != 0)
					{
						++nNum;
					}

					for (int i = 0; i < nNum-1; ++i)
					{
						if (bFullId == false)
						{
							vNames.push_back(string(itSt+1, find(itSt+1, itBeg-1, ' ')));
						}
						else
						{
							vNames.push_back(string(itSt+1, itBeg-1));
						}
						vSeqs.insert(vSeqs.end(), itBeg+i*1848, itBeg+2048+i*1848);
						vLens.push_back(vSeqs.size());
					}
					if (bFullId == false)
					{
						vNames.push_back(string(itSt+1, find(itSt+1, itBeg-1, ' ')));
					}
					else
					{
						vNames.push_back(string(itSt+1, itBeg-1));
					}
					vSeqs.insert(vSeqs.end(), itBeg+(nNum-1)*1848, itEnd);
					vLens.push_back(vSeqs.size());
				}
				else
				{
					if (bFullId == false)
					{
						vNames.push_back(string(itSt+1, find(itSt+1, itBeg-1, ' ')));
					}
					else
					{
						vNames.push_back(string(itSt+1, itBeg-1));
					}
					vSeqs.insert(vSeqs.end(), itBeg, itEnd);
					vLens.push_back(vSeqs.size());
				}

				itSt = itEd;
				itBeg = find(itSt, itStop, '\n');
				itEd = find(itBeg, itStop, '>');
			}

			// char to code
			lnTotalAa += Encode(vSeqs, vFreq);

			for (uint i = 0; i < vLens.size()-1; ++i)
			{
				// -1 or no, I need to think about it
				for (uint j = vLens[i]; j < vLens[i+1]-m_unMer; ++j)
				{
					int nIdx = Tran2Ten(vSeqs, j);
					if (-1 != nIdx)
					{
						// the left 21 bits denotes the index, the right 11 bits denotes the starting position of the seed
						vHash[nIdx].push_back((i<<11)|(j-vLens[i]));
					}
				}
			}

			// output information
			//PrintInfo(vHash);
			//PrintHash(vHash);

			VCOMP vComp(m_unTotalIdx, VUSHORT());
			for (uint i = 0; i < vHash.size(); ++i)
			{
				sort(vHash[i].begin(), vHash[i].end(), CompDbObj(vSeqs, vLens, m_unMer));
				for (uint j = 0; j < vHash[i].size(); ++j)
				{
					uint pos1 = vHash[i][j];
					int nIdx1 = pos1>>11;
					int nLen1 = vLens[nIdx1+1] - vLens[nIdx1];
					int nOff1 = (pos1&0x7ff) + m_unMer;
					uchar* p1 = &vSeqs[vLens[nIdx1]] + nOff1;

					int m = nLen1 - nOff1;
					int n = 0;
					ushort nSuff = 0;
					//string sSuff;
					if (m >= 4)
					{
						nSuff |= (m_aCode2Ten[p1[n]]) << 12;
						nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
						nSuff |= (m_aCode2Ten[p1[++n]]) << 4;
						nSuff |= (m_aCode2Ten[p1[++n]]);
					}
					else if (m == 3)
					{
						nSuff |= (m_aCode2Ten[p1[n]]) << 12;
						nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
						nSuff |= (m_aCode2Ten[p1[++n]]) << 4;
						nSuff |= ONEBYTE;
					}
					else if (m == 2)
					{
						nSuff |= (m_aCode2Ten[p1[n]]) << 12;
						nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
						nSuff |= TWOBYTE;
					}
					else if (m == 1)
					{
						nSuff |= (m_aCode2Ten[p1[n]]) << 12;
						nSuff |= THRBYTE;
					}
					else if (m == 0)
					{
						nSuff |= FOUBYTE;
					}
					vComp[i].push_back(nSuff);
				}
			}

			//serialize
			oa << vSeqs;
			oa << vLens;
			oa << vHash;
			oa << vNames;
			oa << vComp;
			
			uint unTotalWord = 0;
			for (uint i = 0; i < vHash.size(); ++i)
			{
				vWordCnts[i] += vHash[i].size();
				unTotalWord += vHash[i].size();
				vHash[i].clear();
			}
			vLens.clear();
			vSeqs.clear();
			vNames.clear();
			++nBlock;
		}

		// move the data at the end of pool to the beginning
		unLeft = distance(itSt, vPool.end());
		POOL::reverse_iterator itLast = find(vPool.rbegin(), vPool.rend(), '>');
		copy(itSt, vPool.end(), vPool.begin());
    }

	//serialize
	for (uint i = 0; i < vFreq.size(); ++i)
	{
		vFreq[i] /= lnTotalAa;
	}

	oaInfo << nBlock;
	oaInfo << lnSeqNum;
	//oaInfo << lnTotalAa;
	oaInfo << lnAaNum;
	oaInfo << vWordCnts;

	sort(vWordCnts.begin(), vWordCnts.end());
	//nth_element(vWordCnts.begin(), vWordCnts.begin()+m_unTotalIdx/2, vWordCnts.end());
	unMedian = vWordCnts[m_unTotalIdx/2];

	oaInfo << unMedian;
	oaInfo << vFreq;

	ifFile.close();
	of.close();
	ofInfo.close();

	return nBlock;
}


int CHashSearch::BuildQHash(istream& input, int nQueryType, map<string,char>& mTransTable, map<char,char>& mComple, Seg* seg, Seg* segsht, vector<uchar>& vSeqs, vector<uint>& vLens, VNAMES& vNames)
{
	char cIdSt = '>';
	char cSeqEd = '>';
	// init m_bSeqType & m_nIdxScl
	m_nQueryType = nQueryType;
	if (1 == m_nQueryType)
	{
		// nt
        m_bSeqType = true;
		m_nIdxScl = 6;
        printf("Queries are nucleotide sequences in fasta format\n");
	}
	else if (2 == m_nQueryType)
	{
		// aa
        m_bSeqType = false;
        printf("Queries are protein sequences\n");
	}
	else if (3 == m_nQueryType)
	{
		// fastq
        m_bSeqType = true;
		m_nIdxScl = 6;
        printf("Queries are nucleotide sequences in fastq format\n");
		cIdSt = '@';
		cSeqEd = '+';
	}

	int nSeqNum = 0;

	POOL vPool(m_unQSize, 0);
	if (!m_sLeft.empty())
	{
		copy(m_sLeft.begin(), m_sLeft.end(), vPool.begin());
	}

    if (input.good())
    {
		input.read(&vPool[0]+m_sLeft.size(), m_unQSize-m_sLeft.size());
		int nRead = input.gcount();
		if (nRead+m_sLeft.size() == 0)
		{
			return 0;
		}

		if (m_sLeft.size()+nRead < m_unQSize)
		{
			// for last block of data, give it a '>' to split the last sequence
			vPool[m_sLeft.size()+nRead] = cIdSt;
		}

		if (m_nQueryType == 0)
		{
			if (vPool[0] == '>')
			{
				m_nQueryType = GuessQueryType(vPool);
				if (1 == m_nQueryType)
				{
					// nt
					m_bSeqType = true;
					m_nIdxScl = 6;
					printf("Queries are nucleotide sequences in fasta format\n");
				}
				else if (2 == m_nQueryType)
				{
					// aa
					m_bSeqType = false;
					printf("Queries are protein sequences\n");
				}
			}
			else if (vPool[0] == '@')
			{
				// fastq
				m_bSeqType = true;
				m_nIdxScl = 6;
				printf("Queries are nucleotide sequences in fastq format\n");
				cIdSt = '@';
				cSeqEd = '+';
				m_nQueryType = 3;
			}
		}

		ITER itStop = vPool.begin();
		if (m_sLeft.size()+nRead < m_unQSize)
		{
			advance(itStop, m_sLeft.size()+nRead+1);
		}
		else
		{
			if (3 == m_nQueryType)
			{
				ITER itTemp = itStop;
				bool bFound = true;
				while (bFound)
				{
					bFound = true;
					for (int i = 0; i < 4; ++i)
					{
						itTemp = find(++itTemp, vPool.end(), '\n');
						if (vPool.end() == itTemp)
						{
							bFound = false;
							break;
						}
					}
					if (vPool.end() != itTemp)
					{
						itStop = ++itTemp;
					}
				}
				m_sLeft.clear();
				m_sLeft.assign(itStop, vPool.end());
				++itStop;
			}
			else
			{
				POOL::reverse_iterator rit = find(vPool.rbegin(), vPool.rend(), cIdSt);
				while ('\n' != *(++rit))
				{
					rit = find(rit, vPool.rend(), '>');
				}
				int n = distance(rit, vPool.rbegin());
				uint unLeft = m_unQSize + n + 1;
				m_sLeft.clear();
				m_sLeft.assign(vPool.begin()+(unLeft-1), vPool.end());
				itStop = vPool.begin() + unLeft;
			}
		}
		vPool.resize(distance(vPool.begin(), itStop));
		itStop = vPool.end();

		vLens.push_back(0);

		ITER itSt = find(vPool.begin(), itStop, cIdSt);
		ITER itBeg = find(itSt, itStop, '\n');
		ITER itEd = find(itBeg, itStop, cSeqEd);
		if (true == m_bSeqType)
		{
			// query is dna
			vector<char> vTran;
			vTran.reserve(1024);
			while (itEd != itStop)
			{
				++itBeg;
				vNames.push_back(/*'>'+*/string(itSt+1, find(itSt+1, itBeg-1, ' ')));
				vector<char> vS(itBeg, itEd-1);
				for (ITER itUpper = vS.begin(); itUpper != vS.end(); ++itUpper)
				{
					if (*itUpper >= 'a' && *itUpper <= 'z')
					{
						*itUpper = *itUpper - 'a' + 'A';
					}
				}

				vS.erase(remove(vS.begin(), vS.end(), '\r'), vS.end());
				vS.erase(remove(vS.begin(), vS.end(), '\n'), vS.end());
				int x = vS.size();
				for (int nFrame = 0; nFrame < 6; ++nFrame)
				{
					vTran.clear();

					if (3 == nFrame)
					{
						// backward
						reverse(vS.begin(), vS.end());
						for (uint nn = 0; nn < vS.size(); ++nn)
						{
							map<char, char>::iterator it = mComple.find(vS[nn]);
							if (it != mComple.end())
							{
								vS[nn] = it->second;
							}
							else
							{
								vS[nn] = 'N';
							}
						}
					}

					// translate
					char* pSt = &vS[0] + nFrame%3;
					char* pEd = pSt + ((x-nFrame%3)/3)*3;

					for (; pSt < pEd; pSt += 3)
					{
						map<string, char>::iterator it = mTransTable.find(string(pSt, 3));
						if (it != mTransTable.end())
						{
							vTran.push_back(it->second);
						}
						else
						{
							vTran.push_back(UNKNOWN_AA);
						}
					}
					vTran.push_back('\0');

					// mark the sequence
					char* pMasked = NULL;
					if (vTran.size()-1 >= 12)
					{
						pMasked = seg -> maskseq(&vTran[0]);
					}
					else
					{
						pMasked = segsht -> maskseq(&vTran[0]);
					}

					for (uint i = 0; i < strlen(pMasked); ++i)
					{
						if ('X' == pMasked[i] || 'x' == pMasked[i])
						{
							//vTran[i] = pMasked[i];
							vTran[i] += 128;
						}
					}

					delete [] pMasked;

					// a char '\0' was added at the end of vTran, so now ignore it
					vSeqs.insert(vSeqs.end(), vTran.begin(), vTran.end()-1);
					vLens.push_back(vSeqs.size());
				}

				itSt = itEd;
				if (3 == m_nQueryType)
				{
					itSt = find(itSt, itStop, '\n');
					++itSt;
					itSt = find(itSt, itStop, '\n');
					++itSt;
				}
				itBeg = find(itSt, itStop, '\n');
				itEd = find(itBeg, itStop, cSeqEd);
			}
		}
		else
		{
			while (itEd != itStop)
			{
				++itBeg;
				vNames.push_back(string(itSt+1, find(itSt+1, itBeg-1, ' ')));

				for (ITER itUpper = itBeg; itUpper != itEd-1; ++itUpper)
				{
					if (*itUpper >= 'a' && *itUpper <= 'z')
					{
						*itUpper = *itUpper - 'a' + 'A';
					}
				}

				int nIter = vSeqs.size();
				vSeqs.insert(vSeqs.end(), itBeg, itEd-1);
				vSeqs.erase(remove(vSeqs.begin()+nIter, vSeqs.end(), '\r'), vSeqs.end());
				vSeqs.erase(remove(vSeqs.begin()+nIter, vSeqs.end(), '\n'), vSeqs.end());
				vLens.push_back(vSeqs.size());

				itSt = itEd;
				itBeg = find(itSt, itStop, '\n');
				itEd = find(itBeg, itStop, cSeqEd);
			}
		}

		Encode(vSeqs);
		
		nSeqNum += vNames.size();
    }

	return nSeqNum;
}


void CHashSearch::Search(string& sDbPre, int nSeqNum, vector<uchar>& vQSeqs, vector<uint>& vQLens, VNAMES& vQNames)
{
	ifstream ifD(sDbPre.c_str());
	ifstream ifDInfo((sDbPre+".info").c_str());
	if (!ifD.good() || !ifDInfo.good())
	{
		ifstream if1((sDbPre+".des").c_str());
		ifstream if2((sDbPre+".des").c_str());
		ifstream if3((sDbPre+".des").c_str());
		ifstream if4((sDbPre+".des").c_str());
		if (if1.good()
			&& if2.good()
			&& if3.good()
			&& if4.good())
		{
			cout << "This database file comes from RAPSearch1." << endl;
			cout << "Please re-index the database using presearch from RAPSearch2" << endl;
		}
		else
		{
			cout << "Can not open the database file" << endl;
			cout << "Please check the file name" << endl;
		}

		if1.close();
		if2.close();
		if3.close();
		if4.close();
		return;
	}
	
	// construct query package
	CQrPckg Query(vQSeqs, vQLens, vQNames);

	archive::binary_iarchive iaD(ifD);
	archive::binary_iarchive iaDInfo(ifDInfo);

	vector<double> vFreq;
	VUINT vWordCnts(m_unTotalIdx, 0);
	int nDbBlockNum = 0;
	long int lnTotalAa;
	uint unMedian;
	long int lnSeqNum;

	iaDInfo >> nDbBlockNum;
	iaDInfo >> lnSeqNum;
	iaDInfo >> lnTotalAa;
	iaDInfo >> vWordCnts;
	iaDInfo >> unMedian;
	iaDInfo >> vFreq;

	m_lnSeqNum = lnSeqNum;
	m_lnTotalAa = lnTotalAa;

	// set BlastStat
	InitAlignPara(m_bSeqType, lnTotalAa, lnSeqNum, m_nThreadNum);

	pool tp(m_nThreadNum);
	cout << "start " << m_nThreadNum << "  threads" << endl;

	m_vOutIdx.assign(nSeqNum, CIndex());

	for (int j = 0; j < nDbBlockNum; ++j)
	{
		// temp file stream
		if (m_ofTemp.is_open())
		{
			m_ofTemp << m_sOutput;
			m_sOutput.clear();
			m_ofTemp.close();

			m_llOutCum = 0;
			ofstream ofOut((m_sOutBase+".tmp"+lexical_cast<string>(j-1)+".idx").c_str());
			archive::binary_oarchive oaOut(ofOut);
			oaOut << m_vOutIdx;
			ofOut.close();
		}
		for (int nn = 0; nn < nSeqNum; ++nn)
		{
			m_vOutIdx[nn].m_llBeg = 0;
			m_vOutIdx[nn].m_nSize = 0;
		}
		m_nSeqBase = 0;

		m_ofTemp.open((m_sOutBase+".tmp"+lexical_cast<string>(j)).c_str());

		// read db file and store info
		MINDEX vDHash(0, VUINT()); // all k-mer of database
		vector<uint> vDLens;
		vector<uchar> vDSeqs;
		VNAMES vDNames;
		VCOMP vComp;

		iaD >> vDSeqs;
		iaD >> vDLens;
		iaD >> vDHash;
		iaD >> vDNames;
		iaD >> vComp;

		// construct query package
		CDbPckg Db(vDHash, vDSeqs, vDLens, vDNames, vComp, vFreq, vWordCnts, unMedian);

		m_unTotalSubj += vDLens.size() - 1;


		//for (int i = 0; i < nQBlockNum; ++i)
		{

			m_unTotalQuery += vQLens.size() - 1;
			
			// generate results
			for (uint k = 0; k < vQLens.size()-1; k+=m_nIdxScl)
			{
				tp.schedule(bind(&CHashSearch::Searching, this, k, Query, Db));
				//Searching(k, Query, Db);
			}

			tp.wait();

			m_nSeqBase += vQNames.size();
		}
		
	}
	tp.wait();

	if (m_ofTemp.is_open())
	{
		m_ofTemp << m_sOutput;
		m_ofTemp.close();
		m_llOutCum = 0;

		ofstream ofOut((m_sOutBase+".tmp"+lexical_cast<string>(nDbBlockNum-1)+".idx").c_str());
		archive::binary_oarchive oaOut(ofOut);
		oaOut << m_vOutIdx;
		ofOut.close();
		m_vOutIdx.clear();
	}

	// merge nDbBlockNum temp results
	MergeRes(nDbBlockNum, vQNames, sDbPre);

	ifD.close();
	ifDInfo.close();
}


void CHashSearch::Process(char* szDBFile, char* szQFile, char* szOFile, int nStdout, bool bEvalue, bool bLogE, double dThr, int nMaxOut, int nMaxM8, int nQueryType, bool bPrintEmpty, bool bGapExt, bool bAcc, bool bHssp, int nMinLen, bool bXml, uint unDSize, uint unQSize, uint unMer)
{
	m_bEvalue = bEvalue;
	m_bLogE = bLogE;
	if (m_bEvalue == true)
	{
		m_pComptor = new CompEval();
	}
	else
	{
		m_pComptor = new CompBits();
	}
	m_dThr = dThr;
	if (m_bLogE == false)
	{
		//m_dThr = log(m_dThr);
                m_dThr = log(m_dThr) / log(10);
                //needs to be log_10, YY, Sep 2016
	}
	if (nMaxOut == -1)
	{
		m_nMaxOut = LLONG_MAX;
	}
	else
	{
		m_nMaxOut = abs(nMaxOut);
	}
	if (nMaxM8 == -1)
	{
		m_nMaxM8 = LLONG_MAX;
	}
	else
	{
		m_nMaxM8 = abs(nMaxM8);
	}
	m_bPrintEmpty = bPrintEmpty;
	m_bGapExt = bGapExt;
	m_bAcc = bAcc;
	m_bHssp = bHssp;
	m_nMinLen = nMinLen;
	m_bXml = bXml;

	m_unMer = unMer;
	m_unDSize = unDSize;
	m_unQSize = unQSize;
	m_unTotalIdx = lexical_cast<uint>(pow(10.0, int(m_unMer)));

	m_bFast = true;
	if (true == m_bFast)
	{
		m_unMutSeedLen = 10;
		m_vMutation.push_back(lexical_cast<uint>(pow(10.0, int(m_unMer-4-1))));
		m_vMutation.push_back(lexical_cast<uint>(pow(10.0, int(m_unMer-5-1))));
		m_vMutation.push_back(lexical_cast<uint>(pow(10.0, int(m_unMer-3-1))));
		if (m_unMer > 6)
		{
			m_vMutation.push_back(lexical_cast<uint>(pow(10.0, int(m_unMer-6-1))));
		}
	}
	else
	{
		m_unMutSeedLen = 9;
		m_vMutation.push_back(lexical_cast<uint>(pow(10.0, int(m_unMer-3-1))));
		m_vMutation.push_back(lexical_cast<uint>(pow(10.0, int(m_unMer-5-1))));
	}

	m_sQFile = szQFile;
	m_sDFile = szDBFile;
	m_nStdout = nStdout;
	if (szOFile != NULL)
	{
		m_sOutBase.assign(szOFile);

		string sDel = m_sOutBase + ".aln";
		ifstream iffTest;

		iffTest.open(sDel.c_str());
		if (iffTest.good())
		{
			iffTest.close();
			remove(sDel.c_str());
		}

		sDel = m_sOutBase + ".m8";
		iffTest.open(sDel.c_str());
		if (iffTest.good())
		{
			iffTest.close();
			remove(sDel.c_str());
		}
	}
	else
	{
		m_sOutBase = "";
	}

	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	m_sStartTime.assign(asctime(timeinfo));

	ifstream fIn;
	if (m_sQFile != "stdin")
	{
		fIn.open(m_sQFile.c_str());
		if (!fIn.good())
		{
			fIn.close();
			cout << "can not open the file: " << m_sQFile << endl;
			exit(1);
		}
	}
	istream& input = (m_sQFile!="stdin") ? fIn : cin;

	map<string, char> mTransTable;
	map<char, char> mComple;
	Seg* seg = NULL;
	Seg* segsht = NULL;
	//if (true == m_bSeqType)
	{
		// if need, construct paras for translation from nt to aa
		for (int i = 0; i < TOTCODON; ++i)
		{
			const char* p = nt[i];
			mTransTable[string(p, 3)] = aa[i];
		}

		mComple['A'] = 'T';
		mComple['T'] = 'A';
		mComple['a'] = 't';
		mComple['t'] = 'a';
		mComple['C'] = 'G';
		mComple['G'] = 'C';
		mComple['c'] = 'g';
		mComple['g'] = 'c';
		mComple['U'] = 'A';
		mComple['u'] = 'a';
		
		seg = new Seg(LGERSEED);
		segsht = new Seg(DEFSEED);
	}
	
	vector<uchar> vQSeqs;
	vector<uint> vQLens;
	VNAMES vQNames;
	int nSeqNum = 0;
	while ((nSeqNum=BuildQHash(input, nQueryType, mTransTable, mComple, seg, segsht, vQSeqs, vQLens, vQNames)) > 0)
	{
		string sDbOut(szDBFile);
		Search(sDbOut, nSeqNum, vQSeqs, vQLens, vQNames);
		vQSeqs.clear();
		vQLens.clear();
		vQNames.clear();
	}

	if (m_sQFile != "stdin")
	{
		((ifstream&)input).close();
	}

	if (seg != NULL)
	{
		delete seg;
	}
	if (segsht != NULL)
	{
		delete segsht;
	}
}


void CHashSearch::Process(char* szDBFile, char* szDbHash, bool bFullId, int nSplitNum, uint unMer)
{
	m_unMer = unMer;
	m_unTotalIdx = lexical_cast<uint>(pow(10.0, int(m_unMer)));

	string sDbOut(szDbHash);
	BuildDHash(szDBFile, sDbOut, nSplitNum, bFullId);
}


// rewrite this part
void CHashSearch::Searching(int k, CQrPckg& Query, CDbPckg& Db)
{
	//cout << k << endl;
	//cout << "my id:\t" << this_thread::get_id() << endl;
	//int nTreadID = m_mThreadID[this_thread::get_id()];
	//BlastStat* pBlastSig = m_vpBlastSig[m_mThreadID[this_thread::get_id()]];
	int nTreadID = -1;
	muMonitor.lock();
	for (uint i = 0; i < m_vBlastPt.size(); ++i)
	{
		if (-1 == m_vBlastPt[i])
		{
			nTreadID = i;
			m_vBlastPt[i] = 1;
			break;
		}
	}
	muMonitor.unlock();

	int nFoundHit = 0;
	//using namespace boost::chrono;
	//thread_clock::time_point start = thread_clock::now();

	MRESULT mRes;
	for (int nStep = 0; nStep < m_nIdxScl; ++nStep)
	{
		int nQrIdx = k + nStep;
		// the index of frame 0
		int nQDnaIdx = nQrIdx / m_nIdxScl * m_nIdxScl;

		// original length of query
		uint unQLen = Query.m_vLens[nQrIdx+1] - Query.m_vLens[nQrIdx];
		if (unQLen < m_unMer)
		{
			continue;
		}

		int nQOriLen = unQLen;
		if (true == m_bSeqType)
		{
			nQOriLen = 3 * (Query.m_vLens[nQDnaIdx+1]-Query.m_vLens[nQDnaIdx]);
			for (int n = 1; n < 3; ++n)
			{
				if (Query.m_vLens[nQDnaIdx+n+1]-Query.m_vLens[nQDnaIdx+n] == Query.m_vLens[nQDnaIdx+n]-Query.m_vLens[nQDnaIdx+n-1])
				{
					++nQOriLen;
				}
				else
				{
					break;
				}
			}
		}

		// set up BlastStat
		if (true == m_bSeqType)
		{
			// the index of the first seq considering the direction
			int n = nQrIdx / 3 * 3;
			m_vpBlastSig[nTreadID]->blastComputeLengthAdjustmentComp(Query.m_vLens[n+1]-Query.m_vLens[n]);
		}
		else
		{
			m_vpBlastSig[nTreadID]->blastComputeLengthAdjustmentComp(unQLen);
		}

		uchar* pQ = &Query.m_vSeqs[0] + Query.m_vLens[nQrIdx];
		CAlnPckg QrAln(pQ, unQLen, 0);
		
		// build invalid index position
		vector<char> vValid(unQLen, 0);
		for (uint xx = 0; xx < unQLen; ++xx)
		{
			vValid[xx] = m_aCode2Ten[pQ[xx]];
			if (0 != (m_uSeg&pQ[xx]))
			{
				pQ[xx] &= ~m_uSeg;
				vValid[xx] = m_uMask;
			}
		}

		// for consistence with swift
		uint unPrvSdLen = 6;

		for (uint i = 0; i < unQLen - m_unMer; ++i)
		{
			uint unCnt = 0;
		//thread_clock::time_point st = thread_clock::now();
			// pick seed length
			uint unQSeedBeg = QrAln.m_unSeedBeg = i;
			int nSeed = Tran2Ten(QrAln, vValid);
			if (-1 == nSeed)
			{
				continue;
			}
			uint unLocalSeed = 0; 
			uint unIdx = 0;
			if (m_bAcc == false)
			{
				uint unIncr = 0; 
				double	dFold = 0.0;
				int nLeft = unQLen - unQSeedBeg - m_unMer;
				uint unRng = nLeft+1 >= 3 ? 3 : nLeft+1;
				//uint unFreq = Db.m_vHash[nSeed].size();
				uint unFreq = Db.m_vWordCnts[nSeed];
				if(unFreq <= Db.m_unMedian) 
				{
					unLocalSeed = m_unMer;
				}
				else
				{
					double 	dExpFreq = unFreq;
					for(unIncr = 1; unIncr < unRng; unIncr ++)
					{
						if((unIdx = vValid[unQSeedBeg+m_unMer+unIncr-1]) != m_uMask)
						{
							dFold = Db.m_vFreq[unIdx];
						}
						else 
						{
							//dFold = 1.0 / strlen(murphy10r);
							break;
						}
						dExpFreq *= dFold;
						if(dExpFreq <= Db.m_unMedian) 
						{
							break;
						}
					}
					unLocalSeed = m_unMer + unIncr;
				}
				
				// if there is a unacceptable char, give up this seed
				if (m_uMask == unIdx)
				{
					continue;
				}

				if (unLocalSeed < unPrvSdLen - 1)
				{
					unLocalSeed = unPrvSdLen - 1;
				}

				// no enough letters
				if (unQSeedBeg+unLocalSeed > unQLen)
				{
					continue;
				}
			}
			else
			{
				unLocalSeed = 10;
				// no enough letters
				if (unQSeedBeg+unLocalSeed > unQLen)
				{
					continue;
				}
				for (uint i = m_unMer; i < unLocalSeed; ++i)
				{
					if((unIdx = vValid[unQSeedBeg+i]) == m_uMask)
					{
						break;
					}
				}
				if (m_uMask == unIdx)
				{
					continue;
				}
			}

			vector<uchar> vExtra(pQ+unQSeedBeg+m_unMer, pQ+unQSeedBeg+unLocalSeed);
			for (uint idx = 0; idx < vExtra.size(); ++idx)
			{
				vExtra[idx] = m_aCode2Ten[vExtra[idx]];
			}

			if (!Db.m_vHash[nSeed].empty())
			{
				int nCnt = ExtendSeq2Set(nSeed, unLocalSeed, vExtra,
						nQrIdx, QrAln, nQOriLen, vValid,
						Db.m_vHash[nSeed], Db,
						Query.m_vNames, Db.m_vNames,
						mRes, nTreadID);

				if (nCnt > 0)
				{
					unPrvSdLen = unLocalSeed;
				}
				else
				{
					unPrvSdLen = m_unMer;
				}
			}

			if (m_bAcc == true)
			{
				continue;
			}
			
			// mutation, pos: 4, 5, 3, (6)
			// check whether or not the length is enough
			if (unQLen < unQSeedBeg+m_unMutSeedLen)
			{
				continue;
			}
			// check non-aa char
			for (uint j = unQSeedBeg+unLocalSeed; j < unQSeedBeg+m_unMutSeedLen; ++j)
			{
				if((unIdx = vValid[j]) == m_uMask)
				{
					break;
				}
			}
			if (m_uMask == unIdx)
			{
				continue;
			}

			vExtra.assign(pQ+unQSeedBeg+m_unMer, pQ+unQSeedBeg+m_unMutSeedLen);
			for (uint idx = 0; idx < vExtra.size(); ++idx)
			{
				vExtra[idx] = m_aCode2Ten[vExtra[idx]];
			}

			for (uint m = 0; m < m_vMutation.size(); ++m)
			{
				int nVal = (nSeed/m_vMutation[m]) % 10;
				int nMutIdx = nSeed - nVal*m_vMutation[m];
				for (int n = 0; n < 10; ++n)
				{
					if (nMutIdx == nSeed)
					{
						nMutIdx += m_vMutation[m];
						continue;
					}

					if (Db.m_vHash[nMutIdx].empty())
					{
						nMutIdx += m_vMutation[m];
						continue;
					}

					int nCnt = ExtendSeq2Set(nMutIdx, m_unMutSeedLen, vExtra,
							nQrIdx, QrAln, nQOriLen, vValid,
							Db.m_vHash[nMutIdx], Db,
							Query.m_vNames, Db.m_vNames,
							mRes, nTreadID);
					nFoundHit += nCnt;
					unCnt += nCnt;

					nMutIdx += m_vMutation[m];
				}
			}

			/*********************************************************/
			if (6 == m_unMer && m_unMutSeedLen > m_unMer)
			{
				// mutation pos 6
				// think it as a mutation at pos 5, then do set intersection with current seed set
				//int nBase = (nQrIdx % 100000) * 10;
				int nNextNum = m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg+m_unMer]];
				for (int i = 0; i < 10; ++i)
				{
					if (i == nNextNum)
					{
						// if the mutation is equal to the original next char
						continue;
					}

					// change the 6th position
					vExtra[0] = i;
					int nCnt = ExtendSeq2Set(nSeed, m_unMutSeedLen, vExtra,
							nQrIdx, QrAln, nQOriLen, vValid,
							Db.m_vHash[nSeed], Db,
							Query.m_vNames, Db.m_vNames,
							mRes, nTreadID);
					nFoundHit += nCnt;
					unCnt += nCnt;
				}
			}
			//[>*******************************************************<]
		//thread_clock::time_point ed = thread_clock::now();
		//cout << nStep << "\t" << i << "\t" << unCnt << "\tduration:\t" << duration_cast<microseconds>(ed-st).count() << " ms" << endl;
		}
	}
	PrintRes(mRes, nTreadID, Query, Db);

	muMonitor.lock();
	m_vBlastPt[nTreadID] = -1;
	muMonitor.unlock();
}


struct CompSeed
{
	CompSeed(CDbPckg& Db, uint unMer, uchar* aCode2Ten) : m_Db(Db), m_unMer(unMer), m_aCode2Ten(aCode2Ten) {}
	bool operator() (const uint& unPos, const vector<uchar>& vExtra)
	{
		uint unIdx = unPos >> 11;
		uint unDSeedBeg = unPos & 0x000007FF;
		uint unDLen = m_Db.m_vLens[unIdx+1] - m_Db.m_vLens[unIdx];
		uchar* pD = &m_Db.m_vSeqs[0] + m_Db.m_vLens[unIdx];
		int nDOff = unDLen - unDSeedBeg - m_unMer;

		bool bLess = false;
		int nLeast = vExtra.size();
		int nDiff = nDOff>=nLeast ? nLeast : nDOff;

		if (0 == nDiff)
		{
			bLess = true;
		}
		else
		{
			uint i = m_unMer;
			for (; i < m_unMer+nDiff; ++i)
			{
				if (m_aCode2Ten[pD[unDSeedBeg+i]] != vExtra[i-m_unMer])
				{
					bLess = (m_aCode2Ten[pD[unDSeedBeg+i]] < vExtra[i-m_unMer]);
					break;
				}
			}
			// if they are the same for nLeast letters
			if (i==m_unMer+nDiff && bLess==false && m_aCode2Ten[pD[unDSeedBeg+i-1]] == vExtra[i-1-m_unMer])
			{
				if (nDiff < nLeast)
				{
					bLess = true;
				}
				else
				{
					bLess = false;
				}
			}
		}
		return bLess;
	}

	bool operator() (const vector<uchar>& vExtra, uint& unPos)
	{
		uint unIdx = unPos >> 11;
		uint unDSeedBeg = unPos & 0x000007FF;
		uint unDLen = m_Db.m_vLens[unIdx+1] - m_Db.m_vLens[unIdx];
		uchar* pD = &m_Db.m_vSeqs[0] + m_Db.m_vLens[unIdx];
		int nDOff = unDLen - unDSeedBeg - m_unMer;

		bool bLess = false;
		int nLeast = vExtra.size();
		int nDiff = nDOff>=nLeast ? nLeast : nDOff;

		if (0 == nDiff)
		{
			bLess = false;
		}
		else
		{
			uint i = m_unMer;
			for (; i < m_unMer + nDiff; ++i)
			{
				if (m_aCode2Ten[pD[unDSeedBeg+i]] != vExtra[i-m_unMer])
				{
					bLess = (vExtra[i-m_unMer] < m_aCode2Ten[pD[unDSeedBeg+i]]);
					break;
				}
			}
			if (i==m_unMer+nDiff && bLess==false && m_aCode2Ten[pD[unDSeedBeg+i-1]]==vExtra[i-1-m_unMer])
			{
				if (nDiff > nLeast)
				{
					bLess = true;
				}
				else
				{
					bLess = false;
				}
			}
		}
		return bLess;
	}

	CDbPckg& m_Db;
	uint m_unMer;
	uchar* m_aCode2Ten;
};

struct CompShortLow
{
	bool operator() (const ushort& s1, const ushort& s2)
	{
		int nLen1 = 4;
		if ((s1&ONEBYTE) == ONEBYTE)
		{
			--nLen1;
		}
		if ((s1&TWOBYTE) == TWOBYTE)
		{
			--nLen1;
		}
		if ((s1&THRBYTE) == THRBYTE)
		{
			--nLen1;
		}
		if ((s1&FOUBYTE) == FOUBYTE)
		{
			--nLen1;
		}

		int nLen2 = 4;
		if ((s2&ONEBYTE) == ONEBYTE)
		{
			--nLen2;
		}
		if ((s2&TWOBYTE) == TWOBYTE)
		{
			--nLen2;
		}
		if ((s2&THRBYTE) == THRBYTE)
		{
			--nLen2;
		}
		if ((s2&FOUBYTE) == FOUBYTE)
		{
			--nLen2;
		}

		int nLen = nLen1<nLen2 ? nLen1 : nLen2;
		if (0 == nLen)
		{
			return nLen1<nLen2;
		}
		bool b = ((s1>>((4-nLen)<<2)) == (s2>>((4-nLen)<<2)));
		if (true == b)
		{
			return nLen1<nLen2;
		}
		else
		{
			return ((s1>>((4-nLen)<<2)) < (s2>>((4-nLen)<<2)));
		}
	}
};

struct CompShortUp
{
	bool operator() (const ushort& s1, const ushort& s2)
	{
		int nLen1 = 4;
		if ((s1&ONEBYTE) == ONEBYTE)
		{
			--nLen1;
		}
		if ((s1&TWOBYTE) == TWOBYTE)
		{
			--nLen1;
		}
		if ((s1&THRBYTE) == THRBYTE)
		{
			--nLen1;
		}
		if ((s1&FOUBYTE) == FOUBYTE)
		{
			--nLen1;
		}

		int nLen2 = 4;
		if ((s2&ONEBYTE) == ONEBYTE)
		{
			--nLen2;
		}
		if ((s2&TWOBYTE) == TWOBYTE)
		{
			--nLen2;
		}
		if ((s2&THRBYTE) == THRBYTE)
		{
			--nLen2;
		}
		if ((s2&FOUBYTE) == FOUBYTE)
		{
			--nLen2;
		}

		int nLen = nLen1<nLen2 ? nLen1 : nLen2;
		if (0 == nLen)
		{
			return nLen1<nLen2;
		}
		bool b = ((s1>>((4-nLen)<<2)) == (s2>>((4-nLen)<<2)));
		if (true == b)
		{
			return false;
			//return nLen1<nLen2;
		}
		else
		{
			return ((s1>>((4-nLen)<<2)) < (s2>>((4-nLen)<<2)));
		}
	}
};

int CHashSearch::ExtendSeq2Set(int nSeed, uint unLocalSeedLen, vector<uchar>& vExtra,
		int nQSeqIdx, CAlnPckg& QrAln, int nQOriLen, vector<char>& vValid,
		VUINT& vDSet, CDbPckg& Db,
		VNAMES& vQNames, VNAMES& vDNames,
		MRESULT& mRes, int nTreadID)
{
	//using namespace boost::chrono;
	//thread_clock::time_point st = thread_clock::now();
	// find a proper range for the comparisons
	int nSt = 0;
	int nEd = 0;
	if (unLocalSeedLen > m_unMer)
	{
		//VUINT::iterator itSd = lower_bound(vDSet.begin(), vDSet.end(), vExtra, CompSeed(Db, m_unMer, m_aCode2Ten));
		//nSt = itSd - vDSet.begin();

		ushort nExtra = 0;
		for (uint i = 0; i < vExtra.size(); ++i)
		{
			nExtra |= (vExtra[i]) << (12-(i<<2));
		}
		for (int i = vExtra.size(); i < 4; ++i)
		{
			nExtra |= (ONEBYTE) << (12-(i<<2));
		}

		VUSHORT::iterator itShort = lower_bound(Db.m_vComp[nSeed].begin(), Db.m_vComp[nSeed].end(), nExtra, CompShortLow());
		nSt = itShort - Db.m_vComp[nSeed].begin();

		// check if nSt is real hit, if not, break
		{
			if (vDSet.size() == nSt)
			{
				return 0;
			}

			ushort s1 = Db.m_vComp[nSeed][nSt];
			ushort s2 = nExtra;
			int nLen1 = 4;
			if ((s1&ONEBYTE) == ONEBYTE)
			{
				--nLen1;
			}
			if ((s1&TWOBYTE) == TWOBYTE)
			{
				--nLen1;
			}
			if ((s1&THRBYTE) == THRBYTE)
			{
				--nLen1;
			}
			if ((s1&FOUBYTE) == FOUBYTE)
			{
				--nLen1;
			}

			int nLen2 = 4;
			if ((s2&ONEBYTE) == ONEBYTE)
			{
				--nLen2;
			}
			if ((s2&TWOBYTE) == TWOBYTE)
			{
				--nLen2;
			}
			if ((s2&THRBYTE) == THRBYTE)
			{
				--nLen2;
			}
			if ((s2&FOUBYTE) == FOUBYTE)
			{
				--nLen2;
			}

			int nLen = nLen1<nLen2 ? nLen1 : nLen2;
			if (0 == nLen)
			{
				return 0;
			}
			bool b = ((s1>>((4-nLen)<<2)) == (s2>>((4-nLen)<<2)));
			if (true != b)
			{
				return 0;
			}
		}

		//itSd = upper_bound(vDSet.begin(), vDSet.end(), vExtra, CompSeed(Db, m_unMer, m_aCode2Ten));
		//nEd = itSd - vDSet.begin();

		itShort = upper_bound(Db.m_vComp[nSeed].begin(), Db.m_vComp[nSeed].end(), nExtra, CompShortUp());
		nEd = itShort - Db.m_vComp[nSeed].begin();
	}
	else
	{
		nSt = 0;
		nEd = vDSet.size();
	}

	//st = thread_clock::now();
	//sequence extension
	STAlnmnt stAlnmnt;
	for (int j = nSt; j < nEd; ++j)
	{
		if (vDNames[vDSet[j]>>11] == "7719.ENSCINP00000006706")
		{
			int zya = 0;
		}
				
		uint unDLen, unDSeedBeg;
		uchar* pD = GetSeq(Db.m_vSeqs, Db.m_vLens, Db.m_vNames, vDSet[j], unDLen, unDSeedBeg);
		
		// no enough letters
		if (unDLen < unDSeedBeg + unLocalSeedLen)
		{
			continue;
		}
		// have the same previous char, so don't extend them at this time
		if (4 != vExtra.size()	 // if this is a mutation case, do not allow to ignore it
			&& m_bAcc==false 
			&& 0 != QrAln.m_unSeedBeg
			&& 0 != unDSeedBeg
			&& m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg-1]] == m_aCode2Ten[pD[unDSeedBeg-1]]
			&& m_uMask != vValid[QrAln.m_unSeedBeg-1]
			)
		{
			continue;
		}

		CAlnPckg DbAln(pD, unDLen, unDSeedBeg);

		ResetResult(stAlnmnt);
		// for debug

		uint unLocalCopy = unLocalSeedLen;
		uchar* pQAlign = QrAln.m_pSeq + QrAln.m_unSeedBeg;
		uchar* pDAlign =DbAln.m_pSeq + DbAln.m_unSeedBeg;

		int ii = 0;
		for (; ii < unLocalCopy; ++ii)
		{
			stAlnmnt.nScore += m_aSubMatrix[pQAlign[ii]][pDAlign[ii]];
			if (pQAlign[ii] == pDAlign[ii])
			{
				++stAlnmnt.nMatch;
			}
		}
		// forward maximal extension
		uint unTempLen = QrAln.m_unLen-QrAln.m_unSeedBeg < DbAln.m_unLen-DbAln.m_unSeedBeg ? QrAln.m_unLen-QrAln.m_unSeedBeg : DbAln.m_unLen-DbAln.m_unSeedBeg;
		while (ii < unTempLen && m_aCode2Ten[pQAlign[ii]] == m_aCode2Ten[pDAlign[ii]])
		{
			++unLocalCopy;
			stAlnmnt.nScore += m_aSubMatrix[pQAlign[ii]][pDAlign[ii]];
			if (pQAlign[ii] == pDAlign[ii])
			{
				++stAlnmnt.nMatch;
			}
			++ii;
		}
		// backward maximal extension
		uint unQSeed = QrAln.m_unSeedBeg;
		int nRange = QrAln.m_unSeedBeg < DbAln.m_unSeedBeg ? QrAln.m_unSeedBeg : DbAln.m_unSeedBeg;
		nRange = -nRange;
		ii = -1;
		while (ii >= nRange && m_aCode2Ten[pQAlign[ii]] == m_aCode2Ten[pDAlign[ii]])
		{
			++unLocalCopy;
			stAlnmnt.nScore += m_aSubMatrix[pQAlign[ii]][pDAlign[ii]];
			if (pQAlign[ii] == pDAlign[ii])
			{
				++stAlnmnt.nMatch;
			}
			--QrAln.m_unSeedBeg;
			--DbAln.m_unSeedBeg;
			--ii;
		}

		if (stAlnmnt.nScore >= UngapExtSCut && stAlnmnt.nMatch >= MinMatch4Exp)
		{
			AlignSeqs(nSeed, QrAln, DbAln, unLocalCopy, stAlnmnt, nTreadID);

			int nDSeqIdx = vDSet[j] >> 11;
			CalRes(nQSeqIdx, QrAln.m_pSeq, nQOriLen, QrAln.m_unSeedBeg, nDSeqIdx, DbAln.m_pSeq, DbAln.m_unSeedBeg, Db, unLocalCopy, stAlnmnt, mRes, nTreadID);
		}

		if (unQSeed != QrAln.m_unSeedBeg)
		{
			QrAln.m_unSeedBeg = unQSeed;
		}
	}

	return (nEd-nSt);
}


bool CHashSearch::AlignSeqs(int nSeed, CAlnPckg& QrAln, CAlnPckg& DbAln,  uint& unSeedLen, STAlnmnt& stAlnmnt, int nTreadID)
{
	uchar* pQAlign = QrAln.m_pSeq + QrAln.m_unSeedBeg;
	uchar* pDAlign =DbAln.m_pSeq + DbAln.m_unSeedBeg;
	int ext_f = 0;
	int ext_b = 0;
	int match_f = 0;
	int match_b = 0;

	int nScore0 = stAlnmnt.nScore;

	int nQLeft = 0;
	int nDLeft = 0;
	nQLeft = QrAln.m_unLen-QrAln.m_unSeedBeg-unSeedLen;
	nDLeft = DbAln.m_unLen-DbAln.m_unSeedBeg-unSeedLen;
	//if (nQLeft > 0 && nDLeft > 0)
	{
		pQAlign = QrAln.m_pSeq + QrAln.m_unSeedBeg + unSeedLen;
		pDAlign = DbAln.m_pSeq + DbAln.m_unSeedBeg + unSeedLen;
		stAlnmnt.nScore += AlignFwd(pQAlign, pDAlign, nQLeft, nDLeft, &ext_f, &match_f, nScore0);
		stAlnmnt.nMatch += match_f;
		stAlnmnt.nQFwd += ext_f;
		stAlnmnt.nDFwd += ext_f;
	}

	nQLeft = QrAln.m_unSeedBeg-1;
	nDLeft = DbAln.m_unSeedBeg-1;
	//if (nQLeft > 0 && nDLeft > 0)
	{
		pQAlign = QrAln.m_pSeq;
		pDAlign = DbAln.m_pSeq;
		stAlnmnt.nScore += AlignBwd(pQAlign, pDAlign, nQLeft, nDLeft, &ext_b, &match_b, nScore0);
		stAlnmnt.nMatch += match_b;
		stAlnmnt.nQBwd += ext_b;
		stAlnmnt.nDBwd += ext_b;
	}

	int hsplen = unSeedLen + ext_f + ext_b;

	stAlnmnt.vMode.push_back('s');
	stAlnmnt.vLen.push_back(hsplen);

	if (stAlnmnt.nScore < GapExtSCut)
	{
		//cout << "less than cut off" << endl;
		return true;
	}

	if (m_bAcc == false && true == m_bGapExt)
	{
		++m_unGapExt;
		//cout << "1 hit ..." << endl;
		// vector for alignment path
		vector<char> vMode;
		vector<short> vLen;

		// forward gapped alignment
		int exta = 0;
		int extb = 0;
		int gap = 0;
		int nQOff = QrAln.m_unSeedBeg + unSeedLen + stAlnmnt.nQFwd;
		int nDOff = DbAln.m_unSeedBeg + unSeedLen + stAlnmnt.nDFwd;
		nQLeft = QrAln.m_unLen-nQOff;
		nDLeft = DbAln.m_unLen-nDOff;
		if (nQLeft > 2 && nDLeft > 2)
		{
			pQAlign = QrAln.m_pSeq + nQOff;
			pDAlign = DbAln.m_pSeq + nDOff;
			int n1 = AlignGapped(pQAlign, pDAlign, nQLeft, nDLeft, &exta, &extb, &match_f, &gap, vMode, vLen, nTreadID);
			if (n1 > 0)
			{
				stAlnmnt.nScore += n1;
				stAlnmnt.nMatch += match_f;
				stAlnmnt.nQFwd += exta;
				stAlnmnt.nDFwd += extb;
				stAlnmnt.vMode.insert(stAlnmnt.vMode.end(), vMode.rbegin(), vMode.rend());
				stAlnmnt.vLen.insert(stAlnmnt.vLen.end(), vLen.rbegin(), vLen.rend());
			}
		}

		// backward gapped alignment
		nQOff = QrAln.m_unSeedBeg - stAlnmnt.nQBwd;
		nDOff = DbAln.m_unSeedBeg - stAlnmnt.nDBwd;
		if (nQOff > 2 && nDOff > 2)
		{
			vector<uchar> vQ;
			vQ.reserve(nQOff);
			for (int i = nQOff-1; i >= 0; --i)
			{
				vQ.push_back(QrAln.m_pSeq[i]);
			}
			pQAlign = &vQ[0];
			vector<uchar> vD;
			vD.reserve(nDOff);
			for (int i = nDOff-1; i >= 0; --i)
			{
				vD.push_back(DbAln.m_pSeq[i]);
			}
			pDAlign = &vD[0];
			int n2 = AlignGapped(pQAlign, pDAlign, nQOff, nDOff, &exta, &extb, &match_b, &gap, vMode, vLen, nTreadID);
			if (n2 > 0)
			{
				stAlnmnt.nScore += n2;
				stAlnmnt.nMatch += match_b;
				stAlnmnt.nQBwd += exta;
				stAlnmnt.nDBwd += extb;
				stAlnmnt.vMode.insert(stAlnmnt.vMode.begin(), vMode.begin(), vMode.end());
				stAlnmnt.vLen.insert(stAlnmnt.vLen.begin(), vLen.begin(), vLen.end());
			}
		}
	}

	return true;
}


int CHashSearch::AlignFwd(uchar *queryseq, uchar *dataseq, uint len_queryseq, uint len_dataseq, int *extl, int *match, int score0)
{
	int   i, j, l, s, maxs, ma;

	i = j = l = 0;
    ma = 0;
    maxs = s = score0;
    *extl = 0;
    *match = 0;
    while(i < len_queryseq && j < len_dataseq && s >= MINSCORE && s >= maxs - UngapExtDrop)
    {
		// uncompleted
        s += m_aSubMatrix[queryseq[i]][dataseq[j]];
        if(queryseq[i] == dataseq[j])
        {
            ma ++;
        }
        l ++;
        if(s > maxs)
        {
            maxs = s;
            *extl = l;
            *match = ma;
        }
        i ++;
        j ++;
    }
    return maxs - score0;
}


int CHashSearch::AlignBwd(uchar *queryseq, uchar *dataseq, int pos1, int pos2, int *extl, int *match, int score0)
{
    int   i, j, l, s, maxs, ma;

    i = pos1;
    j = pos2;
    l = 0;
    ma = 0;
    maxs = s = score0;
    *match = *extl = 0;
    while(i >= 0 && j >= 0 && s >= MINSCORE && s >= maxs - UngapExtDrop)
    {
        //	Skip stop codons
		//	uncompleted
        s += m_aSubMatrix[queryseq[i]][dataseq[j]];
        if(queryseq[i] == dataseq[j])
        {
            ma ++;
        }
        l ++;
        if(s > maxs)
        {
            maxs = s;
            *extl = l;
            *match = ma;
        }
        i --;
        j --;
    }
    return maxs - score0;
}


int CHashSearch::AlignGapped(uchar *seq1, uchar *seq2, int M, int N, int *ext1, int *ext2, int *match_len, int *gap, vector<char>& vMode, vector<short>& vLen, int nTreadID)
{
    int	i, j;
    int	t, s, e, c, d, wa;
    int	*CC = new int[N + 1]; //note N + 1
    int	*DD = new int[N + 1];
    int	g = GapIni;
    int	h = GapExt;
    int	m = g + h; //gap-create + gap-extend
    int	maxs, E1, E2, match;
    char	trace_e, trace_d;
    maxs = E1 = E2 = match = 0;

    //forward-phase
    CC[0] = 0;
    DD[0] = -g;
    t = -g;

    int	bb = 1; //band_begin
    int	be = int((GapExtDrop - GapIni) / GapExt);
    int	bb_pre, be_pre;
    //these two parameters will be adjusted during the alignment based on the dropoff score

	vector<vector<char> >& trace = m_vTrace[nTreadID];
	vector<vector<char> >& etrace = m_vETrace[nTreadID];
	vector<vector<char> >& dtrace = m_vDTrace[nTreadID];

	// the aligning sequences may be longer than 4096
	bool bModify = false;
	int nMemory = trace.size();
	if (trace.size()-1 < M)
	{
		int nSz = M+1;
		bModify = true;
		trace.clear();
		etrace.clear();
		dtrace.clear();
		trace.assign(nSz, vector<char>(nSz));
		etrace.assign(nSz, vector<char>(nSz));
		dtrace.assign(nSz, vector<char>(nSz));
	}

    trace[0][0] = '0';
    for(j = 1; j <= N && j <= be; j ++)
    {
        CC[j] = t = t - h; //j - 1 ? or j; when j is used, check score is not the same as alignment score
        DD[j] = CC[j] - g;
        if(j == 1)
        {
            trace[0][j] = etrace[0][j] = 'E';
        }
        else
        {
            trace[0][j] = etrace[0][j] = 'e';
        }
        dtrace[0][j] = 'D';
    } //global-alignment, with terminal penalty

    MaxGap = 100;
    for(i = 1; i <= M; i ++)
    {
        bb_pre = bb;
        be_pre = be;
        if(be <= bb) break; //band shrinks to zero
        s = CC[bb - 1];
        if(i == 1)
        {
            trace[i][bb - 1] = dtrace[i][bb - 1] = 'D';
            etrace[i][bb - 1] = 'E';
        }
        else
        {
            trace[i][bb - 1] = dtrace[i][bb - 1] = 'd';
            etrace[i][bb - 1] = 'e';
        }
        if(DD[bb - 1] - h > CC[bb - 1] - m)
        {
            c = DD[bb - 1] - h;
        }
        else
        {
            c = CC[bb - 1] - m;
        }
        CC[bb - 1] = DD[bb - 1] = c; //update it with current row
        e = c - g;
        for(j = bb; j <= be && j <= N; j ++)
        {
            trace_e = 'e'; //insertion extension
            if ((c =   c   - m) >= (e =   e   - h))
            {
                e = c;
                trace_e = 'E';  //new insertion
            }//insertion
            trace_d = 'd'; //deletion extension
            if ((c = CC[j] - m) >= (d = DD[j] - h))
            {
                d = c;
                trace_d = 'D'; //new deletion
            }//deletion
            //here   CC[j]==CC[i-1][j]   DD[j]==DD[i-1][j]

            wa = m_aSubMatrix[seq1[i - 1]][seq2[j - 1]];
            //sij[i - 1][j - 1]; //note i - 1, j - 1
            c = s + wa; //s==CC[i-1][j-1], substitution
            trace[i][j] = 's'; //substitution

            if (e > c)
            {
                c = e;
                trace[i][j] = trace_e;
            }
            if (d > c)
            {
                c = d;
                trace[i][j] = trace_d;
            }
            etrace[i][j] = trace_e;
            dtrace[i][j] = trace_d;
            s = CC[j]; //important for next replace
            CC[j] = c; //CC[i][j]
            DD[j] = d; //DD[i][j]
            if(c > maxs)
            {
                E1 = i;
                E2 = j;
                maxs = c;
            } //local -C
            else if(c < maxs - GapExtDrop && j > E2)   //score drops too much, stop filling this row, note j > E2
            {
                be = j;
                break;
            }
        }
        //after band_e, only allows insertion
        if(be < be_pre) continue;
        for(j = be + 1; j <= N; j ++)
        {
            trace_e = 'e'; //insertion extension
            if ((c =   c   - m) > (e =   e   - h))
            {
                e = c;
                trace_e = 'E';  //new insertion
            }//insertion
            c = e;
            trace[i][j] = trace_e;
            etrace[i][j] = trace_e;

            s = CC[j]; //important for next replace
            CC[j] = c; //CC[i][j]
            DD[j] = c - g;
            if(c > maxs)
            {
                E1 = i;
                E2 = j;
                maxs = c;
            } //local -C
            else if(c < maxs - GapExtDrop)   //score drops too much, stop filling this row
            {
                be = j;
                break;
            }
        }
        //now infer new bb (starting from E2 going backward)
        for(j = E2; j >= bb; j --)
        {
            if(CC[j] < maxs - GapExtDrop)
            {
                bb = j;
                break;
            }
        }
    }

    *ext1 = E1;
    *ext2 = E2;

    delete[] CC;
    delete[] DD;

    //get alignment
    *match_len = 0;
    *gap = 0;

    if(maxs <= 0) return maxs;

	
    if(trace[E1][E2] != 's')
    {
        printf("E1 %d E2 %d, Not end with substitution %c\n", E1, E2, trace[E1][E2]);
        exit(1);
    }
	
    char	mod = trace[E1][E2];
    i = E1;
    j = E2;
    vMode.clear();
    vLen.clear();
    while(mod != '0' && (!(i == 0 && j == 0)))
    {
        if (vMode.empty() || toupper(mod) != toupper(vMode.back()))
        {
			vMode.push_back(mod);
			vLen.push_back(0);
        }
		++vLen.back();
		
        if(mod == 's')
        {
            if(seq1[i - 1] == seq2[j - 1]) *match_len += 1;
            i -= 1;
            j -= 1;
            mod = trace[i][j];
        }
        else if(mod == 'D' || mod == 'd')
        {
            i -= 1;
            if (mod == 'D')	mod = trace[i][j];
            else 	mod = dtrace[i][j];
            *gap += 1;
        }
        else
        {
            j -= 1;
            if (mod == 'E')	mod = trace[i][j];
            else	mod = etrace[i][j];
            *gap += 1;
        }
		if (i<0 || j<0)
		{
			cout << "This is a bug!" << endl;
			for (int m = 0; m < M; ++m)
			{
				cout << m_aCode2Char[seq1[m]];
			}
			cout << endl;
			for (int n = 0; n < N; ++n)
			{
				cout << m_aCode2Char[seq2[n]];
			}
			cout << endl;
			break;
		}
    }

	// reset the size of the buffer
	if (bModify == true)
	{
		trace.clear();
		etrace.clear();
		dtrace.clear();
		trace.assign(nMemory, vector<char>(nMemory));
		etrace.assign(nMemory, vector<char>(nMemory));
		dtrace.assign(nMemory, vector<char>(nMemory));
	}

    return maxs;
}


void CHashSearch::CalRes(int nQIdx, uchar* pQ, int nQOriLen, uint unQSeedBeg, int nDIdx, uchar* pD, uint unDSeedBeg, CDbPckg& Db, uint unLocalSeedLen, STAlnmnt& stAlnmnt, MRESULT& mRes, int nTreadID)
{
	double dEValue = 0.0;
	dEValue = m_vpBlastSig[nTreadID]->rawScore2ExpectLog(stAlnmnt.nScore);
	double dBits = m_vpBlastSig[nTreadID]->rawScore2Bit(stAlnmnt.nScore);
	
	int nTotGap = 0;
	int nGapOpen = 0;
	int nTotAlnLen = 0;

	for (uint i = 0; i < stAlnmnt.vMode.size(); ++i)
	{
		nTotAlnLen += stAlnmnt.vLen[i];
		if ('s' != stAlnmnt.vMode[i])
		{
			++nGapOpen;
			nTotGap += stAlnmnt.vLen[i];
		}
	}

	// evalue criteria
	if (m_bHssp == false && !(stAlnmnt.nScore>SUMHSP_MINRAWSCORE || (m_bEvalue==true && dEValue<=m_dThr) || (m_bEvalue==false && dBits>=m_dThr)))
	{
		return;
	}
	// hssp criteria
	else if (m_bHssp == true && (nTotAlnLen < m_nMinLen || stAlnmnt.nMatch < m_vCriteria[nTotAlnLen]))
	{
		return;
	}

	// compute frame
	//cout << nQIdx << endl;
	int nQSt = 0;
	int nQEd = 0;
	if (m_bSeqType == true)
	{
		if (nQIdx % m_nIdxScl < 3)
		{
			nQSt = 3 * (unQSeedBeg-stAlnmnt.nQBwd) + nQIdx%m_nIdxScl + 1;
			nQEd = 3 * (unQSeedBeg+unLocalSeedLen+stAlnmnt.nQFwd) + nQIdx%m_nIdxScl;
		}
		else
		{
			int nFrame = nQIdx % m_nIdxScl - 3;
			nQSt = nQOriLen - (unQSeedBeg-stAlnmnt.nQBwd)*3 - nFrame;
			nQEd = nQSt - (stAlnmnt.nQBwd+unLocalSeedLen+stAlnmnt.nQFwd)*3 + 1;
		}
	}
	else
	{
		nQSt = unQSeedBeg-stAlnmnt.nQBwd + 1;
		nQEd = unQSeedBeg+unLocalSeedLen+stAlnmnt.nQFwd;
	}

	// print aligned sequences
	uint nAllc = nTotAlnLen>unLocalSeedLen?nTotAlnLen:unLocalSeedLen;
	VUCHAR vQ;
	vQ.reserve(nAllc);
	VUCHAR vD;
	vD.reserve(nAllc);

	uchar* pQAligned = pQ + unQSeedBeg - stAlnmnt.nQBwd;
	uchar* pDAligned = pD + unDSeedBeg - stAlnmnt.nDBwd;

	if (0 == stAlnmnt.vMode.size())
	{
		// only one hits
		vQ.insert(vQ.end(), pQAligned, pQAligned+unLocalSeedLen);
		vD.insert(vD.end(), pDAligned, pDAligned+unLocalSeedLen);
	}
	else if (1 == stAlnmnt.vMode.size())
	{
		// only one hits
		vQ.insert(vQ.end(), pQAligned, pQAligned+stAlnmnt.vLen[0]);
		vD.insert(vD.end(), pDAligned, pDAligned+stAlnmnt.vLen[0]);
	}
	else
	{
		for (uint i = 0; i < stAlnmnt.vMode.size(); ++i)
		{
			char cMode = stAlnmnt.vMode[i];
			if ('s' == cMode)
			{
				vQ.insert(vQ.end(), pQAligned, pQAligned+stAlnmnt.vLen[i]);
				pQAligned += stAlnmnt.vLen[i];
				vD.insert(vD.end(), pDAligned, pDAligned+stAlnmnt.vLen[i]);
				pDAligned += stAlnmnt.vLen[i];
			}
			else if ('D' == cMode || 'd' == cMode)
			{
				vQ.insert(vQ.end(), pQAligned, pQAligned+stAlnmnt.vLen[i]);
				pQAligned += stAlnmnt.vLen[i];
				vD.insert(vD.end(), stAlnmnt.vLen[i], '-');
			}
			else if ('E' == cMode || 'e' == cMode)
			{
				vQ.insert(vQ.end(), stAlnmnt.vLen[i], '-');
				vD.insert(vD.end(), pDAligned, pDAligned+stAlnmnt.vLen[i]);
				pDAligned += stAlnmnt.vLen[i];
			}
		}
	}
		
	string sQ;
	string sD;
	Decode(vQ, sQ);
	Decode(vD, sD);

	string sInfo;
	for (uint i = 0; i < vQ.size(); ++i)
	{
		if (vQ[i] == vD[i])
		{
			sInfo += sQ[i];
		}
		else if (m_aSubMatrix[vQ[i]][vD[i]] > 0)
		{
			sInfo += '+';
		}
		else
		{
			sInfo += ' ';
		}
	}

	MRESULT::iterator it = mRes.lower_bound(pair<int, int>(nQIdx/m_nIdxScl, nDIdx));
	/****************************************************************/
	// for sum evalue, comment this
	// note: here, the hits are stored according to it's real query index, not 1->6 frame query index
	// store all results
	if (mRes.end() != it 
			&& (*it).first.first==nQIdx/m_nIdxScl
			&& (*it).first.second==nDIdx
			&& (*it).second.nFrame==nQIdx%m_nIdxScl
			&& (*it).second.nQSt==unQSeedBeg-stAlnmnt.nQBwd
			&& (*it).second.nDSt==unDSeedBeg-stAlnmnt.nDBwd
			&& (*it).second.nQEd==unQSeedBeg+unLocalSeedLen+stAlnmnt.nQFwd-1
			&& (*it).second.nDEd==unDSeedBeg+unLocalSeedLen+stAlnmnt.nDFwd-1)
	{
		CHitUnit& st = (*it).second;
		if (st.dEValue > dEValue)
		{
			st.nScore = stAlnmnt.nScore;
			st.dBits = dBits;
			st.dEValue = dEValue;
			st.dIdent = stAlnmnt.nMatch*100.0/nTotAlnLen;
			st.nAlnLen = nTotAlnLen;
			st.nMismatch = nTotAlnLen-stAlnmnt.nMatch-nTotGap;
			st.nGapOpen = nGapOpen;
			st.nQBeg = nQSt;
			st.nQEnd = nQEd;
			st.sQ = sQ;
			st.sInfo = sInfo;
			st.sD = sD;
		}
	}
	else
	/****************************************************************/
	{
		MRESULT::iterator itTmp = mRes.insert(it, MRESULT::value_type(pair<int, int>(nQIdx/m_nIdxScl, nDIdx), CHitUnit()));
		CHitUnit& st = (*itTmp).second;
		st.nQrLen = nQOriLen;
		st.nDbIdx =	nDIdx; 
		st.nDbLen = Db.m_vLens[nDIdx+1] - Db.m_vLens[nDIdx];
		st.nScore = stAlnmnt.nScore;
		st.dBits = dBits;
		st.dEValue = dEValue;
		st.dIdent = stAlnmnt.nMatch*100.0/nTotAlnLen;
		st.nAlnLen = nTotAlnLen;
		st.nMismatch = nTotAlnLen-stAlnmnt.nMatch-nTotGap;
		st.nGapOpen = nGapOpen;
		st.nFrame = nQIdx%m_nIdxScl;
		st.nQSt = unQSeedBeg-stAlnmnt.nQBwd;
		st.nQEd = unQSeedBeg+unLocalSeedLen+stAlnmnt.nQFwd-1;
		st.nQBeg = nQSt;
		st.nQEnd = nQEd;
		st.nDSt = unDSeedBeg - stAlnmnt.nDBwd;
		st.nDEd = unDSeedBeg + unLocalSeedLen + stAlnmnt.nDFwd - 1;
		st.sQ = sQ;
		st.sInfo = sInfo;
		st.sD = sD;
		//mRes.insert(it, MRESULT::value_type(pair<int, int>(nQIdx/m_nIdxScl, nDIdx), st));
	}
}


struct SetCompObj
{
	bool operator() (const uint& p1, const uint& p2) const
	{
		return ((p1>>11) < (p2>>11)) && ((p1&0x7ff));
	}
}mySetComp;


void CHashSearch::PrintRes(MRESULT& mRes, int nTreadID, CQrPckg& Query, CDbPckg& Db)
{
	if (mRes.empty())
	{
		return;
	}

	MIT it = mRes.begin();
	int nQrIdx = (*it).first.first;
	MRESULT::iterator itFind = mRes.end();
	vector<CHitUnit> vTemp;
	vTemp.reserve(distance(it, itFind));
	
	// for sum evalue, comment this
	int nDIdx = it->first.second;
	vTemp.push_back(it->second);
	int nSt = 0; 
	MRESULT::iterator itTemp = it;
	++itTemp;
	for (; itTemp != itFind; ++itTemp)
	{
		if (itTemp->first.second != nDIdx)
		{
			if( vTemp.size() - nSt > 1)
			{
				int nLen = Db.m_vLens[nDIdx+1]-Db.m_vLens[nDIdx];
				SumEvalue(vTemp, nSt, vTemp.size(), nLen, nTreadID);
			}

			nDIdx = itTemp->first.second;
			nSt = vTemp.size();
		}
		vTemp.push_back(itTemp->second);
	}
	// process the last one
	if( vTemp.size() - nSt > 1)
	{
		int nLen = Db.m_vLens[nDIdx+1]-Db.m_vLens[nDIdx];
		SumEvalue(vTemp, nSt, vTemp.size(), nLen, nTreadID);
	}

	if (0 == vTemp.size())
	{
		return;
	}

	uint nMax = max(m_nMaxOut, m_nMaxM8);
	nMax = min((uint)vTemp.size(), nMax);
	vector<CHitUnit>::iterator itPrint = vTemp.begin()+nMax;
	partial_sort(vTemp.begin(), itPrint, vTemp.end(), ComptorWrapper(m_pComptor));

	int nBegStrAligned = 10;
	vector<CHitUnit>::iterator itSt = vTemp.begin();
	for (; itSt != itPrint; ++itSt)
	{
		CHitUnit& st= *itSt;
		if ((m_bEvalue==true && st.dEValue>m_dThr) || (m_bEvalue==false && st.dBits<m_dThr))
		{
			break;
		}
		// note: here, the hits are stored according to it's real query index, not 1->6 frame query index

		st.sInfo.insert(0, nBegStrAligned+1, ' ');

		string sQNum = lexical_cast<string>(st.nQBeg);
		st.sQ = string(nBegStrAligned-sQNum.size(), ' ') + sQNum + " " + st.sQ + " " + lexical_cast<string>(st.nQEnd);

		++st.nDSt;
		++st.nDEd;
		int nFac = 0;
		while (0 <= (st.nDbIdx-nFac-1) && Db.m_vNames[st.nDbIdx] == Db.m_vNames[st.nDbIdx-nFac-1])
		{
			++nFac;
		}
		st.nDSt = 1848*nFac+st.nDSt;
		st.nDEd = 1848*nFac+st.nDEd;
		string sDNum = lexical_cast<string>(st.nDSt);
		st.sD = string(nBegStrAligned-sDNum.size(), ' ') + sDNum + " " + st.sD + " " + lexical_cast<string>(st.nDEd);

		st.sQName = Query.m_vNames[nQrIdx];
		st.sDName = Db.m_vNames[st.nDbIdx];
	}

	vTemp.resize(itSt-vTemp.begin());
	if (vTemp.size() != 0)
	{
		// remove possible redundancy in overlapped region
		// there should be only one redundant hit for each overlapped region
		for (uint i = 1; i < vTemp.size(); ++i)
		{
			if (vTemp[i].nScore==vTemp[i-1].nScore
					&& vTemp[i].sDName == vTemp[i-1].sDName
					&& vTemp[i].sQName == vTemp[i-1].sQName
					&& vTemp[i].nDSt == vTemp[i-1].nDSt
					&& vTemp[i].nDEd == vTemp[i-1].nDEd
					&& vTemp[i].nQBeg == vTemp[i-1].nQBeg
					&& vTemp[i].nQEnd == vTemp[i-1].nQEnd)
			{
				vTemp[i].dEValue = 1000000;
			}
		}

		uint nf = 0;
		uint nl = vTemp.size();
		while (nf < nl)
		{
			if (vTemp[nf].dEValue == 1000000)
			{
				--nl;
				swap(vTemp[nf], vTemp[nl]);
			}
			++nf;
		}
		vTemp.resize(nl);

		stringstream sOutput;
		archive::binary_oarchive oa(sOutput);
		oa << vTemp;

		muMonitor.lock();
		long long llBeg = m_llOutCum + m_sOutput.size();
		m_sOutput += sOutput.str();
		int nSize = m_llOutCum + m_sOutput.size() - llBeg;
		m_vOutIdx[m_nSeqBase+nQrIdx].m_llBeg = llBeg;
		m_vOutIdx[m_nSeqBase+nQrIdx].m_nSize = nSize;
		if (m_sOutput.size() > 100000000)
		{
			m_ofTemp << m_sOutput;
			m_llOutCum += m_sOutput.size();
			m_sOutput.clear();
		}
		muMonitor.unlock();
	}

	mRes.clear();
}


void CHashSearch::SumEvalue(vector<CHitUnit>& v, int nSt, int nEd, int nLen, int nTreadID)
{
	typedef vector<CHitUnit>::iterator STIT;
	STIT itSt = v.begin() + nSt;
	STIT itEd = v.begin() + nEd;
	// sort by nFrame
	sort(itSt, itEd, CompFrame());
	CHitUnit st;
	st.nFrame = 3;
	STIT itDir = lower_bound(itSt, itEd, st, CompFrame());
	// if there are more than one hit in one direction
	int nDisPos = distance(itSt, itDir);
	int nDisNeg = distance(itDir, itEd);
	if (nDisPos > 1 || nDisNeg > 1)
	{
		vector<CHitUnit> vRes;
		STIT itStart = itSt;
		STIT itEnd = itDir;
		for (int i = 0; i < 2; ++i)
		{
			if (distance(itStart, itEnd) == 0)
			{
				itStart = itEnd;
				itEnd = itEd;
				continue;
			}
			else if (distance(itStart, itEnd) == 1)
			{
				if ((m_bEvalue==true && itStart->dEValue<=m_dThr) || (m_bEvalue==false && itStart->dBits>=m_dThr))
				{
					vRes.push_back(*itStart);
				}
				itStart = itEnd;
				itEnd = itEd;
				continue;
			}
			// sort by score and start position of query
			sort(itStart, itEnd, CompQSt());
			stable_sort(itStart, itEnd, ComptorWrapper(m_pComptor));
			// check overlap and logevalue
			vector<CHitUnit> vNew;
			vNew.push_back(*itStart);
			for (STIT itTemp = itStart+1; itTemp!=itEnd; ++itTemp)
			{
				int nHalfLen = (itTemp->nQEd - itTemp->nQSt + 1) >> 1;
				int nOverlap = SUMHSP_OVERLAP < nHalfLen ? SUMHSP_OVERLAP : nHalfLen;
				if (itTemp->dEValue >= SUMHSP_MINEVALUE && itTemp->nScore <= SUMHSP_MINRAWSCORE)
				{
					continue;
				}
				bool bNonOvlp = true;
				for (STIT itIso = vNew.begin(); itIso != vNew.end(); ++itIso)
				{
					if ((itTemp->nQSt <= itIso->nQEd - nOverlap
							&& itTemp->nQEd >= itIso->nQSt + nOverlap)
						|| (itIso->nQSt <= itTemp->nQEd - nOverlap
							&& itIso->nQEd >= itTemp->nQSt + nOverlap))
					{
						bNonOvlp = false;
						break;
					}
				}
				if (true == bNonOvlp)
				{
					vNew.push_back(*itTemp);
				}
			}
			if (vNew.size() == 1)
			{
				if ((m_bEvalue==true && vNew[0].dEValue<=m_dThr) || (m_bEvalue==false && vNew[0].dBits>=m_dThr))
				{
					vRes.push_back(vNew[0]);
				}
				//continue;
			}
			else
			{
				// calculate the sum of evalue
				double aRawScore[DEFAULT_SCORE_TOP];
				int nNo = 0;
				for (; nNo < 5 && nNo < vNew.size(); ++nNo)
				{
					aRawScore[nNo] = vNew[nNo].nScore;
				}

				if (m_bEvalue == true)
				{
					double dTmp = m_vpBlastSig[nTreadID]->sumScore2Expect(nNo, aRawScore, nLen);
					double dSumEvalue = -10000.00;
					if (0 != dTmp)
					{
						dSumEvalue = log(dTmp) / LOG10;
					}
					// modify the logevalue
					if (dSumEvalue < m_dThr)
					{
						for (uint i = 0; i < vNew.size(); ++i)
						{
							vNew[i].dEValue = dSumEvalue;
						}
						vRes.insert(vRes.end(), vNew.begin(), vNew.end());
					}
				}
				else
				{
					double dTmpScore = m_vpBlastSig[nTreadID]->sumScore(nNo, aRawScore, nLen);
					double dTmpBits = m_vpBlastSig[nTreadID]->rawScore2Bit(dTmpScore);
					// modify the logevalue
					if (dTmpBits >= m_dThr)
					{
						for (uint i = 0; i < vNew.size(); ++i)
						{
							vNew[i].dBits = dTmpBits;
						}
						vRes.insert(vRes.end(), vNew.begin(), vNew.end());
					}
				}
			}
			itStart = itEnd;
			itEnd = itEd;
		}
		// replace
		if (!vRes.empty())
		{
			v.erase(itSt, itEd);
			v.insert(v.begin()+nSt, vRes.begin(), vRes.end());
		}
	}
}


void CHashSearch::GuessTotSeq(const char* szDBFile, long int& lnSeqNum, long int& lnAaNum)
{
	lnSeqNum = 0;
	lnAaNum = 0;
	ifstream fIn(szDBFile);
	string s;
	while (fIn.good())
	{
		getline(fIn, s);
		if (s[0] == '>')
		{
			lnSeqNum += 1;
		}
		else
		{
			lnAaNum += s.size();
		}
	}
	fIn.close();
}



void CHashSearch::MergeRes(int nDbBlockNum, VNAMES& vQNames, string& sDbPre)
{
	ostream* poAln = NULL;
	ostream* poM8 = NULL;

	if (m_nStdout == 2)
	{
		poAln = &cout;
	}
	else if (!m_sOutBase.empty() && m_nMaxOut != 0)
	{
		poAln = new ofstream((m_sOutBase+".aln").c_str(),ios_base::out|ios_base::app);
		if (!poAln->good())
		{
			((ofstream*)poAln)->close();
			delete poAln;
			cout << "can not open the file: " << m_sOutBase+".aln" << endl;
			exit(1);
		}
	}

	if (m_nStdout == 1)
	{
		poM8 = &cout;
	}
	else if (!m_sOutBase.empty() && m_nMaxM8 != 0)
	{
		poM8 = new ofstream((m_sOutBase+".m8").c_str(),ios_base::out|ios_base::app);
		if (!poM8->good())
		{
			((ofstream*)poM8)->close();
			delete poM8;
			cout << "can not open the file: " << m_sOutBase+".m8" << endl;
			exit(1);
		}
	}

	if (poM8 && !m_sStartTime.empty())
	{
		(*poM8) << "# RAPSearch\n# Job submitted: "
			<< m_sStartTime
			<< "# Query : " << m_sQFile << "\n"
			<< "# Subject : " << m_sDFile << "\n";
		if (m_bLogE == true)
		{
			(*poM8) << "# Fields: Query\tSubject\tidentity\taln-len\tmismatch\tgap-openings\tq.start\tq.end\ts.start\ts.end\tlog(e-value)\tbit-score\n";
		}
		else
		{
			(*poM8) << "# Fields: Query\tSubject\tidentity\taln-len\tmismatch\tgap-openings\tq.start\tq.end\ts.start\ts.end\te-value\tbit-score\n";
		}

		m_sStartTime = "";
	}
	
	if (m_bXml)
	{
		m_ofXml.open((m_sOutBase+".xml").c_str(),ios_base::out|ios_base::app);
		PrintXmlBegin(sDbPre);
	}

	long long unMax = max(m_nMaxOut, m_nMaxM8);
	vector<CHitUnit> v;
	vector<CMergeUnit*> vMergeUnit;

	for (int i = 0; i < nDbBlockNum; ++i)
	{
		string sName = m_sOutBase+".tmp"+lexical_cast<string>(i);
		CMergeUnit* p = new CMergeUnit(sName.c_str());
		vMergeUnit.push_back(p);
	}
	
	int nLastIdx = 0;
	for (int i = 0; i < nDbBlockNum; ++i)
	{
		nLastIdx = max(nLastIdx, vMergeUnit[i]->GetLast());
	}

	/**************************************************************/
	for (int i = 0; i < nLastIdx; ++i)
	{
		for (int j = 0; j < nDbBlockNum; ++j)
		{
			vMergeUnit[j]->Update(i, v);
		}

		/***********************************************************/
		if (0 == v.size())
		{
			if (poAln && true == m_bPrintEmpty)
			{
				(*poAln) << vQNames[i] << "\tNO HIT" << "\n\n";
			}
			continue;
		}
		/***********************************************************/

		uint n = min(unMax, (long long)v.size());
		partial_sort(v.begin(), v.begin()+n, v.end(), ComptorWrapper(m_pComptor));
		v.resize(n);

		if (poAln)
		{
			PrintAln(v, *poAln);
		}

		if (poM8)
		{
			PrintM8(v, *poM8);
		}

		if (m_bXml)
		{
			PrintXml(v, i+1);
		}

		v.clear();
	}

	for (int i = 0; i < nDbBlockNum; ++i)
	{
		delete vMergeUnit[i];
	}

	if (poAln && m_nStdout != 2)
	{
		((ofstream*)poAln)->close();
		delete poAln;
		poAln = NULL;
	}

	if (poM8 && m_nStdout != 1)
	{
		((ofstream*)poM8)->close();
		delete poM8;
		poM8 = NULL;
	}
	if (m_bXml)
	{
		PrintXmlEnd();
		m_ofXml.close();
	}
}


//----------------------------------------------------------------------
int CHashSearch::GuessQueryType(POOL& vPool)
{
	// read some reads
	vector<uchar> seq;
	ITER itStop = vPool.end();
	ITER itSt = find(vPool.begin(), itStop, '>');
	ITER itEd = find(itSt+1, itStop, '>');
	while (itEd!=itStop && seq.size()<10000)
	{
		ITER itBeg = find(itSt, itEd, '\n');
		while (itBeg == itEd)
		{
			itEd = find(itEd+1, itStop, '>');
			itBeg = find(itBeg, itEd, '\n');
		}
		++itBeg;
		int nIter = seq.size();
		seq.insert(seq.end(), itBeg, itEd-1);
		seq.erase(remove(seq.begin()+nIter, seq.end(), '\n'), seq.end());

		itSt = itEd;
		itEd = find(itSt+1, itStop, '>');
	}

    char	nuc[] = "ATCGatcgUu";
    int	i, j;
    int	add = 0;
    for(i = 0; i < seq.size(); i ++)
    {
        for(j = 0; j < 10; j ++)
        {
            if(seq[i] == nuc[j]) break;
        }
        if(j < 10) add += 1;
    }

    if(add > seq.size() * 0.95)
    {
		return 1;
    }
    else
    {
		return 2;
    }
}


void CHashSearch::PrintAln(vector<CHitUnit>& v, ostream& of)
{
	int nPrint = min((long long)v.size(), m_nMaxOut);
	for (int i = 0; i < nPrint; ++i)
	{
		CHitUnit& c = v[i];

		of << c.sQName
			<< " vs "
			<< c.sDName
			<< " bits="	<< c.dBits;
		if (m_bLogE == true)
		{
			of << " log(E-value)="	<< c.dEValue;
		}
		else
		{
			of << " E-value="	<< pow(10,c.dEValue);
		}
		of << " identity=" << c.dIdent << "%"
			<< " aln-len="	<< c.nAlnLen
			<< " mismatch="	<< c.nMismatch
			<< " gap-openings="	<< c.nGapOpen
			<< " nFrame=" << c.nFrame
			<< "\n"
			<< "Query:\t" << c.sQ << "\n"
			<< "      \t" << c.sInfo << "\n"
			<< "Sbjct:\t" << c.sD << "\n"
			<< "\n";
	}
}


void CHashSearch::PrintM8(vector<CHitUnit>& v, ostream& of)
{
	int nPrint = min((long long)v.size(), m_nMaxM8);
	for (int i = 0; i < nPrint; ++i)
	{
		CHitUnit& c = v[i];

		of << c.sQName
			<< "\t" << c.sDName
			<< setprecision(1) << setiosflags(ios::fixed)
			<< "\t" << c.dIdent
			<< "\t"	<< c.nAlnLen
			<< "\t"	<< c.nMismatch
			<< "\t"	<< c.nGapOpen 
			<< "\t" << c.nQBeg
			<< "\t" << c.nQEnd
			<< "\t" << c.nDSt
			<< "\t" << c.nDEd;
		if (m_bLogE == true)
		{
			of << setprecision(1) << setiosflags(ios::fixed)
				<< "\t"	<< c.dEValue;
		}
		else
		{
			c.dEValue = pow(10, c.dEValue);
			if (c.dEValue < 0.01)
			{
				of << setprecision(2) << setiosflags(ios::scientific) << setiosflags(ios::fixed)
					<< "\t"	<< c.dEValue;
				of << resetiosflags(ios::scientific);
			}
			else if (c.dEValue < 10.0)
			{
				of << setprecision(2) << setiosflags(ios::fixed)
					<< "\t"	<< c.dEValue;
			}
			else
			{
				of << setprecision(0) << setiosflags(ios::fixed)
					<< "\t"	<< c.dEValue;
			}
		}
		of << setprecision(1) << setiosflags(ios::fixed)
			<< "\t"	<< c.dBits 
			<< "\n";
	}
}


template<class T>
void CHashSearch::PrintXmlLine(char* sTag, T s)
{
	m_ofXml << string(m_unXmlSp, ' ') << "<" << sTag << ">" << s << "</" << sTag << ">" << "\n";
}


void CHashSearch::PrintXmlTag(char* sTag)
{
	m_ofXml << string(m_unXmlSp, ' ') << "<" << sTag << ">" << "\n";
	m_unXmlSp += 2;
}


void CHashSearch::PrintXmlTagR(char* sTag)
{
	m_unXmlSp -= 2;
	m_ofXml << string(m_unXmlSp, ' ') << "</" << sTag << ">" << "\n";
}


void CHashSearch::PrintXmlBegin(string& sDbPre)
{
	m_ofXml << "<?xml version=\"1.0\"?>" << "\n";
	PrintXmlTag("Output");
	PrintXmlLine("Output_program", "RAPSearch");
	PrintXmlLine("Output_version", "RAPSearch2");
	PrintXmlLine("Output_reference", "Yongan Zhao, Haixu Tang and Yuzhen Ye. RAPSearch2: a fast and memory-efficient protein similarity search tool for next generation sequencing data. Bioinformatics 2012, 28 (1): 125-126");
	PrintXmlLine("Output_db", sDbPre);
	PrintXmlTag("Output_param");
	PrintXmlTag("Parameters");
	PrintXmlLine("Parameters_matrix", "BLOSUM62");
	if (m_bEvalue == true)
	{
		if (m_bLogE == true)
		{
			PrintXmlLine("Parameters_log-expect_evalue", lexical_cast<string>(m_dThr));
		}
		else
		{
			PrintXmlLine("Parameters_expect_evalue", lexical_cast<string>(pow(10,m_dThr)));
		}
	}
	else
	{
		PrintXmlLine("Parameters_bits-expect", lexical_cast<string>(m_dThr));
	}
	PrintXmlLine("Parameters_gap-open", "11");
	PrintXmlLine("Parameters_gap-extend", "1");
	PrintXmlLine("Parameters_filter", "T");
	PrintXmlTagR("Parameters");
	PrintXmlTagR("Output_param");
	PrintXmlTag("Output_iterations");
}


void CHashSearch::PrintXml(vector<CHitUnit>& v, int nIdx)
{
	PrintXmlTag("Iteration");
	PrintXmlLine("Iteration_iter-num", lexical_cast<string>(m_unXmlCnt++));
	//PrintXmlLine("Iteration_query-ID", "lcl|"+lexical_cast<string>(nIdx));
	PrintXmlLine("Iteration_query-def", v[0].sQName);
	PrintXmlLine("Iteration_query-len", v[0].nQrLen);
	PrintXmlTag("Iteration_hits");

	int nPrint = min((long long)v.size(), m_nMaxOut);
	for (int i = 0; i < nPrint; ++i)
	{
		CHitUnit& c = v[i];

		PrintXmlTag("Hit");
		PrintXmlLine("Hit_num", lexical_cast<string>(i+1));
		//PrintXmlLine("Hit_id", "gnl|"+lexical_cast<string>(c.nDbIdx));
		PrintXmlLine("Hit_def", c.sDName);
		//PrintXmlLine("Hit_accession", c.nDbIdx);
		PrintXmlLine("Hit_len", c.nDbLen);
		PrintXmlTag("Hit_hsps");
		PrintXmlTag("Hsp");
		PrintXmlLine("Hsp_num", 1);
		PrintXmlLine("Hsp_bit-score", c.dBits);
		PrintXmlLine("Hsp_score", c.nScore);
		if (m_bLogE == true)
		{
			PrintXmlLine("Hsp_log-evalue", c.dEValue);
		}
		else
		{
			PrintXmlLine("Hsp_evalue", pow(10,c.dEValue));
		}
		PrintXmlLine("Hsp_query-from", c.nQBeg);
		PrintXmlLine("Hsp_query-to", c.nQEnd);
		PrintXmlLine("Hsp_hit-from", c.nDSt);
		PrintXmlLine("Hsp_hit-to", c.nDEd);
		PrintXmlLine("Hsp_query-frame", c.nFrame);
		int nPos = 0;
		int nIdt = 0;
		for (uint j = 0; j < c.sInfo.size(); ++j)
		{
			if (c.sInfo[j] != ' ')
			{
				nPos += 1;
				if (c.sInfo[j] != '+')
				{
					nIdt += 1;
				}
			}
		}
		PrintXmlLine("Hsp_identity", nIdt);
		PrintXmlLine("Hsp_positive", nPos);
		PrintXmlLine("Hsp_align-len", c.nAlnLen);
		size_t n1 = c.sQ.find_first_not_of(" 0123456789");
		size_t n2 = c.sQ.find_last_not_of(" 0123456789");
		PrintXmlLine("Hsp_qseq", c.sQ.substr(n1, n2-n1+1));
		n1 = c.sD.find_first_not_of(" 0123456789");
		n2 = c.sD.find_last_not_of(" 0123456789");
		PrintXmlLine("Hsp_hseq", c.sD.substr(n1, n2-n1+1));
		n1 = c.sInfo.find_first_not_of(" ");
		n2 = c.sInfo.find_last_not_of(" ");
		PrintXmlLine("Hsp_midline", c.sInfo.substr(n1, n2-n1+1));
		PrintXmlTagR("Hsp");
		PrintXmlTagR("Hit_hsps");
		PrintXmlTagR("Hit");
	}

	PrintXmlTagR("Iteration_hits");

	PrintXmlTag("Iteration_stat");
	PrintXmlTag("Statistics");
	PrintXmlLine("Statistics_db-num", m_lnSeqNum);
	PrintXmlLine("Statistics_db-len", m_lnTotalAa);
	PrintXmlLine("Statistics_hsp-len", 0);
	PrintXmlLine("Statistics_eff-space", 0);
	PrintXmlLine("Statistics_kappa", 0.041);
	PrintXmlLine("Statistics_lambda", 0.267);
	PrintXmlLine("Statistics_entropy", 0.14);
	PrintXmlTagR("Statistics");
	PrintXmlTagR("Iteration_stat");

	PrintXmlTagR("Iteration");
}


void CHashSearch::PrintXmlEnd()
{
	PrintXmlTagR("Output_iterations");
	PrintXmlTagR("Output");
}
