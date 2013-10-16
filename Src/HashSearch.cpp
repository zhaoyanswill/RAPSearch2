
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
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/chrono/thread_clock.hpp>
#include "threadpool.hpp"
#include "weight.h"
#include "aa.h"
#include "n2a.h"
#include "Seg.h"
#include "sortUnit.h"
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
	//m_uMask = 238;
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
		for (int j = 0; j < strlen(p); ++j) 
		{
			uint unIdx = (i << 4) + j+1;
			m_aChar2Code[p[j]] = unIdx;
			// add lower-case character
			m_aChar2Code[p[j]+32] = unIdx;
			m_aCode2Char[unIdx] = p[j];
			m_aCode2Ten[unIdx] = i;
		}
	}

	// build substitute matrix for compressed char
	fill_n((int*)m_aSubMatrix, 256*256, -5);
	for (int i = 0; i < strlen(aAAAlph); ++i)
	{
		for (int j = 0; j < strlen(aAAAlph); ++j)
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

	m_sOutBase = "";
	m_sOutput = "";
	m_sOutput.reserve(50000000);
	m_sM8 = "";
	m_sM8.reserve(50000000);
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


int CHashSearch::BuildDHash(const char* szFile, string& sOutFile, int nSplitNum)
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
	int nFileNum = (lnFileSize-1)/lnBlockSize + 1;

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
	uint unWrited = 0;
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
				int nLen = distance(itBeg, itEd-1);
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
						vNames.push_back(string(itSt+1, find(itSt+1, itBeg-1, ' ')));
						vSeqs.insert(vSeqs.end(), itBeg+i*1848, itBeg+2048+i*1848);
						vSeqs.erase(remove(vSeqs.begin()+vLens.back(), vSeqs.end(), '\r'), vSeqs.end());
						vSeqs.erase(remove(vSeqs.begin()+vLens.back(), vSeqs.end(), '\n'), vSeqs.end());
						vLens.push_back(vSeqs.size());
					}
					vNames.push_back(string(itSt+1, find(itSt+1, itBeg-1, ' ')));
					vSeqs.insert(vSeqs.end(), itBeg+(nNum-1)*1848, itEd-1);
					vSeqs.erase(remove(vSeqs.begin()+vLens.back(), vSeqs.end(), '\n'), vSeqs.end());
					vLens.push_back(vSeqs.size());
				}
				else
				{
					vNames.push_back(string(itSt+1, find(itSt+1, itBeg-1, ' ')));
					vSeqs.insert(vSeqs.end(), itBeg, itEd-1);
					vSeqs.erase(remove(vSeqs.begin()+vLens.back(), vSeqs.end(), '\r'), vSeqs.end());
					vSeqs.erase(remove(vSeqs.begin()+vLens.back(), vSeqs.end(), '\n'), vSeqs.end());
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
			for (int i = 0; i < vHash.size(); ++i)
			{
				sort(vHash[i].begin(), vHash[i].end(), CompDbObj(vSeqs, vLens, m_unMer));
				for (int j = 0; j < vHash[i].size(); ++j)
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
						//nSuff |= (p1[n]>>4) << 12;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
						//nSuff |= (p1[++n]>>4) << 8;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= (m_aCode2Ten[p1[++n]]) << 4;
						//nSuff |= (p1[++n]>>4) << 4;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= (m_aCode2Ten[p1[++n]]);
						//nSuff |= p1[++n] >> 4;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
					}
					else if (m == 3)
					{
						nSuff |= (m_aCode2Ten[p1[n]]) << 12;
						//nSuff |= (p1[n]>>4) << 12;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
						//nSuff |= (p1[++n]>>4) << 8;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= (m_aCode2Ten[p1[++n]]) << 4;
						//nSuff |= (p1[++n]>>4) << 4;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= ONEBYTE;
					}
					else if (m == 2)
					{
						nSuff |= (m_aCode2Ten[p1[n]]) << 12;
						//nSuff |= (p1[n]>>4) << 12;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= (m_aCode2Ten[p1[++n]]) << 8;
						//nSuff |= (p1[++n]>>4) << 8;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= TWOBYTE;
					}
					else if (m == 1)
					{
						nSuff |= (m_aCode2Ten[p1[n]]) << 12;
						//nSuff |= (p1[n]>>4) << 12;
						//sSuff += bitset<4>(p1[n]>>4).to_string();
						nSuff |= THRBYTE;
					}
					else if (m == 0)
					{
						nSuff |= FOUBYTE;
					}
					//for (int l = 0; l < (m<4?m:4); ++l)
					//{
						//cout << bitset<8>(p1[l]) << " ";
					//}
					//cout << endl;
					//for (int l = 0; l < (m<4?m:4); ++l)
					//{
						//cout << bitset<4>(m_aCode2Ten[p1[l]]) << " ";
					//}
					//cout << endl;
					//cout << "\t" << bitset<16>(nSuff) << endl;
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
			for (int i = 0; i < vHash.size(); ++i)
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
	//oa << vWordCnts[m_unTotalIdx/2];
	for (int i = 0; i < vFreq.size(); ++i)
	{
		vFreq[i] /= lnTotalAa;
	}

	oaInfo << nBlock;
	oaInfo << lnSeqNum;
	oaInfo << lnTotalAa;
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


int CHashSearch::BuildQHash(const char* szFile, string& sOutFile, int nQueryType)
{
	char cIdSt = '>';
	char cSeqEd = '>';
	// init m_bSeqType & m_nIdxScl
	if (0 == nQueryType)
	{
		GuessQueryType(szFile);
	}
	else if (1 == nQueryType)
	{
		// nt
        m_bSeqType = true;
		m_nIdxScl = 6;
        printf("Queries are nucleotide sequences in fasta format\n");
	}
	else if (2 == nQueryType)
	{
		// aa
        m_bSeqType = false;
        printf("Queries are protein sequences\n");
	}
	else if (3 == nQueryType)
	{
		// fastq
        m_bSeqType = true;
		m_nIdxScl = 6;
        printf("Queries are nucleotide sequences in fastq format\n");
		cIdSt = '@';
		cSeqEd = '+';
	}

    ifstream ifFile(szFile);
    if (!ifFile.good())
    {
        ifFile.close();
        cout << "can not open the file: " << szFile << endl;
        return -1;
    }

	ofstream of(sOutFile.c_str());
	if (!of.good())
	{
		cout << "can not write files..." << endl;
		return -1;
	}
	archive::binary_oarchive oa(of);

	ifFile.seekg(0, ios::end);
	long int lnSize = ifFile.tellg();
	ifFile.seekg(0, ios::beg);
	int nBlock = lnSize / m_unQSize;
	if (0 != lnSize % m_unQSize)
	{
		++nBlock;
	}

	oa << nBlock;

	map<string, char> mTransTable;
	map<char, char> mComple;
	Seg* seg = NULL;
	Seg* segsht = NULL;
	if (true == m_bSeqType)
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
	
	vector<uint> vLens;
	vector<uchar> vSeqs;
	VNAMES vNames;
	uint unWrited = 0;
	int nSeqNum = 0;

	POOL vPool(m_unQSize, 0);
	uint unLeft = 0;
    while (ifFile.good())
    {
		ifFile.read(&vPool[unLeft], m_unQSize-unLeft);
		int nRead = ifFile.gcount();
		ITER itStop = vPool.begin();
		if (unLeft+nRead < m_unQSize)
		{
			// for last block of data, give it a '>' to split the last sequence
			vPool[unLeft+nRead] = cIdSt;
			advance(itStop, unLeft+nRead+1);
		}
		else
		{
			if (3 == nQueryType)
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
			}
			else
			{
				advance(itStop, unLeft+nRead);
			}
		}

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
				vNames.push_back('>'+string(itSt+1, find(itSt+1, itBeg-1, ' ')));
				vector<char> vS(itBeg, itEd-1);
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
						for (int nn = 0; nn < vS.size(); ++nn)
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
					//cout << &vTran[0] << endl;

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

					for (int i = 0; i < strlen(pMasked); ++i)
					{
						if ('X' == pMasked[i] || 'x' == pMasked[i])
						{
							vTran[i] = pMasked[i];
						}
					}

					delete [] pMasked;

					// a char '\0' was added at the end of vTran, so now ignore it
					vSeqs.insert(vSeqs.end(), vTran.begin(), vTran.end()-1);
					vLens.push_back(vSeqs.size());
				}

				itSt = itEd;
				if (3 == nQueryType)
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
				vNames.push_back(string(itSt, find(itSt+1, itBeg-1, ' ')));
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
		
		//serialize
		oa << vSeqs;
		oa << vLens;
		oa << vNames;

		nSeqNum += vNames.size();
		
		vSeqs.clear();
		vLens.clear();
		vNames.clear();

		unLeft = distance(itSt, vPool.end());
		POOL::reverse_iterator itLast = find(vPool.rbegin(), vPool.rend(), cIdSt);
		copy(itSt, vPool.end(), vPool.begin());
    }

	if (true == m_bSeqType)
	{
		delete seg;
		delete segsht;
	}

	ifFile.close();
	of.close();

	return nSeqNum;
	//return nBlock;
}


void CHashSearch::PrintHash(MINDEX& v)
{
	cout << v.size() << endl;
	vector<int> vCnt(100000, 0);
	uint unTotal = 0;
	uint unNot0 = 0;

	for (int i = 0; i < v.size(); ++i)
	{
		++vCnt[v[i].size()];
		unTotal += v[i].size();
		if (v[i].size() != 0)
		{
			++unNot0;
		}
	}

	cout << "*******************" << endl;
	cout << "total #:\t" << unTotal << endl;
	cout << "not 0:\t" << unNot0 << endl;
	//copy(vCnt.begin(), vCnt.end(), ostream_iterator<int>(cout, "\n"));
	for_each(vCnt.begin(), vCnt.end(), cout<<lambda::_1<<"\n");
}


void CHashSearch::PrintInfo(MINDEX& v)
{
	cout << v.size() << endl;
	vector<int> vCnt(1000000, 0);
	uint unTotal = 0;
	uint unNot0 = 0;

	for (int i = 0; i < v.size(); ++i)
	{
		vCnt[i] = v[i].size();
	}

	cout << "*******************" << endl;
	unTotal = accumulate(vCnt.begin(), vCnt.end(), 0);
/*
	sort(vCnt.begin(), vCnt.end(), greater<int>());
	for (int i = 0; i < vCnt.size(); ++i)
	{
		if (vCnt[i] == 0)
		{
			break;
		}

		unTotal += vCnt[i];
		//cout << vCnt[i] << endl;
	}
*/
	cout << "total #:\t" << unTotal << endl;
}

void CHashSearch::Search(string& sDbPre, int& nDbBlockNum, string& sQPre, int& nQBlockNum)
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

	archive::binary_iarchive iaD(ifD);
	archive::binary_iarchive iaDInfo(ifDInfo);

	vector<double> vFreq;
	VUINT vWordCnts(m_unTotalIdx, 0);
	long int lnTotalAa;
	uint unMedian;
	long int lnSeqNum;

	iaDInfo >> nDbBlockNum;
	iaDInfo >> lnSeqNum;
	iaDInfo >> lnTotalAa;
	iaDInfo >> vWordCnts;
	iaDInfo >> unMedian;
	iaDInfo >> vFreq;

	// set BlastStat
	InitAlignPara(m_bSeqType, lnTotalAa, lnSeqNum, m_nThreadNum);

	pool tp(m_nThreadNum);
	cout << "start " << m_nThreadNum << "  threads" << endl;

	int nSeqNum = nQBlockNum;
	m_vOutIdx.assign(nSeqNum, CIndex());
	m_vM8Idx.assign(nSeqNum, CIndex());

	for (int j = 0; j < nDbBlockNum; ++j)
	{
		// temp file stream
		if (m_ofAln.is_open())
		{
			m_ofAln << m_sOutput;
			m_sOutput.clear();
			m_ofAln.close();

			m_llOutCum = 0;
			ofstream ofOut((m_sOutBase+".aln"+".tmp"+lexical_cast<string>(j-1)+".idx").c_str());
			archive::binary_oarchive oaOut(ofOut);
			oaOut << m_vOutIdx;
			ofOut.close();
		}
		if (m_ofM8.is_open())
		{
			m_ofM8 << m_sM8;
			m_sM8.clear();
			m_ofM8.close();

			m_llM8Cum = 0;
			ofstream ofM8((m_sOutBase+".m8"+".tmp"+lexical_cast<string>(j-1)+".idx").c_str());
			archive::binary_oarchive oaM8(ofM8);
			oaM8 << m_vM8Idx;
			ofM8.close();
		}
		for (int nn = 0; nn < nSeqNum; ++nn)
		{
			m_vOutIdx[nn].m_llBeg = 0;
			m_vOutIdx[nn].m_nSize = 0;
			m_vM8Idx[nn].m_llBeg = 0;
			m_vM8Idx[nn].m_nSize = 0;
		}
		m_nSeqBase = 0;

		m_ofAln.open((m_sOutBase+".aln"+".tmp"+lexical_cast<string>(j)).c_str());
		m_ofM8.open((m_sOutBase+".m8"+".tmp"+lexical_cast<string>(j)).c_str());

		// read db file and store info
		MINDEX vDHash(m_unTotalIdx, VUINT()); // all k-mer of database
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

		ifstream ifQ(sQPre.c_str());
		if (!ifQ.good())
		{
			cout << "can not open query index file" << endl;
			return;
		}
		archive::binary_iarchive iaQ(ifQ);
		iaQ >> nQBlockNum;

		vector<uint> vQLens;
		vector<uchar> vQSeqs;
		VNAMES vQNames;

		for (int i = 0; i < nQBlockNum; ++i)
		{
			// read q file and store info
			vQSeqs.clear();
			vQLens.clear();
			vQNames.clear();
			iaQ >> vQSeqs;
			iaQ >> vQLens;
			iaQ >> vQNames;
			// construct query package
			CQrPckg Query(vQSeqs, vQLens, vQNames);

			m_unTotalQuery += vQLens.size() - 1;
			
			// generate results
			for (int k = 0; k < vQLens.size()-1; k+=m_nIdxScl)
			{
				tp.schedule(bind(&CHashSearch::Searching, this, k, Query, Db));
				//Searching(k, Query, Db);
			}

			tp.wait();

			m_nSeqBase += vQNames.size();
		}
		
		ifQ.close();
	}
	tp.wait();

	if (m_ofAln.is_open())
	{
		m_ofAln << m_sOutput;
		m_ofAln.close();
		m_llOutCum = 0;

		ofstream ofOut((m_sOutBase+".aln"+".tmp"+lexical_cast<string>(nDbBlockNum-1)+".idx").c_str());
		archive::binary_oarchive oaOut(ofOut);
		oaOut << m_vOutIdx;
		ofOut.close();
		m_vOutIdx.clear();
	}
	if (m_ofM8.is_open())
	{
		m_ofM8 << m_sM8;
		m_ofM8.close();
		m_llM8Cum = 0;

		ofstream ofM8((m_sOutBase+".m8"+".tmp"+lexical_cast<string>(nDbBlockNum-1)+".idx").c_str());
		archive::binary_oarchive oaM8(ofM8);
		oaM8 << m_vM8Idx;
		ofM8.close();
		m_vM8Idx.clear();
	}


	// merge nDbBlockNum temp results
	MergeRes(nDbBlockNum, sQPre);

	ifD.close();
	ifDInfo.close();
}


void CHashSearch::Process(char* szDBFile, char* szQFile, char* szOFile, double dLogEvaThr, int nMaxOut, int nMaxM8, int nQueryType, bool bPrintEmpty, bool bGapExt, bool bAcc, bool bHssp, int nMinLen, uint unDSize, uint unQSize, uint unMer)
{
	m_dLogEvaThr = dLogEvaThr;
	if (nMaxOut == 0)
	{
		m_nMaxOut = LLONG_MAX;
	}
	else
	{
		m_nMaxOut = nMaxOut;
	}
	if (nMaxM8 == 0)
	{
		m_nMaxM8 = LLONG_MAX;
	}
	else
	{
		m_nMaxM8 = nMaxM8;
	}
	m_bPrintEmpty = bPrintEmpty;
	m_bGapExt = bGapExt;
	m_bAcc = bAcc;
	m_bHssp = bHssp;
	m_nMinLen = nMinLen;

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

	m_sOutBase.assign(szOFile);

	m_ofM8.open((m_sOutBase+".m8").c_str());
	m_ofM8 << "# RAPSearch\n# Job submitted: ";
	time_t rawtime;
	struct tm* timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	m_ofM8 << asctime(timeinfo);
	m_ofM8 << "# Query : " << szQFile << endl;
	m_ofM8 << "# Subject : " << szDBFile << endl;
	m_ofM8 << "# Fields: Query\tSubject\tidentity\taln-len\tmismatch\tgap-openings\tq.start\tq.end\ts.start\ts.end\tlog(e-value)\tbit-score\n";
	m_ofM8.close();

	// bitset for build query hash
	int nDbBlockNum = 1;
	//string sTemp(szDBFile);
	//string sDbOut = "db/dbfile_" + sTemp.substr(sTemp.rfind('/')+1);
	//string sDbOut = "/home/yongzhao/Downloads/HashSearch/db/dbfile" + lexical_cast<string>(m_unMer);
	//nDbBlockNum = BuildDHash(szDBFile, sDbOut);

	string sTemp = szQFile;
	string sQOut  = sTemp + ".tmp";
	
	// test if another rapserach is using the query file
	// if yes, create another temp file
	int i = 0;
	ifstream ifExitTest((sQOut+lexical_cast<string>(i)).c_str());
	while (!ifExitTest.fail())
	{
		ifExitTest.close();
		++i;
		ifExitTest.open((sQOut+lexical_cast<string>(i)).c_str());
	}
	sQOut += lexical_cast<string>(i);
	//string sQOut  = "qr/qrfile_" + sTemp.substr(sTemp.rfind('/')+1);
	//string sQOut  = "qr/qrfile_" + sTemp.substr(sTemp.rfind('/')+1);
	
	//using namespace boost::chrono;
	//thread_clock::time_point st = thread_clock::now();
	int nQBlockNum = BuildQHash(szQFile, sQOut, nQueryType);
	//thread_clock::time_point ed = thread_clock::now();
	//cout << "\t" << "\t" << "\tbuilding query time:\t" << duration_cast<minutes>(ed-st).count() << " ms" << endl;

	string sDbOut(szDBFile);
	Search(sDbOut, nDbBlockNum, sQOut, nQBlockNum);
	remove(sQOut.c_str());
	//cout << "total gap extension:\t" << m_unGapExt << endl;
}


void CHashSearch::Process(char* szDBFile, char* szDbHash, int nSplitNum, uint unMer)
{
	m_unMer = unMer;
	//m_unDSize = unDSize;
	m_unTotalIdx = lexical_cast<uint>(pow(10.0, int(m_unMer)));

	string sDbOut(szDbHash);
	BuildDHash(szDBFile, sDbOut, nSplitNum);
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
	for (int i = 0; i < m_vBlastPt.size(); ++i)
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

		// for consistence with swift
		uint unPrvSdLen = 6;

		for (int i = 0; i < unQLen - m_unMer; ++i)
		{
			uint unCnt = 0;
		//thread_clock::time_point st = thread_clock::now();
			// pick seed length
			uint unQSeedBeg = QrAln.m_unSeedBeg = i;
			int nSeed = Tran2Ten(QrAln);
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
						if((unIdx = m_aCode2Ten[pQ[unQSeedBeg+m_unMer+unIncr-1]]) != m_uMask)
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
				for (int i = m_unMer; i < unLocalSeed; ++i)
				{
					if((unIdx = m_aCode2Ten[pQ[unQSeedBeg+i]]) == m_uMask)
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
			for (int idx = 0; idx < vExtra.size(); ++idx)
			{
				vExtra[idx] = m_aCode2Ten[vExtra[idx]];
			}

			if (!Db.m_vHash[nSeed].empty())
			{
				int nCnt = ExtendSeq2Set(nSeed, unLocalSeed, vExtra,
						nQrIdx, QrAln, nQOriLen,
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
			for (int j = unQSeedBeg+unLocalSeed; j < unQSeedBeg+m_unMutSeedLen; ++j)
			{
				if((unIdx = m_aCode2Ten[pQ[j]]) == m_uMask)
				{
					break;
				}
			}
			if (m_uMask == unIdx)
			{
				continue;
			}

			vExtra.assign(pQ+unQSeedBeg+m_unMer, pQ+unQSeedBeg+m_unMutSeedLen);
			for (int idx = 0; idx < vExtra.size(); ++idx)
			{
				vExtra[idx] = m_aCode2Ten[vExtra[idx]];
			}

			for (int m = 0; m < m_vMutation.size(); ++m)
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
							nQrIdx, QrAln, nQOriLen,
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
							nQrIdx, QrAln, nQOriLen,
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
			int i = m_unMer;
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
			int i = m_unMer;
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
		int nQSeqIdx, CAlnPckg& QrAln, int nQOriLen,
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
		for (int i = 0; i < vExtra.size(); ++i)
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
		uint unDLen, unDSeedBeg;
		uchar* pD = GetSeq(Db.m_vSeqs, Db.m_vLens, Db.m_vNames, vDSet[j], unDLen, unDSeedBeg);
		
		// no enough letters
		if (unDLen < unDSeedBeg + unLocalSeedLen)
		{
			continue;
		}
		// have the same previous char, so don't extend them at this time
		if (m_bAcc==false && 0 != QrAln.m_unSeedBeg && m_uMask != QrAln.m_pSeq[QrAln.m_unSeedBeg-1]
			&& 0 != unDSeedBeg && m_uMask != pD[unDSeedBeg-1]
			&& m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg-1]] == m_aCode2Ten[pD[unDSeedBeg-1]]
			&& 4 != vExtra.size()) // if this is a mutation case, do not allow to ignore it
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
			CalRes(nQSeqIdx, QrAln.m_pSeq, nQOriLen, QrAln.m_unSeedBeg, nDSeqIdx, DbAln.m_pSeq, DbAln.m_unSeedBeg, vDNames, unLocalCopy, stAlnmnt, cout, mRes, nTreadID);
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
	int   i, j, k, l, s, maxs, ma;

	i = j = l = 0;
    ma = 0;
    maxs = s = score0;
    *extl = 0;
    *match = 0;
    while(i < len_queryseq && j < len_dataseq && s >= MINSCORE && s >= maxs - UngapExtDrop)
    {
		// uncompleted
        //if(atypical(queryseq[i]) && atypical(dataseq[j]))	break; //to-be-tested!
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
    int   i, j, k, l, s, maxs, ma;

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
        //if(atypical(queryseq[i]) && atypical(dataseq[j]))	break; //to-be-tested!
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
    int	maxs_row = -1000;
    int	badscore = -1000;
    int	g = GapIni;
    int	h = GapExt;
    int	m = g + h; //gap-create + gap-extend
    int	maxs, E1, E2, match;
    int	ib;
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
        //if(j == 1)
        //{
        //    trace[0][j] = etrace[0][j] = 'E';
        //    dtrace[0][j] = 'D';
        //}
        //else
        //{
        //    trace[0][j] = etrace[0][j] = 'e';
        //    dtrace[0][j] = 'd';
        //}
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


void CHashSearch::CalRes(int nQIdx, uchar* pQ, int nQOriLen, uint unQSeedBeg/*, string& sQName*/, int nDIdx, uchar* pD, uint unDSeedBeg, VNAMES& vDNames, uint unLocalSeedLen, STAlnmnt& stAlnmnt, ostream& out, MRESULT& mRes, int nTreadID)
{
	double dEValue = m_vpBlastSig[nTreadID]->rawScore2ExpectLog(stAlnmnt.nScore);
	if (dEValue > 0)
	{
		dEValue = floor(dEValue*100+0.5) / 100;
	}
	else
	{
		dEValue = floor(dEValue*100-0.5) / 100;
	}
	double dBits = m_vpBlastSig[nTreadID]->rawScore2Bit(stAlnmnt.nScore);
	dBits = floor(dBits*100+0.5) / 100;
	
	int nTotGap = 0;
	int nGapOpen = 0;
	int nTotAlnLen = 0;

	for (int i = 0; i < stAlnmnt.vMode.size(); ++i)
	{
		nTotAlnLen += stAlnmnt.vLen[i];
		if ('s' != stAlnmnt.vMode[i])
		{
			++nGapOpen;
			nTotGap += stAlnmnt.vLen[i];
		}
	}

	// evalue criteria
	if (m_bHssp == false && !(stAlnmnt.nScore>SUMHSP_MINRAWSCORE || dEValue<m_dLogEvaThr))
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
	int nAllc = nTotAlnLen>unLocalSeedLen?nTotAlnLen:unLocalSeedLen;
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
		for (int i = 0; i < stAlnmnt.vMode.size(); ++i)
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
	for (int i = 0; i < vQ.size(); ++i)
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
		STResult& st = (*it).second;
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
			//st.sSeedPair = sSeedPair;
			st.sQ = sQ;
			st.sInfo = sInfo;
			st.sD = sD;
		}
	}
	else
	/****************************************************************/
	{
		MRESULT::iterator itTmp = mRes.insert(it, MRESULT::value_type(pair<int, int>(nQIdx/m_nIdxScl, nDIdx), STResult()));
		STResult& st = (*itTmp).second;
		st.nDbIdx =	nDIdx; 
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


void CHashSearch::FindSeeds(int nSeed, uint unIncr, CAlnPckg& QrAln, MINDEX& vDHash, VUINT& vRes, bool bMutation)
{
	// there are no enough letters
	if (QrAln.m_unLen < QrAln.m_unSeedBeg+m_unMer+unIncr)
	{
		return;
	}
	// if there is non-aa letter in the least length, return
	for (int i = QrAln.m_unSeedBeg+m_unMer; i < QrAln.m_unSeedBeg+m_unMer+unIncr; ++i)
	{
		if (m_uMask == QrAln.m_pSeq[i])
		{
			return;
		}
	}

	// # represented by the previous letter. If it is not aa, # should be m_uMask
	int nPrevNum = m_uMask;
	if (0 != QrAln.m_unSeedBeg)
	{
		nPrevNum = m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg-1]];
	}

	// index represented by 6mer starting from the previous letter. If it is valid, it should be -1
	int nOrgSeed = nSeed;
	int nPrevSeed = -1;
	int nScale = m_unTotalIdx / 10;
	if (m_uMask != nPrevNum)
	{
		nPrevSeed = nPrevNum*nScale + nOrgSeed/10;
	}

	VUINT vTemp;
	vTemp.reserve(vDHash[nSeed].size());
	vRes.reserve(vDHash[nSeed].size());
	if (bMutation == true || 0 == QrAln.m_unSeedBeg || -1 == nPrevSeed || nPrevSeed == nSeed || 0 == vDHash[nPrevSeed].size())
	{
		// the basic seed set is db.hash[nSeed]
		vRes.assign(vDHash[nSeed].begin(), vDHash[nSeed].end());
		
		if (vRes.empty())
		{
			return;
		}
	}
	else
	{
		// the basic seed set is db.hash[nSeed]-db.hash[nPrevSeed]
		/*
		set_difference(vDHash[nSeed].begin(), vDHash[nSeed].end(),
				vDHash[nPrevSeed].begin(), vDHash[nPrevSeed].end(),
				back_inserter(vRes),
				mySetComp);
		*/
		SetDiff(vDHash[nSeed], vDHash[nPrevSeed], -1, vRes);
		
		if (vRes.empty())
		{
			return;
		}
	}

	for (int i = QrAln.m_unSeedBeg+m_unMer; i < QrAln.m_unSeedBeg+m_unMer+unIncr; ++i)
	{
		int nNextNum = m_aCode2Ten[QrAln.m_pSeq[i]];
		// nSeed includes the mutation, so do not get prevNum from the sequence
		nSeed = (nSeed%nScale)*10 + nNextNum;
		if (vDHash[nSeed].empty())
		{
			vRes.clear();
			return;
		}

		/*
		set_intersection(vRes.begin(), vRes.end(),
				vDHash[nSeed].begin(), vDHash[nSeed].end(),
				back_inserter(vTemp),
				mySetComp);
		*/
		SetInter(vRes, vDHash[nSeed], i-QrAln.m_unSeedBeg-m_unMer+1, vTemp);
		
		vTemp.swap(vRes);
		vTemp.clear();
		if (vRes.empty())
		{
			return;
		}
	}

}

struct CompFrameObj
{
	bool operator() (const STResult& st1, const STResult& st2) const
	{
		return st1.nFrame < st2.nFrame;
	}
}compFrame;

struct CompEvalueObj
{
	bool operator() (const STResult& st1, const STResult& st2) const
	{
		return st1.dEValue < st2.dEValue;
	}
}compEvalue;

struct CompQSt
{
	bool operator() (const STResult& st1, const STResult& st2) const
	{
		return st1.nQSt < st2.nQSt;
	}
}compQSt;

struct CompObj
{
	bool operator() (const STResult& st1, const STResult& st2) const
	{
		return st1.dEValue < st2.dEValue;
	}
}myComp;


void CHashSearch::PrintRes(MRESULT& mRes, int nTreadID, CQrPckg& Query, CDbPckg& Db)
{
	if (mRes.empty())
	{
		return;
	}

	MIT it = mRes.begin();
	int nQrIdx = (*it).first.first;
	MRESULT::iterator itFind = mRes.end();
	vector<STResult> vTemp;
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

	sort(vTemp.begin(), vTemp.end(), myComp);

	vector<STResult>::iterator itPrint;
	if (vTemp.size() > m_nMaxOut)
	{
		itPrint = vTemp.begin() + m_nMaxOut;
	}
	else
	{
		itPrint = vTemp.end();
	}

	stringstream sOutput;
	stringstream sM8;
	int nBegStrAligned = 6;
	for (vector<STResult>::iterator itSt = vTemp.begin(); itSt != itPrint; ++itSt)
	{
		STResult& st= *itSt;
		if (st.dEValue >= m_dLogEvaThr)
		{
			break;
		}
		// note: here, the hits are stored according to it's real query index, not 1->6 frame query index

		st.sInfo.insert(0, 7, ' ');

		string sQNum = lexical_cast<string>(st.nQBeg);
		st.sQ = string(nBegStrAligned-sQNum.size(), ' ') + sQNum + " " + st.sQ + " " + lexical_cast<string>(st.nQEnd);

		++st.nDSt;
		++st.nDEd;
		int nFac = 0;
		while (0 <= (st.nDbIdx-nFac-1) && Db.m_vNames[st.nDbIdx] == Db.m_vNames[st.nDbIdx-nFac-1])
		{
			++nFac;
		}
		string sDNum = lexical_cast<string>(1848*nFac+st.nDSt);
		st.sD = string(nBegStrAligned-sDNum.size(), ' ') + sDNum + " " + st.sD + " " + lexical_cast<string>(1848*nFac+st.nDEd);

		sOutput << Query.m_vNames[nQrIdx];
		sOutput << " vs ";
		sOutput << Db.m_vNames[st.nDbIdx];
		//sOutput << " score=" << st.nScore;
		sOutput << " bits="	<< st.dBits;
		sOutput << " log(E-value)="	<< st.dEValue;
		sOutput << " identity=" << st.dIdent << "%";
		sOutput << " aln-len="	<< st.nAlnLen;
		sOutput << " mismatch="	<< st.nMismatch;
		sOutput << " gap-openings="	<< st.nGapOpen;
		sOutput << " nFrame=" << st.nFrame;
		sOutput << "\n";
		//sOutput << "Seed Pair:\t" << st.sSeedPair << "\n";
		sOutput << "Query:\t" << st.sQ << "\n";
		sOutput << "      \t" << st.sInfo << "\n";
		sOutput << "Sbjct:\t" << st.sD << "\n";
		sOutput << "\n";

		sM8 << Query.m_vNames[nQrIdx] 
			<< "\t" << Db.m_vNames[st.nDbIdx]
			<< "\t" << st.dIdent
			<< "\t"	<< st.nAlnLen
			<< "\t"	<< st.nMismatch
			<< "\t"	<< st.nGapOpen 
			<< "\t" << st.nQBeg
			<< "\t" << st.nQEnd
			<< "\t" << st.nDSt
			<< "\t" << st.nDEd
			<< "\t"	<< st.dEValue 
			<< "\t"	<< st.dBits 
			<< endl;
	}

	vector<STResult>::iterator itM8;
	if (vTemp.size() > m_nMaxM8)
	{
		itM8 = vTemp.begin() + m_nMaxM8;
	}
	else
	{
		itM8 = vTemp.end();
	}

	for (vector<STResult>::iterator itSt = itPrint; itSt != itM8; ++itSt)
	{
		STResult& st= *itSt;
		if (st.dEValue >= m_dLogEvaThr)
		{
			break;
		}

		sM8 << Query.m_vNames[nQrIdx] 
			<< "\t" << Db.m_vNames[st.nDbIdx]
			<< "\t" << st.dIdent
			<< "\t"	<< st.nAlnLen
			<< "\t"	<< st.nMismatch
			<< "\t"	<< st.nGapOpen 
			<< "\t" << st.nQBeg
			<< "\t" << st.nQEnd
			<< "\t" << st.nDSt
			<< "\t" << st.nDEd
			<< "\t"	<< st.dEValue 
			<< "\t"	<< st.dBits 
			<< endl;
	}

	//mutex::scoped_lock lock(muMonitor);
	muMonitor.lock();

	long long llBeg = m_llOutCum + m_sOutput.size();
	m_sOutput += sOutput.str();
	int nSize = sOutput.str().size();
	m_vOutIdx[m_nSeqBase+nQrIdx].m_llBeg = llBeg;
	m_vOutIdx[m_nSeqBase+nQrIdx].m_nSize = nSize;
	//m_vOutIdx.push_back(CIndex(nQrIdx, llBeg, nSize));
	if (m_sOutput.size() > 50000000)
	{
		m_ofAln << m_sOutput;
		m_llOutCum += m_sOutput.size();
		m_sOutput.clear();
	}

	llBeg = m_llM8Cum + m_sM8.size();
	m_sM8 += sM8.str();
	nSize = sM8.str().size();
	m_vM8Idx[m_nSeqBase+nQrIdx].m_llBeg = llBeg;
	m_vM8Idx[m_nSeqBase+nQrIdx].m_nSize = nSize;
	//m_vM8Idx.push_back(CIndex(nQrIdx, llBeg, nSize));
	if (m_sM8.size() > 5000000)
	{
		m_ofM8 << m_sM8;
		m_llM8Cum += m_sM8.size();
		m_sM8.clear();
	}

	muMonitor.unlock();

	mRes.clear();
}


void CHashSearch::SumEvalue(vector<STResult>& v, int nSt, int nEd, int nLen, int nTreadID)
{
	typedef vector<STResult>::iterator STIT;
	STIT itSt = v.begin() + nSt;
	STIT itEd = v.begin() + nEd;
	// sort by nFrame
	sort(itSt, itEd, compFrame);
	STResult st;
	st.nFrame = 3;
	STIT itDir = lower_bound(itSt, itEd, st, compFrame);
	// if there are more than one hit in one direction
	int nDisPos = distance(itSt, itDir);
	int nDisNeg = distance(itDir, itEd);
	if (nDisPos > 1 || nDisNeg > 1)
	{
		vector<STResult> vRes;
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
				if (itStart->dEValue < m_dLogEvaThr)
				{
					vRes.push_back(*itStart);
				}
				itStart = itEnd;
				itEnd = itEd;
				continue;
			}
			// sort by score and start position of query
			sort(itStart, itEnd, compQSt);
			stable_sort(itStart, itEnd, compEvalue);
			// check overlap and logevalue
			double dMaxEvalue = 1000;
			vector<STResult> vNew;
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
				if (vNew[0].dEValue < m_dLogEvaThr)
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
				//double dSumEvalue = log(m_vpBlastSig[nTreadID]->sumScore2Expect(nNo, aRawScore, nLen)) / LOG10;
				double dTmp = m_vpBlastSig[nTreadID]->sumScore2Expect(nNo, aRawScore, nLen);
				double dSumEvalue = -10000.00;
				if (0 != dTmp)
				{
					dSumEvalue = log(dTmp) / LOG10;
				}
				// modify the logevalue
				if (dSumEvalue < m_dLogEvaThr)
				{
					for (int i = 0; i < vNew.size(); ++i)
					{
						vNew[i].dEValue = dSumEvalue;
					}
					vRes.insert(vRes.end(), vNew.begin(), vNew.end());
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
	FILE* pf = fopen(szDBFile, "rt");
    if (!pf)
    {
        fprintf(stderr, "Sequence file: %s not found\n", szDBFile);
        exit(-1);
    }
	// bug? some sequences are longer than 1000
    char	str[1001];
    while(!feof(pf))
    {
        fgets(str, 1000, pf);
        if(str[0] == '>') 
		{
			lnSeqNum += 1;
		}
        else
        {
            lnAaNum += strlen(str);
        }
    }
	fclose(pf);
}



void CHashSearch::MergeRes(int nDbBlockNum, string& sQPre)
{
	m_ofAln.open((m_sOutBase+".aln").c_str());
	m_ofM8.open((m_sOutBase+".m8").c_str(), ios::ate|ios::app);

	int nSize = ((1<<31)-1) / nDbBlockNum;

	vector<CSortUnit> v;
	vector<int> vLast(nDbBlockNum);
	vector<CMergeUnit*> vMergeUnit;

	for (int i = 0; i < nDbBlockNum; ++i)
	{
		string sName = m_sOutBase+".aln"+".tmp"+lexical_cast<string>(i);
		CMergeUnit* p = new CMergeUnit(sName.c_str(), nSize, true);
		vMergeUnit.push_back(p);
	}
	
	int nLastIdx = 0;
	for (int i = 0; i < nDbBlockNum; ++i)
	{
		nLastIdx = max(nLastIdx, vMergeUnit[i]->GetLast());
	}

	/**************************************************************/
	VNAMES vQAllNames;
	if (true == m_bPrintEmpty)
	{
		ifstream ifQ(sQPre.c_str());
		if (!ifQ.good())
		{
			cout << "can not open query index file" << endl;
			return;
		}
		archive::binary_iarchive iaQ(ifQ);
		int nQBlockNum;
		iaQ >> nQBlockNum;

		vector<uint> vQLens;
		vector<uchar> vQSeqs;
		VNAMES vQNames;

		for (int i = 0; i < nQBlockNum; ++i)
		{
			// read q file and store info
			vQSeqs.clear();
			vQLens.clear();
			vQNames.clear();
			iaQ >> vQSeqs;
			iaQ >> vQLens;
			iaQ >> vQNames;

			vQAllNames.insert(vQAllNames.end(), vQNames.begin(), vQNames.end());
		}
		
		ifQ.close();
	}
	/**************************************************************/
	for (int i = 0; i < nLastIdx; ++i)
	{
		for (int j = 0; j < nDbBlockNum; ++j)
		{
			vMergeUnit[j]->Update(i, v);
		}

		sort(v.begin(), v.end());

		int n = min(m_nMaxOut, (long long)v.size());
		/***********************************************************/
		if (true == m_bPrintEmpty && 0 == n)
		{
			m_ofAln << vQAllNames[i] << "\tNO HIT" << endl << endl;
		}
		/***********************************************************/
		for (int j = 0; j < n; ++j)
		{
			m_ofAln << v[j].m_sHit;
		}

		v.clear();
	}

	for (int i = 0; i < nDbBlockNum; ++i)
	{
		delete vMergeUnit[i];

		string sName = m_sOutBase+".m8"+".tmp"+lexical_cast<string>(i);
		CMergeUnit* p = new CMergeUnit(sName.c_str(), nSize, false);
		vMergeUnit[i] = p;
	}

	nLastIdx = 0;
	for (int i = 0; i < nDbBlockNum; ++i)
	{
		nLastIdx = max(nLastIdx, vMergeUnit[i]->GetLast());
	}

	for (int i = 0; i < nLastIdx; ++i)
	{
		for (int j = 0; j < nDbBlockNum; ++j)
		{
			vMergeUnit[j]->Update(i, v);
		}

		sort(v.begin(), v.end());

		int n = min(m_nMaxM8, (long long)v.size());
		for (int j = 0; j < n; ++j)
		{
			m_ofM8 << v[j].m_sHit;
		}

		v.clear();
	}

	for (int i = 0; i < nDbBlockNum; ++i)
	{
		delete vMergeUnit[i];
	}

	m_ofAln.close();
	m_ofM8.close();
}


/*
void CHashSearch::MergeRes(int nDbBlockNum)
{
	m_ofAln.open((m_sOutBase+".aln").c_str());
	m_ofM8.open((m_sOutBase+".m8").c_str());

	int nSize = ((1<<31)-1) / nDbBlockNum;

	vector<CSortUnit> v;
	vector<string> vLast(nDbBlockNum);
	vector<CMergeUnit*> vMergeUnit;

	for (int i = 0; i < nDbBlockNum; ++i)
	{
		string sName = m_sOutBase+".aln"+".tmp"+lexical_cast<string>(i);
		CMergeUnit* p = new CMergeUnit(sName.c_str(), nSize, true);
		vMergeUnit.push_back(p);
	}
	
	for (int i = 0; i < nDbBlockNum; ++i)
	{
		vLast[i] = vMergeUnit[i]->Update(v);
	}

	while (v.size() > 0)
	{
		string sMin = *min_element(vLast.begin(), vLast.end());

		sort(v.begin(), v.end());

		string s = v[0].m_sQr;
		uint unLeft = 0;
		int nCnt = 0;
		int i = 0;
		for (i = 0; i < v.size(); ++i)
		{
			if (s == v[i].m_sQr && nCnt < 100)
			{
				m_ofAln << v[i].m_sHit;
				++nCnt;
			}
			else if (s != v[i].m_sQr)
			{
				s = v[i].m_sQr;
				if (s < sMin)
				{
					nCnt = 1;
					m_ofAln << v[i].m_sHit;
				}
				else
				{
					unLeft = distance(v.begin()+i, v.end());
					break;
				}
			}

		}

		copy(v.begin()+i, v.end(), v.begin());
		v.erase(v.begin()+unLeft, v.end());

		for (int i = 0; i < nDbBlockNum; ++i)
		{
			vLast[i] = vMergeUnit[i]->Update(v);
		}
	}

	v.clear();

	for (int i = 0; i < nDbBlockNum; ++i)
	{
		delete vMergeUnit[i];

		string sName = m_sOutBase+".m8"+".tmp"+lexical_cast<string>(i);
		CMergeUnit* p = new CMergeUnit(sName.c_str(), nSize, false);
		vMergeUnit[i] = p;
	}
	
	for (int i = 0; i < nDbBlockNum; ++i)
	{
		vLast[i] = vMergeUnit[i]->Update(v);
	}

	while (v.size() > 0)
	{
		string sMin = *min_element(vLast.begin(), vLast.end());

		sort(v.begin(), v.end());

		string s = v[0].m_sQr;
		uint unLeft = 0;
		int nCnt = 0;
		int i = 0;
		for (i = 0; i < v.size(); ++i)
		{
			if (s == v[i].m_sQr && nCnt < 500)
			{
				m_ofM8 << v[i].m_sHit;
				++nCnt;
			}
			else if (s != v[i].m_sQr)
			{
				s = v[i].m_sQr;
				if (s < sMin)
				{
					nCnt = 1;
					m_ofM8 << v[i].m_sHit;
				}
				else
				{
					unLeft = distance(v.begin()+i, v.end());
					break;
				}
			}

		}

		copy(v.begin()+i, v.end(), v.begin());
		v.erase(v.begin()+unLeft, v.end());

		for (int i = 0; i < nDbBlockNum; ++i)
		{
			vLast[i] = vMergeUnit[i]->Update(v);
		}
	}

	for (int i = 0; i < nDbBlockNum; ++i)
	{
		delete vMergeUnit[i];
	}

	m_ofAln.close();
	m_ofM8.close();
}
*/

//----------------------------------------------------------------------
void CHashSearch::GuessQueryType(const char* szFile)
{
	FILE* fp = fopen(szFile, "rt");
    char	seq[10001], str[1001];
    while (!feof(fp))
    {
        str[0] = 0;
        fgets(str, 1000, fp);
		// add by yongzhao, get rid off '\n'
		str[strlen(str)-1] = 0;
        if(str[0] == '>' || str[0] == '\n') continue;
        else strcat(seq, str);
        if(strlen(seq) >= 9000) break;
    }
    rewind(fp);
    char	nuc[] = "ATCGatcgUu";
    int	i, j;
    int	add = 0;
    for(i = 0; i < strlen(seq); i ++)
    {
        for(j = 0; j < 10; j ++)
        {
            if(seq[i] == nuc[j]) break;
        }
        if(j < 10) add += 1;
    }

    if(add > strlen(seq) * 0.95)
    {
        m_bSeqType = true;
        printf("Queries are nucleotide sequences\n");
    }
    else
    {
        m_bSeqType = false;
        printf("Queries are protein sequences\n");
    }

	// used for convert index from with-fram to non-frame
	if (m_bSeqType == true)
	{
		m_nIdxScl = 6;
	}

	fclose(fp);
}
