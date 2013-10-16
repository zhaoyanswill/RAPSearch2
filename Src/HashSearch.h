#ifndef __HASHSEARCH_H__
#define __HASHSEARCH_H_

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <bitset>
#include <algorithm>
#include <boost/thread/thread.hpp>
#include "BlastStat.h"
#include "paras.h"
#include "cindex.h"
using namespace std;


typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

typedef vector<char> POOL;
typedef POOL::iterator ITER;

typedef vector<uint> VUINT;
typedef vector<uchar> VUCHAR;
typedef vector<ushort> VUSHORT;
typedef vector<VUINT> MINDEX;
typedef vector<string> VNAMES;
typedef vector<VUSHORT> VCOMP;


// query package including sequences, length, and names
class CQrPckg
{
public:
	CQrPckg(VUCHAR& vSeqs, VUINT& vLens, VNAMES& vNames)
		: m_vSeqs(vSeqs), m_vLens(vLens), m_vNames(vNames) {}

public:
	VUCHAR& m_vSeqs;
	VUINT& m_vLens;
	VNAMES& m_vNames;
};


// query package including sequences, length, names, hash table, and statistical information
class CDbPckg
{
public:
	CDbPckg(MINDEX& vHash, VUCHAR& vSeqs, VUINT& vLens, VNAMES& vNames, VCOMP& vComp, vector<double>& vFreq, VUINT& vWordCnts, uint& unMedian)
		: m_vHash(vHash), m_vSeqs(vSeqs), m_vLens(vLens), m_vNames(vNames), m_vComp(vComp), m_vFreq(vFreq), m_vWordCnts(vWordCnts), m_unMedian(unMedian) {}

public:
	MINDEX& m_vHash;
	VUCHAR& m_vSeqs;
	VUINT& m_vLens;
	VNAMES& m_vNames;
	VCOMP& m_vComp;
	vector<double>& m_vFreq;
	VUINT& m_vWordCnts;
	uint m_unMedian;
};


// alignment package including sequence, length and the starting position of the seed
class CAlnPckg
{
public:
	CAlnPckg(uchar* p, uint unLen, uint unSeedBeg)
		: m_pSeq(p), m_unLen(unLen), m_unSeedBeg(unSeedBeg) {}

public:
	uchar* m_pSeq;
	uint m_unLen;
	uint m_unSeedBeg;
};


// structure including the information of a hit
typedef struct STALNMNT
{
	int nScore;
	int nMatch;
	int nQFwd;
	int nDFwd;
	int nQBwd;
	int nDBwd;
	vector<char> vMode;
	vector<int> vLen;
}STAlnmnt;


// structure including the result information of a hit
class STResult
{
public:
	int nDbIdx;
	int nScore;
	double dBits;
	double dEValue;
	double dIdent;
	int nAlnLen;
	int nMismatch;
	int nGapOpen;
	int nFrame;
	int nQSt;
	int nQEd;
	int nQBeg;
	int nQEnd;
	int nDSt;
	int nDEd;
	string sQ;
	string sInfo;
	string sD;
};


/* map for store hit results */
typedef pair<int, int> PIDX;
struct pairComp
{
	bool operator() (const PIDX& lhs, const PIDX& rhs) const
	{
		if (lhs.first != rhs.first)
		{
			return lhs.first < rhs.first;
		}
		else
		{
			return lhs.second < rhs.second;
		}
	}
};
typedef multimap<PIDX, STResult, pairComp> MRESULT;
typedef MRESULT::iterator MIT;

/* the class for indexing, searching */
class CHashSearch
{
public:
	CHashSearch(int nThreadNum);
	~CHashSearch(void)
	{
		if (NULL != m_pBlastSig)
		{
			delete m_pBlastSig;
		}

		//m_ofAln << m_sOutput;
		//m_ofM8 << m_sM8;

		//m_ofAln.close();
		//m_ofM8.close();
	}

	// do the protein database search
	void Process(char* szDBFile, char* szQFile, char* szOFile, double dLogEvaThr, int nMaxOut, int nMaxM8, int nQueryTypeq, bool bPrintEmpty, bool bGapExt, bool bAcc, bool bHssp, int nMinLen, uint unDSize = 300000000, uint unQSize = 50000000, uint unMer = 6);
	// indexing the database
	void Process(char* szDBFile, char* szDbHash, int nSplitNum = 0, uint unMer = 6);

private:
	// read file and build the hash table
	int BuildDHash(const char* szDbFile, string& sOutFile, int nSplitNum); 
	// read file and build the hash table
	int BuildQHash(const char* szQrFile, string& sOutFile, int nQueryType); 
	// probe that what is the type of query
	void GuessQueryType(const char* szFile);
	// char -> compressed code
	void Encode(const string& sIn, vector<uchar>& v); 
	// char -> compressed code
	void Encode(VUCHAR& v); 
	// char -> compressed code & count DB
	long int Encode(VUCHAR& v, vector<double>& vFreq); 
	// compressed code -> char
	void Decode(const vector<uchar>& v, string& sOut); 
	// compressed code -> char
	void Decode(vector<uchar>& v); 
	// char -> 10-base index
	int Tran2Ten(const vector<uchar>& v, uint nBeg); 
	// char -> 10-base index
	int Tran2Ten(CAlnPckg& QrAln); 

	// search all sequences in database
	void Search(string& sDbPre, int& nDbBlockNum, string& sQPre, int& nQBlockNum);
	// search for each search in database
	void Searching(int nQrIdx, CQrPckg& Query, CDbPckg& Db); 
	void ExtendSetPair(int nLen, CQrPckg& Query, CDbPckg& Db); 
	// search for each seed in a entry of database
	int  ExtendSeq2Set(int nSeed, uint unSeedLen, vector<uchar>& vExtra,
			int nQSeqIdx, CAlnPckg& QrAln, int nQOriLen,
			VUINT& vDSet, CDbPckg& Db,
			VNAMES& vQNames, VNAMES& vDNmaes,
			MRESULT& mRes, int nTreadID);
	void FindSeeds(int nSeed, uint unIncr, CAlnPckg& QrAln, MINDEX& vDHash, VUINT& vRes, bool bMutation);
	// set intersection
	void SetInter(vector<uint>& v1, vector<uint>& v2, int nDiff, vector<uint>& v);
	// set difference
	void SetDiff(vector<uint>& v1, vector<uint>& v2, int nDiff, vector<uint>& v);

	// align two sequences
	bool AlignSeqs(int nSeed, CAlnPckg& QrAln, CAlnPckg& DbAln,  uint& unSeedLen, STAlnmnt& stAlnmnt, int nTreadID);
	// forward ungapped alignment
	int AlignFwd(uchar *queryseq, uchar *dataseq, uint len_queryseq, uint len_dataseq, int *extl, int *match, int score0);
	// backward ungapped alignment
	int AlignBwd(uchar *queryseq, uchar *dataseq, int pos1, int pos2, int *extl, int *match, int score0);
	// gapped alignment
	int AlignGapped(uchar *seq1, uchar *seq2, int M, int N, int *ext1, int *ext2, int *match_len, int *gap, vector<char>& vMode, vector<short>& vLen, int nTreadID);
	// retrieve a sequence according to a seed
	uchar* GetSeq(VUCHAR& vSeqs, VUINT& vLen, VNAMES& vNames, uint& unPos, uint& unLen, uint& unSeedBeg);
	
	// reset the struct
	void ResetResult(STAlnmnt& stAlnmnt);

	// calculate e-value and generate the subsequence 
	void CalRes(int nQIdx, uchar* pQ, int nQLen, uint unQSeedBeg/*, string& sQName*/, int nDIdx, uchar* pD, uint unDSeedBeg, VNAMES& vDNames, uint unLocalSeedLen, STAlnmnt& stAlnmnt, ostream& out, MRESULT& mRes, int nTreadID);
	// calculate e-values of multi hits
	void SumEvalue(vector<STResult>& v, int nSt, int nEd, int nLen, int nTreadID);
	// output the result
	void PrintRes(MRESULT& mRes, int nTreadID, CQrPckg& Query, CDbPckg& Db);

	// test 
	void PrintHash(MINDEX& v); 
	// test 
	void PrintInfo(MINDEX& v); 
	template <class type>
	void PrintNum(string sPre, type val);
	bool CompPrev(uchar* p, uint unBeg, uint un);

	// revise the size of database according to swift
	void GuessTotSeq(const char* szFile, long int& lnSeqNum, long int& lnAaNum);

	// merge the result files
	void MergeRes(int nDbBlockNum, string& sQPre);

	// init alignment parameters
	void InitAlignPara(bool bType, long int lnSLen, int nSNum, int nThreadNum);



private:
	uint m_unMer;
	uint m_unDSize;
	uint m_unQSize;
	bool m_bSeqType;
	// if dna, 6; if protein, 1
	int m_nIdxScl;
	uint m_unTotalIdx;

	uchar m_uMask;
	uchar m_aChar2Code[256]; // char -> compressed code
	uchar m_aCode2Char[256]; // compressed code -> char
	uchar m_aCode2Ten[256]; // char -> 10-base

	int m_aSubMatrix[256][256];

	MRESULT m_mRes;

	/* for alignment and calculation */
	double m_dLogEvaThr;
	int GapIni;
	int GapExt;
	int MaxGap;

	BlastStat* m_pBlastSig;
	double	GapExtSCutBits;
	double	GapExtSCut;
	double	UngapExtDropBits;
	double	UngapExtDrop;
	double	GapExtDropBits;
	double	GapExtDrop;
	double	UngapExtSCut;
	int MinMatch4Exp;

	vector<vector<vector<char> > > m_vTrace;
	vector<vector<vector<char> > > m_vETrace;
	vector<vector<vector<char> > > m_vDTrace;

	bool m_bFast;
	uint m_unMutSeedLen;
	VUINT m_vMutation;

	uint m_unTotalSeeds;
	uint m_unTotalQuery;
	uint m_unTotalSubj;

	/* for output */
	long long m_nMaxOut;
	long long m_nMaxM8;
	bool m_bPrintEmpty;
	bool m_bGapExt;
	ofstream m_ofAln;
	ofstream m_ofM8;
	string m_sOutBase;
	// store the ouput
	string m_sOutput;
	string m_sM8;
	vector<CIndex> m_vOutIdx;
	vector<CIndex> m_vM8Idx;
	long long m_llOutCum;
	long long m_llM8Cum;
	int m_nSeqBase;

	// for multithread
	int m_nThreadNum;
	vector<BlastStat*> m_vpBlastSig;
	vector<int> m_vBlastPt; // -1 means available, 1 means used

	// for test on gap extension
	uint m_unGapExt;
	bool m_bAcc;
	bool m_bHssp;
	int m_nMinLen;

	// hssp criteria
	vector<int> m_vCriteria;
};


inline void CHashSearch::InitAlignPara(bool bType, long int lnSLen, int nSNum, int nThreadNUm)
{
	/*
	if (NULL != m_pBlastSig)
	{
		delete m_pBlastSig;
	}
	*/
	for (int i = 0; i < nThreadNUm; ++i)
	{
		BlastStat* pBlastSig = NULL;
		/*
		if (bType)
		{
		*/
		pBlastSig = new BlastStat(1, lnSLen, nSNum);
		/*
		}
		else
		{
			pBlastSig = new BlastStat(0, lnSLen, nSNum);
		}
		*/
		m_vpBlastSig.push_back(pBlastSig);
	}

	GapIni = GAPINI;
	GapExt = GAPEXT;
	MaxGap = MAXGAP;
	// set up cutoff values
	// fix it!!!
	GapExtSCutBits = 25; //the dropoff (in bits) to invoke gapped alignment
	GapExtSCut = m_vpBlastSig[0]->Bits2RawScoreUngapped(GapExtSCutBits);

	UngapExtDropBits = 7; //in bits
	GapExtDropBits = 15; //in bits
	if (m_bAcc == false)
	{
		if(bType)
		{
			UngapExtSCut = 12; //blastx default for -f option
		}
		else
		{
			UngapExtSCut = 11; //blastp default
		}
	}
	else
	{
		if(bType)
		{
			UngapExtSCut = 37; //blastx default for -f option
		}
		else
		{
			UngapExtSCut = 11; //blastp default
		}
	}

	UngapExtDrop = m_vpBlastSig[0]->Bits2RawScoreUngapped(UngapExtDropBits);
	GapExtDrop = m_vpBlastSig[0]->Bits2RawScoreGapped(GapExtDropBits);

	MinMatch4Exp = 4;
}


inline int CHashSearch::Tran2Ten(const vector<uchar>& v, uint nBeg)
{
	if (nBeg >= v.size())
	{
		return -1;
	}
	uint un = 0;
	for (uint i = 0; i < m_unMer; ++i)
	{
		if (m_uMask == m_aCode2Ten[v[nBeg+i]])
		{
			return -1;
		}
		un = un*10 + m_aCode2Ten[v[nBeg+i]];
	}
	return un;
}


inline int CHashSearch::Tran2Ten(CAlnPckg& QrAln)
{
	if (QrAln.m_unSeedBeg >= QrAln.m_unLen-m_unMer+1)
	{
		return -1;
	}
	uint un = 0;
	for (uint i = 0; i < m_unMer; ++i)
	{
		if (m_uMask == m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg+i]])
		{
			return -1;
		}
		un = un*10 + m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg+i]];
	}
	return un;
}


template <class type>
inline void CHashSearch::PrintNum(string sPre, type val)
{
	cout << sPre << ":\t" << val << endl;
}


inline void CHashSearch::Encode(const string& sIn, vector<uchar>& v)
{
	v.reserve(sIn.size());

	for (int i = 0; i < sIn.size(); ++i)
	{
		v.insert(v.end(), m_aChar2Code[sIn[i]]);
	}
}


inline void CHashSearch::Decode(const vector<uchar>& v, string& sOut)
{
	sOut.reserve(v.size());
	for (int i = 0; i < v.size(); ++i)
	{
		sOut += m_aCode2Char[v[i]];
		//if ('-' == v[i])
		//{
			//sOut += '-';
		//}
		//else
		//{
			//sOut += m_aCode2Char[v[i]];
		//}
	}
}


inline void CHashSearch::Encode(vector<uchar>& v)
{
	//cout << v.size() << endl;
	for (int i = 0; i < v.size(); ++i)
	{
		v[i] = m_aChar2Code[v[i]];
	}
}


inline long int CHashSearch::Encode(vector<uchar>& v, vector<double>& vFreq)
{
	long int lnTotalAa = 0;
	//cout << v.size() << endl;
	for (int i = 0; i < v.size(); ++i)
	{
		v[i] = m_aChar2Code[v[i]];
		if (m_uMask != v[i])
		{
			++vFreq[m_aCode2Ten[v[i]]];
			++lnTotalAa;
		}
	}

	return lnTotalAa;
}


inline void CHashSearch::Decode(vector<uchar>& v)
{
	//cout << v.size() << endl;
	for (int i = 0; i < v.size(); ++i)
	{
		v[i] = m_aCode2Char[v[i]];
	}
}


inline uchar* CHashSearch::GetSeq(VUCHAR& vSeqs, VUINT& vLens, VNAMES& vNames, uint& unPos, uint& unLen, uint& unSeedBeg)
{
	uint unIdx = unPos >> 11;
	unSeedBeg = unPos & 0x000007FF;
	unLen = vLens[unIdx+1] - vLens[unIdx];
	return &vSeqs[0] + vLens[unIdx];
}


inline void CHashSearch::ResetResult(STAlnmnt& stAlnmnt)
{
	stAlnmnt.nScore = 0;
	stAlnmnt.nMatch = 0;
	stAlnmnt.nQFwd = 0;
	stAlnmnt.nDFwd = 0;
	stAlnmnt.nQBwd = 0;
	stAlnmnt.nDBwd = 0;
	stAlnmnt.vMode.clear();
	stAlnmnt.vLen.clear();
}


inline void CHashSearch::SetInter(vector<uint>& v1, vector<uint>& v2, int nDiff, vector<uint>& v)
{
	vector<uint>::iterator it1 = v1.begin();
	vector<uint>::iterator it2 = v2.begin();

	while (it1 != v1.end() && it2 != v2.end())
	{
		if (((*it1)>>11)<((*it2)>>11) || ((((*it1)>>11)==((*it2)>>11))&&(((*it1)&0x7ff)<(((*it2)&0x7ff)-nDiff))))
		{
			++it1;
		}
		else if (((*it1)>>11)>((*it2)>>11) || ((((*it1)>>11)==((*it2)>>11))&&(((*it1)&0x7ff)>(((*it2)&0x7ff)-nDiff))))
		{
			++it2;
		}
		else
		{
			v.push_back(*it1);
			++it1;
			++it2;
		}
	}
}


inline void CHashSearch::SetDiff(vector<uint>& v1, vector<uint>& v2, int nDiff, vector<uint>& v)
{
	vector<uint>::iterator it1 = v1.begin();
	vector<uint>::iterator it2 = v2.begin();

	while (it1 != v1.end() && it2 != v2.end())
	{
		if (((*it1)>>11)<((*it2)>>11) || ((((*it1)>>11)==((*it2)>>11))&&(((*it1)&0x7ff)<(((*it2)&0x7ff)-nDiff))))
		{
			v.push_back(*it1);
			++it1;
		}
		else if (((*it1)>>11)>((*it2)>>11) || ((((*it1)>>11)==((*it2)>>11))&&(((*it1)&0x7ff)>(((*it2)&0x7ff)-nDiff))))
		{
			++it2;
		}
		else
		{
			++it1;
			++it2;
		}
	}

	while (it1 != v1.end())
	{
		v.push_back(*it1);
		++it1;
	}
}


inline bool CHashSearch::CompPrev(uchar* p, uint unBeg, uint un)
{
	if (unBeg < un)
	{
		return false;
	}

	for (uint i = unBeg-un; i < unBeg; ++i)
	{
		if (p[i] != p[i+un])
		{
			return false;
		}
	}

	return true;
}


#endif // __HASHSEARCH_H_
