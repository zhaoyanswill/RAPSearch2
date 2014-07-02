#ifndef __HITUNIT_H__
#define __HITUNIT_H__

#include <string>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>



class CHitUnit
{
public:
	friend bool operator< (const CHitUnit& c1, const CHitUnit& c2)
	{
		return (c1.dEValue<c2.dEValue);
	}

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned int version)
	{
		ar & nQrLen;
		ar & nDbIdx;
		ar & nDbLen;
		ar & nScore;
		ar & dBits;
		ar & dEValue;
		ar & dIdent;
		ar & nAlnLen;
		ar & nMismatch;
		ar & nGapOpen;
		ar & nFrame;
		ar & nQSt;
		ar & nQEd;
		ar & nQBeg;
		ar & nQEnd;
		ar & nDSt;
		ar & nDEd;
		ar & sQName;
		ar & sDName;
		ar & sQ;
		ar & sInfo;
		ar & sD;
	}

public:
	int nQrLen;
	int nDbIdx;
	int nDbLen;
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
	std::string sQName;
	std::string sDName;
	std::string sQ;
	std::string sInfo;
	std::string sD;

};


struct HitComptor
{
	virtual bool operator() (const CHitUnit& st1, const CHitUnit& st2) const = 0;
};

struct CompFrame : HitComptor
{
	virtual bool operator() (const CHitUnit& st1, const CHitUnit& st2) const
	{
		return st1.nFrame < st2.nFrame;
	}
};


struct CompQSt : HitComptor
{
	virtual bool operator() (const CHitUnit& st1, const CHitUnit& st2) const
	{
		return st1.nQSt < st2.nQSt;
	}
};


struct CompEval : HitComptor
{
	virtual bool operator() (const CHitUnit& st1, const CHitUnit& st2) const
	{
		return st1.dEValue < st2.dEValue;
	}
};


struct CompBits : HitComptor
{
	virtual bool operator() (const CHitUnit& st1, const CHitUnit& st2) const
	{
		return st1.dBits > st2.dBits;
	}
};


struct ComptorWrapper
{
	HitComptor* _p;
	ComptorWrapper(HitComptor* p) : _p(p) {}

	bool operator() (const CHitUnit& st1, const CHitUnit& st2) const
	{
		return (*_p)(st1, st2);
	}
};



#endif
