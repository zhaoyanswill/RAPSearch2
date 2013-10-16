#ifndef __CINDEX_H__
#define __CINDEX_H__

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>



class CIndex
{
public:
	CIndex() 
	{
		//m_nID = -1;
		m_llBeg = 0;
		m_nSize = 0;
	}

	/*
	CIndex(int nQrIdx, long long llBeg, int nSize)
		: m_nID(nQrIdx), m_llBeg(llBeg), m_nSize(nSize) {}

	friend bool operator< (const CIndex& c1, const CIndex& c2)
	{
		return (c1.m_nID < c2.m_nID);
	}
	*/

private:
	friend class boost::serialization::access;
	template<typename Archive>
	void serialize(Archive& ar, const unsigned int version)
	{
		//ar & m_nID;
		ar & m_llBeg;
		ar & m_nSize;
	}

public:
	//int m_nID;
	long long m_llBeg;
	int m_nSize;
};

#endif
