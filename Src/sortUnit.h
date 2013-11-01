#ifndef __SORTUNIT_H__
#define __SORTUNIT_H__

#include <string>



class CSortUnit
{
public:
	CSortUnit(double d, std::vector<char>::iterator it1, std::vector<char>::iterator it2)
		: m_dEValue(d), m_sHit(it1,it2) {};

	friend bool operator< (const CSortUnit& c1, const CSortUnit& c2)
	{
		return (c1.m_dEValue<c2.m_dEValue);
	}

public:
	double m_dEValue;
	std::string m_sHit;

};

#endif
