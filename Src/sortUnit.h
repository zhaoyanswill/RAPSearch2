#ifndef __SORTUNIT_H__
#define __SORTUNIT_H__

#include <string>
using namespace std;



class CSortUnit
{
public:
	CSortUnit(string s, double d, string sHit)
		:m_sQr(s), m_dEValue(d), m_sHit(sHit) {};

	friend bool operator< (const CSortUnit& c1, const CSortUnit& c2)
	{
		if (c1.m_sQr != c2.m_sQr)
		{
			return (c1.m_sQr.compare(c2.m_sQr) < 0);
		}
		else
		{
			return (c1.m_dEValue<c2.m_dEValue);
		}
	}

public:
	string m_sQr;
	double m_dEValue;
	string m_sHit;

};

#endif
