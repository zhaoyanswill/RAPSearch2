#ifndef __SORTKEY_H__
#define __SORTKEY_H__

#include <string>
using namespace std;



class CSortKey
{
public:
	CSortKey(string s, double d)
		:m_sQr(s), m_dEValue(d) {};

	friend bool operator< (const CSortKey& c1, const CSortKey& c2) const
	{
		if (c1.m_sQr != c2.m_sQr)
		{
			return (c1.m_sQr<c2.m_sQr);
		}
		else
		{
			return (c1.m_dEValue<c2.m_dEValue);
		}
	}

private:
	string m_sQr;
	double m_dEValue;

}

#endif
