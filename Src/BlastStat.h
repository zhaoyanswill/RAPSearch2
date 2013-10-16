//it is the main class for ExpSearch
//reading in the query & database
//make suffix tree of databases & do search

#ifndef __EVALUECAL_H__
#define __EVALUECAL_H__

#include <ctype.h>
#include <map>
#include <math.h>

// blast statistics parameters for **BLOSOM62** and gapped/ungapped alignment
// refer to David W. Mount's bioinformatics book
// data from blastkar.cpp
//composition-based LAMBDA???
#define DEFAULT_L				0.318
#define DEFAULT_K				0.134
#define DEFAULT_H				0.401
#define DEFAULT_ALPHADLAMBDA			2.492397
#define	DEFAULT_BETA				-3.2
#define DEFAULT_GAPPED_L			0.267
#define DEFAULT_GAPPED_K			0.0410
#define DEFAULT_GAPPED_H			0.140
#define DEFAULT_GAPPED_ALPHADLAMBDA		7.116105
#define	DEFAULT_GAPPED_BETA			-30
#define BLAST_GAP_DECAY_RATE 			0.5   	/**< Gap decay rate for ungapped search */
#define BLAST_GAP_DECAY_RATE_GAPPED 		0.1	/**< Gap decay rate for gapped search */


#define	Log2					log(2)
#define	Log10					log(10)
#define	DEFAULT_G				50
#define	LogG					log(DEFAULT_G)
#define	DEFAULT_GAP_DECAY			0.1
#define MAX(m, n)				(m > n)?(m):(n)
//lamda, k, h (h is the average nats/aligned pair)
//G: gap_size

class BlastStat
{
private:
    int	DBSeqNum;
    double	DBLen;
    double	EDBLen;
    double	LogEDBLen;
    int	QueryLen;
    double	EQueryLen;
    double	LogEQueryLen;
    int	ifGapped;
    int	PreAdjustMax;
    int	*AdjustLen;

    double	L;
    double	K;
    double	LogK;
    double	H; //average nats/aligned pair
    double	alpha_d_lambda;
    double	beta;
    double	ExpectedHSPLength; //the length of an HSP that has an expect of 1
    double	gap_decay_rate;

    void	BlastGapDecayCorr(double *e, int nsegs);

public:
    BlastStat(int ifgapped, long int dbseqres, int dbseqnum);
    ~BlastStat(void);
    void	SetPar(int ifgapped);
    void	SetDBInfo(void);
    double	calEffectiveLen(double len);
    double	rawScore2Bit(double rawscore);
    double 	rawScore2Expect(double rawscore);
    double	rawScore2ExpectLog(double rawscore);
    double	sumScore2Expect(int tot, double* score, int subjctlen);
    double	sumScore(int tot, double *score, int subjctlen);
    double	sumScore2Expect(int tot, double sumS, int subjctlen);
    double	fac(int r);
    void	blastComputeLengthAdjustmentComp(int query_length);
    int	blastComputeLengthAdjustment(int query_length, int *length_adjustment);
    int	blastComputeLengthAdjustmentSimple(int query_length, int db_length, int db_num_seqs);
    void	setEffectiveLen(int length_adjustment, int query_length);
    double	Bits2RawScoreUngapped(double bit);
    double	Bits2RawScoreGapped(double bit);
};

#endif
