#include "BlastStat.h"

//---------------------------------------------------
BlastStat::BlastStat(int ifgapped, long int dbseqres, int dbseqnum)
{
    SetPar(ifgapped);
    DBLen = double(dbseqres);
    DBSeqNum = dbseqnum;
    SetDBInfo();
}

void BlastStat::SetPar(int ifgapped)
{
    ifGapped = ifgapped;
    if(ifgapped)
    {
        L = DEFAULT_GAPPED_L;
        K = DEFAULT_GAPPED_K;
        H = DEFAULT_GAPPED_H;
        alpha_d_lambda = DEFAULT_GAPPED_ALPHADLAMBDA;
        beta = DEFAULT_GAPPED_BETA;
        gap_decay_rate = BLAST_GAP_DECAY_RATE_GAPPED;
    }
    else
    {
        L = DEFAULT_L;
        K = DEFAULT_K;
        H = DEFAULT_H;
        alpha_d_lambda = DEFAULT_ALPHADLAMBDA;
        beta = DEFAULT_BETA;
        gap_decay_rate = BLAST_GAP_DECAY_RATE;
    }
    LogK = log(K);
    QueryLen = 0;
    EDBLen = EQueryLen = 0.0;
}

void BlastStat::SetDBInfo()
{
    PreAdjustMax = 1000;
    //precomputed adjustmented length
    AdjustLen = new int[PreAdjustMax];
    int     converge;
    int     i;
    for(i = 0; i < PreAdjustMax; i ++)
    {
        if (i <= 10)    AdjustLen[i] = 0;
        else
        {
            converge = blastComputeLengthAdjustment(i, &AdjustLen[i]);
        }
    }
}

BlastStat::~BlastStat(void)
{
    delete[] AdjustLen;
}

double BlastStat::calEffectiveLen(double len)
{
    double	efflen = len - ExpectedHSPLength;
    if (efflen < 1/ K) return 1 / K;
    else return efflen;
}

//---------------------------------------------------
double BlastStat::rawScore2Bit(double rawscore)
{
    double bit = (L * rawscore - LogK) / Log2;
    //printf("rawscore=%.0f bit=%.2f\n", rawscore, bit);
    return bit;
}

//---------------------------------------------------
double BlastStat::Bits2RawScoreUngapped(double bit)
{
    //rawscore = (bit * Log2 + LogK) / L
    double score = (bit * Log2 + log(DEFAULT_K)) / DEFAULT_L;
    return score;
}

//---------------------------------------------------
double BlastStat::Bits2RawScoreGapped(double bit)
{
    double score = (bit * Log2 + log(DEFAULT_GAPPED_K)) / DEFAULT_GAPPED_L;
    return score;
}

//---------------------------------------------------
double BlastStat::rawScore2Expect(double rawscore)
{
    double 	e = K * EDBLen * EQueryLen * exp(-1 * L * rawscore);

    BlastGapDecayCorr(&e, 1); //to-be-confirmed

    return e;
}

void	BlastStat::BlastGapDecayCorr(double *e, int nsegs)
{
    double	gap_decay_divisor = (1.0 - gap_decay_rate) * pow(gap_decay_rate, nsegs - 1);
    *e /= gap_decay_divisor;
}

//---------------------------------------------------
double BlastStat::rawScore2ExpectLog(double rawscore)
{
    //double loge = (LogK + LogEDBLen + LogEQueryLen -  L * rawscore) / Log10;
    //the above formula does not incoroporate decayrate!

    //double	loge = log(rawScore2Expect(rawscore)) / Log10;
	double d = rawScore2Expect(rawscore);
	if (0 == d)
	{
		return -10000.00;
	}

	double	loge = log(d) / Log10;
	return loge;

}

//---------------------------------------------------
double BlastStat::sumScore2Expect(int tot, double *score, int subjctlen)
{
    double	sumS = sumScore(tot, score, subjctlen);
    double	e =sumScore2Expect(tot, sumS, subjctlen);
    //BlastGapDecayCorr(&e, tot); //to-be-confirmed
    return e;
}

//---------------------------------------------------
double BlastStat::sumScore(int tot, double *score, int subjctlen)
{
    double totalscore = 0;
    for(int i = 0; i < tot; i ++)
    {
        totalscore += score[i];
    }
    double	ESubjctLen = calEffectiveLen(subjctlen);
    double	lgkmn = log(K * EQueryLen * ESubjctLen);
    double 	n_score = L * totalscore - lgkmn - (tot - 1) * (LogK + 2 * LogG) - log(fac(tot));
    //printf("sumscore=%f bits=%f\n", n_score, n_score/Log2);
    return n_score;
}

//---------------------------------------------------
double BlastStat::sumScore2Expect(int tot, double sumS, int subjctlen)
{
    //sumScore to P-value
    double 	sumP = exp(-1 * sumS) * pow(sumS, tot - 1) / (fac(tot) * fac(tot - 1));
    //printf("sumP=%e\n", sumP);

    //correcting for multiple tests
    double 	sumP_corrected = sumP / (pow(DEFAULT_GAP_DECAY, tot - 1) * (1 - DEFAULT_GAP_DECAY));
    //printf("sumP_corrected=%e\n", sumP_corrected);

    //correcting for database size (Effective DB Length, actual query length)
    double	expect = (EDBLen / subjctlen) * sumP_corrected;

    return expect;
}

//---------------------------------------------------
double BlastStat::fac(int r)
{
    int	n = 1;
    for(int i = r; i > 1; i --)
    {
        n *= i;
    }
    return n;
}

//ref: old way for calculating adjustment length
int BlastStat::blastComputeLengthAdjustmentSimple(int query_length, int db_length, int db_num_seqs)
{
    double	len = log(K * query_length * db_length) / H;
    //printf("expected HSP length (simple) = %.3f\n", len);
    //getchar();
}

//from blast_stat.c

/**
 * Computes the adjustment to the lengths of the query and database sequences
 * that is used to compensate for edge effects when computing evalues.
 *
 * The length adjustment is an integer-valued approximation to the fixed
 * point of the function
 *
 *    f(ell) = beta +
 *               (alpha/lambda) * (log K + log((m - ell)*(n - N ell)))
 *
 * where m is the query length n is the length of the database and N is the
 * number of sequences in the database. The values beta, alpha, lambda and
 * K are statistical, Karlin-Altschul parameters.
 *
 * The value of the length adjustment computed by this routine, A,
 * will always be an integer smaller than the fixed point of
 * f(ell). Usually, it will be the largest such integer.  However, the
 * computed length adjustment, A, will also be so small that
 *
 *    K * (m - A) * (n - N * A) > MAX(m,n).
 *
 * Moreover, an iterative method is used to compute A, and under
 * unusual circumstances the iterative method may not converge.
 *
 * @param K      the statistical parameter K
 * @param logK   the natural logarithm of K
 * @param alpha_d_lambda    the ratio of the statistical parameters
 *                          alpha and lambda (for ungapped alignments, the
 *                          value 1/H should be used)
 * @param beta              the statistical parameter beta (for ungapped
 *                          alignments, beta == 0)
 * @param query_length      the length of the query sequence
 * @param db_length         the length of the database
 * @param db_num_seqs       the number of sequences in the database
 * @param length_adjustment the computed value of the length adjustment [out]
 *
 * @return   0 if length_adjustment is known to be the largest integer less
 *           than the fixed point of f(ell); 1 otherwise.
Int4
BLAST_ComputeLengthAdjustment(double K,
                             double logK,
                             double alpha_d_lambda,
                             double beta,
                             Int4 query_length,
                             Int8 db_length,
                             Int4 db_num_seqs,
                             Int4 * length_adjustment)
*/

void BlastStat::blastComputeLengthAdjustmentComp(int query_length)
{
    if(query_length < PreAdjustMax)
    {
        setEffectiveLen(AdjustLen[query_length], query_length);
    }
    else
    {
        int	length_adjustment;
        blastComputeLengthAdjustment(query_length, &length_adjustment);
    }
}

int BlastStat::blastComputeLengthAdjustment(int query_length, int *length_adjustment)
{
    int 	i;                     /* iteration index */
    int	kMaxIterations = 20;     /* maximum allowed iterations */
    double m = (double) query_length;
    //double n = (double) DBLen;
    double n = DBLen;
    double N = (double) DBSeqNum;

    double ell;                 /* A double value of the length adjustment */
    double ss;                  /* effective size of the search space */
    double ell_min = 0, ell_max;   /* At each iteration i,
                                         * ell_min <= ell <= ell_max. */
    int converged    = 0;       /* True if the iteration converged */
    double ell_next = 0;        /* Value the variable ell takes at iteration
                                 * i + 1 */
    double	logK = LogK;

//printf("K=%.3f alpha_d_lambda=%f beta=%f query_length=%d db_length=%d db_num_seq=%d\n", K, alpha_d_lambda, beta, query_length, DBLen, DBSeqNum);

    /* Choose ell_max to be the largest nonnegative value that satisfies
    *
    *    K * (m - ell) * (n - N * ell) > MAX(m,n)
    *
    * Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c )) */
    {
        /* scope of a, mb, and c, the coefficients in the quadratic formula * (the variable mb is -b) */
        double a  = N;
        double mb = m * N + n;
        double	maxmn;
        if(m > n) maxmn = m;
        else maxmn = n;
        //double c  = n * m - MAX(m, n) / K;
        double	c = n * m - maxmn / K;

        if(c < 0)
        {
            *length_adjustment = 0;
            return 1;
        }
        else
        {
            ell_max = 2 * c / (mb + sqrt(mb * mb - 4 * a * c));
        }
    } /* end scope of a, mb and c */

    for(i = 1; i <= kMaxIterations; i++)        /* for all iteration indices */
    {
        double ell_bar;         /* proposed next value of ell */
        ell      = ell_next;
        ss       = (m - ell) * (n - N * ell);
        ell_bar  = alpha_d_lambda * (logK + log(ss)) + beta;
        if(ell_bar >= ell)   /* ell is no bigger than the true fixed point */
        {
            ell_min = ell;
            if(ell_bar - ell_min <= 1.0)
            {
                converged = 1;
                break;
            }
            if(ell_min == ell_max)   /* There are no more points to check */
            {
                break;
            }
        }
        else   /* else ell is greater than the true fixed point */
        {
            ell_max = ell;
        }
        if(ell_min <= ell_bar && ell_bar <= ell_max)
        {
            /* ell_bar is in range. Accept it */
            ell_next = ell_bar;
        }
        else   /* else ell_bar is not in range. Reject it */
        {
            ell_next = (i == 1) ? ell_max : (ell_min + ell_max) / 2;
        }
    } /* end for all iteration indices */
    if(converged)   /* the iteration converged */
    {
        /* If ell_fixed is the (unknown) true fixed point, then we
        * wish to set (*length_adjustment) to floor(ell_fixed).  We
        * assume that floor(ell_min) = floor(ell_fixed) */

        *length_adjustment = (int) ell_min;

        /* But verify that ceil(ell_min) != floor(ell_fixed) */
        ell = ceil(ell_min);
        if( ell <= ell_max )
        {
            ss = (m - ell) * (n - N * ell);
            if(alpha_d_lambda * (logK + log(ss)) + beta >= ell)
            {
                /* ceil(ell_min) == floor(ell_fixed) */
                *length_adjustment = (int) ell;
            }
        }
    }
    else   /* else the iteration did not converge. */
    {
        /* Use the best value seen so far */
        *length_adjustment = (int) ell_min;
    }

//printf("blast_stat: length %d\n", *length_adjustment);

    setEffectiveLen(*length_adjustment, query_length);

    return converged ? 0 : 1;
}

void BlastStat::setEffectiveLen(int length_adjustment, int query_length)
{
    ExpectedHSPLength = length_adjustment;
    QueryLen = query_length;
    EQueryLen = query_length - ExpectedHSPLength;
    LogEQueryLen = log(EQueryLen);
    EDBLen = DBLen - DBSeqNum * ExpectedHSPLength;
    LogEDBLen = log(EDBLen);
    //printf("ExpectedHSPLength = %.0f QueryLen = %d Effective-query-len = %.0f DBLen = %d EDBLen = %.0f\n",
    //	ExpectedHSPLength, QueryLen, EQueryLen, DBLen, EDBLen);
}
