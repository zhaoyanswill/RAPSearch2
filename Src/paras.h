#ifndef __PARAS_H__
#define __PARAS_H__

// define all of parameters

const int ALPHLEN = 32;
const char aAAAlph[32] = "ARNDCQEGHILKMFPSTWYV.ZXxs";
const int NEGSCORE = -5;

const int GAPINI = 11;
const int GAPEXT = 1;
const int MAXGAP = 3;
const int MINSCORE = -20;

const int SUMHSP_OVERLAP = 10;
const int SUMHSP_MINEVALUE = 1;
const int SUMHSP_MINRAWSCORE = 30;
const double DEFAULT_LOGEVALTHRESH = 1.0;
const int DEFAULT_SCORE_TOP = 5;
const double LOG10 = 2.3025850929940459;

const int DEFSEED = 8;
const int LGERSEED = 12;


#endif
