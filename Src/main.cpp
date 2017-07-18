#include "HashSearch.h"
#include <iostream>
#include "paras.h"
using namespace std;


#define OPTION_QUERY "q"
#define OPTION_SUBJECT "d"
#define OPTION_OUTPUT "o"
#define OPTION_STDOUT "u"
#define OPTION_EVAL "e"
#define OPTION_THREADNUM	"z"
#define OPTION_MAXHIT "v"
#define OPTION_MAXALN "b"
#define OPTION_QUERYTYPE "t"
#define OPTION_PRINTEMPTY "p"
#define OPTION_GAPEXT "g"
#define OPTION_ACCELERATE "a"
#define OPTION_HSSP "w"
#define OPTION_MINLEN "l"
#define OPTION_XML "x"
#define OPTION_BITS "i"
#define OPTION_LOGE "s"
#define	Program	"rapsearch"
#define	Version "2.24"

void printUsage(char *error);

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        printUsage("");
        return 1;
    }

    char* szDbHash = NULL;
    char* szQrFile = NULL;
    char* szOutFile = NULL;
	int nThreadNum = 0;
	double dLogEvalueCut = DEFAULT_LOGEVALTHRESH;
	int	nMaxHitPer = 500;
	int	nMaxAlnPer = 100;
	int nQueryType = 0;	//0:unknown, 1:nt, 2:aa
	bool bPrintEmpty = false;
	bool bGapExt = true;
	bool bAcc = false;
	bool bHssp = false;
	int nMinLen = 0;
	int bXml = false;
	int nStdout = 0;
	double dBitsCut = 0.0;
	bool bBitsCut = false;
	bool bLogE = true;

    int	i;
    for (i = 0; i < argc; i ++)
    {
        if(argv[i][0] != '-') continue;
        if(argv[i][1] == OPTION_QUERY[0] && argc > i + 1)  // string : query file (FASTA format)
        {
            szQrFile = argv[i+1];
        }
        else if(argv[i][1] == OPTION_SUBJECT[0] && argc > i + 1) // string : subject file (base name only)
        {
            szDbHash = argv[i+1];
        }
        else if(argv[i][1] == OPTION_OUTPUT[0] && argc > i + 1) // string : output file name
        {
            szOutFile = argv[i+1];
        }
		else if(argv[i][1] == OPTION_THREADNUM[0] && argc > i + 1) // string : suffix array file (binary file)
		{
			sscanf(argv[i+1], "%d", &nThreadNum);
		}
		else if(argv[i][1] == OPTION_EVAL[0] && argc > i + 1) // evalue 
		{
			sscanf(argv[i+1], "%lf", &dLogEvalueCut);
		}
		else if(argv[i][1] == OPTION_MAXHIT[0] && argc > i + 1) // maxhit 
		{
			sscanf(argv[i+1], "%d", &nMaxHitPer);
		}
		else if(argv[i][1] == OPTION_MAXALN[0] && argc > i + 1) // maxaln 
		{
			sscanf(argv[i+1], "%d", &nMaxAlnPer);
		}
		else if(argv[i][1] == OPTION_QUERYTYPE[0] && argc > i + 1) // query type 
		{
			char cType = 'n';
			sscanf(argv[i+1], "%c", &cType);
			if ('n' == cType || 'N' == cType)
			{
				nQueryType = 1;
			}
			else if ('a' == cType || 'A' == cType)
			{
				nQueryType = 2;
			}
			else if ('q' == cType || 'Q' == cType)
			{
				nQueryType = 3;
			}
			else
			{
				nQueryType = 0;
			}
		}
		else if(argv[i][1] == OPTION_PRINTEMPTY[0]) // print query without hits 
		{
			char cPrint = 'f';
			sscanf(argv[i+1], "%c", &cPrint);
			if ('t' == cPrint || 'T' == cPrint)
			{
				bPrintEmpty = true;
			}
			else
			{
				bPrintEmpty = false;
			}
		}
		else if(argv[i][1] == OPTION_GAPEXT[0]) // gap extension 
		{
			char cGapExt = 't';
			sscanf(argv[i+1], "%c", &cGapExt);
			if ('f' == cGapExt || 'F' == cGapExt)
			{
				bGapExt = false;
			}
			else
			{
				bGapExt = true;
			}
		}
		else if(argv[i][1] == OPTION_ACCELERATE[0]) // accelerate mode 
		{
			char cAcc = 'f';
			sscanf(argv[i+1], "%c", &cAcc);
			if ('f' == cAcc || 'F' == cAcc)
			{
				bAcc = false;
			}
			else
			{
				bAcc = true;
			}
		}
		else if(argv[i][1] == OPTION_HSSP[0]) // HSSP criteria 
		{
			char cHssp = 'f';
			sscanf(argv[i+1], "%c", &cHssp);
			if ('f' == cHssp || 'F' == cHssp)
			{
				bHssp = false;
			}
			else
			{
				bHssp = true;
			}
		}
		else if(argv[i][1] == OPTION_MINLEN[0] && argc > i + 1) // min alignment length 
		{
			sscanf(argv[i+1], "%d", &nMinLen);
		}
		else if(argv[i][1] == OPTION_XML[0]) // HSSP criteria 
		{
			char cXml = 'f';
			sscanf(argv[i+1], "%c", &cXml);
			if ('f' == cXml || 'F' == cXml)
			{
				bXml = false;
			}
			else
			{
				bXml = true;
			}
		}
        else if(argv[i][1] == OPTION_STDOUT[0] && argc > i + 1) // string : which result outputed into stdout
        {
			sscanf(argv[i+1], "%d", &nStdout);
        }
		else if(argv[i][1] == OPTION_BITS[0] && argc > i + 1) // bits 
		{
			sscanf(argv[i+1], "%lf", &dBitsCut);
			bBitsCut = true;
		}
		else if(argv[i][1] == OPTION_LOGE[0]) // full id
		{
			char cLogE = 't';
			sscanf(argv[i+1], "%c", &cLogE);
			if ('f' == cLogE || 'F' == cLogE)
			{
				bLogE = false;
			}
			else
			{
				bLogE = true;
			}
		}
    }

    if (!szQrFile)	printUsage("Error: No query");
    if (!szDbHash)	printUsage("Error: No database");
    if (!szOutFile && nStdout == 0)	printUsage("Error: No output");

	bool bEvalue = true;
	double dThr = dLogEvalueCut;
	if (bLogE == false && dThr == DEFAULT_LOGEVALTHRESH)
	{
		dThr = 10.0;
	}
	if (bBitsCut == true)
	{
		bEvalue = false;
		dThr = dBitsCut;
	}

    time_t jobstart = time(NULL);

    printf("QueryFileName %s\n", szQrFile);

	CHashSearch hs(nThreadNum);
	hs.Process(szDbHash, szQrFile, szOutFile, nStdout, bEvalue, bLogE, dThr, nMaxAlnPer, nMaxHitPer, nQueryType, bPrintEmpty, bGapExt, bAcc, bHssp, nMinLen, bXml);

    printf(">>>Main END\n");
    time_t jobfinished = time(NULL);
    double	timeused = difftime(jobfinished, jobstart);
    printf("Time used: %d min (%.1f hours)\n", int(timeused / 60), timeused / 3600.0);

    return 0;
}


void printUsage(char *error)
{
    fprintf(stderr, "%s\n", error);
    fprintf(stderr,
            "%s v%s: Fast protein similarity search tool for short reads\n"
            "-------------------------------------------------------------------------\n"
            " Options: \n"
            "\t-" OPTION_QUERY        " string : query file or stdin (FASTA or FASTQ format)\n"
            "\t-" OPTION_SUBJECT      " string : database file (**base name only**, with full path)\n"
            "\t-" OPTION_OUTPUT       " string : output file name\n"
            "\t-" OPTION_STDOUT       " int    : stream one result through stdout [1: m8 result, 2: aln result, default: don't stream any result through stdout]\n"
            "\t-" OPTION_THREADNUM    " int    : number of threads [default: %d]\n"
	    "\t-" OPTION_LOGE       " char   : report log10(E-value) or E-value for each hit [t/T: log10(E-value), the default; f/F: E-value]\n"
            "\t-" OPTION_EVAL         " double : E-value threshold, given in the format of log10(E-value), or E-value (when -s is set to f) [default: %.1f/10.0]. \n"
            "\t-" OPTION_BITS         " double : threshold of bit score [default: %.1f]. It is the alternative option to limit the hits to report.\n"
            "\t-" OPTION_MINLEN         " int    : threshold of minimal alignment length [default: %d]\n"
	    "\t-" OPTION_MAXHIT       " int    : number of database sequences to show one-line descriptions [default: %d]. If it's -1, all results will be shown.\n"
	    "\t-" OPTION_MAXALN       " int    : number of database sequence to show alignments [default: %d]. If it's -1, all results will be shown.\n"
	    "\t-" OPTION_QUERYTYPE       " char   : type of query sequences [u/U:unknown, n/N:nucleotide, a/A:amino acid, q/Q:fastq, default: %s]\n"
	    "\t-" OPTION_PRINTEMPTY       " char   : output ALL/MATCHED query reads into the alignment file [t/T: output all query reads, f/F: output matched reads, default: %s]\n"
	    "\t-" OPTION_GAPEXT       " char   : apply gap extension [t/T: yes, f/F: no, default: %s]\n"
	    "\t-" OPTION_ACCELERATE       " char   : use fast mode (10~30 fold) [t/T: yes, f/F: no, default: %s]\n"
	    "\t-" OPTION_HSSP       " char   : apply HSSP criterion instead of E-value criterion [t/T: HSSP, f/F: E-value criteria, default: %s]\n"
     	    "\t-" OPTION_XML       " char   : print hits in xml format [t/T: yes, f/F: no, default: %s]\n"
            "-------------------------------------------------------------------------\n"
            "example 1> %s -" OPTION_QUERY " query.fa -" OPTION_SUBJECT " nr -" OPTION_OUTPUT " output_file\n"
            "example 2> %s -" OPTION_QUERY " query.fa -" OPTION_SUBJECT " nr -" OPTION_OUTPUT " output_file -" OPTION_BITS " 40 -" OPTION_MINLEN " 25\n"
	    "   this setting only reports the hits with bit score >= 40 and alignment length >= 25\n"
            "example 3> %s -" OPTION_QUERY " query.fa -" OPTION_SUBJECT " nr -" OPTION_OUTPUT " output_file -" OPTION_EVAL " -5\n"
	    "   this setting only reports hits with log(E-value) <= -5 (i.e., E-value <= 1e-5)\n"
            "example 4> %s -" OPTION_QUERY " query.fa -" OPTION_SUBJECT " nr -" OPTION_OUTPUT " output_file -" OPTION_EVAL " 1e-5 -" OPTION_LOGE " f\n"
	    "   this setting only reports the hits with E-value <= 1e-5\n"
	    "the difference between example 3 & 4 is that the former lists log(E-value) while the latter lists E-value for each hit in the output file\n"
            ,
            Program, Version, 1, 1.0, 0.0, 0, 500, 100, "u", "f", "t", "f", "f", "f", Program, Program, Program, Program
           );
    exit(-1);
}
