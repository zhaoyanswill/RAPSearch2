#include "HashSearch.h"

void printUsage(char *error);

int main(int argc, char *argv[])
{
	char* szDbSeq = NULL;
	char* szDbHash = NULL;
	int	nSplitNum = 0;

	if(argc < 5) printUsage(""); 

	int	i, k;
	for (i = 0; i < argc; i ++) {
		if(argv[i][0] != '-') continue;
		if(argv[i][1] == 'd' && argc > i + 1) // string : database file (FASTA format)
		{
			szDbSeq = argv[i+1];
		}
		else if(argv[i][1] == 'n' && argc > i + 1) // string : suffix array file (binary file)
		{
			szDbHash = argv[i+1];
		}
		else if(argv[i][1] == 's' && argc > i + 1) // string : suffix array file (binary file)
		{
			sscanf(argv[i+1], "%d", &nSplitNum);
		}
	}

	if (!szDbSeq)	printUsage("Error: No database\n");
	if (!szDbHash)	printUsage("Error: No hash file\n");

	printf("now building hash file\n");
	CHashSearch hashSearch(0);
	hashSearch.Process(szDbSeq, szDbHash, nSplitNum);

	printf("hash file saved to file %s\n", szDbHash);

	printf(">>>Main END\n");

	return 0;
}

void printUsage(char *info)
{
	printf("%s", info);
	printf("Usage: prerapsearch -d database -n swift-file(base name only)\n");
	printf("       optional parameter:\n");
	printf("          -s splits-num (splits the database into specified number of files)\n"); 
	printf("Example: prerapsearch -d nr -n nr\n");
	exit(0);
}
