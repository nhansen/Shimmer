/**************************************************************************
#
#  This software/database is "United States Government Work" under the terms of
#  the United States Copyright Act.  It was written as part of the authors'
#  official duties for the United States Government and thus cannot be
#  copyrighted.  This software/database is freely available to the public for
#  use without a copyright notice.  Restrictions cannot be placed on its present
#  or future use.
#
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and data, the National Human Genome Research
#  Institute (NHGRI) and the U.S. Government does not and cannot warrant the
#  performance or results that may be obtained by using this software or data.
#  NHGRI and the U.S.  Government disclaims all warranties as to performance,
#  merchantability or fitness for any particular purpose.
# 
#  In any work or product derived from this material, proper attribution of the
#  authors as the source of the software or data should be made, using "NHGRI
#  Genome Technology Branch" as the citation.
#
**************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "printcomp.h"

extern Params *parameters;

main(int argc, char **argv)
{
	int i,j,k;
	int baseLength;
	FILE *fil;
	long chunkSize;
	long position;
	long lastPosition;
	char *charPtr;
	char *chrPtr;
	char cmd[1000000];
	char tmpcmd[1000000];
	char mapQual, baseQual;
	int baseQualNo;
        char nextChar;
        char refBase;
        int fieldNo;
	int noBases[10];
	int indelLength;
	int posLength;
	int refID; /* store current reference base as an integer */
	int sampleShift;
	char posBases[3000000];
	char posString[3000000];
	char chromosome[1000];
	char* baseArray = "ATGCN";
	int maxBases;
	int maxBase;
	int nextMaxBases;
	int nextMaxBase;
	int baseCount;

        get_params(argc, argv);

        if (!parameters->bam1 || !parameters->bam2 || !parameters->fasta) {
            fprintf(stderr, "Usage: printCompCounts -fasta <reference fasta> -bam1 <normal bam file> -bam2 <tumor bam file>\n");
            exit(1);
        }

	/* specify reference fasta file, then chromosome, then normal bam file, then mutated bam file */
        chunkSize = 4096;

	/* base command which is always run */
	sprintf(cmd, "samtools mpileup -f %s -BQ0 -d100000 -A %s %s", parameters->fasta, parameters->bam1, parameters->bam2);

        if (strcmp(parameters->region, "") != 0) {
		sprintf(tmpcmd, "%s -r %s", &cmd, parameters->region);
		sprintf(cmd, "%s", &tmpcmd);
	}

	if (parameters->mapqual != 0) {
		sprintf(tmpcmd, "%s -q %i", &cmd, parameters->mapqual);
		sprintf(cmd, "%s", &tmpcmd);
	}

        /* check to be sure we can open mpileup */
	if ((fil = popen(cmd,"r")) == NULL) {
		printf("error printCompCounts: unable to run mpileup: %s\n",cmd);
		exit(1);
	}

	fieldNo = 0;
	noBases[0] = 0;
	noBases[1] = 0;
	noBases[2] = 0;
	noBases[3] = 0;
	noBases[4] = 0;
	noBases[5] = 0;
	noBases[6] = 0;
	noBases[7] = 0;
	noBases[8] = 0;
	noBases[9] = 0;
	baseLength = 0;
	while((nextChar = fgetc( fil )) != EOF) {
		if (nextChar == '\t') {
			fieldNo++;
		}
		if (fieldNo == 0) { /* populate chromosome with the current chromosome */
			chrPtr = chromosome;
			*chrPtr = nextChar;
			chrPtr++;
			while ((nextChar = fgetc( fil )) != '\t') {
				*chrPtr = nextChar;
				chrPtr++;
			}
			*chrPtr = '\0';
			fieldNo++;
		}
		if (fieldNo == 1) { /* populate position with current position */
			posLength = 0;
			while ((nextChar = fgetc( fil )) != '\t') {
				posString[posLength++] = nextChar;
			}
			posString[posLength] = '\0';
			fieldNo = 2;
			position = atoi(posString);
		}
		if (fieldNo == 2) { /* record reference base */
			nextChar = fgetc ( fil );
			refBase = nextChar;
			switch( refBase ) {
				case 'A' :
				case 'a' :
					refID = 0;
					break;
				case 'T' :
				case 't' :
					refID = 1;
					break;
				case 'G' :
				case 'g' :
					refID = 2;
					break;
				case 'C' :
				case 'c' :
					refID = 3;
					break;
			}
			while ((nextChar = fgetc( fil )) != '\t') {
			}
			fieldNo = 3;
		}
		if ((fieldNo == 4) || (fieldNo == 7)) { /* record read bases */
			sampleShift = 0;
			if (fieldNo == 7) {
				sampleShift += 5; /* will store counts for second sample in indices 5-9 */
			}

			baseLength = 0;
			while ((nextChar = fgetc( fil )) != '\t') {
				posBases[baseLength++] = nextChar;
			}
			posBases[baseLength] = '\0';
			fieldNo++;
			/* now read quality scores as needed */
			for (i=0; i<baseLength; i++) {
				if (posBases[i] == '^') { /* new read */
					mapQual = posBases[++i];
					/* printf("New read qual %c\n", mapQual); */
					/* might store stacked quals here */
				}
				else if (posBases[i] == '$') { /* end of read */
					/* shift the stacked quals */
				}
				else if ((posBases[i] == '+') ||
						(posBases[i] == '-')) { /* indel */
					indelLength = atoi(&posBases[++i]);
					i += indelLength;
					/* for (j=0; j<indelLength; j++) {
						putchar(posBases[++i]);
					} */
					while (indelLength/10 > 0) {
						indelLength /= 10;
						i++;
					}
					/* printf("\n"); */
				}
				else if (posBases[i] == '*') { /* placeholder */
				}
				else { /* a base */
					baseQual = fgetc( fil );
					baseQualNo = baseQual - '!';
					/* printf("New base qual %c score %i\n", baseQual, baseQualNo); */

                                        /* can calculate map quality of read from stacked quals here (stored cyclically?) */

					if((parameters->minqual == 0) || (baseQualNo >= parameters->minqual)) {
						switch( posBases[i] ) {
							case 'A' :
							case 'a' :
								noBases[0+sampleShift]++;
								break;
							case 'T' :
							case 't' :
								noBases[1+sampleShift]++;
								break;
							case 'G' :
							case 'g' :
								noBases[2+sampleShift]++;
								break;
							case 'C' :
							case 'c' :
								noBases[3+sampleShift]++;
								break;
							case '.' :
							case ',' :
								noBases[refID+sampleShift]++;
								break;
						}
					}
				}
			}
		}
		if (nextChar == '\n') { /* end of line */
			/* determine most frequent base, then next most frequent base */
			maxBases = 0;
			maxBase = 4;
			nextMaxBases = 0;
			nextMaxBase = 4;
			for (j=0; j<4; j++) {
				baseCount = noBases[j] + noBases[j+sampleShift]; /* how many of this base? */
				if (baseCount > maxBases) { /* move things down the ladder */
					nextMaxBase = maxBase;
					nextMaxBases = maxBases;
					maxBase = j;
					maxBases = baseCount;
				}
				else if (baseCount > nextMaxBases) {
					nextMaxBase = j;
					nextMaxBases = baseCount;
				}
			}
			/* if ((nextMaxBases > 1) || (position < lastPosition) ||  */
                                   /* (position - lastPosition > 1000)) {  */
                                   /* want to print for testing */
				printf("%s\t%ld\t%c\t%c\t%i\t%i\t%c\t%i\t%i\n", chromosome, position, refBase, baseArray[maxBase], noBases[maxBase], 
					noBases[maxBase+sampleShift], baseArray[nextMaxBase], noBases[nextMaxBase], noBases[nextMaxBase+sampleShift]);
			/* } */
			fieldNo=0;
			noBases[0] = 0;
			noBases[1] = 0;
			noBases[2] = 0;
			noBases[3] = 0;
			noBases[4] = 0;
			noBases[5] = 0;
			noBases[6] = 0;
			noBases[7] = 0;
			noBases[8] = 0;
			noBases[9] = 0;
 			lastPosition = position;
		}
	}
        /* free_params(); */
}
