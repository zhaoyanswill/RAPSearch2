#include "Seg.h"
#include "lnfac.h"

unsigned char   aaflag[128];

int state_cmp(const void *s1, const void *s2)
{
    return ( *(int*)s2 - *(int*)s1 );
    //return *s2 - *s1;
}

//---------------------------------
Seg::Seg(void)
{
    initialize(12);
}
Seg::Seg(int window_inp)
{
    initialize(window_inp);
}

void Seg::initialize(int window_inp)
{
    window = window_inp;
    locut = 2.2;
    hicut = 2.5;

    downset = 0;
    upset = 1;

    hilenmin = 0;
    overlaps = FALSE;
    hionly = FALSE;
    loonly = FALSE;
    entinfo = TRUE;
    singleseq = FALSE;
    prettyseq = FALSE;
    prettytree = TRUE;
    charline = 60;
    maxtrim = 100;

    getwin_init();

    entropy_init(window);

}

Seg::~Seg(void)
{
    free(entray);
}

void Seg::getwin_init(void)
{
    char    *cp, *cp0;
    int             i;
    char    c;

    for (i = 0; i < sizeof(aaindex)/sizeof(aaindex[0]); ++i)
    {
        aaindex[i] = 20;
        aaflag[i] = TRUE;
    }

    for (cp = cp0 = "ACDEFGHIKLMNPQRSTVWY"; (c = *cp) != '\0'; ++cp)
    {
        i = cp - cp0;
        aaindex[c] = i;
        aaindex[tolower(c)] = i;
        aachar[i] = tolower(c);
        aaflag[c] = FALSE;
        aaflag[tolower(c)] = FALSE;
    }
    return;
}

//---------------------------------
char* Seg::maskseq(const char *fastaseq)
{
    Segment *segs;
    Sequen *seq = openseq(fastaseq);

    segs = (Segment *) NULL;
    segseq(seq, &segs, 0);
    mergesegs(seq, segs);
    char *masked = singreport(seq, segs);

    freesegs(segs);
    closeseq(seq);

    return masked;
}

/*---------------------------------------------------------------(segment)---*/

void Seg::segseq(Sequen *seq, Segment **segs, int offset)
{
    Segment *seg, *leftsegs;
    Sequen *leftseq;
    int first, last, lowlim;
    int loi, hii, i;
    int leftend, rightend, lend, rend;
    double *H;
    H = seqent(seq);
    if (H==NULL) return;

    first = downset;
    last = seq->length - upset;
    lowlim = first;

    for (i=first; i<=last; i++)
    {
        if (H[i]<=locut && H[i]!=-1)
        {
            loi = findlo(i, lowlim, H);
            hii = findhi(i, last, H);

            leftend = loi - downset;
            rightend = hii + upset - 1;
            //rightend = hii + upset;

            trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend);

            if (i+upset-1<leftend)     /* check for trigger window in left trim */
            {
                lend = loi - downset;
                rend = leftend - 1;

                leftseq = openwin(seq, lend, rend-lend+1);
                leftsegs = (Segment *) NULL;
                segseq(leftseq, &leftsegs, offset+lend);
                if (leftsegs!=NULL)
                {
                    if (*segs==NULL) *segs = leftsegs;
                    else appendseg(*segs, leftsegs);
                }
                closewin(leftseq);

            }

            seg = (Segment *) malloc(sizeof(Segment));
            seg->begin = leftend + offset;
            seg->end = rightend + offset;
            seg->next = (Segment *) NULL;

            if (*segs==NULL) *segs = seg;
            else appendseg(*segs, seg);

            i = min(hii, rightend+downset);
            lowlim = i + 1;
        }
    }

    free(H);
    return;
}

int Seg::min(int a, int b)
{
    if (a<b)
    {
        return(a);
    }
    else
    {
        return(b);
    }
}

Sequen* Seg::openwin(Sequen *parent, int start, int length)
{
    Sequen *win;

    if (start<0 || length<0 || start+length>parent->length)
    {
        return((Sequen *) NULL);
    }

    win = (Sequen *) malloc(sizeof(Sequen));

    win->start = start;
    win->length = length;
    win->seq = parent->seq + start;
    win->punctuation = FALSE;

    win->parent = parent;

    win->entropy = -2.;
    win->composition = (int *) NULL;
    win->state = (int *) NULL;

    stateon(win);

    return win;
}

void Seg::stateon(Sequen *win)
{
    register int	aa, nel, c;

    if (win->composition == NULL)
        compon(win);

    win->state = (int *) malloc(21*sizeof(win->state[0]));

    for (aa = nel = 0; aa < 20; ++aa)
    {
        if ((c = win->composition[aa]) == 0)
            continue;
        win->state[nel++] = c;
    }
    for (aa = nel; aa < 21; ++aa)
        win->state[aa] = 0;

    qsort(win->state, nel, sizeof(int), state_cmp);

    win->entropy = entropy_cal(win->state);

    return;
}

void Seg::compon(Sequen *win)
{
    register int	*comp;
    register int	aa;
    register char	*seq, *seqmax;

    win->composition = comp = (int *) calloc(20*sizeof(*comp), 1);
    seq = win->seq;
    seqmax = seq + win->length;

    while (seq < seqmax)
    {
        aa = *seq++;
        if (!aaflag[aa])
            comp[aaindex[aa]]++;
    }

    return;
}

void Seg::enton(Sequen *win)
{
    if (win->state==NULL)
    {
        stateon(win);
    }

    win->entropy = entropy_cal(win->state);

    return;
}

/*----------------------------------------------------------------(seqent)---*/

double* Seg::seqent(Sequen *seq)
{
    Sequen *win;
    double *H;
    int i, first, last;

    if (window>seq->length)
    {
        return((double *) NULL);
    }

    H = (double *) malloc(seq->length*sizeof(double));

    for (i=0; i<seq->length; i++)
    {
        H[i] = -1.;
    }

    win = openwin(seq, 0, window);
    enton(win);

    first = downset;
    last = seq->length - upset;

    for (i=first; i<=last; i++)
    {
        if (seq->punctuation && hasdash(win))
        {
            H[i] = -1;
            shiftwin1(win);
            continue;
        }
        H[i] = win->entropy;
        shiftwin1(win);
    }

    closewin(win);
    return(H);
}

Bool Seg::shiftwin1(Sequen *win)
{
    register int	j, length;
    register int	*comp;

    length = win->length;
    comp = win->composition;

    if ((++win->start + length) > win->parent->length)
    {
        --win->start;
        return FALSE;
    }

    if (!aaflag[j = win->seq[0]])
        decrementsv(win->state, comp[aaindex[j]]--);

    j = win->seq[length];
    ++win->seq;

    if (!aaflag[j])
        incrementsv(win->state, comp[aaindex[j]]++);

    if (win->entropy > -2.)
        win->entropy = entropy_cal(win->state);

    return TRUE;
}

/*---------------------------------------------------------------(hasdash)---*/

Bool Seg::hasdash(Sequen *win)
{
    register char	*seq, *seqmax;

    seq = win->seq;
    seqmax = seq + win->length;

    while (seq < seqmax)
    {
        if (*seq++ == '-')
            return TRUE;
    }
    return FALSE;
}

/*----------------------------------------------------------------(findlo)---*/

int Seg::findlo(int i, int limit, double *H)
{
    int j;

    for (j=i; j>=limit; j--)
    {
        if (H[j]==-1) break;
        if (H[j]>hicut) break;
    }

    return(j+1);
}

/*----------------------------------------------------------------(findhi)---*/

int Seg::findhi(int i, int limit, double *H)
{
    int j;
    for (j=i; j<=limit; j++)
    {
        if (H[j]==-1) break;
        if (H[j]>hicut) break;
    }

    return(j-1);
}

/*------------------------------------------------------------------(trim)---*/

void Seg::trim(Sequen *seq, int *leftend, int *rightend)
{
    Sequen *win;
    double prob, minprob;
    int shift, len, i;
    int lend, rend;
    int minlen;

    /* fprintf(stderr, "%d %d\n", *leftend, *rightend);  */

    lend = 0;
    rend = seq->length - 1;
    minlen = 1;
    if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;

    minprob = 1.;
    for (len=seq->length; len>minlen; len--)
    {
        win = openwin(seq, 0, len);
        i = 0;

        shift = TRUE;
        while (shift)
        {
            prob = getprob(win->state, len);
            if (prob<minprob)
            {
                minprob = prob;
                lend = i;
                rend = len + i - 1;
            }
            shift = shiftwin1(win);
            i++;
        }
        closewin(win);
    }

    /* fprintf(stderr, "%d-%d ", *leftend, *rightend);  */

    *leftend = *leftend + lend;
    *rightend = *rightend - (seq->length - rend - 1);

    /* fprintf(stderr, "%d-%d\n", *leftend, *rightend);  */

    closewin(seq);
    return;
}

/*---------------------------------------------------------------(getprob)---*/

#define LN20	2.9957322735539909
double Seg::getprob(int *sv, int total)
{
    double ans, totseq;

    totseq = ((double) total) * LN20;

    ans = lnass(sv) + lnperm(sv, total) - totseq;

    return(ans);
}

/*----------------------------------------------------------------(lnperm)---*/

double Seg::lnperm(int *sv, int tot)
{
    double ans;
    int i;

    ans = lnfac[tot];

    for (i=0; sv[i]!=0; i++)
    {
        ans -= lnfac[sv[i]];
    }

    return(ans);
}

/*-----------------------------------------------------------------(lnass)---*/

//double Seg::lnass(register int *sv)
double Seg::lnass(int *sv)
{
    double	ans;
    register int	svi, svim1;
    register int	thisclass, total;
    register int    i;

    ans = lnfac[20];
    if (sv[0] == 0)
        return ans;

    total = 20;
    thisclass = 1;
    svim1 = sv[0];
    for (i=0;; svim1 = svi)
    {
        if (++i==20)
        {
            ans -= lnfac[thisclass];
            break;
        }
        else if ((svi = *++sv) == svim1)
        {
            thisclass++;
            continue;
        }
        else
        {
            total -= thisclass;
            ans -= lnfac[thisclass];
            if (svi == 0)
            {
                ans -= lnfac[total];
                break;
            }
            else
            {
                thisclass = 1;
                continue;
            }
        }
    }

    return ans;
}

/*-------------------------------------------------------------(mergesegs)---*/

void Seg::mergesegs(Sequen *seq, Segment *segs)
{
    Segment *seg, *nextseg;
    int len;

    if (overlaps) return;
    if (segs==NULL) return;

    if (segs->begin<hilenmin) segs->begin = 0;

    seg = segs;
    nextseg = seg->next;

    while (nextseg!=NULL)
    {
        if (seg->end>=nextseg->begin)                  /* overlapping segments */
        {
            seg->end = nextseg->end;
            seg->next = nextseg->next;
            free(nextseg);
            nextseg = seg->next;
            continue;
        }
        len = nextseg->begin - seg->end - 1;
        if (len<hilenmin)                             /* short hient segment */
        {
            seg->end = nextseg->end;
            seg->next = nextseg->next;
            free(nextseg);
            nextseg = seg->next;
            continue;
        }
        seg = nextseg;
        nextseg = seg->next;
    }

    len = seq->length - seg->end - 1;
    if (len<hilenmin) seg->end = seq->length - 1;

    return;
}

/*------------------------------------------------------------(singreport)---*/

char* Seg::singreport(Sequen *seq, Segment *segs)
{
    char	*proseq, *proseqmax;
    Segment	*seg;
    int	begin, end, i, ctr;
    char	*maskedseq = new char[seq -> length + 1];
    maskedseq[0] = 0;
    int	add = 0;

    proseq = seq->seq;
    proseqmax = proseq + seq->length;
    upper(proseq, seq->length);

    for (seg=segs; seg!=NULL; seg=seg->next)
    {
        begin = seg->begin;
        end = seg->end;
        memset(proseq + begin, 'x', end - begin +1);
        //memset(proseq + begin, 'X', end - begin +1);
    }

    for (i=0, ctr=0; proseq < proseqmax; ++i, ++ctr, ++proseq)
    {
        maskedseq[add ++] = *proseq;
    }

    maskedseq[add] = 0;

    return maskedseq;
}


/*-------------------------------------------------------------(appendseg)---*/

void Seg::appendseg(Segment *segs, Segment *seg)
{
    Segment *temp;

    temp = segs;
    while (1)
    {
        if (temp->next==NULL)
        {
            temp->next = seg;
            break;
        }
        else
        {
            temp = temp->next;
        }
    }

    return;
}

/*--------------------------------------------------------------(freesegs)---*/

void Seg::freesegs(Segment *segs)
{
    Segment *temp;

    while (segs!=NULL)
    {
        temp = segs->next;
        free(segs);
        segs = temp;
    }
}

void Seg::closewin(Sequen *win)
{
    if (win==NULL) return;

    if (win->state!=NULL)       free(win->state);
    if (win->composition != NULL) free(win->composition);
    //if (win->seq != NULL)	free(win->seq); //seq is a pointer!

    free(win);
    return;
}


#define LN2	0.69314718055994530941723212145818

void Seg::entropy_init(int window)
{
    int	i;
    double	x, xw;

    entray = (double *)malloc((window+1) * sizeof(*entray));
    xw = window;
    for (i = 1; i <= window; ++i)
    {
        x = i / xw;
        entray[i] = -x * log(x) / LN2;
    }
    thewindow = window;
}

double Seg::entropy_cal(register int *sv)
{
    int	*sv0 = sv;
    register double	ent;
    register int	i, total;
    register int	*svmax;
    register double	xtotrecip, xsv;

    for (total = 0; (i = *sv) != 0; ++sv)
        total += i;
    svmax = sv;
    ent = 0.0;
    if (total == thewindow)
    {
        for (sv = sv0; sv < svmax; )
        {
            ent += entray[*sv++];
        }
        return ent;
    }
    if (total == 0)
        return 0.;

    xtotrecip = 1./(double)total;
    for (sv = sv0; sv < svmax; )
    {
        xsv = *sv++;
        ent += xsv * log(xsv * xtotrecip);
    }
    return -ent * xtotrecip / LN2;
}

Sequen* Seg::openseq(const char *fastaseq)
{
    Sequen* win = (Sequen *) malloc(sizeof(Sequen));
    int	length = strlen(fastaseq);
    win->length = length;
    win->seq = (char *) malloc(sizeof(char)*length + 1);
    //memcpy(win->seq, (parent->seq)+start, length);
    memcpy(win->seq, fastaseq, length);
    win->seq[length] = '\0';

    win->parent = (Sequen *) NULL;
    win->state = NULL;
    win->composition = NULL;

    stateon(win);
    return win;
}

void Seg::closeseq(Sequen *seq)
{
    if (seq==NULL) return;

    if (seq->state!=NULL)       free(seq->state);
    if (seq->seq != NULL)	free(seq->seq);
    if (seq->composition != NULL) free(seq->composition);
    free(seq);

    return;
}

void Seg::decrementsv(register int *sv, register int thisclass)
{
    register int	svi;

    while ((svi = *sv++) != 0)
    {
        if (svi == thisclass && *sv < thisclass)
        {
            sv[-1] = svi - 1;
            break;
        }
    }
}

/*-----------------------------------------------------------(incrementsv)---*/

void Seg::incrementsv(register int *sv, int thisclass)
{
    for (;;)
    {
        if (*sv++ == thisclass)
        {
            sv[-1]++;
            break;
        }
    }
}

void Seg::upper(register char *string, size_t len)
{
    register char   *stringmax, c;

    for (stringmax = string + len; string < stringmax; ++string)
        if (islower(c = *string))
            *string = toupper(c);
}
