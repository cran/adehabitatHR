#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


/* ***********************************************************************
 *                                                                       *
 *                          Declaration of functions                     *
 *                                                                       *
 * ********************************************************************* */



/* Functions coming from the package ade4 */

void vecpermut (double *A, int *num, double *B);
double alea (void);
void aleapermutvec (double *a);
void trirapideintswap (int *v, int i, int j);
void trirapideint (int *x , int *num, int gauche, int droite);
void sqrvec (double *v1);
void getpermutation (int *numero, int repet);
void prodmatABC (double **a, double **b, double **c);
void prodmatAtAB (double **a, double **b);
void prodmatAtBC (double **a, double **b, double **c);
void prodmatAAtB (double **a, double **b);
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut);
void taballoc (double ***tab, int l1, int c1);
void vecalloc (double **vec, int n);
void vecintalloc (int **vec, int n);
void freetab (double **tab);
void freevec (double *vec);
void freeintvec (int *vec);
void matcentrage (double **A, double *poili, char *typ);
void matmodiffc (double **tab, double *poili);
void matmodifcp (double **tab, double *poili);
void matmodifcs (double **tab, double *poili);
void matmodifcn (double **tab, double *poili);
void matmodifcm (double **tab, double *poili);
void DiagobgComp (int n0, double **w, double *d, int *rang);





/* Functions from the package adehabitat */
void epa(double *X, double *Y, double *xl, double *yl, 
	 double *val, double *fen);
void kernelhr(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *fen, double *xlo, double *ylo);
void CVmise(int *nloc, double *xlo, double *ylo,
	    double *hvec, double *CV, int *nhteste);
void calcvolume(double *grille, int *ncolgri, int *nliggri, double *cellsize);
void integrno(double *XG, double *X1, double *X2, 
	      double *T, double *sig1,
	      double *sig2, double *alpha, double *res);
void kernelbb(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *sig1, double *sig2, 
	      double *xlo, double *ylo, double *Tr, int *controlbox, 
	      int *nalpha);
void trouveclustmin(double **xy, int *clust, int *lo1, int *lo2,
		    int *lo3, double *dist);
void trouveclustminr(double *xyr, int *nr, int *clustr, int *lo1, int *lo2,
		     int *lo3, double *dist);
void nndistclust(double **xy, double *xyp, double *dist);
void parclust(double **xy, int *clust, int *noclust, 
	      int *noloc, double *dist);
void trouveminclust(double **xy, int *liclust, int *clust, 
		    int *noclust, int *noloc, double *dist);
void choisnvclust(double **xy, int *liclust, int *clust, int *ordre);
void clusterhr(double **xy, int *facso, int *nolocso, int *cluso);
void longfacclust(double **xy, int *len2);
void longfacclustr(double *xyr, int *nr, int *len2);
void clusterhrr(double *xyr, int *nr, int *facsor, 
		int *nolocsor, int *clusor, int *len);
void kcprcirc(double **xyd, double *h, double *x, double t, 
	      double *val);
void kcprlin(double **xyd, double *h, double *x, double t, 
	     double *val);
void kernelkcr(double *xydr, double *tcalcr, int *nlr, double *gridr,
	       double *xgri, double *ygri, int *nliggri, int *ncolgri, 
	       double *hr, int *circularr);
void findmaxgrid(double *grille, int *nlig, int *ncol);






/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of ADE-4                    *****
 *********               --------------------                    *****
 *********************************************************************
 *********************************************************************
 */



/**************************/
double alea (void)
{
    double w;
    w = ((double) rand())/ (double)RAND_MAX;
    return (w);
}

/*************************/
void aleapermutvec (double *a)
{
    /* Randomly permutes the elements of a vector a
       Manly p. 42 The vector is modified
       from Knuth 1981 p. 139 */
    int lig, i,j, k;
    double z;
    
    lig = a[0];
    for (i=1; i<=lig-1; i++) {
	j=lig-i+1;
	k = (int) (j*alea()+1);
	/* k = (int) (j*genrand()+1); */
	if (k>j) k=j;
	z = a[j];
	a[j]=a[k];
	a[k] = z;
    }
}


/*******************/	
void vecpermut (double *A, int *num, double *B)
{
/*---------------------------------------
 * A is a vector with n elements
 * B is a vector with n elements
 * num is a random permutation of the n first integers
 * B contains in output the permuted elements of A
 * ---------------------------------------*/
    
    int lig, lig1, lig2, i, k;
    
    lig = A[0];
    lig1 = B[0];
    lig2 = num[0];
    
    
    if ( (lig!=lig1) || (lig!=lig2) ) {
	/* err_message ("Illegal parameters (vecpermut)");
	   closelisting(); */
    }
    
    for (i=1; i<=lig; i++) {
	k=num[i];
	B[i] = A[k];
    }
}

/********* Centring accrding to row weights poili **********/	
void matcentrage (double **A, double *poili, char *typ)
{
    
    if (strcmp (typ,"nc") == 0) {
	return;
    } else if (strcmp (typ,"cm") == 0) {
	matmodifcm (A, poili);
	return;
    } else if (strcmp (typ,"cn") == 0) {
	matmodifcn (A, poili);
	return;
    } else if (strcmp (typ,"cp") == 0) {
	matmodifcp (A, poili);
	return;
    } else if (strcmp (typ,"cs") == 0) {
	matmodifcs (A, poili);
	return;
    } else if (strcmp (typ,"fc") == 0) {
	matmodiffc (A, poili);
	return;
    } else if (strcmp (typ,"fl") == 0) {
	matmodifcm (A, poili);
	return;
    }
}

/*********************/
void matmodifcm (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a complete disjonctive table with n rows and m columns
 * poili is a vector with n components
 * The process returns tab centred by column
 * with weighting poili (sum=1)
 * centring type multple correspondances
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    for (i=1;i<=l1;i++) tab[i][j] = 0;
	} else {
	    
	    for (i=1;i<=l1;i++) {
		z = tab[i][j]/x - 1.0;
		tab[i][j] = z;
	    }
	}
    }
    freevec (poimoda);
}

/*********************************************************/
void matmodifcn (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows and p columns
 * poili is a vector with n components
 * the function returns tab normed by column
 * with the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid, x, z, y, v2;
    int 			i, j, l1, c1;
    double		*moy, *var;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    
    vecalloc(&moy, c1);
    vecalloc(&var, c1);
    
    
/*--------------------------------------------------
 * centred and normed table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    for (i=1;i<=l1;i++) {
	poid=poili[i];
	for (j=1;j<=c1;j++) {
	    x = tab[i][j] - moy[j];
	    var[j] = var[j] + poid * x * x;
	}
    }
    
    for (j=1;j<=c1;j++) {
	v2 = var[j];
	if (v2<=0) v2 = 1;
	v2 = sqrt(v2);
	var[j] = v2;
    }
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	y = var[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    z = z / y;
	    tab[j][i] = z;
	}
    }
    
    freevec(moy);
    freevec(var);
    
}

/*********************************************************/
void matmodifcs (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows, p columns
 * poili is a vector with n components
 * The function returns tab standardised by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
	double		x,poid, z, y, v2;
	int 			i, j, l1, c1;
	double		*var;
	
	l1 = tab[0][0];
	c1 = tab[1][0];
	vecalloc(&var, c1);
	

/*--------------------------------------------------
 * calculation of the standardised table
 --------------------------------------------------*/
	
	for (i=1;i<=l1;i++) {
	    poid=poili[i];
	    for (j=1;j<=c1;j++) {
		x = tab[i][j];
		var[j] = var[j] + poid * x * x;
	    }
	}
	
	for (j=1;j<=c1;j++) {
	    v2 = var[j];
	    if (v2<=0) v2 = 1;
	    v2 = sqrt(v2);
	    var[j] = v2;
	}
	
	for (i=1;i<=c1;i++) {
	    y = var[i];
	    for (j=1;j<=l1;j++) {
		z = tab[j][i];
		z = z / y;
		tab[j][i] = z;
	    }
	}
	freevec(var);
}


/**********/
void matmodifcp (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and p colonnes
 * poili is a vector with n components
 * The function returns tab centred by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, c1;
    double		*moy, x, z;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    vecalloc(&moy, c1);
    
    
/*--------------------------------------------------
 * Centred table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    tab[j][i] = z;
	}
    }
    freevec(moy);
}

/*********************/
void matmodiffc (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and m columns
 * of number >=0
 * poili is a vector with n components
 * The function returns tab doubly centred
 * for the weighting poili (sum=1)
 * centring type simple correspondance analysis
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	x = 0;
	for (j=1;j<=m1;j++) {
	    x = x + tab[i][j];
	}
	if (x!=0) {
	    for (j=1;j<=m1;j++) {
		tab[i][j] = tab[i][j]/x;
	    }
	}	
    }
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    /* err_message("column has a nul weight (matmodiffc)"); */
	}
	
	for (i=1;i<=l1;i++) {
	    z = tab[i][j]/x - 1.0;
	    tab[i][j] = z;
	}
    }
    freevec (poimoda);
}









/*****************/
void getpermutation (int *numero, int repet)
/*----------------------
 * affects a random permutation of the first n integers
 * in an integer vector of length n
 * First vecintalloc is needed
 * *numero is a vector of integer
 * repet is an integer which can take any arbitrary value
 * used in the seed of the pseudo-random number generation process
 * if it is increased in repeated calls (e.g. simulation), it is ensured that
 * two calls returns different results (seed=clock+repet)
 ------------------------*/
{
    int i, n, seed;
    int *alea;
    
    n=numero[0];
    vecintalloc (&alea,n);
    
    /*-------------
     * numbering in numero
     -----------*/
    for (i=1;i<=n;i++) {
	numero[i]=i;
    }
    
    /*-------------
     * affects random numbers in alea
     ----------------*/
    seed = clock();
    seed = seed + repet;
    srand(seed);
    for (i=1;i<=n;i++) {
	alea[i]=rand();
    }
    
    trirapideint (alea , numero, 1, n);
    freeintvec (alea);
}

/*****************************************/
/* Sorting: used in getpermutation */

void trirapideint (int *x , int *num, int gauche, int droite)
{
    int j, dernier, milieu, t;
    
    if ( (droite-gauche)<=0) return;
    
    milieu = (gauche+droite)/2;
    trirapideintswap (x, gauche, milieu);
    trirapideintswap (num, gauche, milieu);
    
    t=x[gauche];
    dernier=gauche;
    for (j = gauche+1; j<=droite; j++) {
	if (x[j] < t) {
	    dernier = dernier + 1;
	    trirapideintswap (x, dernier, j);	
	    trirapideintswap (num, dernier, j);
	}
    }
    trirapideintswap (x, gauche, dernier);
    trirapideintswap (num, gauche, dernier);
    
    trirapideint (x, num, gauche, dernier-1);
    trirapideint (x, num, dernier+1, droite);
    
}

/**************************************/
/* Sorting: used in trirapideint */

void trirapideintswap (int *v, int i, int j)
{
    int provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}

/***********************************************************************/
void sqrvec (double *v1)
/*--------------------------------------------------
 * Square root of the elements of a vector
 --------------------------------------------------*/
{
    int i, c1;
    double v2;
    
    c1 = v1[0];
    
    for (i=1;i<=c1;i++) {
	v2 = v1[i];
	/* if (v2 < 0.0) err_message("Error: Square root of negative number (sqrvec)"); */
	v2 = sqrt(v2);
	v1[i] = v2;
    }
}

/***********************************************************************/
void DiagobgComp (int n0, double **w, double *d, int *rang)
/*--------------------------------------------------
 * Eigenstructure of a matrix. See
 * T. FOUCART Analyse factorielle de tableaux multiples,
 * Masson, Paris 1984,185p., p. 62. D'apr?s VPROP et TRIDI,
 * de LEBART et coll.
 --------------------------------------------------*/
{
    double			*s;
    double			a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
    double			dble;
    int				ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
    
    vecalloc(&s, n0);
    a = 0.000000001;
    ni = 100;
    if (n0 == 1) {
	d[1] = w[1][1];
	w[1][1] = 1.0;
	*rang = 1;
	freevec (s);
	return;
    }
    
    for (i2=2;i2<=n0;i2++) {
	
	b=0.0;
	c=0.0;
	i=n0-i2+2;
	k=i-1;
	if (k < 2) goto Et1;
	for (l=1;l<=k;l++) {
	    c = c + fabs((double) w[i][l]);
	}
	if (c != 0.0) goto Et2;
	
    Et1:	s[i] = w[i][k];
	goto Etc;
	
    Et2:	for (l=1;l<=k;l++) {
	x = w[i][l] / c;
	w[i][l] = x;
	b = b + x * x;
    }
	xp = w[i][k];
	ix = 1;
	if (xp < 0.0) ix = -1;
		
/*		q = -sqrt(b) * ix; */
	dble = b;
	dble = -sqrt(dble);
	q = dble * ix;
	
	s[i] = c * q;
	b = b - xp * q;
	w[i][k] = xp - q;
	xp = 0;
	for (m=1;m<=k;m++) {
	    w[m][i] = w[i][m] / b / c;
	    q = 0;
	    for (l=1;l<=m;l++) {
		q = q + w[m][l] * w[i][l];
	    }
	    m1 = m + 1;
	    if (k < m1) goto Et3;
	    for (l=m1;l<=k;l++) {
		q = q + w[l][m] * w[i][l];
	    }
	    
	Et3:		s[m] = q / b;
	    xp = xp + s[m] * w[i][m];
	}
	bp = xp * 0.5 / b;
	for (m=1;m<=k;m++) {
	    xp = w[i][m];
	    q = s[m] - bp * xp;
	    s[m] = q;
	    for (l=1;l<=m;l++) {
		w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
	    }
	}
	for (l=1;l<=k;l++) {
	    w[i][l] = c * w[i][l];
	}
	
    Etc:	d[i] = b;
    } /* for (i2=2;i2<n0;i2++) */
    
    s[1] = 0.0;
    d[1] = 0.0;
    
    for (i=1;i<=n0;i++) {
	
	k = i - 1;
	if (d[i] == 0.0) goto Et4;
	for (m=1;m<=k;m++) {
	    q = 0.0;
	    for (l=1;l<=k;l++) {
		q = q + w[i][l] * w[l][m];
	    }
	    for (l=1;l<=k;l++) {
		w[l][m] = w[l][m] - q * w[l][i];
	    }
	}
	
    Et4:	d[i] = w[i][i];
	w[i][i] = 1.0;
	if (k < 1) goto Et5;
	for (m=1;m<=k;m++) {
	    w[i][m] = 0.0;
	    w[m][i] = 0.0;
	}
	
    Et5:;
    }
    
    for (i=2;i<=n0;i++) {
	s[i-1] = s[i];
    }
    s[n0] = 0.0;
    
    for (k=1;k<=n0;k++) {
	
	m = 0;
	
    Et6: 	for (j=k;j<=n0;j++) {
	if (j == n0) goto Et7;
	ab = fabs((double) s[j]);
	ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
	if (ab < ep) goto Et7;
    }
	
    Et7: 	isnou = 1;
	h = d[k];
	if (j == k) goto Eta;
	if (m < ni) goto Etd;
	
	/* err_message("Error: can't compute matrix eigenvalues"); */
	
    Etd:	m = m + 1;
	q = (d[k+1]-h) * 0.5 / s[k];
	
/*		t = sqrt(q * q + 1.0); */
	dble = q * q + 1.0;
	dble = sqrt(dble);
	t = dble;
	
	if (q < 0.0) isnou = -1;
	q = d[j] - h + s[k] / (q + t * isnou);
	u = 1.0;
	v = 1.0;
	h = 0.0;
	jk = j-k;
	for (ijk=1;ijk<=jk;ijk++) {
	    i = j - ijk;
	    xp = u * s[i];
	    b = v * s[i];
	    if (fabs((double) xp) < fabs((double) q)) goto Et8;
	    u = xp / q;
	    
/*			t = sqrt(u * u + 1); */
	    dble = u * u + 1.0;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = q * t;
	    v = 1 / t;
	    u = u * v;
	    goto Et9;
	    
	Et8:		v = q / xp;
	    
/*			t = sqrt(1 + v * v); */
	    dble = 1.0 + v * v;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = t * xp;
	    u = 1 / t;
	    v = v * u;
	    
	Et9:
	    q = d[i+1] - h;
	    t = (d[i] - q) * u + 2.0 * v * b;
	    h = u * t;
	    d[i+1] = q + h;
	    q = v * t - b;
	    for (l=1;l<=n0;l++) {
		xp = w[l][i+1];
		w[l][i+1] = u * w[l][i] + v * xp;
		w[l][i] = v * w[l][i] - u * xp;
	    }
	}
	d[k] = d[k] - h;
	s[k] = q;
	s[j] = 0.0;
	
	goto Et6;
	
    Eta:;
    } /* for (k=1;k<=n0;k++) */
    
    for (ij=2;ij<=n0;ij++) {
	
	i = ij - 1;
	l = i;
	h = d[i];
	for (m=ij;m<=n0;m++) {
	    if (d[m] >= h) {
		l = m;
		h = d[m];
	    }
	}
	if (l == i) {
	    goto Etb;
	} else {
	    d[l] = d[i];
	    d[i] = h;
	}
	for (m=1;m<=n0;m++) {
	    h = w[m][i];
	    w[m][i] = w[m][l];
	    w[m][l] = h;
	}
	
    Etb:;
    } /* for (ij=2;ij<=n0;ij++) */
    
    /* final:; */
    *rang = 0;
    for (i=1;i<=n0;i++) {
	/*
	  if (d[i] / d[1] < 0.00001) d[i] = 0.0;
	  if (d[i] != 0.0) *rang = *rang + 1;
	*/
	if (d[i] > 0.0) *rang = *rang + 1;
    }
    freevec(s);
} /* DiagoCompbg */







/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
* Matrix product AB
--------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (i=1;i<=lig;i++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (j=1;j<=col;j++) {
		s = s + a[i][j] * b[j][k];
	    }
	    c[i][k] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
* Matrix product AtA
--------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=j;k<=col;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][k] * a[i][j];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtBC (double **a, double **b, double **c)
/*--------------------------------------------------
 * Matrix product AtB
 --------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][j] * b[i][k];
	    }
	    c[j][k] = s;
	}		
    }
}


/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
 * Matrix product B = AAt
 --------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=lig;j++) {
	for (k=j;k<=lig;k++) {
	    s = 0;
	    for (i=1;i<=col;i++) {
		s = s + a[j][i] * a[k][i];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/*******************/
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut)
/*--------------------------------------------------
 * Produit matriciel AtB
 * les lignes de B sont permutÚes par la permutation permut
 --------------------------------------------------*/
{
    int j, k, i, i0, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		i0 = permut[i];
		s = s + a[i][j] * b[i0][k];
	    }
	    c[j][k] = s;
	}		
    }
}

/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
 * Dynamic Memory Allocation for a table (l1, c1)
 --------------------------------------------------*/
{
    int i, j;
    
    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
	for (i=0;i<=l1;i++) {
	    if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
		return;
		for (j=0;j<i;j++) {
		    free(*(*tab+j));
		}
	    }
	}
    }
    
    **(*tab) = l1;
    **(*tab+1) = c1;
}

/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
 * Memory Allocation for a vector of length n
 --------------------------------------------------*/
{
    if ( (*vec = (double *) calloc(n+1, sizeof(double))) != 0) {
	**vec = n;
	return;
    } else {
	return;
    }
}

/*****************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
 * Memory allocation for an integer vector of length  n
 --------------------------------------------------*/
{
    if ( (*vec = (int *) calloc(n+1, sizeof(int))) != NULL) {
	**vec = n;
	return;
    } else {
	return;
    }
}

/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
 * Free memory for a table
 --------------------------------------------------*/
{
    int 	i, n;
    
    n = *(*(tab));
    for (i=0;i<=n;i++) {
	free((char *) *(tab+i) );
    }
    free((char *) tab);
}

/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
 * Free memory for a vector
 --------------------------------------------------*/
{
    free((char *) vec);	
}

/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
* Free memory for an integer  vector
--------------------------------------------------*/
{
    
    free((char *) vec);
    
}














/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of adehabitat               *****
 *********               -------------------------               *****
 *********************************************************************
 *********************************************************************
 */




/* ****************************************************************
   *                                                              *
   *   epa: bivariate normal kernel                               *
   *                                                              *
   **************************************************************** */

int selectptsbo(double *xl, double *yl, double *box, 
		int *indcons)
{
    int i,nl,cons,k;
    nl = (int) xl[0];
    
    k=0;
    for (i=1; i<=nl; i++) {
	cons = 0;
	if (xl[i] < box[1]) {
	    if (xl[i] > box[2]) {
		if (yl[i] < box[3]) {
		    if (yl[i] > box[4]) {
			cons = 1;
		    }
		}
	    }
	}
	if (cons == 1) {
	    k++;
	    indcons[k] = i;
	}
    }
    return(k);
}

void epa(double *X, double *Y, double *xl, double *yl, 
	 double *val, double *fen)
{
    /* Declaration of local variables */
    int k, nl, *indcons, ncons;
    double *xy, kx, di2, h, *box;
    
    /* Bases */
    nl = (int) xl[0];
    vecalloc(&xy, 2);
    vecintalloc(&indcons, nl);
    vecalloc(&box, 4);
    *val = 0;
    h = *fen;
    kx = 0;
    
    /* Keep only the points no further than 4*fen of the current pixel */
    box[1] = *X + (4 * h);
    box[2] = *X - (4 * h);
    box[3] = *Y + (4 * h);
    box[4] = *Y - (4 * h);
    ncons = selectptsbo(xl, yl, box, indcons);
        
    
    /* The bivariate normal kernel */
    if (ncons>0) {
	for (k=1; k<=ncons; k++) {
	    xy[1] = (xl[indcons[k]] - *X);
	    xy[2] = (yl[indcons[k]] - *Y);
	    di2 = xy[1]*xy[1] + xy[2]*xy[2];
	    
	    kx = exp(-di2/(2*h*h));
	    *val = *val + kx;
	}
	*val = *val * (1/(((double) nl)*h*h*2*3.14159265359));
    } else {
	*val=0;
    }
    freevec(xy);
    freeintvec(indcons);
    freevec(box);
}




/* ****************************************************************
   *                                                              *
   *             Kernel home range                                *
   *                                                              *
   **************************************************************** */


void kernelhr(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *fen, double *xlo, double *ylo)
{
    /* Declaration of local variables */
    int i, j, k, ncg, nlg, nlo;
    double **gri, *xg, *yg, *xl, *yl, X, Y, tmp;
    
    /* Memory Allocation */
    ncg = *ncolgri;
    nlg = *nliggri;
    nlo = *nloc;
    tmp = 0;
    
    taballoc(&gri,nlg, ncg);
    vecalloc(&xg, nlg);
    vecalloc(&yg, ncg);
    vecalloc(&xl, nlo);
    vecalloc(&yl, nlo);
    
    /* R objects -> C objects */
  
    for (i=1; i<=nlo; i++) {
	xl[i] = xlo[i-1];
	yl[i] = ylo[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xg[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	yg[i] = ygri[i-1];
    }
    
    /* loop on the grid */
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    X = xg[i];
	    Y = yg[j];
	    epa(&X, &Y, xl, yl, &tmp, fen);
	    gri[i][j] = tmp;
	}
    }
    
    /* C objects -> R objects */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    grille[k] = gri[i][j];
	    k++;
	}
    }

    /* Memory Free */
    freetab(gri);
    freevec(xg);
    freevec(yg);
    freevec(xl);
    freevec(yl);
}




/* ****************************************************************
   *                                                              *
   *   Epanechnikov estimation thanks to a kernel                 *
   *                                                              *
   **************************************************************** */

void epanechnikov(double *Xo, double *Yo, double *xg, double *yg, 
		  double *fen, double **grille, int nlo)
{
    /* Declaration of local variables */
    
    int i, j, ncg, nlg, imin, imax, jmin, jmax;
    double X, Y, h, *xgb, *ygb, tmp;
    
    nlg = xg[0];
    ncg = yg[0];
    h = *fen;
    X = *Xo;
    Y = *Yo;
    vecalloc(&xgb, nlg);
    vecalloc(&ygb, ncg);
    imin=0;
    jmin=0;
    imax=0;
    jmax=0;
    
    /* Computes again the values xg and yg */
    for (i=1; i<=nlg; i++) {
	xgb[i] = abs(xg[i]-X);
	if (xgb[i] < h) {
	    if (imin == 0) {
		imin = i;
	    }
	}
	if (xgb[i] > h) {
	    if (imin != 0) {
		imax = i;
	    }
	}
    }
    for (i=1; i<=ncg; i++) {
	ygb[i] = abs(yg[i]-Y);
	if (ygb[i] < h) {
	    if (jmin == 0) {
		jmin = i;
	    }
	}
	if (ygb[i] > h) {
	    if (jmin != 0) {
		jmax = i;
	    }
	}
    }
    
    for (i=imin; i<=imax; i++) {
	for (j=jmin; j<=jmax; j++) {
	    tmp = ( (xgb[i] / h) * (xgb[i] / h) ) + ( (ygb[j] / h) * (ygb[j] / h) );
	    if (tmp < 1) {
		grille[i][j] = grille[i][j] + 
		    2 * (1 - tmp) / (3.14159265359 * nlo *  h * h);
	    }
	}
    }
    
    freevec(xgb);
    freevec(ygb);
}


/* For R */

void kernepan(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *fen, double *xlo, double *ylo)
{
    /* Declaration */
    int i, j, k, ncg, nlg, nlo;
    double **gri, *xg, *yg, *xl, *yl, X, Y, tmp;
    
    /* Memory Allocation */
    ncg = *ncolgri;
    nlg = *nliggri;
    nlo = *nloc;
    tmp = 0;
    
    taballoc(&gri,nlg, ncg);
    vecalloc(&xg, nlg);
    vecalloc(&yg, ncg);
    vecalloc(&xl, nlo);
    vecalloc(&yl, nlo);
  
    /* R to C */
    
    for (i=1; i<=nlo; i++) {
	xl[i] = xlo[i-1];
	yl[i] = ylo[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xg[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	yg[i] = ygri[i-1];
    }
  
    /* Loop on the relocations */
    for (i=1; i<=nlo; i++) {
	X = xl[i];
	Y = yl[i];
	epanechnikov(&X, &Y, xg, yg, fen, gri, nlo);
    }
    
    /* C to R */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    grille[k] = gri[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(gri);
    freevec(xg);
    freevec(yg);
    freevec(xl);
    freevec(yl);
}



/* ****************************************************************
   *                                                              *
   *   Find Minimum LSCV                                          *
   *                                                              *
   **************************************************************** */


double L(double smooth, int nlo, int ndist, double *dists) 
{
  int ii;
  double resL,n;  

  n = (double) nlo;
  resL = 0.;
  
  for(ii=0; ii < ndist; ii++){ 
      resL+= (exp(-pow(dists[ii],2)/(4. * pow(smooth,2)))) - (4. * (exp(-pow(dists[ii],2)/(2. * pow(smooth,2.)))));
  }
  
  resL = 1./(3.14159265359 * pow(smooth,2.) * n) + (2*resL -3*n)/(3.14159265359 * 4. * pow(smooth,2.) * pow(n, 2.));

  return(resL);
}




double euclidean_distance(double x1, double y1, double x2, double y2)
{
    double out = 0.0;
    out = pow((x2-x1), 2) + pow((y2-y1), 2);
    return sqrt(out);
}


double comdi(double *x, double *y, double *dists, 
	     int n)
{
    int ii,jj,kk;
    int nn;

    nn = n*(n-1)/2;
    kk=0;
    
    for(ii=1; ii <= n-1; ii++){
	for(jj=ii+1; jj<=n; jj++){
	    dists[kk] = euclidean_distance(x[ii], y[ii], x[jj],y[jj]);
	    kk++;
	}
    }
    return (kk);
}



void CVmise(int *nloc, double *xlo, double *ylo,
	    double *hvec, double *CV, int *nhteste)
{
    /* Declaration */
    int i, nlo, nh, ndist;
    double *xl, *yl, h, di2, *dists;
    
    /* Allocation de mémoire */
    nlo = *nloc;
    nh = *nhteste;
    
    vecalloc(&xl, nlo);
    vecalloc(&yl, nlo);
    vecalloc(&dists, (nlo-1)*nlo);
    
    /* R to C */
    for (i=1; i<=nlo; i++) {
	xl[i] = xlo[i-1];
	yl[i] = ylo[i-1];
    }
    
    /* Compute the distances */
    ndist=comdi(xl, yl, dists, nlo);
    
    /* Loop on the window of h */
    for (i=1; i<=nh; i++) {
	h = hvec[i-1];
	CV[i-1]=L(h, nlo, ndist, dists);
    }
	
    /* Free Memory */
    freevec(dists);
    freevec(xl);
    freevec(yl);
}

  



/* ****************************************************************
   *                                                              *
   *        Computation of the volume under the UD                *
   *                                                              *
   **************************************************************** */


void calcvolume(double *grille, int *ncolgri, int *nliggri, double *cellsize)
{
    /* Declaration */
    int i, j, k, nl, nc;
    double cs, **gri;
    
    /* Memory Allocation */
    nl = *nliggri;
    nc = *ncolgri;
    cs = *cellsize;
    
    taballoc(&gri, nl, nc);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    gri[i][j] = grille[k];
	    k++;
	}
    }
    
    /* Volume of the grid */
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    gri[i][j] = gri[i][j]*cs*cs;
	}
    }
    
    /* C to R */
    k = 0;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    grille[k] = gri[i][j];
	    k++;
	}
    }
    
    /* Free Memory */
    freetab(gri);
}







/* *********************************************************************
 *                                                                     *
 *                   Brownian bridge kernel                            *
 *                                                                     *
 ***********************************************************************/


/* Function normal 2D for brownian bridge */

void norm2d(double x1, double y1, double moyx, double moyy,
	    double var, double *res)
{
    double cste;
    cste = (1 / (2.0 * 3.141592653589793238 * var));
    cste = cste * exp( (-1.0 / (2.0 * var)) * (((x1 - moyx) * (x1 - moyx))+((y1 - moyy) * (y1 - moyy))));
    *res = cste;
}


double maxh(double sig1, double sig2, double *alpha, double maxt)
{
    int na,i;
    double res, tmp, a;
    
    res = 0;
    na = alpha[0];
    
    for (i = 1; i <= na; i++) {
	a = alpha[i];
	tmp = (maxt * a * (1 - a) * sig1) + 
	    ((pow(a,2) + pow((1-a), 2)) * sig2);
	if (tmp > res)
	    res = tmp;
    }
    return(sqrt(res));
}


double maxdt(double *T)
{
    int i,nt;
    double res, tmp;

    res = 0;
    nt = T[0];

    for (i = 2; i <= nt; i++) {
	if ((T[i]-T[i-1]) > res)
	    res = (T[i]-T[i-1]);
    }
    return(res);
}

/* keeps all the steps for which at least one relocation is 
   available in the box */
int consdanslabox(double *Xg, double **xy, 
		  int nl, int *indcons, double maxvh, int controlbox)
{
    int i,k,cons;
    double tmp1, tmp2, a, b;
    /* On a besoin d'une boucle sur les pas */
    k=0;
    
    for (i = 1; i<nl; i++) {
	
	cons = 0;
	
	if (xy[i][1] > (Xg[1] - (controlbox * maxvh)) ) {
	    if (xy[i][1] < (Xg[1] + (controlbox * maxvh)) ) {
		if (xy[i][2] > (Xg[2] - (controlbox * maxvh)) ) {
		    if (xy[i][2] < (Xg[2] + (controlbox * maxvh)) ) {
			cons = 1;
		    }
		}
	    }
	}
	if (xy[i+1][1] > (Xg[1] - (controlbox * maxvh)) ) {
	    if (xy[i+1][1] < (Xg[1] + (controlbox * maxvh)) ) {
		if (xy[i+1][2] > (Xg[2] - (controlbox * maxvh)) ) {
		    if (xy[i+1][2] < (Xg[2] + (controlbox * maxvh)) ) {
			cons = 1;
		    }
		}
	    }
	}
	
	if (cons == 0) {
	    a = (xy[i+1][2] - xy[i][2]) / (xy[i+1][1] - xy[i][1]);
	    b = xy[i+1][2] - a * xy[i+1][1];
	    tmp1 = a * (Xg[1] - (controlbox * maxvh)) + b;
	    tmp2 = a * (Xg[1] + (controlbox * maxvh)) + b;
	    
	    if (tmp1 <= (Xg[2] + (controlbox * maxvh))) {
		if (tmp1 >= (Xg[2] - (controlbox * maxvh))) {
		    cons = 1;
		}
	    }
	    
	    if (tmp2 <= (Xg[2] + (controlbox * maxvh))) {
		if (tmp2 >= (Xg[2] - (controlbox * maxvh))) {
		    cons = 1;
		}
	    }
	}
	
	
	if (cons==1) {
	    k++;
	    indcons[k]=i;
	}
	
    }
    
    return(k);
}



/* Integral of norm2d on alpha */
void integrno(double *XG, double *X1, double *X2, 
	      double *T, double *sig1,
	      double *sig2, double *alpha, double *res)
{
    /* Declaration */
    int i, na;
    double *val, tmp, *XX, var, nx1, ny1, nx2, ny2, ny, moyx, moyy, al;
    
    /* Memory allocation */
    na = alpha[0];
    vecalloc(&val, na);
    vecalloc(&XX, 2);
    
    XX[1] = X2[1] - X1[1];
    XX[2] = X2[2] - X1[2];
    *res = 0;
    
    
    /* loop for the computation of the value */
    for (i = 1; i<= na; i++) {
	al = alpha[i];
	
	var = (*T) * al * (1.0 - al) * (*sig1);
	var = var + (((al * al) + ((1.0 - al) * (1.0 - al))) * (*sig2));
	
	moyx = X1[1] + al * XX[1];
	moyy = X1[2] + al * XX[2];
	
	norm2d(XG[1], XG[2], moyx, moyy, var, &tmp);
	
	val[i] = tmp;
    }
    
    /* loop for the computation of the integral */
    for (i = 2; i<= na; i++) {
	nx1 = alpha[i-1];
	ny1 = val[i-1];
	nx2 = alpha[i];
	ny2 = val[i];
	ny = ny1;
	if (ny2 <= ny1)
	    ny = ny2;
	*res = *res + (nx2 - nx1) * (ny + (abs(ny2 - ny1) / 2));
    }
    
    /* Free memory */
    freevec(val);
    freevec(XX);
}




/* Computes UD at a node of the grid */
void udbbnoeud(double *XG, double **XY, double *T, double *sig1,
	       double *sig2, double *alpha, double *res, int ncons, 
	       int *indcons)
{
    /* Declaration */
    int i, nlo;
    double *Xtmp1, *Xtmp2, dt, poids, dttot, tmp;
    
    /* Memory allocation */
    vecalloc(&Xtmp1, 2);
    vecalloc(&Xtmp2, 2);
    nlo = XY[0][0];
    dttot = T[nlo] - T[1];
    *res = 0;
    
    /* for each step */
    for (i = 1; i <= ncons; i++) {
	
	/* Computes weights and time lags */
	dt = T[indcons[i]+1] - T[indcons[i]];
	poids = dt / dttot;
	
	/* Output of the relocation values at i, and use of the function integrno */
	Xtmp1[1] = XY[indcons[i]][1];
	Xtmp1[2] = XY[indcons[i]][2];
	Xtmp2[1] = XY[indcons[i]+1][1];
	Xtmp2[2] = XY[indcons[i]+1][2];
	
	integrno(XG, Xtmp1, Xtmp2, &dt, sig1, sig2, alpha, &tmp);
	*res = *res + (poids * tmp);
    }
    freevec(Xtmp1);
    freevec(Xtmp2);
}



/* Main Function */
void kernelbb(double *grille, double *xgri, double *ygri, int *ncolgri,
	      int *nliggri, int *nloc, double *sig1, double *sig2, 
	      double *xlo, double *ylo, double *Tr, int *controlbox, 
	      int *nalpha)
{
    /* Declaration */
    int i, j, k, ncg, nlg, nlo, *indcons, ncons;
    double **gri, *xg, *yg, **XY, tmp, *alpha, *Xgr, *T, maxt,maxvh, res, vol;
    
    /* Memory Allocation */
    ncg = *ncolgri;
    nlg = *nliggri;
    nlo = *nloc;
    tmp = 0;
    
    taballoc(&gri,nlg, ncg);
    taballoc(&XY, nlo, 2);
    vecalloc(&xg, nlg);
    vecalloc(&T, nlo);
    vecalloc(&yg, ncg);
    vecalloc(&Xgr, 2);
    vecalloc(&alpha, *nalpha);
    vecintalloc(&indcons, nlo);
    
    /* R to C */    
    for (i=1; i<=nlo; i++) {
	XY[i][1] = xlo[i-1];
	XY[i][2] = ylo[i-1];
	T[i] = Tr[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xg[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	yg[i] = ygri[i-1];
    }
    
    /* Build the vector alpha */
    alpha[1] = 0;
    for (i = 2; i <= *nalpha; i++) {
	alpha[i] = ((double) i) / ((double) *nalpha);
    }
    
    /* Maximum dt and sigma for the normal distribution*/
    maxt = maxdt(T);
    maxvh = maxh(*sig1, *sig2, alpha, maxt);
    
	
    /* Loop on the grid */
    vol = 0;
    res = xg[2] - xg[1];
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    Xgr[1] = xg[i];
	    Xgr[2] = yg[j];
	    /*
	      ncons = consdanslabox(Xgr, XY, nlo, indcons, maxvh, *controlbox);
	    */
	    ncons = nlo-1;
	    for (k = 1; k < nlo; k++)
		indcons[k]=k;
	    udbbnoeud(Xgr, XY, T, sig1, sig2, alpha, &tmp, ncons, indcons);
	    gri[i][j] = tmp;
	    vol+=tmp;
	}
    }
    
    /* Standardization of the volume */

    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    gri[i][j] =  gri[i][j] / (vol * pow(res,2));
	}
    }
  
    /* C to R */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    grille[k] = gri[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(gri);
    freevec(xg);
    freevec(yg);
    freevec(T);
    freetab(XY);
    freevec(Xgr);
    freevec(alpha);
    freeintvec(indcons);
}


/* *********************************************************************
   Maximisation of the likelihood for the Brownian bridge
   
   *********************************************************************/
void CVL(double *xyr, double *Tr, 
	 int *nloc, double *Lr, double *sigma, int *nsig, double *sigma2)
{
    int i, j, k, nlo, ns, r;
    double **xy, *T,s2,ai,*mui,sigmai,res;
    
    nlo = *nloc;
    ns = *nsig;
    s2 = *sigma2;

    taballoc(&xy, nlo, 2);
    vecalloc(&T, nlo);
    vecalloc(&mui, 2);
    
    /* C to R */
    k = 0;
    for (i=1; i <= nlo; i++) {
	for (j = 1; j <= 2; j++) {
	    xy[i][j] = xyr[k];
	    k++;
	}
	T[i] = Tr[i-1];
    }
    
    /* Indices of odd locations */
    for (r = 1; r <= ns; r++) {
	Lr[r-1] = 0;
	k=1;
	for (i=1; i < nlo; i++) {
	    if (k == 2) {
		ai = (T[i] - T[i-1])/(T[i+1] - T[i-1]);
		
		mui[1] = xy[i-1][1] + ai * (xy[i+1][1] - xy[i-1][1]);
		mui[2] = xy[i-1][2] + ai * (xy[i+1][2] - xy[i-1][2]);
		
		sigmai = ((T[i+1]-T[i-1]) * ai * (1-ai) * sigma[r-1]) + (pow((1 - ai),2) * (*sigma2)) + (pow(ai,2) * (*sigma2));
		
		norm2d(xy[i][1], xy[i][2], 
		       mui[1], mui[2], sigmai, &res);
		Lr[r-1] = Lr[r-1] + log(res);
		k=1;
	    }
	    k++;
	}
    }

    /* Free memory */
    freetab(xy);
    freevec(T);
    freevec(mui);
}







/* *********************************************************************
 *                                                                     *
 *              Home range by Clustering (Kenward et al. 2001)         *
 *                                                                     *
 ***********************************************************************/


/* finds the cluster with the minimum average distance between the 3 points
   not assigned to a cluster */

void trouveclustmin(double **xy, int *clust, int *lo1, int *lo2,
		    int *lo3, double *dist)
{
    /* Declaration */
    int i, j, k, m, npas, nr, *indice;
    double **xy2, di1, di2, di3, ditmp;
    
    /* Memory allocation */
    nr = (int) xy[0][0];
    npas = 0;
    di1 = 0;
    di2 = 0;
    di3 = 0;
    ditmp = 0;
    
    /* Number of non assigned points */
    for (i = 1; i <= nr; i++) {
	if (clust[i] == 0) {
	    npas++;
	}
    }
    taballoc(&xy2, npas, 2);
    vecintalloc(&indice, npas);
    
    /* The non assigned points are stored in xy2 */
    k = 1;
    for (i = 1; i <= nr; i++) {
	if (clust[i] == 0) {
	    xy2[k][1] = xy[i][1];
	    xy2[k][2] = xy[i][2];
	    indice[k] = i;
	    k++;
	}
    }
    
    /* Computes the distane between the relocations */
    *dist = 0;
    m=0;
    for (i = 1; i <= (npas-2); i++) {
	for (j = (i+1); j <= (npas-1); j++) {
	    for (k = (j+1); k <= npas; k++) {
		di1 = sqrt((xy2[i][1] - xy2[j][1]) * (xy2[i][1] - xy2[j][1]) + 
			   (xy2[i][2] - xy2[j][2]) * (xy2[i][2] - xy2[j][2]) );
		di2 = sqrt((xy2[i][1] - xy2[k][1]) * (xy2[i][1] - xy2[k][1]) + 
			   (xy2[i][2] - xy2[k][2]) * (xy2[i][2] - xy2[k][2]));
		di3 = sqrt((xy2[k][1] - xy2[j][1]) * (xy2[k][1] - xy2[j][1]) + 
			   (xy2[k][2] - xy2[j][2]) * (xy2[k][2] - xy2[j][2]));
		/* average distance */
		ditmp = (di1 + di2 + di3) / 3;
		/* minimum distance */
		if ((m == 0) || (ditmp < *dist)) {
		    *dist = ditmp;
		    *lo1 = indice[i];
		    *lo2 = indice[j];
		    *lo3 = indice[k];
		}
		m = 1;
	    }
	}
    }
    /* free memory */
    freeintvec(indice);
    freetab(xy2);
}


/* For external call from within R */
void trouveclustminr(double *xyr, int *nr, int *clustr, int *lo1, int *lo2,
		     int *lo3, double *dist)
{
    /* Declaration */
    double **xy;
    int i, j, k, *clust;

    /* Memory allocation */
    taballoc(&xy, *nr, 2);
    vecintalloc(&clust, *nr);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= *nr; i++) {
	for (j = 1; j <= 2; j++) {
	    xy[i][j] = xyr[k];
	    k++;
	}
    }
    for (i = 1; i <= *nr; i++) {
	clust[i] = clustr[i-1];
    }
    
    /* main function */
    trouveclustmin(xy, clust, lo1, lo2, lo3, dist);
    
    /* Free memory */
    freetab(xy);
    freeintvec(clust);
}


/* Finds the distance between a cluster of points and the nearest point*/
void nndistclust(double **xy, double *xyp, double *dist)
{
    /* Declaration */
    int n, i, m;
    double di;

    m = 0;
    di =0;
    n = (int) xy[0][0];
    *dist = 0;
    
    /* finds the distance and the corresponding point */
    for (i = 1; i <= n; i++) {
	di = sqrt( (xy[i][1] - xyp[1]) * (xy[i][1] - xyp[1]) + 
		   (xy[i][2] - xyp[2]) * (xy[i][2] - xyp[2]) );
	if ( (di < *dist) || (m == 0) ) {
	    *dist = di;
	}
	m = 1;
    }
}


/* The function nndistclust is applied for all available clusters */
void parclust(double **xy, int *clust, int *noclust, 
	      int *noloc, double *dist)
{
    /* Declaration */
    int i, k, m, nr2, nr, nocl;
    double **xy2, *xyp, di, di2;
    
    /* Memory allocation */
    nocl = *noclust;
    nr = xy[0][0];
    nr2 = 0;
    
    /* The number of available clusters */
    for (i = 1; i <= nr; i++) {
	if (clust[i] == nocl) {
	    nr2++;
	}
    }

    taballoc(&xy2, nr2, 2);
    vecalloc(&xyp, 2);
    
    /* stores the non assigned points in xy2 */
    k = 1;
    for (i = 1; i <= nr; i++) {
	if (clust[i] == nocl) {
	    xy2[k][1] = xy[i][1];
	    xy2[k][2] = xy[i][2];
	    k++;
	}
    }
    
    /* Finds the minimum distance between a point and a cluster, 
       performed for all clusters */
    di = 0;
    di2 = 0;
    m = 0;
    *dist = 0;
    for (i = 1; i <= nr; i++) {
	if (clust[i] != nocl) {
	    xyp[1] = xy[i][1];
	    xyp[2] = xy[i][2];
	    nndistclust(xy2, xyp, &di);
	    if ( (di < *dist) || (m == 0) ) {
		*dist = di;
		*noloc = i;
	    }
	    m = 1;
	}
    }

    /* Free memory */
    freetab(xy2);
    freevec(xyp);
}


/* The function trouveminclust identifies the cluster for which the nearest 
   point is the closest */
void trouveminclust(double **xy, int *liclust, int *clust, 
		    int *noclust, int *noloc, double *dist)
{
    /* Declaration */
    int i, nr, nc, m, labclust, nolo;
    double di;
    
    nr = (int) xy[0][0];
    nc = 0;
    di = 0;
    labclust = 0;
    nolo = 0;
    
    /* Assigned clusters */
    for (i = 1; i <= nr; i++) {
	if (liclust[i] > 0) {
	    nc++;
	}
    }
    
    /* finds the minimum distance between a cluster and its nearest point 
       (the cluster name and the point ID are searched) */
    m = 0;
    *dist = 0;
    for (i = 1; i <= nc; i++) {
	labclust = liclust[i];
	parclust(xy, clust, &labclust, &nolo, &di);
	if ( (m == 0) || (di < *dist) ) {
	    *dist = di;
	    *noloc = nolo;
	    *noclust = labclust;
	}
	m = 1;
    }
}


/* What should be done: create a new cluster or add a relocation 
   to an existing one ? */

void choisnvclust(double **xy, int *liclust, int *clust, int *ordre)
{
    /* Declaration */
    int i, k, nr, noloat, cluat, nolo1, nolo2, nolo3, maxclust;
    int maxindiceclust, clu1, *liclub, nz;
    double dmoyclust, dminloc;
    
    /* Memory allocation */
    nz = 0;
    i = 0;
    k = 0;
    nr = (int) xy[0][0];
    maxclust = 0;
    maxindiceclust = 0;
    nolo1 = 0;
    nolo2 = 0;
    nolo3 = 0;
    noloat = 0;
    cluat = 0;
    clu1 = 0;
    vecintalloc(&liclub, nr);
    
    /* finds the max label for the cluster */
    for (i = 1; i <= nr; i++) {
	if (clust[i] != 0) {
	    if (clust[i] > maxclust) {
		maxclust = clust[i];
	    }
	    if (liclust[i] != 0) {
		maxindiceclust = i;
	    }
	}
    }
    
    /* Finds the min distance between 3 relocations */
    trouveminclust(xy, liclust, clust, &cluat, &noloat, &dminloc);
    
    /* Computes the average distance between the locs of the smaller cluster */
    /* First, one verifies that there is at least Three non assigned locs */
    dmoyclust = dminloc +1;
    for (i = 1; i <= nr; i++) {
	if (clust[i] == 0) {
	    nz++;
	}
    }
    if (nz > 3) {
	dmoyclust = 0;
	trouveclustmin(xy, clust, &nolo1, &nolo2, &nolo3, &dmoyclust);
    }
    
    /* First case: A new cluster independent from the others */
    if (dmoyclust < dminloc) {
	ordre[nolo1] = 1;
	ordre[nolo2] = 1;
	ordre[nolo3] = 1;
	
	clust[nolo1] = maxclust + 1;
	clust[nolo2] = maxclust + 1;
	clust[nolo3] = maxclust + 1;
	
	liclust[maxindiceclust+1] = maxclust + 1;
	
    } else {
	/* Second case: one loc is added to a cluster */
	
	/* Case 2.1: the loc does not belong to one cluster */
	if (clust[noloat] == 0) {
	    ordre[noloat] = 1;
	    clust[noloat] = cluat;
	} else {
	    
	    /* Case 2.2: the loc belong to one cluster: fusion */
	    clu1 = clust[noloat];
	    for (i = 1; i <= nr; i++) {
		if (clust[i] == clu1) {
		    clust[i] = cluat;
		    ordre[i] = 1;
		}
		if (liclust[i] == clu1) {
		    liclust[i] = 0;
		}
	    }
	    /* and cleaning of liclust */
	    k = 1;
	    for (i = 1; i <= nr; i++) {
		if (liclust[i] != 0) {
		    liclub[k] = liclust[i];
		    k++;
		}
	    }
	    for (i = 1; i <= nr; i++) {
		liclust[i] = liclub[i];
	    }
	}
    }
    freeintvec(liclub);
}


/* The main function for home range computation */

void clusterhr(double **xy, int *facso, int *nolocso, int *cluso)
{
    /* Declaration */
    int i, nr, lo1, lo2, lo3, *clust, len, con, *ordre, *liclust, courant;
    double di;

    /* Memory allocation */
    courant = 1;
    nr = (int) xy[0][0];
    vecintalloc(&clust, nr);
    vecintalloc(&ordre, nr);
    vecintalloc(&liclust, nr);
    lo1 = 0;
    lo2 = 0;
    lo3 = 0;
    di = 0;
    con = 1;
    len = 0;
    
    /* Begin: Search for the first cluster */
    trouveclustmin(xy, clust, &lo1, &lo2,
		   &lo3, &di);
    
    clust[lo1] = 1;
    clust[lo2] = 1;
    clust[lo3] = 1;
    liclust[1] = 1;
    len = 3;
    
    /* We store it in the output */
    cluso[1] = 1;
    cluso[2] = 1;
    cluso[3] = 1;
    nolocso[1] = lo1;
    nolocso[2] = lo2;
    nolocso[3] = lo3;
    facso[1] = 1;
    facso[2] = 1;
    facso[3] = 1;
    
    /* Then repeat until all relocations belong to the same cluster */
    while (con == 1) {
	courant++;
	
	for (i = 1; i <= nr; i++) {
	    ordre[i] = 0;
	}
	
	choisnvclust(xy, liclust, clust, ordre);
	
	for (i = 1; i <= nr; i++) {
	    if (ordre[i] != 0) {
		len++;
		cluso[len] = clust[i];
		nolocso[len] = i;
		facso[len] = courant;
	    }
	}
	
	con = 0;
	for (i = 2; i <= nr; i++) {
	    if (clust[i] != clust[1])
		con = 1;
	}
	if (con == 0) {
	    con = 0;
	}
    }
    
    /* Free memory */
    freeintvec(clust);
    freeintvec(ordre);
    freeintvec(liclust);
}



/* Finds the length of the output for the table containing the home range */

void longfacclust(double **xy, int *len2)
{
    /* Declaration */
    int i, nr, lo1, lo2, lo3, *clust, len, con, *ordre, *liclust, courant;
    double di;
    
    /* Memory allocation */
    courant = 1;
    nr = (int) xy[0][0];
    vecintalloc(&clust, nr);
    vecintalloc(&ordre, nr);
    vecintalloc(&liclust, nr);
    lo1 = 0;
    lo2 = 0;
    lo3 = 0;
    di = 0;
    con = 1;
    len = 0;
    
    /* Begin: search for the first cluster */
    trouveclustmin(xy, clust, &lo1, &lo2,
		   &lo3, &di);
    clust[lo1] = 1;
    clust[lo2] = 1;
    clust[lo3] = 1;
    liclust[1] = 1;
    len = 3;
    
    /* Counts the number of rows needed for the table, which will contain the results */
    while (con == 1) {
	courant++;
	
	for (i = 1; i <= nr; i++) {
	    ordre[i] = 0;
	}
	
	choisnvclust(xy, liclust, clust, ordre);
	
	for (i = 1; i <= nr; i++) {
	    if (ordre[i] != 0) {
		len++;
	    }
	}
	con = 0;
	for (i = 2; i <= nr; i++) {
	    if (clust[i] != clust[1])
		con = 1;
	}
	if (con == 0) {
	    con = 0;
	}
    }

    *len2 = len;

    /* Free memory */
    freeintvec(clust);
    freeintvec(ordre);
    freeintvec(liclust);
}


/* For external call from within R */
void longfacclustr(double *xyr, int *nr, int *len2)
{
    /* Declaration */
    double **xy;
    int i, j, k;
    
    /* Memory allocation */
    taballoc(&xy, *nr, 2);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= *nr; i++) {
	for (j = 1; j <= 2; j++) {
	    xy[i][j] = xyr[k];
	    k++;
	}    
    }
    
    /* Main function */
    longfacclust(xy, len2);

    /* Free memory */
    freetab(xy);
}



/* For external call of clusterhrr from within R */

void clusterhrr(double *xyr, int *nr, int *facsor, 
		int *nolocsor, int *clusor, int *len)
{
    /* Declaration */
    double **xy;
    int i, j, k, *facso, *nolocso, *cluso;
    
    /* Memory allocation */
    taballoc(&xy, *nr, 2);
    vecintalloc(&facso, *len);
    vecintalloc(&nolocso, *len);
    vecintalloc(&cluso, *len);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= *nr; i++) {
	for (j = 1; j <= 2; j++) {
	    xy[i][j] = xyr[k];
	    k++;
	}
    }
    
    /* Main function */
    clusterhr(xy, facso, nolocso, cluso);
    
    /* C to R */
    for (i = 1; i <= *len; i++) {
	facsor[i-1] = facso[i];
	nolocsor[i-1] = nolocso[i];
	clusor[i-1] = cluso[i];
    }
    
    /* Free memory */
    freetab(xy);
    freeintvec(facso);
    freeintvec(nolocso);
    freeintvec(cluso);
}







/* *********************************************************************
 *                                                                     *
 *           Kernel in time and space (Keating and Cherry, 2008)       *
 *                                                                     *
 ***********************************************************************/



void kcprcirc(double **xyd, double *h, double *x, double t, 
	      double *val)
{
    int i, j, k, nlo;
    double tmp, tmp2, vi, som;
    
    nlo = xyd[0][0];
    som=0;
    
    for (i=1; i<=nlo; i++) {
	
	tmp2 = 1;
	
	/* spatial coordinates */
	for (j=1; j<=2; j++) {
	    vi = (x[j] - xyd[i][j]) / h[j];
	    if (fabs(vi) < 1.0) {
		tmp = (((double) 15)/((double) 16)) * ((1- (vi * vi)) * (1- (vi * vi)));
		tmp2 = tmp2 * tmp;
	    } else {
		tmp2 = 0.0;
	    }
	}
	
	/* time */
	vi = (t - xyd[i][3]);
	tmp2 = tmp2 * (1-h[3]) * (1 - h[3] * h[3]) / 
	    (2 * 3.14153265359 * ( 1 +  (h[3] * h[3]) - (2 * h[3] * cos(vi) )));
	som = som + tmp2;
    }
    *val = ( 1/( ((double) nlo) * h[1] * h[2] * h[3])) * som;
}

void kcprlin(double **xyd, double *h, double *x, double t, 
	      double *val)
{
    int i, j, k, nlo;
    double tmp, tmp2, vi, som;
    
    nlo = xyd[0][0];
    
    som = 0;
    
    for (i=1; i<=nlo; i++) {
	
	tmp2 = 1;
	
	/* spatial coordinates and time */
	for (j=1; j<=2; j++) {
	    vi = (x[j] - xyd[i][j]) / h[j];
	    if (fabs(vi) < 1.0) {
		tmp = (((double) 15)/((double) 16)) * ((1- (vi * vi)) * (1- (vi * vi)));
		tmp2 = tmp2 * tmp;
	    } else {
		tmp2 = 0.0;
	    }
	}
	
	/* time */
	vi = (t - xyd[i][3]) / h[3];
	if (fabs(vi) < 1.0) {
	    tmp2 = tmp2 * (((double) 15)/((double) 16)) * ((1- (vi * vi)) * (1- (vi * vi)));
	} else {
	    tmp2 = 0.0;
	}
	som = som + tmp2;
    }
    *val = ( 1/( ((double) nlo) * h[1] * h[2] * h[3])) * som;
}


void kernelkcr(double *xydr, double *tcalcr, int *nlr, double *gridr,
	       double *xgri, double *ygri, int *nliggri, int *ncolgri, 
	       double *hr, int *circularr)
{
    /* Declaration of local variables */
    int i, j, k, ncg, nlg, nlo, circular;
    double **gri, **xyd, *xx, *xg, *yg, tmp, *h, tca;
    
    /* Memory Allocation */
    ncg = *ncolgri;
    nlg = *nliggri;
    nlo = *nlr;
    tca = *tcalcr;
    circular = *circularr;
    
    taballoc(&gri, nlg, ncg);
    taballoc(&xyd, nlo, 3);
    vecalloc(&xg, nlg);
    vecalloc(&yg, ncg);
    vecalloc(&h, 3);
    vecalloc(&xx, 2);
    
    /* R objects -> C objects */

    for (i=1; i<=3; i++) {
	h[i] = hr[i-1];
    }
    
    for (i=1; i<=nlg; i++) {
	xg[i] = xgri[i-1];
    }
    
    for (i=1; i<=ncg; i++) {
	yg[i] = ygri[i-1];
    }
    
    k = 0;
    for (i=1; i<=nlo; i++) {
	for (j=1; j<=3; j++) {
	    xyd[i][j] = xydr[k];
	    k++;
	}
    }
    
    /* loop on the grid */
    if (circular==1) {
	for (i=1; i<=nlg; i++) {
	    for (j=1; j<=ncg; j++) {
		xx[1] = xg[i];
		xx[2] = yg[j];
		kcprcirc(xyd, h, xx, tca, &tmp);
		gri[i][j] = tmp;
	    }
	}
    } else {
	for (i=1; i<=nlg; i++) {
	    for (j=1; j<=ncg; j++) {
		xx[1] = xg[i];
		xx[2] = yg[j];
		kcprlin(xyd, h, xx, tca, &tmp);
		gri[i][j] = tmp;
	    }
	}
    }
    
    /* C objects -> R objects */
    k = 0;
    for (i=1; i<=nlg; i++) {
	for (j=1; j<=ncg; j++) {
	    gridr[k] = gri[i][j];
	    k++;
	}
    }

    /* Memory Free */
    freetab(gri);
    freetab(xyd);
    freevec(xg);
    freevec(yg);
    freevec(h);
    freevec(xx);
}











/* *********************************************************************
 *                                                                     *
 *                      find local maxima/minima on a map              *
 *                                                                     *
 ***********************************************************************/


void findmaxgrid(double *grille, int *nlig, int *ncol)
{
    /* declaration */
    int i,j,k,nl,nc,sto;
    double **x, **grille2, r1,r2,r3,r4,r5,r6,r7,r8;
    
    nl = *nlig;
    nc = *ncol;
    sto=0;
    
    /* Memory alloocation */
    taballoc(&x,nl,nc);
    taballoc(&grille2,nl,nc);
    
    /* R to C */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    x[i][j]=grille[k];
	    k++;
	}
    }
    
    /* Loop */
    for (i = 2; i <= (nl-1); i++) {
	for (j=2; j<= (nc-1); j++) {
	    
	    r1=x[i-1][j-1] - x[i][j];
	    r2=x[i][j-1] - x[i][j];
	    r3=x[i+1][j-1] - x[i][j];
	    r4=x[i+1][j] - x[i][j];
	    r5=x[i+1][j+1] - x[i][j];
	    r6=x[i][j+1] - x[i][j];
	    r7=x[i-1][j+1] - x[i][j];
	    r8=x[i-1][j] - x[i][j];
	    
	    sto=1;
	    
	    if (r1 < -0.000000000001) {
		if (r2 < -0.000000000001) {
		    if (r3 < -0.000000000001) {
			if (r4 < -0.000000000001) {
			    if (r5 < -0.000000000001) {
				if (r6 < -0.000000000001) {
				    if (r7 < -0.000000000001) {
					if (r8 < -0.000000000001) {
					    sto=0;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	    
	    if (sto==0)
		grille2[i][j] = 1;
	}
    }
        
    
    /* C to R */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    grille[k]=grille2[i][j];
	    k++;
	}
    }
    
    /* Memory */
    freetab(x);
    freetab(grille2);
    
}




/* ***********************************************************************
 *                                                                       *
 *                          BRB                                          *
 *                                                                       *
 * ********************************************************************* */



int HBT(double xt, double yt, SEXP hab, SEXP nrow, SEXP cs, double xll2, 
	double yll2)
{
    int hh, nl, nc;
    nl = (int) ftrunc(((xt - xll2)/REAL(cs)[0]) + REAL(cs)[0]*0.000001); 
    nc = (int) ftrunc(((yt - yll2)/REAL(cs)[0]) + REAL(cs)[0]*0.000001);
    hh = INTEGER(hab)[ nl + (nc * (INTEGER(nrow)[0])) ];
    return(hh);
}


int HBTl(SEXP xl, SEXP yl, SEXP PAtmp, SEXP hab, SEXP nrow, SEXP cs, double xll2, 
	double yll2, int k, int i)
{
    double t1, xt, yt;
    int n, j, hh, th;
    SEXP habp;
    
    PROTECT(habp = allocVector(INTSXP, k+1));
    

    t1 = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
    n = (int) round(t1);
    if (n <1)
	n = 1;
    for (j = 0; j < k+1; j++) {
	INTEGER(habp)[j] = 0;
    }
    
    /* identify the habitat at each step */
    for (j = 0; j <= n; j++) {
	xt = REAL(xl)[i] + (((double) j)/((double) n)) * (REAL(xl)[i+1] - REAL(xl)[i]);
	yt = REAL(yl)[i] + (((double) j)/((double) n)) * (REAL(yl)[i+1] - REAL(yl)[i]);
	hh = HBT(xt, yt, hab, nrow, cs, xll2, yll2);
	if (hh != NA_INTEGER) {
	    INTEGER(habp)[hh]++;
	} else {
	    INTEGER(habp)[k]++;
	}
    }
    hh=0;
    for (j = 0; j < k+1; j ++) {
	if (INTEGER(habp)[j] == (n+1)) {
	    th = j;
	    hh++;
	}
    }
    
    if (hh > 0) {
	UNPROTECT(1);
	return(th);
    } else {
	UNPROTECT(1);
	return(NA_INTEGER);
    }
}


/* df contient x,y,date en posix */
SEXP fillsegments(SEXP df, SEXP Tmaxr, SEXP taur, SEXP hminr, SEXP D, SEXP Lminr, 
		  SEXP b, SEXP hab, SEXP xll, SEXP yll, SEXP cs, SEXP nrowc, SEXP PA)
{
    int nrow, nnr, ni, i, m, k, nh, h, lp;
    double dt, dta, Tmax, tau, h2min, Lmin, dist, hmin, xll2, yll2, hm;
    SEXP x, y, date, resux, resuy, resuh, dfso, hmax, h2max, PAtmp, PA2;
    
    /* Le nombre de lignes de ce data.frame */
    nrow = length(VECTOR_ELT(df,0));
    nnr = 0;
    nh = length(D);
    if (nh > 1) {
	xll2 = REAL(xll)[0] - REAL(cs)[0]/2.0;
	yll2 = REAL(yll)[0] - REAL(cs)[0]/2.0;
    }
    
    Tmax = REAL(Tmaxr)[0];
    hmin = REAL(hminr)[0];
    h2min = R_pow(hmin, 2.0);
    
    PROTECT(hmax = allocVector(REALSXP, nh+1));
    PROTECT(h2max = allocVector(REALSXP, nh+1));

    REAL(hmax)[0] = sqrt(((1.0 - REAL(b)[0]) * h2min) + (REAL(D)[0] * Tmax / 2.0));
    hm = REAL(hmax)[0];
    if (nh > 1) {
	for (i = 1; i < nh; i++) {
	    REAL(hmax)[i] = sqrt(((1.0 - REAL(b)[0]) * h2min) + (REAL(D)[i] * Tmax / 2.0));
	    if (hm < REAL(hmax)[i])
		hm = REAL(hmax)[i];
	}
	REAL(hmax)[nh] = hm;
    }
    Lmin = REAL(Lminr)[0];
    tau = REAL(taur)[0];
    REAL(h2max)[0] = R_pow(REAL(hmax)[0], 2.0);
    if (nh > 1) {
	for (i = 1; i <= nh; i++) {
	    REAL(h2max)[i] = R_pow(REAL(hmax)[i], 2.0);
	}
    }

    /* Get the coordinates and date */
    PROTECT(x = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(df,1), REALSXP));
    PROTECT(date = coerceVector(VECTOR_ELT(df,2), REALSXP));
    lp = length(PA);
    PROTECT(PAtmp = allocVector(REALSXP, nrow));
    PROTECT(PA2 = coerceVector(PA, REALSXP));
    if (lp > 1) {
	REAL(PAtmp)[0] = 0.0;
	for (i = 1; i < nrow; i++) {
	    REAL(PAtmp)[i] = REAL(PAtmp)[i-1] + (REAL(PA2)[i-1] * (REAL(date)[i] - REAL(date)[i-1]));
	}	
    } else {
	for (i = 0; i < nrow; i++) {
	    REAL(PAtmp)[i] = REAL(date)[i];
	}
    }
    
    
    /* for each segment, calculates the number of points to add */
    for (i = 0; i < (nrow-1); i++) {
	dt = REAL(date)[i+1] - REAL(date)[i];
	dta = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
	dist = pythag(REAL(x)[i+1] - REAL(x)[i], REAL(y)[i+1] - REAL(y)[i]);
	if ((dt < Tmax)&&(dist > Lmin)) {
	    nnr = nnr + (int) round(dta/tau);
	}
    }
    nnr++;
    
    /* prepares the vector of output */
    PROTECT(resux = allocVector(REALSXP, nnr));
    PROTECT(resuy = allocVector(REALSXP, nnr));
    PROTECT(resuh = allocVector(REALSXP, nnr));
    
    /* and finds the coordinates of the points */
    k=0;
    for (i = 0; i < (nrow-1); i++) {

	dt = REAL(date)[i+1] - REAL(date)[i];
	dta = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
	dist = pythag(REAL(x)[i+1] - REAL(x)[i], REAL(y)[i+1] - REAL(y)[i]);

	if ((dt < Tmax)&&(dist>Lmin)&&(dta>0.0000001)) {
	    ni = (int) round(dta/tau);
	    for (m = 0; m < ni; m++) {
		REAL(resux)[k] = REAL(x)[i]+ 
			((double) m) * (REAL(x)[i+1] - REAL(x)[i])/((double) ni);
		REAL(resuy)[k] = REAL(y)[i]+ 
		    ((double) m) * (REAL(y)[i+1] - REAL(y)[i])/((double) ni);
		if (nh < 2) {
		    REAL(resuh)[k] = sqrt(h2min + 4.0*(((double) m)/((double) ni))*
					  (1 - ((double) m)/((double) ni))*(REAL(h2max)[0] - h2min)*dta/Tmax);
		} else {
		    h = HBT(REAL(resux)[k], REAL(resuy)[k], hab, nrowc, cs, xll2, 
			    yll2);
		    if (h == NA_INTEGER) {
			REAL(resuh)[k] = sqrt(h2min + 4.0*(((double) m)/((double) ni))*
					      (1 - ((double) m)/((double) ni))*
					      (REAL(h2max)[nh] - h2min)*dta/Tmax);
			
		    } else {
			REAL(resuh)[k] = sqrt(h2min + 4.0*(((double) m)/((double) ni))*
					      (1 - ((double) m)/((double) ni))*
					      (REAL(h2max)[h] - h2min)*dta/Tmax);
			
		    }
		}
		k++;
	    }
	}
    }
    
    REAL(resux)[k] = REAL(x)[nrow-1];
    REAL(resuy)[k] = REAL(y)[nrow-1];
    REAL(resuh)[k] = sqrt(h2min);
    
    PROTECT(dfso = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(dfso, 0, resux);
    SET_VECTOR_ELT(dfso, 1, resuy);
    SET_VECTOR_ELT(dfso, 2, resuh);

    UNPROTECT(11);

    return(dfso);
}


/* On calcule maintenant, sur la base d'une grille passée, l'estimation kernel */
SEXP mkde(SEXP xyh, SEXP grid)
{
    
    int n, nl, i, j;
    SEXP x, y, h, dens, xg, yg, gridso;
    double xmin, ymin, xmax, ymax, hmax, dist;
    
    /* on ajuste alors le noyau */
    n = length(VECTOR_ELT(grid,0));
    nl = length(VECTOR_ELT(xyh,0));


    PROTECT(x = coerceVector(VECTOR_ELT(xyh,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyh,1), REALSXP));
    PROTECT(h = coerceVector(VECTOR_ELT(xyh,2), REALSXP));
    PROTECT(xg = coerceVector(VECTOR_ELT(grid,0), REALSXP));
    PROTECT(yg = coerceVector(VECTOR_ELT(grid,1), REALSXP));
    PROTECT(dens = allocVector(REALSXP, n));
    
    
    xmin = REAL(x)[0];
    ymin = REAL(y)[0];
    xmax = REAL(x)[0];
    ymax = REAL(y)[0];
    hmax = REAL(h)[0];
    for (j = 1; j < nl; j++) {
	if (REAL(x)[j] < xmin)
	    xmin = REAL(x)[j];
	if (REAL(x)[j] > xmax)
	    xmax = REAL(x)[j];
	if (REAL(y)[j] < ymin)
	    ymin = REAL(y)[j];
	if (REAL(y)[j] > ymax)
	    ymax = REAL(y)[j];
	if (REAL(h)[j] > hmax)
	    hmax = REAL(h)[j];
    }
    hmax = hmax * 3.0;
    
    for (i = 0; i < n; i++) {
	R_CheckUserInterrupt();
	REAL(dens)[i] = 0.0;
	if ((xmin - REAL(xg)[i] < hmax)&&
	    (ymin - REAL(yg)[i] < hmax)&&
	    (REAL(xg)[i] - xmax < hmax)&&
	    (REAL(yg)[i] - ymax < hmax)) {

	    for (j = 0; j < nl; j++) {
		dist= pythag(REAL(x)[j] -REAL(xg)[i], REAL(y)[j] -REAL(yg)[i]);
		if (dist < 3.0*REAL(h)[j]) {
		    REAL(dens)[i] = REAL(dens)[i] + exp(-(R_pow(dist,2.0))/
							(2.0 * R_pow(REAL(h)[j], 2.0))) / 
			R_pow(REAL(h)[j], 2.0);
		}
	    }
	    REAL(dens)[i] = (1.0/(2.0 * M_PI * ((double) nl)))*(REAL(dens)[i]);
	}
    }

    PROTECT(gridso = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(gridso, 0, xg);
    SET_VECTOR_ELT(gridso, 1, yg);
    SET_VECTOR_ELT(gridso, 2, dens);

    UNPROTECT(7);
    return(gridso);
}




/* bis */
SEXP mkdeb(SEXP xyh, SEXP xll, SEXP yll, SEXP cs, SEXP nrow, SEXP ncol)
{
    
    int nl, nr, nc, nro, nco, i, j, l, c, hmaxdis;
    SEXP x, y, h, dens, xg, yg, gridso;
    double xlo, ylo, hmax, dist, xll2, yll2;
    
    /* on ajuste alors le noyau */
    nl = length(VECTOR_ELT(xyh,0));


    PROTECT(x = coerceVector(VECTOR_ELT(xyh,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyh,1), REALSXP));
    PROTECT(h = coerceVector(VECTOR_ELT(xyh,2), REALSXP));
    PROTECT(xg = allocVector(REALSXP, INTEGER(nrow)[0]*INTEGER(ncol)[0]));
    PROTECT(yg = allocVector(REALSXP, INTEGER(nrow)[0]*INTEGER(ncol)[0]));
    PROTECT(dens = allocVector(REALSXP, INTEGER(nrow)[0]*INTEGER(ncol)[0]));
    nro = INTEGER(nrow)[0];
    nco = INTEGER(ncol)[0];

    for (j = 0; j < nco; j++) {
	for (i = 0; i < nro; i++) {
	    REAL(xg)[ i + (j * (INTEGER(nrow)[0])) ] = REAL(xll)[0] + ((double) i)*REAL(cs)[0];
	    REAL(yg)[ i + (j * (INTEGER(nrow)[0])) ] = REAL(yll)[0] + ((double) j)*REAL(cs)[0];
	}
    }
	    
    for (i = 0; i < INTEGER(nrow)[0]*INTEGER(ncol)[0]; i++) {
	REAL(dens)[i] = 0.0;
    }
    
    hmax = REAL(h)[0];
    for (j = 1; j < nl; j++) {
	if (REAL(h)[j] > hmax)
	    hmax = REAL(h)[j];
    }
    hmax = hmax * 3.0;
    xll2 = REAL(xll)[0] - REAL(cs)[0]/2.0;
    yll2 = REAL(yll)[0] - REAL(cs)[0]/2.0;
    hmaxdis = (int) round(hmax / REAL(cs)[0]);
    
    for (j = 0; j < nl; j++) {
	R_CheckUserInterrupt();
	xlo = REAL(x)[j];
	ylo = REAL(y)[j];
	nr = (int) ftrunc(((xlo - xll2)/REAL(cs)[0]) + REAL(cs)[0]*0.000001); 
	nc = (int) ftrunc(((ylo - yll2)/REAL(cs)[0]) + REAL(cs)[0]*0.000001);
	for (l = (nr-hmaxdis-1); l <(nr+hmaxdis+1); l++) {
	    for (c = (nc-hmaxdis-1); c <(nc+hmaxdis+1); c++) {
		if ((l<nro)&&(l>0)) {
		    if ((c<nco)&&(c>0)) {
			dist= pythag(xlo -REAL(xg)[l + (c * (INTEGER(nrow)[0]))], 
				     ylo -REAL(yg)[l + (c * (INTEGER(nrow)[0]))]);
			REAL(dens)[ l + (c * (INTEGER(nrow)[0])) ] =
			    REAL(dens)[ l + (c * (INTEGER(nrow)[0])) ] +
			    exp(-(R_pow(dist,2.0))/
				(2.0 * R_pow(REAL(h)[j], 2.0))) / 
			    R_pow(REAL(h)[j], 2.0)/(2.0 * M_PI * ((double) nl));
		    }
		}
	    }
	}
    }
	

    PROTECT(gridso = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(gridso, 0, xg);
    SET_VECTOR_ELT(gridso, 1, yg);
    SET_VECTOR_ELT(gridso, 2, dens);

    UNPROTECT(7);
    return(gridso);
}



/* */
SEXP CalculD(SEXP tra, SEXP Tmaxr, SEXP Lmin, SEXP PA)
{
    double Tmax, t1, t2, xt, yt, delta2, D, l1, l2;
    int n, i, Nc, lp;
    SEXP x, y, date, Ds, PA2, PAtmp;
    
    Tmax = REAL(Tmaxr)[0];
    n = length(VECTOR_ELT(tra,0));
    xt = 0.0;
    yt = 0.0;

    /* Get the coordinates and date */
    PROTECT(x = coerceVector(VECTOR_ELT(tra,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(tra,1), REALSXP));
    PROTECT(date = coerceVector(VECTOR_ELT(tra,2), REALSXP));
    lp = length(PA);
    PROTECT(PAtmp = allocVector(REALSXP, n));
    PROTECT(PA2 = coerceVector(PA, REALSXP));
    if (lp > 1) {
	REAL(PAtmp)[0] = 0.0;
	for (i = 1; i < n; i++) {
	    REAL(PAtmp)[i] = REAL(PAtmp)[i-1] + (REAL(PA2)[i-1] * 
						 (REAL(date)[i] - REAL(date)[i-1]));
	}	
    } else {
	for (i = 0; i < n; i++) {
	    REAL(PAtmp)[i] = REAL(date)[i];
	}
    }


    D = 0.0;
    Nc = 0;
    for (i = 0; i < (n-2); i++) {	
	t1 = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
	t2 = REAL(PAtmp)[i+2] - REAL(PAtmp)[i+1];
	l1 = pythag(REAL(x)[i+1] - REAL(x)[i], 
		    REAL(y)[i+1] - REAL(y)[i]);
	l2 = pythag(REAL(x)[i+2] - REAL(x)[i+1],
		    REAL(y)[i+2] - REAL(y)[i+1]);
	if ((REAL(date)[i+2]-REAL(date)[i]) < Tmax) {
	    if (t1 > 0.0000000001) {
		if (t2 > 0.0000000001) {
		    if (t1 < 2.0 * t2) {
			if (t1 > t2/2.0) {
			    if (l1 <  2.0 * l2) {
				if (l1 >  l2/2.0) {
				    if (l1 > REAL(Lmin)[0]) {
					if (l2 > REAL(Lmin)[0]) {
					    xt = REAL(x)[i] + (REAL(x)[i+2] - REAL(x)[i])*(t1/(t1+t2));
					    yt = REAL(y)[i] + (REAL(y)[i+2] - REAL(y)[i])*(t1/(t1+t2));
					    delta2 = R_pow((xt - REAL(x)[i+1]), 2.0) + 
						R_pow((yt - REAL(y)[i+1]), 2.0);
					    D = D + (delta2*((1.0/t1) + (1.0/t2)));
					    Nc++;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    D = D/(4.0 * ((double) Nc));
    PROTECT(Ds = allocVector(REALSXP, 2));
    REAL(Ds)[0] = (double) Nc;
    REAL(Ds)[1] = D;
    UNPROTECT(6);

    return(Ds);    
}








/* Differences with the program of Simon Benhamou: trunc on a integer i stored as a double
   can return i or i-1
 */
SEXP calculDparhab(SEXP df, SEXP hab, SEXP xll, SEXP yll, SEXP cs, SEXP nrow,
		   SEXP Lmin, SEXP nombrehab, SEXP PA, SEXP Tmax)
{
    SEXP xl, yl, tem, typpas, habp, Nc, Dh, dfso, PAtmp, PA2;
    int i, k, nlocs, lp;
    double t1, t2, l1, l2, xt, yt, delta2, tau, xll2, yll2;

    k = INTEGER(nombrehab)[0];
    nlocs = length(VECTOR_ELT(df,0));
    tau = REAL(cs)[0]/100.0;

    PROTECT(xl = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(yl = coerceVector(VECTOR_ELT(df,1), REALSXP));
    PROTECT(tem = coerceVector(VECTOR_ELT(df,2), REALSXP));
    PROTECT(typpas = allocVector(INTSXP, nlocs-1));
    PROTECT(habp = allocVector(INTSXP, k+1));
    lp = length(PA);
    PROTECT(PAtmp = allocVector(REALSXP, nlocs));
    PROTECT(PA2 = coerceVector(PA, REALSXP));

    /* From center to the corner of lower left pixel */
    xll2 = REAL(xll)[0] - (REAL(cs)[0]/2.0);
    yll2 = REAL(yll)[0] - (REAL(cs)[0]/2.0);

    /*  Take into account the proportion of activity time */
    if (lp > 1) {
	REAL(PAtmp)[0] = 0.0;
	for (i = 1; i < nlocs; i++) {
	    REAL(PAtmp)[i] = REAL(PAtmp)[i-1] + (REAL(PA2)[i-1] * (REAL(tem)[i] - REAL(tem)[i-1]));
	}	
    } else {
	for (i = 0; i < nlocs; i++) {
	    REAL(PAtmp)[i] = REAL(tem)[i];
	}
    }
    
    
    /* for each step */
    for (i = 0; i < nlocs-1; i++) {
	INTEGER(typpas)[i] = HBTl(xl, yl, PAtmp, hab, nrow, cs, xll2, yll2, k, i);
    }
    
    
    /* calculates the D coefficient */
    PROTECT(Nc = allocVector(INTSXP, k));
    PROTECT(Dh = allocVector(REALSXP, k));

    for (i = 0; i < k; i++) {
	REAL(Dh)[i] = 0.0;
	INTEGER(Nc)[i] = 0;
    }
    
    for (i = 0; i < (nlocs-2); i++) { 
	if ((INTEGER(typpas)[i+1] != NA_INTEGER)&&
	    (INTEGER(typpas)[i+1] == INTEGER(typpas)[i])) {
	    l2 = pythag(REAL(xl)[i+2] - REAL(xl)[i+1], REAL(yl)[i+2] - REAL(yl)[i+1]);
	    l1 = pythag(REAL(xl)[i+1] - REAL(xl)[i], REAL(yl)[i+1] - REAL(yl)[i]);
	    t2 = REAL(PAtmp)[i+2] - REAL(PAtmp)[i+1];
	    t1 = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
	    if (t1 > 0.0000000001) {
		if (t2 > 0.0000000001) {
		    if ((REAL(tem)[i+2]-REAL(tem)[i]) < REAL(Tmax)[0]) {
			if (t1 < 2.0 * t2) {
			    if (t1 > t2/2.0) {
				if (l1 <  2.0 * l2) {
				    if (l1 >  l2/2.0) {
					if (l1 > REAL(Lmin)[0]) {
					    if (l2 > REAL(Lmin)[0]) {
						xt = REAL(xl)[i] + (REAL(xl)[i+2] - 
								    REAL(xl)[i])*(t1/(t1+t2));
						yt = REAL(yl)[i] + (REAL(yl)[i+2] - 
								    REAL(yl)[i])*(t1/(t1+t2));
						delta2 = R_pow((xt - REAL(xl)[i+1]), 2.0) + 
						    R_pow((yt - REAL(yl)[i+1]), 2.0);
						REAL(Dh)[INTEGER(typpas)[i]] = 
						    REAL(Dh)[INTEGER(typpas)[i]] +
						    (delta2*((1.0/t1) + (1.0/t2)));
						INTEGER(Nc)[INTEGER(typpas)[i]]++;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    
    for (i = 0; i < k; i++) {
	REAL(Dh)[i] = 
	    REAL(Dh)[i] / 
	    (4.0 * ((double) INTEGER(Nc)[i]));	    
    }

    PROTECT(dfso = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dfso, 0, Nc);
    SET_VECTOR_ELT(dfso, 1, Dh);


    UNPROTECT(10);
    return(dfso);    
}

/* **********************************************************


********************************************************** */

/* Calcul D par maximum de vraisemblance */

double calcv(SEXP xl, SEXP yl, SEXP da, double D, SEXP pc)
{
    int n, i,k;
    double vrais, d, T, t;
    
    n = length(xl);
    vrais = 0.0;
    k = 0;
    for (i = 1; i < n-1; i++) {
	if (k == 0) {
	    if (INTEGER(pc)[i] == 1) {
		T = REAL(da)[i+1]-REAL(da)[i-1];
		t = REAL(da)[i]-REAL(da)[i-1];
		d = pythag(REAL(xl)[i] - REAL(xl)[i-1] - (t/T)*(REAL(xl)[i+1] - REAL(xl)[i-1]),
			   REAL(yl)[i] - REAL(yl)[i-1] - (t/T)*(REAL(yl)[i+1] - REAL(yl)[i-1]));
		vrais = vrais + log(T/(4.0*M_PI*D*t*(T-t))) - R_pow(d,2.0)/(4.0*D*t*(T-t)/T);
		k++;
	    } 
	} else {
	    k=0;
	}
    }
    return(vrais);
}


double compteN(SEXP xl, SEXP pc)
{
    int n, i,k, cons;
    
    n = length(xl);
    k = 1;
    cons=0;
    for (i = 1; i < n-1; i++) {
	if (k == 0) {
	    if (INTEGER(pc)[i] == 1) {
		cons++;
		k++;
	    } 
	} else {
	    k=0;
	}
    }
    return((double) cons);
}


SEXP Dmv(SEXP df, SEXP Dr, SEXP pcr)
{
    SEXP xl, yl, da, D, sor, pc;
    double fx1, fx2, fx3, fx4, x1, x2, x3, x4, phi, z, tmp;
    int n, ndr, conv;
    
    ndr = length(Dr);
    PROTECT(D = coerceVector(Dr, REALSXP));
    PROTECT(pc = coerceVector(pcr, INTSXP));
    PROTECT(xl = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(yl = coerceVector(VECTOR_ELT(df,1), REALSXP));
    PROTECT(da = coerceVector(VECTOR_ELT(df,2), REALSXP));
    PROTECT(sor = allocVector(REALSXP, 2));
    
    /* Golden section search */
    x1 = REAL(D)[0];
    x3 = REAL(D)[1];
    phi = (-1.0 + sqrt(5.0))/2.0;
    fx1 = calcv(xl, yl, da, x1, pc);
    fx3 = calcv(xl, yl, da, x3, pc);
    
    conv = 0;
    while (!conv) {
	x2 = x3 - phi*(x3-x1);
	x4 = x1 + phi*(x3-x1);
	fx2 = calcv(xl, yl, da, x2, pc);
	fx4 = calcv(xl, yl, da, x4, pc);

	if (fx2 < fx4) {
	    x1 = x2;
	    fx1 = fx2;
	} else {
	    x3 = x4;
	}
	if (fabs(x3-x1)<0.00000001) {
	    conv = 1;
	    x4 = (x3+x1)/2.0;
	}
    }
    REAL(sor)[0] = compteN(xl, pc);
    REAL(sor)[1] = x4;


    UNPROTECT(6);
    return(sor);
}
