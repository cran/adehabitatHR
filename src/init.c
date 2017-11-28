#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void calcvolume(void *, void *, void *, void *);
extern void clusterhrr(void *, void *, void *, void *, void *, void *);
extern void CVL(void *, void *, void *, void *, void *, void *, void *);
extern void CVmise(void *, void *, void *, void *, void *, void *);
extern void findmaxgrid(void *, void *, void *);
extern void kernelbbc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kernelhr(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kernelkcr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kernepan(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void longfacclustr(void *, void *, void *);

/* .Call calls */
extern SEXP CalculD(SEXP, SEXP, SEXP, SEXP);
extern SEXP calculDparhab(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Dmv(SEXP, SEXP, SEXP);
extern SEXP fillsegments(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP HRresidtime(SEXP, SEXP, SEXP);
extern SEXP mkdeb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nvisits(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"calcvolume",    (DL_FUNC) &calcvolume,     4},
    {"clusterhrr",    (DL_FUNC) &clusterhrr,     6},
    {"CVL",           (DL_FUNC) &CVL,            7},
    {"CVmise",        (DL_FUNC) &CVmise,         6},
    {"findmaxgrid",   (DL_FUNC) &findmaxgrid,    3},
    {"kernelbbc",      (DL_FUNC) &kernelbbc,      12},
    {"kernelhr",      (DL_FUNC) &kernelhr,       9},
    {"kernelkcr",     (DL_FUNC) &kernelkcr,     10},
    {"kernepan",      (DL_FUNC) &kernepan,       9},
    {"longfacclustr", (DL_FUNC) &longfacclustr,  3},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"CalculD",       (DL_FUNC) &CalculD,        4},
    {"calculDparhab", (DL_FUNC) &calculDparhab, 10},
    {"Dmv",           (DL_FUNC) &Dmv,            3},
    {"fillsegments",  (DL_FUNC) &fillsegments,  14},
    {"HRresidtime",   (DL_FUNC) &HRresidtime,    3},
    {"mkdeb",         (DL_FUNC) &mkdeb,          8},
    {"nvisits",       (DL_FUNC) &nvisits,        3},
    {NULL, NULL, 0}
};

void R_init_adehabitatHR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
