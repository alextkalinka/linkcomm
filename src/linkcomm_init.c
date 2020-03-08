#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void cutTreeAt(void *, void *, void *, void *, void *, void *);
extern void edgeDuplicates(void *, void *, void *, void *, void *, void *);
extern void getEdgeSimilarities(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void getEdgeSimilarities_all(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void getJaccards(void *, void *, void *, void *, void *, void *);
extern void getLinkDensities(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void getNumClusters(void *, void *, void *, void *, void *, void *);
extern void getOCGclusters(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hclustLinkComm(void *, void *, void *, void *, void *, void *);
extern void hclustPlotOrder(void *, void *, void *, void *);
extern void readOCG(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"cutTreeAt",               (DL_FUNC) &cutTreeAt,                6},
    {"edgeDuplicates",          (DL_FUNC) &edgeDuplicates,           6},
    {"getEdgeSimilarities",     (DL_FUNC) &getEdgeSimilarities,     12},
    {"getEdgeSimilarities_all", (DL_FUNC) &getEdgeSimilarities_all, 13},
    {"getJaccards",             (DL_FUNC) &getJaccards,              6},
    {"getLinkDensities",        (DL_FUNC) &getLinkDensities,        15},
    {"getNumClusters",          (DL_FUNC) &getNumClusters,           6},
    {"getOCGclusters",          (DL_FUNC) &getOCGclusters,           9},
    {"hclustLinkComm",          (DL_FUNC) &hclustLinkComm,           6},
    {"hclustPlotOrder",         (DL_FUNC) &hclustPlotOrder,          4},
    {"readOCG",                 (DL_FUNC) &readOCG,                  3},
    {NULL, NULL, 0}
};

void R_init_linkcomm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

