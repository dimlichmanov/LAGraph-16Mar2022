//------------------------------------------------------------------------------
// LG_BreadthFirstSearch_vanilla:  BFS using only GraphBLAS API
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

// Contributed by Scott McMillan, derived from examples in the appendix of
// The GraphBLAS C API Specification, v1.3.0

//------------------------------------------------------------------------------

// This is a Basic algorithm (no extra G-> properties are required),
// but it is not user-callable (see LAGr_BreadthFirstSearch instead).

#define LG_FREE_WORK        \
{                           \
    GrB_free (&frontier);   \
}

#define LG_FREE_ALL         \
{                           \
    LG_FREE_WORK ;          \
    GrB_free (&l_parent);   \
    GrB_free (&l_level);    \
}

#include "LG_internal.h"

//****************************************************************************
int LG_BreadthFirstSearch_vanilla
(
    GrB_Vector    *level,
    GrB_Vector    *parent,
    const LAGraph_Graph G,
    GrB_Index      src,
    char          *msg
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    GrB_Vector frontier = NULL;     // the current frontier
    GrB_Vector l_parent = NULL;     // parent vector
    GrB_Vector l_level = NULL;      // level vector

    bool compute_level  = (level != NULL);
    bool compute_parent = (parent != NULL);
    if (compute_level ) (*level ) = NULL;
    if (compute_parent) (*parent) = NULL;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;

    if (!(compute_level || compute_parent))
    {
        // nothing to do
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // get the problem size and properties
    //--------------------------------------------------------------------------

    GrB_Matrix A = G->A ;

    GrB_Index n;
    GRB_TRY( GrB_Matrix_nrows (&n, A) );
    LG_ASSERT_MSG (src < n, GrB_INVALID_INDEX, "invalid source node") ;

    // create a sparse boolean vector frontier, and set frontier(src) = true
    GRB_TRY (GrB_Vector_new(&frontier, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_setElement(frontier, true, src)) ;

    GRB_TRY (GrB_Vector_new(&l_level, int_type, n)) ;

    //--------------------------------------------------------------------------
    // BFS traversal and label the nodes
    //--------------------------------------------------------------------------
    GrB_Index nq = 1 ;          // number of nodes in the current level
    GrB_Index last_nq = 0 ;
    GrB_Index current_level = 1;
    GrB_Index nvals = 1;

    // {!mask} is the set of unvisited nodes
    GrB_Vector mask = l_level ;

    // parent BFS
    do
    {
        // assign levels: l_level<s(frontier)> = current_level
        GRB_TRY( GrB_assign(l_level, frontier, GrB_NULL,
                            current_level, GrB_ALL, n, GrB_DESC_S) );
        ++current_level;
        // frontier = kth level of the BFS
        // mask is l_parent if computing parent, l_level if computing just level
        //WARNING! USING USER-VERSION OF A SEMIRING
        GRB_TRY( GrB_vxm(frontier, mask, GrB_NULL, lablas::LogicalOrAndSemiring<bool>(),
                         frontier, A, GrB_DESC_RSC) );

        // done if frontier is empty
        GRB_TRY( GrB_Vector_nvals(&nvals, frontier) );
    } while (nvals > 0);

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    (*level ) = l_level ;
    return (GrB_SUCCESS) ;
}
