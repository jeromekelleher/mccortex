#ifndef GRAPH_PATH_H_
#define GRAPH_PATH_H_

#include "cortex_types.h"
#include "db_graph.h"
#include "db_node.h"
#include "path_hash.h"
#include "path_file_reader.h"
#include "graph_file_reader.h"
#include "path_set.h"

typedef struct {
  Colour *ctxcols, *ctpcols;
  size_t n;
} GraphPathPairing;

static inline void gp_alloc(GraphPathPairing *gp, size_t n) {
  gp->ctxcols = ctx_malloc(2*n*sizeof(*gp->ctxcols));
  gp->ctpcols = gp->ctxcols+n;
  gp->n = n;
}

static inline void gp_dealloc(GraphPathPairing *gp) {
  ctx_free(gp->ctxcols);
  gp->n = 0; gp->ctxcols = gp->ctpcols = NULL;
}

// De-fragment the linked list used by PathStore
// Uses memory or disk to de-fragment the file
// Warning: The PathStore must not be modified during the running of the
//          function. That means it is NOT thread-safe.
void graph_paths_defragment(dBGraph *db_graph, FILE* tmp_fh);

// Similar to path_file_reader.c:path_file_load_check()
// Check kmer size matches and sample names match
void graphs_paths_compatible(const GraphFileReader *graphs, size_t num_graphs,
                             const PathFileReader *paths, size_t num_paths);

// Returns true if new to colour, false otherwise
// packed points to <PathLen><PackedSeq>
// Returns address of path in PathStore by setting newidx
bool graph_paths_find_or_add_mt(dBNode node, Colour ctpcol,
                                const uint8_t *packed, size_t plen,
                                PathStore *pstore, PathIndex *newidx);

//
// Update PathStore from PathSet which has had entries removed / reordered
//
void graph_paths_update_from_set(PathStore *pstore, const PathSet *set, hkey_t hkey);

//
// Remove all redundant paths
//
// Remove colours from paths with counts < threshold
void graph_paths_clean(dBGraph *db_graph, size_t num_threads, uint8_t threshold);

//
// Integrity checks on graph+paths
//

// col is graph colour
// packed is just <PackedBases>
bool graph_paths_check_valid(dBNode node, size_t ctxcol, const uint8_t *packed,
                             size_t nbases, const dBGraph *db_graph);

bool graph_paths_check_all_paths(const GraphPathPairing *gp,
                                    const dBGraph *db_graph);

bool graph_paths_check_path(hkey_t node, PathIndex pindex,
                              const GraphPathPairing *gp,
                              const dBGraph *db_graph);

// For debugging
void graph_paths_check_counts(const dBGraph *db_graph);

#endif /* GRAPH_PATH_H_ */
