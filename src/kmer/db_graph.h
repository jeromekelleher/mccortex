#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <inttypes.h>
#include "graph_typedef.h"

#define db_graph_node_assigned(graph,hkey) HASH_ENTRY_ASSIGNED((graph)->ht.table[hkey])

void db_graph_alloc(dBGraph *db_graph, size_t kmer_size,
                    size_t num_of_cols, size_t num_edge_cols,
                    uint64_t capacity);

void db_graph_realloc(dBGraph *graph, size_t num_of_cols, size_t num_edge_cols);

void db_graph_dealloc(dBGraph *db_graph);

// Get an oriented bkmer
#define db_graph_oriented_bkmer(graph,hkey,or) \
        db_node_oriented_bkmer(db_node_bkmer(graph,hkey),or,(graph)->kmer_size)

//
// Add to the de bruijn graph
//

// Note: node may alreay exist in the graph
hkey_t db_graph_find_or_add_node(dBGraph *db_graph, BinaryKmer bkey, Colour col);

// In the case of self-loops in palindromes the two edges collapse into one
void db_graph_add_edge(dBGraph *db_graph, Colour colour,
                       hkey_t src_node, hkey_t tgt_node,
                       Orientation src_orient, Orientation tgt_orient);

// For debugging + healthcheck
void db_graph_check_edges(const dBGraph *db_graph,
                          hkey_t src_node, hkey_t tgt_node,
                          Orientation src_orient, Orientation tgt_orient);

//
// Graph Traversal
//

void db_graph_next_node(const dBGraph *db_graph, const BinaryKmer node_bkey,
                        Nucleotide next_nuc, Orientation orient,
                        hkey_t *next_node, Orientation *next_orient);

// edges are forward+reverse, db_graph_next_nodes orients them
// fw_nucs is the nuc you would add when walking forward
size_t db_graph_next_nodes(const dBGraph *db_graph, const BinaryKmer node_bkey,
                           Orientation orient, Edges edges,
                           hkey_t nodes[4], Orientation orients[4],
                           Nucleotide fw_nucs[4]);

//
// Healthcheck
//

void db_graph_healthcheck(const dBGraph *db_graph);

//
// Functions applying to whole graph
//
void db_graph_wipe_colour(dBGraph *db_graph, Colour col);

// Add edges between all kmers with k-1 bases overlapping
void db_graph_add_all_edges(dBGraph *db_graph);

void db_graph_dump_paths_by_kmer(const dBGraph *db_graph);

// Get a random node from the graph
// call seed_random() before any calls to this function please
hkey_t db_graph_rand_node(const dBGraph *db_graph);

#endif /* DB_GRAPH_H_ */
