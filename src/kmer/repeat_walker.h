#ifndef REPEAT_WALKER_H_
#define REPEAT_WALKER_H_

#include "graph_walker.h"

typedef struct
{
  uint64_t *const visited, *const bloom;
  const size_t bloom_nbits, mem_bytes;
  const uint32_t mask;
} RepeatWalker;

// GraphWalker wlk is proposing node and orient as next move
// We determine if it is safe to make the traversal without getting stuck in
// a loop/cycle in the graph
static inline boolean walker_attempt_traverse(RepeatWalker *rpt,
                                              const GraphWalker *wlk,
                                              hkey_t node, Orientation orient,
                                              const BinaryKmer bkmer)
{
  if(!db_node_has_traversed(rpt->visited, node, orient)) {
    db_node_set_traversed(rpt->visited, node, orient);
    return true;
  }
  else {
    uint32_t hash32 = graph_walker_fasthash(wlk, bkmer) & rpt->mask;
    boolean collision = bitset_has(rpt->bloom, hash32);
    bitset_set(rpt->bloom, hash32);
    return !collision;
  }
}

static inline void walker_alloc(RepeatWalker *rpt,
                                size_t hash_capacity, size_t nbits)
{
  assert(nbits > 0 && nbits <= 32);
  size_t visited_words = round_bits_to_words64(hash_capacity*2);
  size_t repeat_words = round_bits_to_words64(1UL<<nbits);
  size_t nbytes = (visited_words + repeat_words) * sizeof(uint64_t);
  uint64_t *mem = calloc2(visited_words+repeat_words, sizeof(uint64_t));
  uint32_t mask = ~(uint32_t)0 >> (32-nbits); // beware: UINT32_MAX if nbits == 0
  RepeatWalker tmp = {.visited = mem, .bloom = mem+visited_words,
                      .bloom_nbits = nbits, .mem_bytes = nbytes, .mask = mask};
  memcpy(rpt, &tmp, sizeof(RepeatWalker));
}

static inline void walker_dealloc(RepeatWalker *rpt)
{
  free(rpt->visited);
}

static inline void walker_clear(RepeatWalker *rpt)
{
  memset(rpt->visited, 0, rpt->mem_bytes);
}

static inline void walker_fast_clear(RepeatWalker *rpt, hkey_t *nodes, size_t n)
{
  size_t i, bmem = round_bits_to_words64(1UL<<rpt->bloom_nbits)*sizeof(uint64_t);
  for(i = 0; i < n; i++) db_node_fast_clear_traversed(rpt->visited, nodes[i]);
  memset(rpt->bloom, 0, bmem);
}

#endif /* REPEAT_WALKER_H_ */
