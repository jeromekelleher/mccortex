#include "global.h"
#include "build_graph.h"
#include "db_graph.h"
#include "db_node.h"
#include "seq_reader.h"
#include "async_read_io.h"
#include "loading_stats.h"
#include "util.h"

#include <pthread.h>
#include "seq_file.h"

typedef struct
{
  pthread_t thread;
  dBGraph *const db_graph;
  MsgPool *pool;
  volatile size_t *rcounter; // counter of entries taken from the pool
  LoadingStats *file_stats; // Array of stats for diff input files
} BuildGraphWorker;

//
// Check for PCR duplicates
//

// Read start (duplicate removal during read loading)
#define db_node_has_read_start_mt(graph,node) \
        bitset_get_mt((volatile uint8_t*)(graph)->readstrt, 2*(node).key+(node).orient)
#define db_node_set_read_start_mt(graph,node) \
        bitset_set_mt((volatile uint8_t*)(graph)->readstrt, 2*(node).key+(node).orient)

// Returns true if start1, start2 set and reads should be added
static bool seq_reads_are_novel(read_t *r1, read_t *r2,
                                uint8_t fq_cutoff1, uint8_t fq_cutoff2,
                                uint8_t hp_cutoff, ReadMateDir matedir,
                                LoadingStats *stats, dBGraph *db_graph)
{
  // Remove SAM/BAM duplicates
  if(r1->from_sam && r1->bam->core.flag & BAM_FDUP &&
     (r2 == NULL || (r2->from_sam && r2->bam->core.flag & BAM_FDUP))) {
    return false;
  }

  seq_reader_orient_mp_FF(r1, r2, matedir);

  const size_t kmer_size = db_graph->kmer_size;
  size_t start1, start2 = 0;
  bool got_kmer1 = false, got_kmer2 = false;
  BinaryKmer bkmer1, bkmer2;
  dBNode node1 = DB_NODE_INIT, node2 = DB_NODE_INIT;

  start1 = seq_contig_start(r1, 0, kmer_size, fq_cutoff1, hp_cutoff);
  got_kmer1 = (start1 < r1->seq.end);

  if(r2) {
    start2 = seq_contig_start(r2, 0, kmer_size, fq_cutoff2, hp_cutoff);
    got_kmer2 = (start2 < r2->seq.end);
  }

  bool found1 = false, found2 = false;

  // Look up first kmer
  if(got_kmer1) {
    bkmer1 = binary_kmer_from_str(r1->seq.b + start1, kmer_size);
    node1 = db_graph_find_or_add_node_mt(db_graph, bkmer1, &found1);
  }

  // Look up second kmer
  if(got_kmer2) {
    bkmer2 = binary_kmer_from_str(r2->seq.b + start2, kmer_size);
    node2 = db_graph_find_or_add_node_mt(db_graph, bkmer2, &found2);
  }

  stats->num_kmers_novel += !found1 + !found2;

  // Each read gives no kmer or a duplicate kmer
  // used find_or_insert so if we have a kmer we have a graph node
  if((!got_kmer1 || db_node_has_read_start_mt(db_graph, node1)) &&
     (!got_kmer2 || db_node_has_read_start_mt(db_graph, node2)))
  {
    return false;
  }

  // Read is novel
  if(got_kmer1) db_node_set_read_start_mt(db_graph, node1);
  if(got_kmer2) db_node_set_read_start_mt(db_graph, node2);

  return true;
}


//
// Add to the de bruijn graph
//

// Threadsafe
// Sequence must be entirely ACGT and len >= kmer_size
// Returns number of novel kmers loaded
size_t build_graph_from_str_mt(dBGraph *db_graph, size_t colour,
                               const char *seq, size_t len)
{
  ctx_assert(len >= db_graph->kmer_size);
  const size_t kmer_size = db_graph->kmer_size;
  BinaryKmer bkmer;
  Nucleotide nuc;
  dBNode prev, curr;
  size_t i, num_novel_kmers = 0;
  size_t edge_col = db_graph->num_edge_cols == 1 ? 0 : colour;
  bool found;

  bkmer = binary_kmer_from_str(seq, kmer_size);
  prev = db_graph_find_or_add_node_mt(db_graph, bkmer, &found);
  db_graph_update_node_mt(db_graph, prev, colour);
  num_novel_kmers += !found;

  for(i = kmer_size; i < len; i++)
  {
    nuc = dna_char_to_nuc(seq[i]);
    bkmer = binary_kmer_left_shift_add(bkmer, kmer_size, nuc);
    curr = db_graph_find_or_add_node_mt(db_graph, bkmer, &found);
    db_graph_update_node_mt(db_graph, curr, colour);
    db_graph_add_edge_mt(db_graph, edge_col, prev, curr);
    num_novel_kmers += !found;
    prev = curr;
  }

  return num_novel_kmers;
}

// Already found a start position
static void load_read(const read_t *r, uint8_t qual_cutoff, uint8_t hp_cutoff,
                      LoadingStats *stats, Colour colour, dBGraph *db_graph)
{
  const size_t kmer_size = db_graph->kmer_size;
  size_t contig_start, contig_end, contig_len;
  size_t num_contigs = 0, search_start = 0, num_novel_kmers;

  while((contig_start = seq_contig_start(r, search_start, kmer_size,
                                         qual_cutoff, hp_cutoff)) < r->seq.end)
  {
    contig_end = seq_contig_end(r, contig_start, kmer_size,
                                qual_cutoff, hp_cutoff, &search_start);

    contig_len = contig_end - contig_start;
    num_novel_kmers = build_graph_from_str_mt(db_graph, colour,
                                              r->seq.b+contig_start, contig_len);

    stats->total_bases_loaded += contig_len;
    stats->num_kmers_loaded += contig_len + 1 - kmer_size;
    stats->num_kmers_novel += num_novel_kmers;
    num_contigs++;
  }

  stats->contigs_loaded += num_contigs;

  if(num_contigs) stats->num_good_reads++;
  else stats->num_bad_reads++;
}

void build_graph_from_reads_mt(read_t *r1, read_t *r2,
                               uint8_t fq_offset1, uint8_t fq_offset2,
                               uint8_t fq_cutoff, uint8_t hp_cutoff,
                               bool remove_pcr_dups, ReadMateDir matedir,
                               LoadingStats *stats, size_t colour,
                               dBGraph *db_graph)
{
  // status("r1: '%s' '%s'", r1->name.b, r1->seq.b);
  // if(r2) status("r2: '%s' '%s'", r2->name.b, r2->seq.b);

  uint8_t fq_cutoff1 = fq_cutoff, fq_cutoff2 = fq_cutoff;

  if(fq_cutoff) {
    fq_cutoff1 += fq_offset1;
    fq_cutoff2 += fq_offset2;
  }

  stats->total_bases_read += r1->seq.end + (r2 ? r2->seq.end : 0);

  if(r2) stats->num_pe_reads += 2;
  else stats->num_se_reads++;

  // printf(">%s %zu\n", r1->name.b, colour);

  if(remove_pcr_dups && !seq_reads_are_novel(r1, r2,
                                             fq_cutoff1, fq_cutoff2, hp_cutoff,
                                             matedir, stats, db_graph))
  {
    if(r2) stats->num_dup_pe_pairs++;
    else stats->num_dup_se_reads++;
  }
  else {
    load_read(r1, fq_cutoff1, hp_cutoff, stats, colour, db_graph);
    if(r2) load_read(r2, fq_cutoff2, hp_cutoff, stats, colour, db_graph);
  }
}

// Print progress every 5M reads
#define REPORT_RATE 5000000

static void build_graph_print_progress(size_t n)
{
  if(n % REPORT_RATE == 0)
  {
    char num_str[100];
    long_to_str(n, num_str);
    status("[BuildGraph] Read %s entries (reads / read pairs)", num_str);
  }
}

// pthread method, loop: reads from pool, add to graph
static void grab_reads_from_pool(void *ptr)
{
  BuildGraphWorker *wrkr = (BuildGraphWorker*)ptr;
  MsgPool *pool = wrkr->pool;
  BuildGraphTask *task;
  AsyncIOData *data;
  int pos;
  read_t *r2;

  while((pos = msgpool_claim_read(pool)) != -1)
  {
    memcpy(&data, msgpool_get_ptr(pool, pos), sizeof(AsyncIOData*));
    task = (BuildGraphTask*)data->ptr;

    r2 = data->r2.name.end == 0 && data->r2.seq.end == 0 ? NULL : &data->r2;

    build_graph_from_reads_mt(&data->r1, r2,
                              data->fq_offset1, data->fq_offset2,
                              task->fq_cutoff, task->hp_cutoff,
                              task->remove_pcr_dups, task->matedir,
                              &wrkr->file_stats[task->idx],
                              task->colour, wrkr->db_graph);

    msgpool_release(pool, pos, MPOOL_EMPTY);

    // Print progress
    size_t n = __sync_fetch_and_add(wrkr->rcounter, 1);
    build_graph_print_progress(n);
  }
}

// One thread used per input file, num_build_threads used to add reads to graph
void build_graph(dBGraph *db_graph, BuildGraphTask *files,
                 size_t num_files, size_t num_build_threads)
{
  ctx_assert(db_graph->bktlocks != NULL);

  size_t i, f;

  AsyncIOData *data = ctx_malloc(MSGPOOLSIZE * sizeof(AsyncIOData));
  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_alloc(&data[i]);

  MsgPool pool;
  msgpool_alloc(&pool, MSGPOOLSIZE, sizeof(AsyncIOData*), USE_MSG_POOL);
  msgpool_iterate(&pool, asynciodata_pool_init, data);

  // Start async io reading
  AsyncIOReadInput *async_tasks = ctx_malloc(num_files * sizeof(AsyncIOReadInput));

  for(f = 0; f < num_files; f++) {
    files[f].idx = f;
    files[f].files.ptr = &files[f];
    memcpy(&async_tasks[f], &files[f].files, sizeof(AsyncIOReadInput));
  }

  BuildGraphWorker *workers = ctx_malloc(num_build_threads * sizeof(BuildGraphWorker));
  size_t rcounter = 0;

  for(i = 0; i < num_build_threads; i++)
  {
    BuildGraphWorker tmp_wrkr = {.db_graph = db_graph, .pool = &pool,
                                 .rcounter = &rcounter};
    tmp_wrkr.file_stats = ctx_calloc(num_files, sizeof(LoadingStats));
    memcpy(&workers[i], &tmp_wrkr, sizeof(BuildGraphWorker));
  }

  // Create a lot of workers to build the graph
  asyncio_run_threads(&pool, async_tasks, num_files, grab_reads_from_pool,
                      workers, num_build_threads, sizeof(BuildGraphWorker));

  // start_build_graph_workers(&pool, db_graph, files, num_files, num_build_threads);

  ctx_free(async_tasks);
  msgpool_dealloc(&pool);

  // Clean up workers one by one...
  for(i = 0; i < num_build_threads; i++) {
    // Merge stats
    for(f = 0; f < num_files; f++)
      loading_stats_merge(&files[f].stats, &workers[i].file_stats[f]);

    // Free memory
    ctx_free(workers[i].file_stats);
  }

  ctx_free(workers);

  // Copy stats into ginfo
  size_t max_col = 0;
  for(f = 0; f < num_files; f++) {
    max_col = MAX2(max_col, files[f].colour);
    graph_info_update_stats(&db_graph->ginfo[files[f].colour], &files[f].stats);
  }

  db_graph->num_of_cols_used = MAX2(db_graph->num_of_cols_used, max_col+1);

  for(i = 0; i < MSGPOOLSIZE; i++) asynciodata_dealloc(&data[i]);
  ctx_free(data);
}


void build_graph_task_print(const BuildGraphTask *task)
{
  const AsyncIOReadInput *io = &task->files;
  char fqOffset[30] = "auto-detect", fqCutoff[30] = "off", hpCutoff[30] = "off";

  if(io->fq_offset > 0) sprintf(fqOffset, "%u", io->fq_offset);
  if(task->fq_cutoff > 0) sprintf(fqCutoff, "%u", task->fq_cutoff);
  if(task->hp_cutoff > 0) sprintf(hpCutoff, "%u", task->hp_cutoff);

  status("[task] %s%s%s; FASTQ offset: %s, threshold: %s; "
         "cut homopolymers: %s; remove PCR duplicates: %s; colour: %zu\n",
         io->file1->path,
         io->file2 ? ", " : "", io->file2 ? io->file2->path : "",
         fqOffset, fqCutoff, hpCutoff, task->remove_pcr_dups ? "yes" : "no",
         task->colour);
}

void build_graph_task_print_stats(const BuildGraphTask *task)
{
  const LoadingStats stats = task->stats;
  const AsyncIOReadInput *io = &task->files;

  status("[task] input: %s%s%s colour: %zu",
         io->file1->path, io->file2 ? ", " : "",
         io->file2 ? io->file2->path : "", task->colour);

  char se_reads_str[50], pe_reads_str[50];
  char good_reads_str[50], bad_reads_str[50];
  char dup_se_reads_str[50], dup_pe_pairs_str[50];
  char bases_read_str[50], bases_loaded_str[50];
  char num_contigs_str[50], num_kmers_loaded_str[50], num_kmers_novel_str[50];

  ulong_to_str(stats.num_se_reads, se_reads_str);
  ulong_to_str(stats.num_pe_reads, pe_reads_str);
  ulong_to_str(stats.num_good_reads, good_reads_str);
  ulong_to_str(stats.num_bad_reads, bad_reads_str);
  ulong_to_str(stats.num_dup_se_reads, dup_se_reads_str);
  ulong_to_str(stats.num_dup_pe_pairs, dup_pe_pairs_str);
  ulong_to_str(stats.total_bases_read, bases_read_str);
  ulong_to_str(stats.total_bases_loaded, bases_loaded_str);
  ulong_to_str(stats.contigs_loaded, num_contigs_str);
  ulong_to_str(stats.num_kmers_loaded, num_kmers_loaded_str);
  ulong_to_str(stats.num_kmers_novel, num_kmers_novel_str);

  status("  SE reads: %s  PE reads: %s", se_reads_str, pe_reads_str);
  status("  good reads: %s  bad reads: %s", good_reads_str, bad_reads_str);
  status("  dup SE reads: %s  dup PE pairs: %s", dup_se_reads_str, dup_pe_pairs_str);
  status("  bases read: %s  bases loaded: %s", bases_read_str, bases_loaded_str);
  status("  num contigs: %s  num kmers: %s novel kmers: %s",
         num_contigs_str, num_kmers_loaded_str, num_kmers_novel_str);
}
