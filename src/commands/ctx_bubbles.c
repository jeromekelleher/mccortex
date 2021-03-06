#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_node.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_paths.h"
#include "graph_walker.h"
#include "bubble_caller.h"
#include "path_file_reader.h"

// Long flanks help us map calls
// increasing allele length can be costly
#define DEFAULT_MAX_FLANK 1000
#define DEFAULT_MAX_ALLELE 300

const char bubbles_usage[] =
"usage: "CMD" bubbles [options] <in.ctx> [in2.ctx ...]\n"
"\n"
"  Find bubbles in the graph, which are potential variants.\n"
"\n"
"  -h, --help              This help message\n"
"  -o, --out <bub.txt.gz>  Output file [required]\n"
"  -m, --memory <mem>      Memory to use\n"
"  -n, --nkmers <kmers>    Number of hash table entries (e.g. 1G ~ 1 billion)\n"
"  -t, --threads <T>       Number of threads to use [default: "QUOTE_VALUE(DEFAULT_NTHREADS)"]\n"
"  -p, --paths <in.ctp>    Load path file (can specify multiple times)\n"
//
"  -H, --haploid <col>     Colour is haploid, can use repeatedly [e.g. ref colour]\n"
"  -a, --max-allele <len>  Max bubble branch length in kmers [default: "QUOTE_VALUE(DEFAULT_MAX_ALLELE)"]\n"
"  -f, --max-flank <len>   Max flank length in kmers [default: "QUOTE_VALUE(DEFAULT_MAX_FLANK)"]\n"
"\n"
"  When loading path files with -p, use offset (e.g. 2:in.ctp) to specify\n"
"  which colour to load the data into.\n"
"\n";

static struct option longopts[] =
{
// General options
  {"help",         no_argument,       NULL, 'h'},
  {"out",          required_argument, NULL, 'o'},
  {"memory",       required_argument, NULL, 'm'},
  {"nkmers",       required_argument, NULL, 'n'},
  {"threads",      required_argument, NULL, 't'},
  {"paths",        required_argument, NULL, 'p'},
// command specific
  {"haploid",      required_argument, NULL, 'H'},
  {"max-allele",   required_argument, NULL, 'a'},
  {"max-flank",    required_argument, NULL, 'f'},
  {NULL, 0, NULL, 0}
};

#include "objbuf_macro.h"
create_objbuf(size_buf, SizeBuffer, size_t);

int ctx_bubbles(int argc, char **argv)
{
  size_t num_of_threads = 0;
  struct MemArgs memargs = MEM_ARGS_INIT;
  const char *out_path = NULL;
  size_t max_allele_len = 0, max_flank_len = 0;

  SizeBuffer haploidbuf;
  size_buf_alloc(&haploidbuf, 8);

  PathFileReader tmp_pfile;
  PathFileBuffer pfilesbuf;
  pfile_buf_alloc(&pfilesbuf, 8);

  // tmp
  size_t tmp_col;

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int c;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o':
        if(out_path != NULL) cmd_print_usage(NULL);
        out_path = optarg;
        break;
      case 'p':
        tmp_pfile = INIT_PATH_READER;
        path_file_open(&tmp_pfile, optarg, true);
        pfile_buf_add(&pfilesbuf, tmp_pfile);
        break;
      case 't':
        if(num_of_threads) die("%s set twice", cmd);
        num_of_threads = cmd_parse_arg_uint32_nonzero(cmd, optarg);
        break;
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case 'H':
        tmp_col = cmd_parse_arg_uint32(cmd, optarg);
        size_buf_add(&haploidbuf, tmp_col);
        break;
      case 'a':
        if(max_allele_len) die("%s set twice", cmd);
        max_allele_len = cmd_parse_arg_uint32_nonzero(cmd, optarg);
        break;
      case 'f':
        if(max_flank_len) die("%s set twice", cmd);
        max_flank_len = cmd_parse_arg_uint32_nonzero(cmd, optarg);
        break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" bubbles -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  // Defaults for unset values
  if(out_path == NULL) out_path = "-";
  if(num_of_threads == 0) num_of_threads = DEFAULT_NTHREADS;
  if(max_allele_len == 0) max_allele_len = DEFAULT_MAX_ALLELE;
  if(max_flank_len == 0) max_flank_len = DEFAULT_MAX_FLANK;

  if(optind >= argc) cmd_print_usage("Require input graph files (.ctx)");

  //
  // Open graph files
  //
  const size_t num_gfiles = argc - optind;
  char **graph_paths = argv + optind;
  ctx_assert(num_gfiles > 0);

  GraphFileReader *gfiles = ctx_calloc(num_gfiles, sizeof(GraphFileReader));
  size_t i, ncols, ctx_max_kmers = 0, ctx_sum_kmers = 0;

  ncols = graph_files_open(graph_paths, gfiles, num_gfiles,
                           &ctx_max_kmers, &ctx_sum_kmers);

  //
  // Open path files
  //
  size_t path_max_usedcols = 0;

  for(i = 0; i < pfilesbuf.len; i++) {
    path_max_usedcols = MAX2(path_max_usedcols,
                             path_file_usedcols(&pfilesbuf.data[i]));
  }

  // Check graph + paths are compatible
  graphs_paths_compatible(gfiles, num_gfiles, pfilesbuf.data, pfilesbuf.len);

  //
  // Check haploid colours are valid
  //
  for(i = 0; i < haploidbuf.len; i++) {
    if(haploidbuf.data[i] >= ncols) {
      cmd_print_usage("-H,--haploid <col> is greater than max colour [%zu > %zu]",
                      haploidbuf.data[i], ncols-1);
    }
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem, path_mem, thread_mem;
  char thread_mem_str[100];

  // edges(1bytes) + kmer_paths(8bytes) + in_colour(1bit/col) +
  // visitedfw/rv(2bits/thread)

  bits_per_kmer = sizeof(Edges)*8 + sizeof(PathIndex)*8 +
                  ncols + 2*num_of_threads;

  kmers_in_hash = cmd_get_kmers_in_hash2(memargs.mem_to_use,
                                         memargs.mem_to_use_set,
                                         memargs.num_kmers,
                                         memargs.num_kmers_set,
                                         bits_per_kmer,
                                         ctx_max_kmers, ctx_sum_kmers,
                                         false, &graph_mem);

  // Thread memory
  thread_mem = roundup_bits2bytes(kmers_in_hash) * 2;
  bytes_to_str(thread_mem * num_of_threads, 1, thread_mem_str);
  status("[memory] (of which threads: %zu x %zu = %s)\n",
          num_of_threads, thread_mem, thread_mem_str);

  // Path Memory
  path_mem = path_files_mem_required(pfilesbuf.data, pfilesbuf.len,
                                     false, false, path_max_usedcols, 0);
  cmd_print_mem(path_mem, "paths");

  size_t total_mem = graph_mem + thread_mem + path_mem;
  cmd_check_mem_limit(memargs.mem_to_use, total_mem);

  //
  // Open output file
  //
  gzFile gzout = futil_gzopen_output(out_path);

  // Allocate memory
  dBGraph db_graph;
  db_graph_alloc(&db_graph, gfiles[0].hdr.kmer_size, ncols, 1, kmers_in_hash);

  // Edges merged into one colour
  db_graph.col_edges = ctx_calloc(db_graph.ht.capacity, sizeof(uint8_t));

  // In colour
  size_t bytes_per_col = roundup_bits2bytes(db_graph.ht.capacity);
  db_graph.node_in_cols = ctx_calloc(bytes_per_col*ncols, sizeof(uint8_t));

  // Paths
  path_store_alloc(&db_graph.pstore, path_mem, false,
                   db_graph.ht.capacity, path_max_usedcols);

  //
  // Load graphs
  //
  LoadingStats stats = LOAD_STATS_INIT_MACRO;

  GraphLoadingPrefs gprefs = {.db_graph = &db_graph,
                              .boolean_covgs = false,
                              .must_exist_in_graph = false,
                              .must_exist_in_edges = NULL,
                              .empty_colours = true};

  for(i = 0; i < num_gfiles; i++) {
    graph_load(&gfiles[i], gprefs, &stats);
    graph_file_close(&gfiles[i]);
    gprefs.empty_colours = false;
  }
  ctx_free(gfiles);

  hash_table_print_stats(&db_graph.ht);

  // Load path files (does nothing if num_fpiles == 0)
  paths_format_merge(pfilesbuf.data, pfilesbuf.len, false, false,
                     num_of_threads, &db_graph);

  for(i = 0; i < pfilesbuf.len; i++) path_file_close(&pfilesbuf.data[i]);

  // Now call variants
  BubbleCallingPrefs call_prefs = {.max_allele_len = max_allele_len,
                                   .max_flank_len = max_flank_len,
                                   .haploid_cols = haploidbuf.data,
                                   .num_haploid = haploidbuf.len};

  invoke_bubble_caller(num_of_threads, call_prefs,
                       gzout, out_path, &db_graph);

  status("  saved to: %s\n", out_path);
  gzclose(gzout);

  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
