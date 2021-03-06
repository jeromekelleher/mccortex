#include "global.h"

#include "commands.h"
#include "util.h"
#include "file_util.h"
#include "db_graph.h"
#include "graph_format.h"
#include "path_format.h"
#include "graph_paths.h"

const char pview_usage[] = ""
"usage: "CMD" pview [options] <in.ctp>\n"
"\n"
"  View and check a paths file.\n"
"\n"
"  -p, --print  Print paths\n"
"  -c, --check  Check path file integrity\n"
"\n";

int print_paths = 0, do_paths_check = 0;

static struct option longopts[] =
{
  {"memory", required_argument, NULL,            'm'},
  {"nkmers", required_argument, NULL,            'n'},
  {"print",  no_argument,       &print_paths,    1},
  {"check",  no_argument,       &do_paths_check, 1},
  {NULL, 0, NULL, 0}
};

int ctx_pview(int argc, char **argv)
{
  struct MemArgs memargs = MEM_ARGS_INIT;

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
      case 'm': cmd_mem_args_set_memory(&memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&memargs, optarg); break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" pview -h` for help. Bad option: %s", argv[optind-1]);
      default:
        die("Programmer fail. Tell Isaac.");
    }
  }

  if(optind+1 != argc) cmd_print_usage("Require one input path file (.ctp)");

  char *input_paths_file = argv[optind];

  // Open paths file
  PathFileReader pfile = INIT_PATH_READER;
  path_file_open(&pfile, input_paths_file, true);

  PathFileHeader *phdr = &pfile.hdr;
  char num_paths_str[100], path_bytes_str[100], kmers_with_paths_str[100];
  ulong_to_str(phdr->num_of_paths, num_paths_str);
  bytes_to_str(phdr->num_path_bytes, 1, path_bytes_str);
  ulong_to_str(phdr->num_kmers_with_paths, kmers_with_paths_str);

  // Print header
  printf("version: %u\n", phdr->version);
  printf("kmer size: %u\n", phdr->kmer_size);
  printf("colours: %u\n", phdr->num_of_cols);
  printf("paths: %s\n", num_paths_str);
  printf("bytes: %s\n", path_bytes_str);
  printf("kmers starting paths: %s\n", kmers_with_paths_str);

  size_t col;
  for(col = 0; col < phdr->num_of_cols; col++) {
    printf(" colour %zu: %s\n", col, phdr->sample_names[col].buff);
  }

  //
  // Decide on memory
  //
  size_t bits_per_kmer, kmers_in_hash, graph_mem;

  bits_per_kmer = sizeof(uint64_t) * 8;
  kmers_in_hash = cmd_get_kmers_in_hash2(memargs.mem_to_use,
                                         memargs.mem_to_use_set,
                                         memargs.num_kmers,
                                         memargs.num_kmers_set,
                                         bits_per_kmer,
                                         phdr->num_kmers_with_paths,
                                         phdr->num_kmers_with_paths,
                                         true, &graph_mem);

  cmd_check_mem_limit(memargs.mem_to_use, graph_mem);

  // Allocate memory
  // db graph is required to store the end position for each kmer list
  dBGraph db_graph;
  db_graph_alloc(&db_graph, phdr->kmer_size, phdr->num_of_cols, 0, kmers_in_hash);

  path_file_set_graph_sample_names(&pfile, &db_graph);

  // Paths
  path_store_alloc(&db_graph.pstore, phdr->num_path_bytes, false,
                   db_graph.ht.capacity, phdr->num_of_cols);

  // Pretend we've read all the kmers in
  db_graph.num_of_cols_used = phdr->num_of_cols;

  // Add kmers as reading
  bool add_kmers = true;

  paths_format_load(&pfile, add_kmers, &db_graph);

  if(print_paths)
    db_graph_dump_paths_by_kmer(&db_graph);

  // Check data store
  if(do_paths_check) {
    status("Checking path store integrity...");
    path_store_integrity_check(&db_graph.pstore);
  }

  path_file_close(&pfile);
  db_graph_dealloc(&db_graph);

  return EXIT_SUCCESS;
}
