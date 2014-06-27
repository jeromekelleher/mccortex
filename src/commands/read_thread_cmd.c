#include "global.h"
#include "cmd.h"
#include "read_thread_cmd.h"
#include "gpath_checks.h"

//
// ctx_thread.c and ctx_correct.c use many of the same command line arguments
// Therefore this file provides the common functionality to both
//

void read_thread_args_alloc(struct ReadThreadCmdArgs *args)
{
  correct_aln_input_buf_alloc(&args->inputs, 16);
  gpfile_buf_alloc(&args->gpfiles, 16);
}

void read_thread_args_dealloc(struct ReadThreadCmdArgs *args)
{
  size_t i;
  for(i = 0; i < args->inputs.len; i++) asyncio_task_close(&args->inputs.data[i].files);
  for(i = 0; i < args->gpfiles.len; i++) gpath_reader_close(&args->gpfiles.data[i]);

  correct_aln_input_buf_dealloc(&args->inputs);
  gpfile_buf_dealloc(&args->gpfiles);
}

void read_thread_args_parse(struct ReadThreadCmdArgs *args,
                            int argc, char **argv,
                            const struct option *longopts, bool correct_cmd)
{
  size_t i;
  CorrectAlnInput task = CORRECT_ALN_INPUT_INIT;
  uint8_t fq_offset = 0;
  size_t dump_seq_n = 0, dump_mp_n = 0; // how many times are -g -G specified
  // PathFileReader tmp_pfile;
  GPathReader tmp_gpfile;

  CorrectAlnInputBuffer *inputs = &args->inputs;

  // Arg parsing
  char cmd[100];
  char shortopts[300];
  cmd_long_opts_to_short(longopts, shortopts, sizeof(shortopts));
  int used = 1, c;
  char *tmp_path;

  // silence error messages from getopt_long
  // opterr = 0;

  while((c = getopt_long_only(argc, argv, shortopts, longopts, NULL)) != -1) {
    cmd_get_longopt_str(longopts, c, cmd, sizeof(cmd));
    switch(c) {
      case 0: /* flag set */ break;
      case 'h': cmd_print_usage(NULL); break;
      case 'o':
        if(args->out_ctp_path != NULL) cmd_print_usage(NULL);
        args->out_ctp_path = optarg;
        break;
      case 'p':
        memset(&tmp_gpfile, 0, sizeof(GPathReader));
        gpath_reader_open(&tmp_gpfile, optarg, true);
        gpfile_buf_add(&args->gpfiles, tmp_gpfile);
        break;
      case 't':
        if(args->num_of_threads != 0) die("%s set twice", cmd);
        args->num_of_threads = cmd_uint32_nonzero(cmd, optarg);
        break;
      case 'm': cmd_mem_args_set_memory(&args->memargs, optarg); break;
      case 'n': cmd_mem_args_set_nkmers(&args->memargs, optarg); break;
      case 'c': args->colour = cmd_uint32(cmd, optarg); break;
      case '1':
      case '2':
      case 'i':
        used = 1;
        correct_aln_input_buf_add(inputs, task);
        asyncio_task_parse(&inputs->data[inputs->len-1].files, c, optarg,
                           fq_offset, correct_cmd ? &tmp_path : NULL);
        if(correct_cmd) inputs->data[inputs->len-1].out_base = tmp_path;
        break;
      case 'f': task.matedir = READPAIR_FR; used = 0; break;
      case 'F': task.matedir = READPAIR_FF; used = 0; break;
      case 'r': task.matedir = READPAIR_RF; used = 0; break;
      case 'R': task.matedir = READPAIR_RR; used = 0; break;
      case 'w': task.crt_params.one_way_gap_traverse = true; used = 0; break;
      case 'W': task.crt_params.one_way_gap_traverse = false; used = 0; break;
      case 'q': fq_offset = cmd_uint8(cmd, optarg); used = 0; break;
      case 'Q': task.fq_cutoff = cmd_uint8(cmd, optarg); used = 0; break;
      case 'H': task.hp_cutoff = cmd_uint8(cmd, optarg); used = 0; break;
      case 'g': task.crt_params.ins_gap_min = cmd_uint32(cmd, optarg); used = 0; break;
      case 'G': task.crt_params.ins_gap_max = cmd_uint32(cmd, optarg); used = 0; break;
      case 'd': task.crt_params.gap_wiggle = cmd_udouble(cmd, optarg); used = 0; break;
      case 'D': task.crt_params.gap_variance = cmd_udouble(cmd, optarg); used = 0; break;
      case 'X': task.crt_params.max_context = cmd_uint32(cmd, optarg); used = 0; break;
      case 'e': task.crt_params.use_end_check = true; used = 0; break;
      case 'E': task.crt_params.use_end_check = false; used = 0; break;
      case 'S': args->dump_seq_sizes = optarg; dump_seq_n++; break;
      case 'M': args->dump_mp_sizes = optarg; dump_mp_n++; break;
      case 'u': args->use_new_paths = true; break;
      case 'x': gen_paths_print_contigs = true; break;
      case 'y': gen_paths_print_paths = true; break;
      case 'z': gen_paths_print_reads = true; break;
      case ':': /* BADARG */
      case '?': /* BADCH getopt_long has already printed error */
        // cmd_print_usage(NULL);
        die("`"CMD" thread -h` for help. Bad option: %s", argv[optind-1]);
      default: abort();
    }
  }

  if(args->num_of_threads == 0) args->num_of_threads = DEFAULT_NTHREADS;

  // Check that optind+1 == argc
  if(optind+1 > argc)
    cmd_print_usage("Expected exactly one graph file");
  else if(optind+1 < argc)
    cmd_print_usage("Expected only one graph file. What is this: '%s'", argv[optind]);

  char *graph_path = argv[optind];
  status("Reading graph: %s", graph_path);

  if(!used) cmd_print_usage("Ignored arguments after last --seq");

  if(dump_seq_n > 1) die("Cannot specify --seq-gaps <out> more than once");
  if(dump_mp_n > 1) die("Cannot specify --mp-gaps <out> more than once");

  //
  // Open graph graph file
  //
  GraphFileReader *gfile = &args->gfile;
  graph_file_open(gfile, graph_path, true);
  file_filter_update_intocol(&gfile->fltr, 0);
  if(!correct_cmd && graph_file_usedcols(gfile) > 1)
    die("Please specify a single colour e.g. %s:0", gfile->fltr.file_path.buff);

  //
  // Open path files
  //
  size_t path_max_usedcols = 0;
  for(i = 0; i < args->gpfiles.len; i++) {
    // file_filter_update_intocol(&args->pfiles.data[i].fltr, 0);
    if(!correct_cmd && file_filter_usedcols(&args->gpfiles.data[i].fltr) > 1) {
      die("Please specify a single colour e.g. %s:0",
          args->gpfiles.data[i].fltr.file_path.buff);
    }
    path_max_usedcols = MAX2(path_max_usedcols,
                             file_filter_usedcols(&args->gpfiles.data[i].fltr));
  }
  args->path_max_usedcols = path_max_usedcols;

  // Check for compatibility between graph files and path files
  graphs_gpaths_compatible(gfile, 1, args->gpfiles.data, args->gpfiles.len);

  // if no paths loaded, set all max_context values to 1, since >1 kmer only
  // useful if can pickup paths
  if(args->gpfiles.len == 0) {
    for(i = 0; i < inputs->len; i++)
      inputs->data[i].crt_params.max_context = 1;
  }

  // Check ins_gap_min < ins_gap_max
  for(i = 0; i < inputs->len; i++)
  {
    CorrectAlnInput *t = &inputs->data[i];
    t->files.ptr = t;
    if(t->crt_params.ins_gap_min > t->crt_params.ins_gap_max) {
      die("--min-ins %u is greater than --max-ins %u",
          t->crt_params.ins_gap_min, t->crt_params.ins_gap_max);
    }
    correct_aln_input_print(&inputs->data[i]);
    args->max_gap_limit = MAX2(args->max_gap_limit, t->crt_params.ins_gap_max);
  }
}
