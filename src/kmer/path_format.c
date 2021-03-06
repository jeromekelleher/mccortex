#include "global.h"
#include "path_format.h"
#include "db_graph.h"
#include "db_node.h"
#include "hash_table.h"
#include "path_store.h"
#include "util.h"
#include "file_util.h"
#include "path_set.h"
#include "graph_paths.h"

// Format:
// -- Header --
// "PATHS"<uint32_t:version><uint32_t:kmersize><uint32_t:num_of_cols>
// <uint64_t:num_of_paths><uint64_t:num_path_bytes><uint64_t:num_kmers_with_paths>
// -- Colours --
// <uint32_t:sname_len><uint8_t x sname_len:sample_name> x num_of_cols
// -- Data --
// <uint8_t x num_path_bytes:path_data>
// <binarykmer x num_kmers_with_paths><uint64_t:path_index>

void paths_header_alloc(PathFileHeader *h, size_t num_of_cols)
{
  size_t i, old_cap = h->capacity;

  if(h->capacity == 0) {
    h->sample_names = ctx_malloc(num_of_cols * sizeof(StrBuf));
  }
  else if(num_of_cols > h->capacity) {
    h->sample_names = ctx_realloc(h->sample_names, num_of_cols * sizeof(StrBuf));
  }

  for(i = old_cap; i < num_of_cols; i++) {
    strbuf_alloc(&h->sample_names[i], 256);
    strbuf_set(&h->sample_names[i], "noname");
  }

  h->capacity = MAX2(old_cap, num_of_cols);
}

void paths_header_dealloc(PathFileHeader *h)
{
  size_t i;
  if(h->capacity > 0) {
    for(i = 0; i < h->capacity; i++) strbuf_dealloc(&h->sample_names[i]);
    ctx_free(h->sample_names);
    h->capacity = 0;
  }
}

// Set path header variables based on PathStore
void paths_header_update(PathFileHeader *header, const PathStore *paths)
{
  header->num_of_cols = (uint32_t)paths->num_of_cols;
  header->num_of_paths = paths->num_of_paths;
  header->num_path_bytes = paths->num_of_bytes;
  header->num_kmers_with_paths = paths->num_kmers_with_paths;
}

// Returns number of bytes read or -1 on error (if fatal is false)
int paths_file_read_header(FILE *fh, PathFileHeader *h,
                            bool fatal, const char *path)
{
  int bytes_read = 0;
  char sig[6] = {0};

  SAFE_READ(fh, sig, 5, "PATHS", path, fatal);
  SAFE_READ(fh, &h->version, sizeof(uint32_t), "version", path, fatal);
  SAFE_READ(fh, &h->kmer_size, sizeof(uint32_t), "kmer_size", path, fatal);
  SAFE_READ(fh, &h->num_of_cols, sizeof(uint32_t), "num_of_cols", path, fatal);
  SAFE_READ(fh, &h->num_of_paths, sizeof(uint64_t), "num_of_paths", path, fatal);
  SAFE_READ(fh, &h->num_path_bytes, sizeof(uint64_t),
            "num_path_bytes", path, fatal);
  SAFE_READ(fh, &h->num_kmers_with_paths, sizeof(uint64_t),
            "num_kmers_with_paths", path, fatal);

  bytes_read += 5 + sizeof(uint32_t)*3 + sizeof(uint64_t)*3;

  if(h->num_of_cols > 10000) die("Large number of colours: %u", h->num_of_cols);

  // paths_header_alloc will only alloc or realloc only if it needs to
  paths_header_alloc(h, h->num_of_cols);

  // Read sample names
  size_t i;
  uint32_t len;
  StrBuf *sbuf;
  for(i = 0; i < h->num_of_cols; i++)
  {
    sbuf = h->sample_names + i;
    SAFE_READ(fh, &len, sizeof(uint32_t), "sample name length", path, fatal);
    strbuf_ensure_capacity(sbuf, len);
    SAFE_READ(fh, sbuf->buff, len, "sample name", path, fatal);
    sbuf->buff[sbuf->len = len] = '\0';
    bytes_read += sizeof(uint32_t) + len;
  }

  // Checks
  if(h->version < 1 || h->version > 1) {
    if(!fatal) return -1;
    die("file version not supported [version: %u; path: %s]", h->version, path);
  }

  if(strncmp(sig, "PATHS", 5) != 0) {
    if(!fatal) return -1;
    die("File is not valid paths file [path: %s]", path);
  }

  if(h->kmer_size % 2 == 0) {
    if(!fatal) return -1;
    die("kmer size is not an odd number [kmer_size: %u; path: %s]\n",
        h->kmer_size, path);
  }

  if(h->kmer_size < 3) {
    if(!fatal) return -1;
    die("kmer size is less than three [kmer_size: %u; path: %s]\n",
        h->kmer_size, path);
  }

  if(h->num_of_cols == 0) {
    if(!fatal) return -1;
    die("number of colours is zero [path: %s]\n", path);
  }

  return bytes_read;
}

size_t paths_get_max_usedcols(PathFileReader *files, size_t num_files)
{
  size_t i, ncols, used_cols = path_file_usedcols(&files[0]);

  for(i = 1; i < num_files; i++) {
    ncols = path_file_usedcols(&files[i]);
    used_cols = MAX2(used_cols, ncols);
  }
  return used_cols;
}


#define pfile_tmp_mem(pf,ps) \
  ((pf).hdr.num_path_bytes + \
   (pf).hdr.num_of_paths * \
    MAX2(((pf).hdr.num_of_cols+7)/8 - (long)(ps)->colset_bytes, 0))

// Used internally only
// Get min tmp memory required to load files
// if we already have paths, return the max, otherwise second max
static size_t path_files_tmp_mem_required(const PathStore *ps,
                                          const PathFileReader *files,
                                          size_t num_files,
                                          bool remove_substr)
{
  bool file0_needs_tmp = (remove_substr || ps->extra_bytes || ps->num_of_paths);
  if(num_files == 0) return 0;
  if(num_files == 1) return file0_needs_tmp ? files[0].hdr.num_path_bytes+1 : 0;

  // We need the size of the first and second largest file path_mem
  size_t i, tmp, s0, s1;

  // make s0 > s1
  s0 = pfile_tmp_mem(files[0],ps);
  s1 = pfile_tmp_mem(files[1],ps);
  if(s1 > s0) { SWAP(s0, s1); }

  for(i = 2; i < num_files; i++) {
    tmp = pfile_tmp_mem(files[i],ps);
    if(tmp > s0) { s1 = s0; s0 = tmp; }
    else if(tmp > s1) { s0 = tmp; }
  }

  // Return largest file size if already have paths
  return file0_needs_tmp ? s0 : s1;
}

#define pfile_mem(pf,cbytes,ebytes) \
  ((pf).hdr.num_path_bytes + \
   (pf).hdr.num_of_paths * \
    (ebytes + (cbytes - ((pf).fltr.ncols+7)/8)))

// Get min memory required to load files. Returns memory required in bytes.
// remove_substr requires extra memory if only loading one file
// (as if loading two files)
size_t path_files_mem_required(const PathFileReader *files, size_t num_files,
                               bool remove_substr, bool use_path_hash,
                               size_t ncols, size_t extra_bytes)
{
  size_t i, tmp, s0, s1;
  if(num_files == 0) return 0;

  // get cbytes
  size_t cbytes = (ncols+7)/8;

  bool file0needs_tmp = (remove_substr || extra_bytes);
  size_t multiplier = (file0needs_tmp ? 2 : 1) * (use_path_hash ? 2 : 1);

  s0 = pfile_mem(files[0],cbytes,extra_bytes);
  if(num_files == 1) return s0 * multiplier;

  // We need the size of the first and second largest file path_mem
  // s0 > s1
  s1 = pfile_mem(files[1],cbytes,extra_bytes);

  if(s1 > s0) { SWAP(s0, s1); }

  for(i = 2; i < num_files; i++) {
    tmp = pfile_mem(files[i],cbytes,extra_bytes);
    if(tmp > s0) { s1 = s0; s0 = tmp; }
    else if(tmp > s1) { s1 = tmp; }
  }

  return (file0needs_tmp ? s0 * 2 : s0 + s1) * (use_path_hash ? 2 : 1);
}

// Print some output
static void paths_loading_print_status(const PathFileReader *file)
{
  const PathFileHeader *hdr = &file->hdr;
  const FileFilter *fltr = &file->fltr;
  char kmers_str[100], paths_str[100], mem_str[100], filesize_str[100];

  ulong_to_str(hdr->num_kmers_with_paths, kmers_str);
  ulong_to_str(hdr->num_of_paths, paths_str);
  bytes_to_str(hdr->num_path_bytes, 1, mem_str);
  bytes_to_str(fltr->file_size, 1, filesize_str);

  file_filter_status(fltr);
  status("  %s paths, %s path-bytes, %s kmers, %s filesize",
         paths_str, mem_str, kmers_str, filesize_str);
}

// Update sample names of the graph using path files
// Only updates colours where sample name has not been set
static void path_files_update_empty_sample_names(const PathFileReader *files,
                                                 size_t num_files,
                                                 dBGraph *db_graph)
{
  size_t i, j, fromcol;
  for(i = 0; i < db_graph->num_of_cols; i++) {
    if(strcmp(db_graph->ginfo[i].sample_name.buff,"undefined") == 0) {
      for(j = 0; j < num_files; j++) {
        if(file_filter_iscolloaded(&files[j].fltr, i)) {
          fromcol = path_file_fromcol(&files[j], i - files[j].fltr.intocol);
          strbuf_set(&db_graph->ginfo[i].sample_name,
                     files[j].hdr.sample_names[fromcol].buff);
          break;
        }
      }
    }
  }
}

// if insert is true, insert missing kmers into the graph
void paths_format_load(PathFileReader *file, bool insert_missing_kmers,
                       dBGraph *db_graph)
{
  const PathFileHeader *hdr = &file->hdr;
  FileFilter *fltr = &file->fltr;
  FILE *fh = fltr->fh;
  const char *path = fltr->file_path.buff;
  PathStore *pstore = &db_graph->pstore;

  // If you want to use a file filter you must use paths_format_merge
  // Check file filter and PathStore are compatible
  ctx_assert(path_store_fltr_compatible(pstore, fltr));

  // Check PathStore has not been used yet
  ctx_assert(pstore->extra_bytes == 0);
  ctx_assert(pstore->next == pstore->store);
  ctx_assert(pstore->num_of_paths == 0 && pstore->num_kmers_with_paths == 0);

  path_files_update_empty_sample_names(file, 1, db_graph);
  path_file_load_check(file, db_graph);

  // Print some output
  paths_loading_print_status(file);

  size_t i;
  BinaryKmer bkmer;
  hkey_t hkey;
  bool found;
  PathIndex pindex;

  // Load paths
  ctx_assert((ptrdiff_t)hdr->num_path_bytes <= pstore->end - pstore->store);
  safe_fread(fh, pstore->store, hdr->num_path_bytes, "pstore->store", path);
  pstore->next = pstore->store + hdr->num_path_bytes;
  pstore->num_of_paths = hdr->num_of_paths;
  pstore->num_kmers_with_paths = hdr->num_kmers_with_paths;
  pstore->num_of_bytes = hdr->num_path_bytes;

  // Load kmer pointers to paths
  for(i = 0; i < hdr->num_kmers_with_paths; i++)
  {
    safe_fread(fh, bkmer.b, sizeof(BinaryKmer), "bkmer", path);

    if(insert_missing_kmers) {
      hkey = hash_table_find_or_insert(&db_graph->ht, bkmer, &found);
    }
    else if((hkey = hash_table_find(&db_graph->ht, bkmer)) == HASH_NOT_FOUND) {
      char kmer_str[MAX_KMER_SIZE+1];
      binary_kmer_to_str(bkmer, db_graph->kmer_size, kmer_str);
      die("Node missing: %s [path: %s]", kmer_str, path);
    }

    safe_fread(fh, &pindex, sizeof(uint64_t), "kmer_index", path);
    if(pindex > hdr->num_path_bytes) {
      die("Path index out of bounds [%zu > %zu]",
          (size_t)pindex, (size_t)hdr->num_path_bytes);
    }

    pstore_set_pindex(pstore, hkey, pindex);
  }

  // Test that this is the end of the file
  uint8_t end;
  if(fread(&end, 1, 1, fh) != 0)
    warn("End of file not reached when loading! [path: %s]", path);
}

void paths_load_colour(PathFileReader *pfile,
                       bool insert_missing_kmers,
                       size_t colour_idx, size_t intocol,
                       dBGraph *db_graph)
{
  ctx_assert(colour_idx < pfile->fltr.ncols);
  ctx_assert(intocol < db_graph->num_of_cols);

  // Copy current values
  FileFilter tmp = pfile->fltr;

  // Set new values
  size_t newcol = pfile->fltr.cols[colour_idx];
  pfile->fltr.cols = &newcol;
  pfile->fltr.ncols = 1;
  file_filter_update_intocol(&pfile->fltr, intocol);

  // Load paths
  paths_format_load(pfile, insert_missing_kmers, db_graph);

  // Restore values
  pfile->fltr = tmp;
}

// pindex is index of last path
static void load_path_set(hkey_t hkey, PathIndex pindex,
                          const PathSet *set, PathStore *ps)
{
  if(set->members.len == 0) return;

  size_t i;
  const PathEntry *entry;

  for(i = 0; i < set->members.len; i++) {
    entry = &set->members.data[i];
    pindex = path_store_add_packed(ps, hkey, pindex, entry->orient, entry->plen,
                                   path_set_colset(entry,set),
                                   path_set_seq(entry,set));
  }

  pstore_set_pindex(ps, hkey, pindex);
}

// set0 is already loaded
// we are considering loading from set1
static void load_linkedlist(hkey_t hkey, PathIndex loadindex,
                            PathFileReader *pfile,
                            PathSet *set0, PathSet *set1,
                            bool rmv_redundant, PathStore *ps)
{
  const FileFilter *fltr = &pfile->fltr;

  // Might not need to use filter
  if(path_store_fltr_compatible(ps,fltr)) fltr = NULL;

  // Create PathSets, set0 from main PathStore, set1 from file (temp store)
  PathIndex pindex = pstore_get_pindex(ps, hkey);
  path_set_init2(set0, ps->colset_bytes, ps->store, pindex, NULL);
  path_set_init2(set1, ps->colset_bytes, ps->tmpstore, loadindex, fltr);

  // Remove redundant paths in paths from the file
  if(rmv_redundant) path_set_slim(set1);

  // Remove elements from set1 (from file) that are already in set0 (loaded)
  path_set_merge(set0, set1, rmv_redundant, ps->store);

  // Store new paths
  load_path_set(hkey, pindex, set1, ps);
}

// Load 1 or more path files; can be called consecutively
// if `rmv_redundant` is true we remove non-informative paths
//  `thread_limit` is the number of threads to use for removing redundant paths
void paths_format_merge(PathFileReader *files, size_t num_files,
                        bool insert_missing_kmers,
                        bool rmv_redundant, size_t thread_limit,
                        dBGraph *db_graph)
{
  if(num_files == 0) return;

  PathStore *pstore = &db_graph->pstore;
  size_t tmp_pmem;

  tmp_pmem = path_files_tmp_mem_required(pstore, files, num_files, rmv_redundant);
  status("[PathFormat] With %zu files, require %zu tmp memory [%zu extra bytes]",
         num_files, tmp_pmem, pstore->extra_bytes);

  // Check number of bytes for colour bitset (path in which cols)
  // This should have been dealt with in the setup of the PathStore
  size_t required_ncols = paths_get_max_usedcols(files, num_files);
  size_t required_nbytes = roundup_bits2bytes(required_ncols);
  ctx_assert(required_ncols <= pstore->num_of_cols);
  ctx_assert(required_nbytes <= pstore->colset_bytes);

  // load files one at a time
  FileFilter *fltr;
  PathFileHeader *hdr;
  FILE *fh;
  const char *path;
  BinaryKmer bkey;
  hkey_t hkey;
  PathIndex tmpindex;
  bool found;
  size_t i, k, first_file = 0;

  // Update sample names of the graph
  path_files_update_empty_sample_names(files, num_files, db_graph);

  for(i = 0; i < num_files; i++)
    path_file_load_check(&files[i], db_graph);

  // Temporary sets used in loading and removing duplicates
  PathSet pset0, pset1;
  path_set_alloc(&pset0);
  path_set_alloc(&pset1);

  // Load first file into main pstore,
  // if no paths loaded yet and no extra bytes padding needed per path
  while(first_file < num_files &&
        pstore->next == pstore->store && pstore->extra_bytes == 0 &&
        path_store_fltr_compatible(pstore, &files[first_file].fltr))
  {
    // Currently no paths loaded
    if(!rmv_redundant)
    {
      paths_format_load(&files[first_file], insert_missing_kmers, db_graph);
      first_file++;
    }
    else if(num_files == 1)
    {
      // Load whole file and remove duplicates
      paths_format_load(&files[first_file], insert_missing_kmers, db_graph);

      // Slim paths store
      graph_paths_clean(db_graph, thread_limit, 0);

      // done
      first_file++;
    }
    else
      break;
  }

  if(tmp_pmem) path_store_setup_tmp(pstore, tmp_pmem);

  for(i = first_file; i < num_files; i++)
  {
    fltr = &files[i].fltr;
    hdr = &files[i].hdr;
    path = fltr->orig_path.buff;
    fh = fltr->fh;

    // Print some output
    paths_loading_print_status(&files[i]);

    ctx_assert(hdr->num_path_bytes == 0 || pstore->tmpstore != NULL);
    ctx_assert(hdr->num_path_bytes <= pstore->tmpsize);
    ctx_assert(!hdr->num_path_bytes == !hdr->num_kmers_with_paths);

    safe_fread(fh, pstore->tmpstore, hdr->num_path_bytes, "paths->store", path);

    // Load kmer pointers to paths
    for(k = 0; k < hdr->num_kmers_with_paths; k++)
    {
      safe_fread(fh, bkey.b, sizeof(BinaryKmer), "bkey", path);

      if(insert_missing_kmers) {
        hkey = hash_table_find_or_insert(&db_graph->ht, bkey, &found);
      }
      else if((hkey = hash_table_find(&db_graph->ht, bkey)) == HASH_NOT_FOUND)
      {
        char kmer_str[MAX_KMER_SIZE+1];
        binary_kmer_to_str(bkey, db_graph->kmer_size, kmer_str);
        die("Node missing: %s [path: %s]", kmer_str, path);
      }

      safe_fread(fh, &tmpindex, sizeof(uint64_t), "kmer_index", path);
      if(tmpindex > hdr->num_path_bytes) {
        die("Path index out of bounds [%zu > %zu]",
            (size_t)tmpindex, (size_t)hdr->num_path_bytes);
      }

      // Merge into currently loaded paths
      load_linkedlist(hkey, tmpindex, &files[i],
                      &pset0, &pset1, rmv_redundant, pstore);
    }

    // Test that this is the end of the file
    uint8_t end;
    if(fread(&end, 1, 1, fh) != 0)
      warn("End of file not reached when loading! [path: %s]", path);
  }

  path_store_print_status(pstore);

  if(tmp_pmem) path_store_release_tmp(pstore);

  path_set_dealloc(&pset0);
  path_set_dealloc(&pset1);
}


//
// Write
//

// returns number of bytes written
size_t paths_format_write_header_core(const PathFileHeader *header, FILE *fout)
{
  size_t mem = fwrite("PATHS", 1, 5, fout) +
               fwrite(&header->version, 1, sizeof(uint32_t), fout) +
               fwrite(&header->kmer_size, 1, sizeof(uint32_t), fout) +
               fwrite(&header->num_of_cols, 1, sizeof(uint32_t), fout) +
               fwrite(&header->num_of_paths, 1, sizeof(uint64_t), fout) +
               fwrite(&header->num_path_bytes, 1, sizeof(uint64_t), fout) +
               fwrite(&header->num_kmers_with_paths, 1, sizeof(uint64_t), fout);

  const size_t expmem = 5 + sizeof(uint32_t)*3 + sizeof(uint64_t)*3;
  if(mem != expmem) die("Couldn't write header core");
  return mem;
}

// returns number of bytes written
size_t paths_format_write_header(const PathFileHeader *header, FILE *fout)
{
  size_t i, bytes, written;
  uint32_t len;
  const StrBuf *buf;

  written = bytes = paths_format_write_header_core(header, fout);

  for(i = 0; i < header->num_of_cols; i++)
  {
    buf = &header->sample_names[i];
    len = (uint32_t)buf->len;
    written += fwrite(&len, 1, sizeof(uint32_t), fout);
    written += fwrite(buf->buff, 1, len, fout);
    bytes += sizeof(uint32_t) + len;
  }

  if(written != bytes) die("Couldn't write header");
  return bytes;
}

static inline void write_optimised_paths(hkey_t hkey, PathIndex *pidx_ptr,
                                         dBGraph *db_graph, FILE *const fout)
{
  const PathStore *pstore = &db_graph->pstore;
  PathIndex pindex, newidx, pidx = *pidx_ptr;
  const uint8_t *path;
  size_t mem;

  ctx_assert(fout != NULL);
  ctx_assert(pstore->kmer_paths_read != NULL);

  pindex = pstore_get_pindex(pstore, hkey);

  // Return if not paths associated with this kmer
  if(pindex == PATH_NULL) return;

  pstore_set_pindex(pstore, hkey, pidx);

  do
  {
    path = pstore->store+pindex;
    mem = packedpath_mem(path, pstore->colset_bytes);
    pindex = packedpath_get_prev(path);

    pidx += mem;
    newidx = (pindex == PATH_NULL ? PATH_NULL : pidx);

    // printf("%zu writing %zu paths bytes next: %zu\n", pidx, mem, newidx);

    if(fwrite(&newidx, 1, sizeof(PathIndex), fout) +
       fwrite(path+sizeof(PathIndex), 1, mem-sizeof(PathIndex), fout) != mem)
    {
      die("Couldn't write to file");
    }
  }
  while(pindex != PATH_NULL);

  *pidx_ptr = pidx;
}

void paths_format_write_optimised_paths_only(dBGraph *db_graph, FILE *fout)
{
  ctx_assert(db_graph->pstore.kmer_paths_read != NULL);
  PathIndex poffset = 0;
  HASH_ITERATE(&db_graph->ht, write_optimised_paths, &poffset, db_graph, fout);
}

static inline void write_kmer_path_indices(hkey_t hkey, const dBGraph *db_graph,
                                           FILE *fout)
{
  const PathStore *pstore = &db_graph->pstore;
  size_t written;
  PathIndex pindex;
  BinaryKmer bkmer;

  if((pindex = pstore_get_pindex(pstore, hkey)) != PATH_NULL)
  {
    bkmer = db_node_get_bkmer(db_graph, hkey);

    written = fwrite(&bkmer, 1, sizeof(BinaryKmer), fout) +
              fwrite(&pindex, 1, sizeof(PathIndex), fout);

    if(written != sizeof(BinaryKmer)+sizeof(PathIndex))
      die("Couldn't write to file");
  }
}

// Corrupts paths so they cannot be used elsewhere
// unless you reload the optimised paths from fout
void paths_format_write_optimised_paths(dBGraph *db_graph, FILE *fout)
{
  ctx_assert(db_graph->pstore.kmer_paths_read != NULL);
  paths_format_write_optimised_paths_only(db_graph, fout);
  HASH_ITERATE(&db_graph->ht, write_kmer_path_indices, db_graph, fout);
}
