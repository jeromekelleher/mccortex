#include "global.h"
#include "file_util.h"

#include <libgen.h> // dirname

//
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
//

// Returns 0 on success, -1 on failure
static int do_mkdir(const char *path, mode_t mode)
{
  struct stat st;

  // mkdir returns zero on success
  //               nonzero on error, setting errno to EEXIST if dir already exists
  // stat returns nonzero if cannot stat file
  if((mkdir(path, mode) != 0 && errno != EEXIST) ||
     (stat(path, &st) != 0 || !S_ISDIR(st.st_mode)))
    return -1;

  return 0;
}

/**
** Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
** Returns 0 on success, -1 on failure
*/
int futil_mkpath(const char *path, mode_t mode)
{
  char *pp;
  char *sp;
  int status;
  char *copypath = strdup(path);

  status = 0;
  pp = copypath;
  while (status == 0 && (sp = strchr(pp, '/')) != 0)
  {
    if (sp != pp)
    {
      /* Neither root nor double slash in path */
      *sp = '\0';
      status = do_mkdir(copypath, mode);
      *sp = '/';
    }
    pp = sp + 1;
  }
  if(status == 0)
    status = do_mkdir(path, mode);
  free(copypath);
  return status;
}

//
//
//

bool futil_file_exists(const char *file)
{
  return (access(file, F_OK) != -1);
}

bool futil_is_file_readable(const char *file)
{
  FILE *fp = fopen(file, "r");
  if(fp == NULL) return false;
  else fclose(fp);
  return true;
}

// Creates file if it can write
bool futil_is_file_writable(const char *file)
{
  FILE *fp = fopen(file, "a");
  if(fp == NULL) return false;
  else fclose(fp);
  return true;
}

// Returns -1 on failure
off_t futil_get_file_size(const char* filepath)
{
  struct stat st;

  if (stat(filepath, &st) == 0)
      return st.st_size;

  warn("Cannot determine size of %s: %s\n", filepath, strerror(errno));

  return -1;
}

// Open and return output file.
// If "-" return stdout, if cannot open die with error message
FILE* futil_open_output(const char *path)
{
  FILE *fout;

  if(strcmp(path, "-") == 0)
    fout = stdout;
  else if(futil_file_exists(path))
    die("Output file already exists: %s", path);
  else
    fout = fopen(path, "w");

  if(fout == NULL)
    die("Cannot open output file: %s", futil_outpath_str(path));

  return fout;
}

// Open and return gzip output file.
// If "-" return stdout, if cannot open die with error message
gzFile futil_gzopen_output(const char *path)
{
  ctx_assert(path != NULL);
  gzFile gzout;

  if(strcmp(path, "-") == 0)
    gzout = gzdopen(fileno(stdout), "w");
  else if(futil_file_exists(path))
    die("Output file already exists: %s", path);
  else
    gzout = gzopen(path, "w");

  if(gzout == NULL)
    die("Cannot open output file: %s", futil_outpath_str(path));

  return gzout;
}

bool futil_generate_filename(const char *base_fmt, StrBuf *str)
{
  int i;

  for(i = 0; i < 10000; i++)
  {
    strbuf_reset(str);
    strbuf_sprintf(str, base_fmt, i);
    struct stat st;
    if(stat(str->buff, &st) != 0) return true;
  }

  return false;
}

// Remember to free the result
void futil_get_strbuf_of_dir_path(const char *path, StrBuf *dir)
{
  char *tmp = strdup(path);
  strbuf_set(dir, dirname(tmp));
  strbuf_append_char(dir, '/');
  ctx_free(tmp);
}

char* futil_get_current_dir(char abspath[PATH_MAX+1])
{
  char cwd[PATH_MAX + 1];
  if(getcwd(cwd, PATH_MAX + 1) != NULL)
    return realpath(cwd, abspath);
  else
    return NULL;
}

// Usage:
//     FILE **tmp_files = futil_create_tmp_files(num_tmp);
// to clear up:
//     for(i = 0; i < num_tmp; i++) fclose(tmp_files[i]);
//     ctx_free(tmp_files);
FILE** futil_create_tmp_files(size_t num_tmp_files)
{
  size_t i;
  StrBuf tmppath;
  strbuf_alloc(&tmppath, 1024);
  FILE **tmp_files = ctx_malloc(num_tmp_files * sizeof(FILE*));

  int r = rand() & ((1<<20)-1);

  for(i = 0; i < num_tmp_files; i++)
  {
    strbuf_reset(&tmppath);
    strbuf_sprintf(&tmppath, "/tmp/cortex.tmp.%i.%zu", r, i);
    if((tmp_files[i] = fopen(tmppath.buff, "r+")) == NULL) {
      die("Cannot write temporary file: %s [%s]", tmppath.buff, strerror(errno));
    }
    unlink(tmppath.buff); // Immediately unlink to hide temp file
  }
  strbuf_dealloc(&tmppath);

  return tmp_files;
}

// Merge temporary files, closes tmp files
void futil_merge_tmp_files(FILE **tmp_files, size_t num_files, FILE *fout)
{
  #define TMP_BUF_SIZE (1<<25) /* 32MB */

  char *data = ctx_malloc(TMP_BUF_SIZE);
  size_t i, len;
  FILE *tmp_file;

  for(i = 0; i < num_files; i++)
  {
    tmp_file = tmp_files[i];
    if(fseek(tmp_file, 0L, SEEK_SET) == -1) die("gzseek error");

    while((len = fread(data, 1, TMP_BUF_SIZE, tmp_file)) > 0)
      if(fwrite(data, 1, len, fout) != len)
        die("write error [%s]", strerror(errno));

    fclose(tmp_file);
  }

  ctx_free(data);

  #undef TMP_BUF_SIZE
}
