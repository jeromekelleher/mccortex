#ifndef VCF_PARSING_H_
#define VCF_PARSING_H_

#include "string_buffer/string_buffer.h"

// VCF field indices
#define VCFCHROM   0
#define VCFPOS     1
#define VCFID      2
#define VCFREF     3
#define VCFALT     4
#define VCFQUAL    5
#define VCFFILTER  6
#define VCFINFO    7
#define VCFFORMAT  8
#define VCFSAMPLES 9

typedef struct
{
  StrBuf *cols;
  // Unpacked into:
  StrBuf *alts;
  size_t num_alts, alts_capacity;
  StrBuf *info;
  size_t num_info, info_capacity;
  StrBuf *lf, *rf; // flanks
} vcf_entry_t;

void vcf_entry_alloc(vcf_entry_t *entry, size_t num_samples);
void vcf_entry_dealloc(vcf_entry_t *entry, size_t num_samples);

void vcf_entry_alt_capacity(vcf_entry_t *entry, size_t num_alts);
void vcf_entry_info_capacity(vcf_entry_t *entry, size_t num_info);

void vcf_entry_cpy(vcf_entry_t *dst, const vcf_entry_t *src, size_t num_samples);
void vcf_entry_parse(StrBuf *line, vcf_entry_t *entry, size_t num_samples);
void vcf_entry_revcmp(vcf_entry_t *entry);
size_t vcf_entry_longest_allele(const vcf_entry_t *entry);
void vcf_entry_print(const vcf_entry_t *entry, FILE *out, size_t num_samples);
void vcf_entry_add_filter(vcf_entry_t *entry, const char *status);

// DP=asdf;TXT="a;b;c";F1;X=4
// ends:  ^           ^  ^   ^
char* vcf_info_tag_end(char *str);
StrBuf* vcf_info_tag_find(vcf_entry_t *entry, const char* str);
void vcf_info_tag_add(vcf_entry_t *entry, const char* fmt, ...)
__attribute__((format(printf, 2, 3)));
void vcf_info_tag_del(vcf_entry_t *entry, const char* tag);


#endif /* VCF_PARSING_H_ */
