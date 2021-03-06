#include "global.h"
#include "dna.h"

const char dna_nuc_to_char_arr[4] = "ACGT";

// 0:Adenine, 1:Cytosine, 2:Guanine, 3:Thymine, 4:N 8:other
const Nucleotide dna_char_to_nuc_arr[256]
  = {8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,0,8,1,8,8,8,2,8,8,8,8,8,8,4,8, // A C G N
     8,8,8,8,3,8,8,8,8,8,8,8,8,8,8,8, // T
     8,0,8,1,8,8,8,2,8,8,8,8,8,8,4,8, // a c g n
     8,8,8,8,3,8,8,8,8,8,8,8,8,8,8,8, // t
     // 128-256
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
     8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8};

const unsigned char dna_complement_char_arr[256]
  = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
      10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
      20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
      30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
      40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
      50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
      60, 61, 62, 63, 64,'T','B','G','D','E', // A C
     'F','C','H','I','J','K','L','M','N','O', // G N
     'P','Q','R','S','A','U','V','W','X','Y', // T
     'Z', 91, 92, 93, 94, 95, 96,'t','b','g', // a c
     'd','e','f','c','h','i','j','k','l','m', // g
     'n','o','p','q','r','s','a','u','v','w', // n t
     'x','y','z',123,124,125,126,127,128,129,
     130,131,132,133,134,135,136,137,138,139,
     140,141,142,143,144,145,146,147,148,149,
     150,151,152,153,154,155,156,157,158,159,
     160,161,162,163,164,165,166,167,168,169,
     170,171,172,173,174,175,176,177,178,179,
     180,181,182,183,184,185,186,187,188,189,
     190,191,192,193,194,195,196,197,198,199,
     200,201,202,203,204,205,206,207,208,209,
     210,211,212,213,214,215,216,217,218,219,
     220,221,222,223,224,225,226,227,228,229,
     230,231,232,233,234,235,236,237,238,239,
     240,241,242,243,244,245,246,247,248,249,
     250,251,252,253,254,255};

// length is the length in number of bases
// the char* should have one MORE base than that allocated, to hold '\0'
char* dna_reverse_complement_str(char *str, size_t length)
{
  ctx_assert(strlen(str) >= length);

  if(length == 0) return str;
  if(length == 1) { str[0] = dna_char_complement(str[0]); return str; }

  size_t i, j;
  for(i = 0, j = length - 1; i <= j; i++, j--)
  {
    char tmp = str[i];
    str[i] = dna_char_complement(str[j]);
    str[j] = dna_char_complement(tmp);
  }

  return str;
}

// Generate a random dna str "ACGT" of length cap-1, terminated with a \0 at
// position cap-1. If cap is 0, does nothing. Useful for testing
char* dna_rand_str(char *str, size_t cap)
{
  const char bases[4] = "ACGT";
  size_t i, r = 0;

  if(cap == 0) return str;
  if(cap == 1) { str[0] = '\0'; return str; }

  for(i = 0; i < cap-1; i++) {
    if((i & 15) == 0) r = (size_t)rand(); // 2 bits per cycle, 32 bits in rand()
    str[i] = bases[r&3];
    r >>= 2;
  }

  str[cap-1] = '\0';

  return str;
}
