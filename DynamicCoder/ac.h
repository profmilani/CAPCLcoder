#ifndef AC_HEADER
#define AC_HEADER

#include <stdio.h>

#define Code_value_bits 16

#define Top_value (((long)1<<Code_value_bits)-1)
#define First_qtr (Top_value/4+1)
#define Half	  (2*First_qtr)
#define Third_qtr (3*First_qtr)
#define Max_frequency 16383

#define error(m)                                           \
do  {                                                      \
  fflush (stdout);                                         \
  fprintf (stderr, "%s:%d: error: ", __FILE__, __LINE__);  \
  fprintf (stderr, m);                                     \
  fprintf (stderr, "\n");                                  \
  exit (1);                                                \
}  while (0)

#define check(b,m)                                         \
do  {                                                      \
  if (b)                                                   \
    error (m);                                             \
}  while (0)

typedef struct {
  FILE *fp;
  long low;
  long high;
  long fbits;
  int buffer;
  int bits_to_go;
  long total_bits;
} ac_encoder;

typedef struct {
  FILE *fp;
  long value;
  long low;
  long high;
  int buffer;
  int bits_to_go;
  int garbage_bits;
} ac_decoder;

typedef struct {
  int nsym;
  int *freq;
  int *cfreq;
  int adapt;
} ac_model;

void ac_encoder_init (ac_encoder *, const char *);
void ac_encoder_done (ac_encoder *);
void ac_decoder_init (ac_decoder *, const char *);
void ac_decoder_done (ac_decoder *);
void ac_model_init (ac_model *, int, int *, int);
void ac_model_done (ac_model *);
long ac_encoder_bits (ac_encoder *);
void ac_encode_symbol (ac_encoder *, ac_model *, int);
int ac_decode_symbol (ac_decoder *, ac_model *);

#endif
