struct nmm_base;
struct nmm_codon;
struct nmm_frame_state;

struct nmm_ccode
{
    char a;
    char b;
    char c;
};

/* Codon */
struct nmm_codon *nmm_codon_create(struct imm_abc const *abc);
int    nmm_codon_set_lprob(struct nmm_codon *codon, struct nmm_ccode const *ccode, double lprob);
double nmm_codon_get_lprob(struct nmm_codon const *codon, struct nmm_ccode const *ccode);
int    nmm_codon_normalize(struct nmm_codon *codon);
void   nmm_codon_destroy(struct nmm_codon *codon);
struct imm_abc const *nmm_codon_get_abc(struct nmm_codon const *codon);

/* Base */
struct nmm_base *     nmm_base_create(struct imm_abc const *abc);
int                   nmm_base_set_lprob(struct nmm_base *base, char nucleotide, double lprob);
double                nmm_base_get_lprob(struct nmm_base const *base, char nucleotide);
int                   nmm_base_normalize(struct nmm_base *base);
void                  nmm_base_destroy(struct nmm_base *base);
struct imm_abc const *nmm_base_get_abc(struct nmm_base const *base);

/* Frame state */
struct nmm_frame_state *nmm_frame_state_create(char const *name, struct nmm_base const *base,
                                               struct nmm_codon const *codon, double epsilon);
double nmm_frame_state_lposterior(struct nmm_frame_state *state, struct nmm_ccode const *ccode,
                                  char const *seq, int seq_len);
double nmm_frame_state_decode(struct nmm_frame_state *state, char const *seq, int seq_len,
                              struct nmm_ccode *ccode);
void   nmm_frame_state_destroy(struct nmm_frame_state *state);
