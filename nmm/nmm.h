struct nmm_base;
struct nmm_codont;
struct nmm_frame_state;

struct nmm_codon
{
    char a;
    char b;
    char c;
};

/* Codon table */
struct nmm_codont *nmm_codont_create(struct imm_abc const *abc);
int    nmm_codont_set_lprob(struct nmm_codont *codont, struct nmm_codon const *codon, double lprob);
double nmm_codont_get_lprob(struct nmm_codont const *codont, struct nmm_codon const *codon);
int    nmm_codont_normalize(struct nmm_codont *codont);
void   nmm_codont_destroy(struct nmm_codont *codont);
struct imm_abc const *nmm_codont_get_abc(struct nmm_codont const *codont);

/* Base */
struct nmm_base *     nmm_base_create(struct imm_abc const *abc);
int                   nmm_base_set_lprob(struct nmm_base *base, char nucleotide, double lprob);
double                nmm_base_get_lprob(struct nmm_base const *base, char nucleotide);
int                   nmm_base_normalize(struct nmm_base *base);
void                  nmm_base_destroy(struct nmm_base *base);
struct imm_abc const *nmm_base_get_abc(struct nmm_base const *base);

/* Frame state */
struct nmm_frame_state *nmm_frame_state_create(char const *name, struct nmm_base const *base,
                                               struct nmm_codont const *codont, double epsilon);
double nmm_frame_state_lposterior(struct nmm_frame_state *state, struct nmm_codon const *codon,
                                  char const *seq, int seq_len);
double nmm_frame_state_decode(struct nmm_frame_state *state, char const *seq, int seq_len,
                              struct nmm_codon *codon);
void   nmm_frame_state_destroy(struct nmm_frame_state *state);
