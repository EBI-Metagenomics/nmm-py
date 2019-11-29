struct nmm_baset;
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

/* Base table */
struct nmm_baset *    nmm_baset_create(struct imm_abc const *abc);
int                   nmm_baset_set_lprob(struct nmm_baset *baset, char nucleotide, double lprob);
double                nmm_baset_get_lprob(struct nmm_baset const *baset, char nucleotide);
int                   nmm_baset_normalize(struct nmm_baset *baset);
void                  nmm_baset_destroy(struct nmm_baset *baset);
struct imm_abc const *nmm_baset_get_abc(struct nmm_baset const *baset);

/* Frame state */
struct nmm_frame_state *nmm_frame_state_create(char const *name, struct nmm_baset const *baset,
                                               struct nmm_codont const *codont, double epsilon);
double nmm_frame_state_lposterior(struct nmm_frame_state *state, struct nmm_codon const *codon,
                                  char const *seq, int seq_len);
double nmm_frame_state_decode(struct nmm_frame_state *state, char const *seq, int seq_len,
                              struct nmm_codon *codon);
void   nmm_frame_state_destroy(struct nmm_frame_state *state);
