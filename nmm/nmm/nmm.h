struct nmm_base;
struct nmm_codon;
struct nmm_baset;
struct nmm_codonp;
struct nmm_codont;
struct nmm_frame_state;

struct nmm_triplet
{
    char a;
    char b;
    char c;
};

/* Base */
struct nmm_base const *nmm_base_create(struct imm_abc const *abc);
void                   nmm_base_destroy(struct nmm_base const *base);
struct imm_abc const * nmm_base_get_abc(struct nmm_base const *base);

/* Base table */
struct nmm_baset const *nmm_baset_create(struct nmm_base const *base, double a, double b, double c,
                                         double d);
double                  nmm_baset_lprob(struct nmm_baset const *baset, char nucleotide);
void                    nmm_baset_destroy(struct nmm_baset const *baset);
struct nmm_base const * nmm_baset_get_base(struct nmm_baset const *baset);

/* Codon */
struct nmm_codon *     nmm_codon_create(struct nmm_base const *base);
void                   nmm_codon_destroy(struct nmm_codon const *codon);
struct nmm_base const *nmm_codon_get_base(struct nmm_codon const *codon);
int                    nmm_codon_set_triplet(struct nmm_codon *codon, struct nmm_triplet triplet);
struct nmm_triplet     nmm_codon_get_triplet(struct nmm_codon const *codon);

/* Codon table */
struct nmm_codont const *nmm_codont_create(struct nmm_codonp const *codonp);
double nmm_codont_lprob(struct nmm_codont const *codont, struct nmm_codon const *codon);
void   nmm_codont_destroy(struct nmm_codont const *codont);
struct nmm_base const *nmm_codont_get_base(struct nmm_codont const *codont);

/* Codon probability */
struct nmm_codonp *nmm_codonp_create(struct nmm_base const *base);
int    nmm_codonp_set_lprob(struct nmm_codonp *codonp, struct nmm_codon const *codon, double lprob);
double nmm_codonp_get_lprob(struct nmm_codonp const *codonp, struct nmm_codon const *codon);
int    nmm_codonp_normalize(struct nmm_codonp *codonp);
void   nmm_codonp_destroy(struct nmm_codonp const *codonp);
struct nmm_base const *nmm_codonp_get_base(struct nmm_codonp const *codonp);

/* Frame state */
struct nmm_frame_state const *nmm_frame_state_create(char const *name, struct nmm_baset const *baset,
                                                     struct nmm_codont const *codont, double epsilon);
double nmm_frame_state_lposterior(struct nmm_frame_state const *state, struct nmm_codon const *codon,
                                  struct imm_seq const *seq);
double nmm_frame_state_decode(struct nmm_frame_state const *state, struct imm_seq const *seq,
                              struct nmm_codon *codon);
void   nmm_frame_state_destroy(struct nmm_frame_state const *state);
