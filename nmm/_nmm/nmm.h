struct nmm_base_abc;
struct nmm_codon;
struct nmm_base_table;
struct nmm_codon_lprob;
struct nmm_codon_table;
struct nmm_frame_state;

struct nmm_triplet
{
    char a;
    char b;
    char c;
};

/* Base abc */
struct nmm_base_abc const *nmm_base_abc_create(struct imm_abc const *abc);
void                       nmm_base_abc_destroy(struct nmm_base_abc const *base);
struct imm_abc const *     nmm_base_abc_cast(struct nmm_base_abc const *base);

/* Base table */
struct nmm_base_table const *nmm_base_table_create(struct nmm_base_abc const *base, double a,
                                                   double b, double c, double d);
double                     nmm_base_table_lprob(struct nmm_base_table const *baset, char nucleotide);
void                       nmm_base_table_destroy(struct nmm_base_table const *baset);
struct nmm_base_abc const *nmm_base_table_get_base(struct nmm_base_table const *baset);

/* Amino abc */
struct nmm_amino_abc const *nmm_amino_abc_create(struct imm_abc const *abc);
void                        nmm_amino_abc_destroy(struct nmm_amino_abc const *amino_abc);
struct imm_abc const *      nmm_amino_abc_cast(struct nmm_amino_abc const *amino_abc);

/* Codon */
struct nmm_codon *         nmm_codon_create(struct nmm_base_abc const *base);
void                       nmm_codon_destroy(struct nmm_codon const *codon);
struct nmm_base_abc const *nmm_codon_get_base(struct nmm_codon const *codon);
int                        nmm_codon_set_triplet(struct nmm_codon *codon, struct nmm_triplet triplet);
struct nmm_triplet         nmm_codon_get_triplet(struct nmm_codon const *codon);

/* Codon table */
struct nmm_codon_table const *nmm_codon_table_create(struct nmm_codon_lprob const *codonp);
double nmm_codon_table_lprob(struct nmm_codon_table const *codont, struct nmm_codon const *codon);
void   nmm_codon_table_destroy(struct nmm_codon_table const *codont);
struct nmm_base_abc const *nmm_codon_table_get_base(struct nmm_codon_table const *codont);

/* Codon probability */
struct nmm_codon_lprob *nmm_codon_lprob_create(struct nmm_base_abc const *base);
int nmm_codon_lprob_set(struct nmm_codon_lprob *codonp, struct nmm_codon const *codon, double lprob);
double nmm_codon_lprob_get(struct nmm_codon_lprob const *codonp, struct nmm_codon const *codon);
int    nmm_codon_lprob_normalize(struct nmm_codon_lprob *codonp);
void   nmm_codon_lprob_destroy(struct nmm_codon_lprob const *codonp);
struct nmm_base_abc const *nmm_codon_lprob_get_base_abc(struct nmm_codon_lprob const *codonp);

/* Frame state */
struct nmm_frame_state const *nmm_frame_state_create(char const *                  name,
                                                     struct nmm_base_table const * baset,
                                                     struct nmm_codon_table const *codont,
                                                     double                        epsilon);
double nmm_frame_state_lposterior(struct nmm_frame_state const *state, struct nmm_codon const *codon,
                                  struct imm_seq const *seq);
double nmm_frame_state_decode(struct nmm_frame_state const *state, struct imm_seq const *seq,
                              struct nmm_codon *codon);
void   nmm_frame_state_destroy(struct nmm_frame_state const *state);
