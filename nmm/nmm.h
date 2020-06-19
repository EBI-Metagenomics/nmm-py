struct imm_abc;
struct imm_dp;
struct imm_hmm;
struct imm_seq;
struct imm_state;
struct nmm_amino_abc;
struct nmm_amino_table;
struct nmm_base_abc;
struct nmm_base_table;
struct nmm_codon;
struct nmm_codon_lprob;
struct nmm_codon_state;
struct nmm_codon_table;
struct nmm_frame_state;
struct nmm_input;
struct nmm_model;
struct nmm_output;

#define NMM_CODON_STATE_TYPE_ID 0x10
#define NMM_FRAME_STATE_TYPE_ID 0x11

struct nmm_triplet
{
    char a;
    char b;
    char c;
};

/* Amino abc */
struct nmm_amino_abc const *nmm_amino_abc_create(char const *symbols, char const any_symbol);
struct nmm_amino_abc const *nmm_amino_abc_derived(struct imm_abc const *abc);
void                        nmm_amino_abc_destroy(struct nmm_amino_abc const *amino_abc);
struct imm_abc const *      nmm_amino_abc_super(struct nmm_amino_abc const *amino_abc);

/* Amino table */
struct nmm_amino_table const *nmm_amino_table_create(struct nmm_amino_abc const *abc,
                                                     double const *              lprobs);
void                          nmm_amino_table_destroy(struct nmm_amino_table const *tbl);
struct nmm_amino_abc const *  nmm_amino_table_abc(struct nmm_amino_table const *tbl);
double nmm_amino_table_lprob(struct nmm_amino_table const *tbl, char const amino);

/* Base abc */
struct nmm_base_abc const *nmm_base_abc_create(char const *symbols, char const any_symbol);
struct nmm_base_abc const *nmm_base_abc_derived(struct imm_abc const *abc);
void                       nmm_base_abc_destroy(struct nmm_base_abc const *base_abc);
struct imm_abc const *     nmm_base_abc_super(struct nmm_base_abc const *base_abc);

/* Base table */
struct nmm_base_abc const *  nmm_base_table_abc(struct nmm_base_table const *baset);
struct nmm_base_table const *nmm_base_table_create(struct nmm_base_abc const *abc, double a, double b,
                                                   double c, double d);
void                         nmm_base_table_destroy(struct nmm_base_table const *baset);
double nmm_base_table_lprob(struct nmm_base_table const *baset, char const base);

/* Codon */
struct nmm_base_abc const *nmm_codon_abc(struct nmm_codon const *codon);
struct nmm_codon *         nmm_codon_create(struct nmm_base_abc const *base);
void                       nmm_codon_destroy(struct nmm_codon const *codon);
struct nmm_triplet         nmm_codon_get_triplet(struct nmm_codon const *codon);
void                       nmm_codon_init(struct nmm_codon *codon, struct nmm_base_abc const *base);
int                        nmm_codon_set_triplet(struct nmm_codon *codon, struct nmm_triplet triplet);

/* Codon probability */
struct nmm_base_abc const *nmm_codon_lprob_abc(struct nmm_codon_lprob const *codonp);
struct nmm_codon_lprob *   nmm_codon_lprob_create(struct nmm_base_abc const *abc);
void                       nmm_codon_lprob_destroy(struct nmm_codon_lprob const *codonp);
double nmm_codon_lprob_get(struct nmm_codon_lprob const *codonp, struct nmm_codon const *codon);
int    nmm_codon_lprob_normalize(struct nmm_codon_lprob *codonp);
int nmm_codon_lprob_set(struct nmm_codon_lprob *codonp, struct nmm_codon const *codon, double lprob);

/* Codon state */
struct nmm_codon_lprob const *nmm_codon_state_codon_lprob(struct nmm_codon_state const *state);
struct nmm_codon_state const *nmm_codon_state_create(char const *                  name,
                                                     struct nmm_codon_lprob const *codonp);
struct nmm_codon_state const *nmm_codon_state_derived(struct imm_state const *state);
void                          nmm_codon_state_destroy(struct nmm_codon_state const *state);
struct imm_state const *      nmm_codon_state_read(FILE *stream, struct nmm_model const *model);
struct imm_state const *      nmm_codon_state_super(struct nmm_codon_state const *state);
int nmm_codon_state_write(struct imm_state const *state, struct nmm_model const *model, FILE *stream);

/* Codon table */
struct nmm_base_abc const *   nmm_codon_table_abc(struct nmm_codon_table const *codont);
struct nmm_codon_table const *nmm_codon_table_create(struct nmm_codon_lprob const *prob);
void                          nmm_codon_table_destroy(struct nmm_codon_table const *codont);
double nmm_codon_table_lprob(struct nmm_codon_table const *codont, struct nmm_codon const *codon);

/* Frame state */
struct nmm_base_table const * nmm_frame_state_base_table(struct nmm_frame_state const *state);
struct nmm_codon_table const *nmm_frame_state_codon_table(struct nmm_frame_state const *state);
struct nmm_frame_state const *nmm_frame_state_create(char const *                  name,
                                                     struct nmm_base_table const * baset,
                                                     struct nmm_codon_table const *codont,
                                                     double                        epsilon);
double nmm_frame_state_decode(struct nmm_frame_state const *state, struct imm_seq const *seq,
                              struct nmm_codon *codon);
struct nmm_frame_state const *nmm_frame_state_derived(struct imm_state const *state);
void                          nmm_frame_state_destroy(struct nmm_frame_state const *state);
double                        nmm_frame_state_epsilon(struct nmm_frame_state const *state);
double nmm_frame_state_lposterior(struct nmm_frame_state const *state, struct nmm_codon const *codon,
                                  struct imm_seq const *seq);
struct imm_state const *nmm_frame_state_read(FILE *stream, struct nmm_model const *model);
struct imm_state const *nmm_frame_state_super(struct nmm_frame_state const *state);
int nmm_frame_state_write(struct imm_state const *state, struct nmm_model const *model, FILE *stream);

/* Input */
int                     nmm_input_close(struct nmm_input *input);
struct nmm_input *      nmm_input_create(char const *filepath);
int                     nmm_input_destroy(struct nmm_input *input);
bool                    nmm_input_eof(struct nmm_input const *input);
int                     nmm_input_fseek(struct nmm_input *input, long offset);
long                    nmm_input_ftell(struct nmm_input *input);
struct nmm_model const *nmm_input_read(struct nmm_input *input);

/* Model */
struct imm_abc const *        nmm_model_abc(struct nmm_model const *model);
struct nmm_base_table const * nmm_model_base_table(struct nmm_model const *model, uint32_t index);
struct nmm_codon_lprob const *nmm_model_codon_lprob(struct nmm_model const *model, uint32_t index);
struct nmm_codon_table const *nmm_model_codon_table(struct nmm_model const *model, uint32_t index);
struct nmm_model const *      nmm_model_create(struct imm_hmm *hmm, struct imm_dp const *dp);
void                          nmm_model_destroy(struct nmm_model const *model);
uint32_t                      nmm_model_nbase_tables(struct nmm_model const *model);
uint32_t                      nmm_model_ncodon_lprobs(struct nmm_model const *model);
uint32_t                      nmm_model_ncodon_tables(struct nmm_model const *model);
struct imm_hmm *              nmm_model_hmm(struct nmm_model const *model);
struct imm_dp const *         nmm_model_dp(struct nmm_model const *model);
struct imm_state const *      nmm_model_state(struct nmm_model const *model, uint32_t i);
uint32_t                      nmm_model_nstates(struct nmm_model const *model);
struct nmm_model const *      nmm_model_read(FILE *stream);
int                           nmm_model_write(struct nmm_model const *io, FILE *stream);

/* Output */
int                nmm_output_close(struct nmm_output *output);
struct nmm_output *nmm_output_create(char const *filepath);
int                nmm_output_destroy(struct nmm_output *output);
int                nmm_output_write(struct nmm_output *output, struct nmm_model const *model);

/* Triplet */
struct nmm_triplet NMM_TRIPLET(char a, char b, char c);
