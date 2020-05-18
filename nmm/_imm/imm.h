struct imm_abc;
struct imm_abc_table;
struct imm_dp;
struct imm_hmm;
struct imm_input;
struct imm_model;
struct imm_mute_state;
struct imm_normal_state;
struct imm_output;
struct imm_path;
struct imm_result;
struct imm_results;
struct imm_seq_table;
struct imm_state;
struct imm_step;
struct imm_table_state;
struct imm_window;

#define IMM_MUTE_STATE_TYPE_ID 0x00
#define IMM_NORMAL_STATE_TYPE_ID 0x01
#define IMM_TABLE_STATE_TYPE_ID 0x02

enum imm_symbol_type
{
    IMM_SYMBOL_UNKNOWN = 0,
    IMM_SYMBOL_NORMAL = 1,
    IMM_SYMBOL_ANY = 2,
};

struct imm_seq
{
    struct imm_abc const *abc;
    char const *          string;
    unsigned              length;
};

struct imm_subseq
{
    struct imm_seq const *super;
    struct imm_seq        seq;
};

/* Alphabet */
char                  imm_abc_any_symbol(struct imm_abc const *abc);
struct imm_abc const *imm_abc_clone(struct imm_abc const *abc);
struct imm_abc const *imm_abc_create(char const *symbols, char any_symbol);
void                  imm_abc_destroy(struct imm_abc const *abc);
bool                  imm_abc_has_symbol(struct imm_abc const *abc, char symbol_id);
uint8_t               imm_abc_length(struct imm_abc const *abc);
struct imm_abc *      imm_abc_read(FILE *stream);
char                  imm_abc_symbol_id(struct imm_abc const *abc, uint8_t symbol_idx);
uint8_t               imm_abc_symbol_idx(struct imm_abc const *abc, char symbol_id);
enum imm_symbol_type  imm_abc_symbol_type(struct imm_abc const *abc, char symbol_id);
char const *          imm_abc_symbols(struct imm_abc const *abc);
uint8_t               imm_abc_type_id(struct imm_abc const *abc);
int                   imm_abc_write(struct imm_abc const *abc, FILE *stream);

/* Alphabet table */
struct imm_abc_table const *imm_abc_table_create(struct imm_abc const *abc, double const *lprobs);
void                        imm_abc_table_destroy(struct imm_abc_table const *abc_table);
struct imm_abc const *      imm_abc_table_abc(struct imm_abc_table const *abc_table);
double                      imm_abc_table_lprob(struct imm_abc_table const *abc_table, char symbol);

/* DP */
void                      imm_dp_destroy(struct imm_dp const *dp);
struct imm_results const *imm_dp_viterbi(struct imm_dp const *dp, struct imm_seq const *seq,
                                         unsigned window_length);

/* HMM */
int imm_hmm_add_state(struct imm_hmm *hmm, struct imm_state const *state, double start_lprob);
struct imm_hmm *     imm_hmm_create(struct imm_abc const *abc);
struct imm_dp const *imm_hmm_create_dp(struct imm_hmm const *hmm, struct imm_state const *end_state);
int                  imm_hmm_del_state(struct imm_hmm *hmm, struct imm_state const *state);
void                 imm_hmm_destroy(struct imm_hmm const *hmm);
double               imm_hmm_get_trans(struct imm_hmm const *hmm, struct imm_state const *src_state,
                                       struct imm_state const *tgt_state);
double               imm_hmm_likelihood(struct imm_hmm const *hmm, struct imm_seq const *seq,
                                        struct imm_path const *path);
int                  imm_hmm_normalize(struct imm_hmm *hmm);
int                  imm_hmm_normalize_start(struct imm_hmm *hmm);
int                  imm_hmm_normalize_trans(struct imm_hmm *hmm, struct imm_state const *src_state);
double               imm_hmm_get_start(struct imm_hmm const *hmm, struct imm_state const *state);
int imm_hmm_set_start(struct imm_hmm *hmm, struct imm_state const *state, double lprob);
int imm_hmm_set_trans(struct imm_hmm *hmm, struct imm_state const *src_state,
                      struct imm_state const *tgt_state, double lprob);

/* Input */
int                     imm_input_close(struct imm_input *input);
struct imm_input *      imm_input_create(char const *filepath);
int                     imm_input_destroy(struct imm_input *input);
bool                    imm_input_eof(struct imm_input const *input);
struct imm_model const *imm_input_read(struct imm_input *input);

/* Model */
struct imm_abc const *  imm_model_abc(struct imm_model const *model);
struct imm_model *      imm_model_create(struct imm_hmm *hmm, struct imm_dp const *dp);
void                    imm_model_destroy(struct imm_model const *model);
struct imm_dp const *   imm_model_dp(struct imm_model const *model);
struct imm_hmm *        imm_model_hmm(struct imm_model const *model);
uint32_t                imm_model_nstates(struct imm_model const *model);
struct imm_model const *imm_model_read(FILE *stream);
struct imm_state const *imm_model_state(struct imm_model const *model, uint32_t i);
int                     imm_model_write(struct imm_model const *model, FILE *stream);

/* Mute state */
struct imm_mute_state const *imm_mute_state_create(char const *name, struct imm_abc const *abc);
struct imm_mute_state const *imm_mute_state_derived(struct imm_state const *state);
void                         imm_mute_state_destroy(struct imm_mute_state const *state);
struct imm_state const *     imm_mute_state_read(FILE *stream, struct imm_abc const *abc);
struct imm_state const *     imm_mute_state_super(struct imm_mute_state const *state);
int imm_mute_state_write(struct imm_state const *state, struct imm_model const *model, FILE *stream);

/* Normal state */
struct imm_normal_state const *imm_normal_state_create(char const *name, struct imm_abc const *abc,
                                                       double const *lprobs);
struct imm_normal_state const *imm_normal_state_derived(struct imm_state const *state);
void                           imm_normal_state_destroy(struct imm_normal_state const *state);
struct imm_state const *       imm_normal_state_read(FILE *stream, struct imm_abc const *abc);
struct imm_state const *       imm_normal_state_super(struct imm_normal_state const *state);
int imm_normal_state_write(struct imm_state const *state, struct imm_model const *model,
                           FILE *stream);

/* Output */
int                imm_output_close(struct imm_output *output);
struct imm_output *imm_output_create(char const *filepath);
int                imm_output_destroy(struct imm_output *output);
int                imm_output_write(struct imm_output *output, struct imm_model const *model);

/* Path */
void                   imm_path_append(struct imm_path *path, struct imm_step *step);
struct imm_path *      imm_path_clone(struct imm_path const *path);
struct imm_path *      imm_path_create(void);
void                   imm_path_destroy(struct imm_path const *path);
struct imm_step const *imm_path_first(struct imm_path const *path);
void                   imm_path_free(struct imm_path const *path);
struct imm_step const *imm_path_next(struct imm_path const *path, struct imm_step const *step);
void                   imm_path_prepend(struct imm_path *path, struct imm_step *step);

/* Result */
void                   imm_result_destroy(struct imm_result const *result);
void                   imm_result_free(struct imm_result const *result);
double                 imm_result_loglik(struct imm_result const *result);
struct imm_path const *imm_result_path(struct imm_result const *result);
struct imm_subseq      imm_result_subseq(struct imm_result const *result);

/* Results */
struct imm_results *     imm_results_create(struct imm_seq const *seq, unsigned nresults);
void                     imm_results_destroy(struct imm_results const *results);
void                     imm_results_free(struct imm_results const *results);
struct imm_result const *imm_results_get(struct imm_results const *results, unsigned idx);
void     imm_results_set(struct imm_results *results, unsigned idx, struct imm_subseq subseq,
                         struct imm_path const *path, double loglik);
unsigned imm_results_size(struct imm_results const *results);

/* Sequence */
struct imm_seq        IMM_SEQ(struct imm_abc const *abc, char const *string, unsigned length);
struct imm_seq const *imm_seq_clone(struct imm_seq const *seq);
struct imm_seq const *imm_seq_create(char const *string, struct imm_abc const *abc);
void                  imm_seq_destroy(struct imm_seq const *seq);
struct imm_abc const *imm_seq_get_abc(struct imm_seq const *seq);
unsigned              imm_seq_length(struct imm_seq const *seq);
char const *          imm_seq_string(struct imm_seq const *seq);

/* Sequence table */
struct imm_abc const *imm_seq_table_abc(struct imm_seq_table const *table);
int imm_seq_table_add(struct imm_seq_table *table, struct imm_seq const *seq, double lprob);
struct imm_seq_table *imm_seq_table_clone(struct imm_seq_table const *table);
struct imm_seq_table *imm_seq_table_create(struct imm_abc const *abc);
void                  imm_seq_table_destroy(struct imm_seq_table const *table);
double  imm_seq_table_lprob(struct imm_seq_table const *table, struct imm_seq const *seq);
uint8_t imm_seq_table_max_seq(struct imm_seq_table const *table);
uint8_t imm_seq_table_min_seq(struct imm_seq_table const *table);
int     imm_seq_table_normalize(struct imm_seq_table *table);

/* State */
struct imm_state const *imm_state_create(char const *name, struct imm_abc const *abc,
                                         struct imm_state_vtable vtable, void *derived);
void                    imm_state_destroy(struct imm_state const *state);
struct imm_abc const *  imm_state_get_abc(struct imm_state const *state);
char const *            imm_state_get_name(struct imm_state const *state);
double                  imm_state_lprob(struct imm_state const *state, struct imm_seq const *seq);
uint8_t                 imm_state_max_seq(struct imm_state const *state);
uint8_t                 imm_state_min_seq(struct imm_state const *state);
uint8_t                 imm_state_type_id(struct imm_state const *state);

/* Step */
struct imm_step *       imm_step_clone(struct imm_step const *step);
struct imm_step *       imm_step_create(struct imm_state const *state, unsigned seq_len);
void                    imm_step_destroy(struct imm_step const *step);
unsigned                imm_step_seq_len(struct imm_step const *step);
struct imm_state const *imm_step_state(struct imm_step const *step);

/* Subsequence */
struct imm_seq const *imm_subseq_cast(struct imm_subseq const *subseq);
struct imm_subseq     imm_subseq_init(struct imm_subseq *subseq, struct imm_seq const *seq,
                                      unsigned start, unsigned length);
unsigned              imm_subseq_length(struct imm_subseq const *subseq);
void                  imm_subseq_set(struct imm_subseq *subseq, unsigned start, unsigned length);
struct imm_subseq     imm_subseq_slice(struct imm_seq const *seq, unsigned start, unsigned length);
unsigned              imm_subseq_start(struct imm_subseq const *subseq);

/* Table state */
struct imm_table_state const *imm_table_state_create(char const *                name,
                                                     struct imm_seq_table const *table);
struct imm_table_state const *imm_table_state_derived(struct imm_state const *state);
void                          imm_table_state_destroy(struct imm_table_state const *state);
struct imm_state const *      imm_table_state_read(FILE *stream, struct imm_abc const *abc);
struct imm_state const *      imm_table_state_super(struct imm_table_state const *state);
int imm_table_state_write(struct imm_state const *state, struct imm_model const *model, FILE *stream);

/* Window */
struct imm_window *imm_window_create(struct imm_seq const *seq, unsigned length);
void               imm_window_destroy(struct imm_window const *window);
struct imm_subseq  imm_window_get(struct imm_window const *window, unsigned index);
unsigned           imm_window_size(struct imm_window const *window);
