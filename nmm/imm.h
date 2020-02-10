struct imm_abc;
struct imm_hmm;
struct imm_mute_state;
struct imm_normal_state;
struct imm_path;
struct imm_seq;
struct imm_seq_table;
struct imm_state;
struct imm_step;
struct imm_table_state;

enum imm_symbol_type
{
    IMM_SYMBOL_UNKNOWN = 0,
    IMM_SYMBOL_NORMAL = 1,
    IMM_SYMBOL_ANY = 2,
};

/* Alphabet */
struct imm_abc const *imm_abc_create(char const *symbols, char any_symbol);
struct imm_abc const *imm_abc_clone(struct imm_abc const *abc);
void                  imm_abc_destroy(struct imm_abc const *abc);
unsigned              imm_abc_length(struct imm_abc const *abc);
char const *          imm_abc_symbols(struct imm_abc const *abc);
bool                  imm_abc_has_symbol(struct imm_abc const *abc, char symbol_id);
int                   imm_abc_symbol_idx(struct imm_abc const *abc, char symbol_id);
char                  imm_abc_symbol_id(struct imm_abc const *abc, unsigned symbol_idx);

/* Sequence */
struct imm_seq const *imm_seq_create(char const *seq, struct imm_abc const *abc);
struct imm_abc const *imm_seq_get_abc(struct imm_seq const *seq);
unsigned              imm_seq_length(struct imm_seq const *seq);
char const *          imm_seq_string(struct imm_seq const *seq);
struct imm_seq const *imm_seq_clone(struct imm_seq const *seq);
void                  imm_seq_destroy(struct imm_seq const *seq);

/* Squence table */
struct imm_seq_table *imm_seq_table_create(struct imm_abc const *abc);
struct imm_seq_table *imm_seq_table_clone(struct imm_seq_table const *table);
void                  imm_seq_table_destroy(struct imm_seq_table const *table);
int    imm_seq_table_add(struct imm_seq_table *table, struct imm_seq const *seq, double lprob);
int    imm_seq_table_normalize(struct imm_seq_table *table);
double imm_seq_table_lprob(struct imm_seq_table const *table, struct imm_seq const *seq);
struct imm_abc const *imm_seq_table_get_abc(struct imm_seq_table const *table);
unsigned              imm_seq_table_min_seq(struct imm_seq_table const *table);
unsigned              imm_seq_table_max_seq(struct imm_seq_table const *table);

/* State */
char const *            imm_state_get_name(struct imm_state const *state);
double                  imm_state_lprob(struct imm_state const *state, struct imm_seq const *seq);
unsigned                imm_state_min_seq(struct imm_state const *state);
unsigned                imm_state_max_seq(struct imm_state const *state);
struct imm_state const *imm_state_cast_c(void const *state);

/* Normal state */
struct imm_normal_state const *imm_normal_state_create(char const *name, struct imm_abc const *abc,
                                                       double const *lprobs);
void                           imm_normal_state_destroy(struct imm_normal_state const *state);

/* Mute state */
struct imm_mute_state const *imm_mute_state_create(char const *name, struct imm_abc const *abc);
void                         imm_mute_state_destroy(struct imm_mute_state const *state);

/* Table state */
struct imm_table_state *imm_table_state_create(char const *name, struct imm_seq_table const *table);
void                    imm_table_state_destroy(struct imm_table_state const *state);

/* HMM */
struct imm_hmm *imm_hmm_create(struct imm_abc const *abc);
void            imm_hmm_destroy(struct imm_hmm *hmm);
int    imm_hmm_add_state(struct imm_hmm *hmm, struct imm_state const *state, double start_lprob);
int    imm_hmm_del_state(struct imm_hmm *hmm, struct imm_state const *state);
int    imm_hmm_set_start(struct imm_hmm *hmm, struct imm_state const *state, double lprob);
int    imm_hmm_set_trans(struct imm_hmm *hmm, struct imm_state const *src_state,
                         struct imm_state const *dst_state, double lprob);
double imm_hmm_get_trans(struct imm_hmm const *hmm, struct imm_state const *src_state,
                         struct imm_state const *dst_state);
double imm_hmm_likelihood(struct imm_hmm const *hmm, struct imm_seq const *seq,
                          struct imm_path const *path);
struct imm_results const *imm_hmm_viterbi(struct imm_hmm const *hmm, struct imm_seq const *seq,
                                          struct imm_state const *end_state, unsigned window_length);
int                       imm_hmm_normalize(struct imm_hmm *hmm);
int                       imm_hmm_normalize_trans(struct imm_hmm *hmm, struct imm_state const *src);

/* Path */
struct imm_path *imm_path_create(void);
void             imm_path_destroy(struct imm_path const *path);
int imm_path_append(struct imm_path *path, struct imm_state const *state, unsigned seq_len);
int imm_path_prepend(struct imm_path *path, struct imm_state const *state, unsigned seq_len);
struct imm_step const *imm_path_first(struct imm_path const *path);
struct imm_step const *imm_path_last(struct imm_path const *path);
struct imm_step const *imm_path_next(struct imm_path const *path, struct imm_step const *step);

/* Step */
struct imm_state const *imm_step_state(struct imm_step const *step);
unsigned                imm_step_seq_len(struct imm_step const *step);
