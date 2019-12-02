struct imm_abc;
struct imm_hmm;
struct imm_mute_state;
struct imm_normal_state;
struct imm_path;
struct imm_state;
struct imm_step;
struct imm_table_state;

/* Alphabet */
struct imm_abc *imm_abc_create(char const *symbols, char any_symbol);
struct imm_abc *imm_abc_clone(struct imm_abc const *abc);
void            imm_abc_destroy(struct imm_abc *abc);
int             imm_abc_length(struct imm_abc const *abc);
char const *    imm_abc_symbols(struct imm_abc const *abc);
int             imm_abc_has_symbol(struct imm_abc const *abc, char symbol_id);
int             imm_abc_symbol_idx(struct imm_abc const *abc, char symbol_id);
char            imm_abc_symbol_id(struct imm_abc const *abc, int symbol_idx);

/* State */
char const *            imm_state_get_name(struct imm_state const *state);
double                  imm_state_lprob(struct imm_state const *state, char const *seq, int seq_len);
int                     imm_state_min_seq(struct imm_state const *state);
int                     imm_state_max_seq(struct imm_state const *state);
struct imm_state const *imm_state_cast_c(void const *state);

/* Normal state */
struct imm_normal_state *imm_normal_state_create(char const *name, struct imm_abc const *abc,
                                                 double const *lprobs);
void                     imm_normal_state_destroy(struct imm_normal_state *state);
int                      imm_normal_state_normalize(struct imm_normal_state *state);

/* Mute state */
struct imm_mute_state *imm_mute_state_create(char const *name, struct imm_abc const *abc);
void                   imm_mute_state_destroy(struct imm_mute_state *state);

/* Table state */
struct imm_table_state *imm_table_state_create(char const *name, struct imm_abc const *abc);
void                    imm_table_state_destroy(struct imm_table_state *state);
void imm_table_state_add(struct imm_table_state *state, char const *seq, double lprob);
int  imm_table_state_normalize(struct imm_table_state *state);

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
double imm_hmm_likelihood(struct imm_hmm const *hmm, char const *seq, struct imm_path const *path);
double imm_hmm_viterbi(struct imm_hmm const *hmm, char const *seq, struct imm_state const *end_state,
                       struct imm_path *path);
int    imm_hmm_normalize(struct imm_hmm *hmm);
int    imm_hmm_normalize_trans(struct imm_hmm *hmm, struct imm_state const *src);

/* Path */
struct imm_path *imm_path_create(void);
void             imm_path_destroy(struct imm_path *path);
int              imm_path_append(struct imm_path *path, struct imm_state const *state, int seq_len);
int              imm_path_prepend(struct imm_path *path, struct imm_state const *state, int seq_len);
struct imm_step const *imm_path_first(struct imm_path const *path);
struct imm_step const *imm_path_last(struct imm_path const *path);
struct imm_step const *imm_path_next(struct imm_path const *path, struct imm_step const *step);

/* Step */
struct imm_state const *imm_step_state(struct imm_step const *step);
int                     imm_step_seq_len(struct imm_step const *step);
