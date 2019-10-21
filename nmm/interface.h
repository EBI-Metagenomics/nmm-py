struct imm_abc;
struct imm_mute_state;
struct imm_normal_state;
struct imm_state;
struct imm_table_state;
struct nmm_base;
struct nmm_codon;
struct nmm_frame_state;

struct imm_abc *imm_abc_create(const char *symbols);
void            imm_abc_destroy(struct imm_abc *abc);
int             imm_abc_length(const struct imm_abc *abc);
int             imm_abc_has_symbol(const struct imm_abc *abc, char symbol_id);
int             imm_abc_symbol_idx(const struct imm_abc *abc, char symbol_id);
char            imm_abc_symbol_id(const struct imm_abc *abc, int symbol_idx);

struct imm_state *      imm_state_create(const char *name, const struct imm_abc *abc,
                                         struct imm_state_funcs funcs, void *impl);
void                    imm_state_destroy(struct imm_state *state);
const char *            imm_state_get_name(const struct imm_state *state);
const struct imm_abc *  imm_state_get_abc(const struct imm_state *state);
double                  imm_state_lprob(const struct imm_state *state, const char *seq, int seq_len);
int                     imm_state_min_seq(const struct imm_state *state);
int                     imm_state_max_seq(const struct imm_state *state);
struct imm_state *      imm_state_cast(void *state);
const struct imm_state *imm_state_cast_c(const void *state);
void *                  imm_state_get_impl(struct imm_state *state);
const void *            imm_state_get_impl_c(const struct imm_state *state);

struct imm_normal_state *imm_normal_state_create(const char *name, const struct imm_abc *abc,
                                                 const double *lprobs);
void                     imm_normal_state_destroy(struct imm_normal_state *state);
int                      imm_normal_state_normalize(struct imm_normal_state *state);

struct imm_mute_state *imm_mute_state_create(const char *name, const struct imm_abc *abc);
void                   imm_mute_state_destroy(struct imm_mute_state *state);

struct imm_table_state *imm_table_state_create(const char *name, const struct imm_abc *abc);
void                    imm_table_state_destroy(struct imm_table_state *state);
void imm_table_state_add(struct imm_table_state *state, const char *seq, double lprob);
int  imm_table_state_normalize(struct imm_table_state *state);

struct nmm_codon *nmm_codon_create(const struct imm_abc *abc);
int               nmm_codon_set_lprob(struct nmm_codon *codon, char a, char b, char c, double lprob);
double            nmm_codon_get_lprob(const struct nmm_codon *codon, char a, char b, char c);
int               nmm_codon_normalize(struct nmm_codon *codon);
void              nmm_codon_destroy(struct nmm_codon *codon);

struct nmm_frame_state *nmm_frame_state_create(const char *name, const struct nmm_base *base,
                                               const struct nmm_codon *codon, double epsilon);
void                    nmm_frame_state_destroy(struct nmm_frame_state *state);

struct nmm_base *nmm_base_create(const struct imm_abc *abc);
int              nmm_base_set_lprob(struct nmm_base *base, char nucleotide, double lprob);
double           nmm_base_get_lprob(const struct nmm_base *base, char nucleotide);
int              nmm_base_normalize(struct nmm_base *base);
void             nmm_base_destroy(struct nmm_base *base);
