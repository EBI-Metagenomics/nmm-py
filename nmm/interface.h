struct imm_abc;
struct imm_state;

struct imm_abc *imm_abc_create(const char *symbols);
void imm_abc_destroy(struct imm_abc *abc);
int imm_abc_length(const struct imm_abc *abc);
int imm_abc_has_symbol(const struct imm_abc *abc, char symbol_id);
int imm_abc_symbol_idx(const struct imm_abc *abc, char symbol_id);
char imm_abc_symbol_id(const struct imm_abc *abc, int symbol_idx);

struct imm_normal_state *imm_normal_state_create(const char *name,
                                                 const struct imm_abc *abc,
                                                 const double *lprobs);
void imm_normal_state_destroy(struct imm_normal_state *state);
int imm_normal_state_normalize(struct imm_normal_state *state);
