struct imm_abc;
struct imm_state;

struct imm_normal_state *imm_normal_state_create(const char *name,
                                                 const struct imm_abc *abc,
                                                 const double *lprobs);
void imm_normal_state_destroy(struct imm_normal_state *state);
int imm_normal_state_normalize(struct imm_normal_state *state);
