#ifndef PANPP_UTILS
#define PANPP_UTILS

#define fm6_comp(a) ((a) >= 1 && (a) <= 4 ? 5 - (a) : (a))

#define fm6_set_intv(e, c, ik)                                                 \
  ((ik).x[0] = (e)->cnt[(int)(c)],                                             \
   (ik).x[2] = (e)->cnt[(int)(c) + 1] - (e)->cnt[(int)(c)],                    \
   (ik).x[1] = (e)->cnt[fm6_comp(c)], (ik).info = 0)

static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
    5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

static const std::vector<std::string> int2char({"$", "A", "C", "G", "T", "N"});

static inline int kputsn(const char *p, int l, kstring_t *s) {
  if (s->l + l + 1 >= s->m) {
    char *tmp;
    s->m = s->l + l + 2;
    kroundup32(s->m);
    if ((tmp = (char *)realloc(s->s, s->m)))
      s->s = tmp;
    else
      return EOF;
  }
  memcpy(s->s + s->l, p, l);
  s->l += l;
  s->s[s->l] = 0;
  return l;
}

static void seq_char2nt6(int l, char *s) {
  for (int i = 0; i < l; ++i)
    s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
}

#endif
