void integrate(mp_real *x,mp_real delta,int poincmapiters);
void f(mp_real *x, mp_real *dx, mp_real params[], mp_real delta, int lblock, int rblock);
void runge4(mp_real *xprev, mp_real *xnext, mp_real params[],mp_real delta, int lblock, int rblock, mp_real step);
int get_next_impact(mp_real *pimpact,mp_real *nimpact,mp_real *accumTime,mp_real params[],mp_real delta,int *lblock, int *rblock);
