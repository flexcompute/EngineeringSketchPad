#ifndef __LIMITERS_H__
#define __LIMITERS_H__

/* ------------- datatypes ------------------------ */
#ifndef __HAVE_LIMTYPE__
#define __HAVE_LIMTYPE__
typedef enum {                /* ..limiter options in */
              NOLIMITER,      /* order of increaseing */
              BJ_LIM,         /*          dissipation */
              VL_LIM,
              SINLIM,
              VA_LIM,
              MINLIM
              } limtype;
#endif /* __HAVE_LIMTYPE__ */


/* --------------  Interface Routines ------------- */

#ifdef __HAVE_GRID_STRUCTURES__  /* ...then safe to prototype */
void limitGradients(const p_tsGrid, const limtype);
#endif /* __HAVE_GRID_STRUCTURES__ */

#endif /* __LIMITERS_H__ */
