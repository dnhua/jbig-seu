#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "jbig.h"

 /* optional export of arithmetic coder functions for test purposes */
#ifdef TEST_CODEC
#define ARITH
#define ARITH_INL
#else
#define ARITH static
#ifdef __GNUC__
#define ARITH_INL static __inline__
#else
#define ARITH_INL static
#endif
#endif

#define MX_MAX 127 /* maximal supported mx offset for \
                    * adaptive template in the encoder */

#define TPB2CX 0x195 /* contexts for TP special pixels */
#define TPB3CX 0x0e5
#define TPDCX 0xc3f

/* marker codes */
#define MARKER_STUFF 0x00
#define MARKER_RESERVE 0x01
#define MARKER_SDNORM 0x02
#define MARKER_SDRST 0x03
#define MARKER_ABORT 0x04
#define MARKER_NEWLEN 0x05
#define MARKER_ATMOVE 0x06
#define MARKER_COMMENT 0x07
#define MARKER_ESC 0xff

/* loop array indices */
#define STRIPE 0
#define LAYER 1
#define PLANE 2

/* special jbg_buf pointers (instead of NULL) */
#define SDE_DONE ((struct jbg_buf *)-1)
#define SDE_TODO ((struct jbg_buf *)0)


/*
 * the following array specifies for each combination of the 3
 * ordering bits, which ii[] variable represents which dimension
 * of s->sde.
 */
static const int iindex[8][3] = {
	{2, 1, 0},    /* no ordering bit set */
	{-1, -1, -1}, /* SMID -> illegal combination */
	{2, 0, 1},    /* ILEAVE */
	{1, 0, 2},    /* SMID + ILEAVE */
	{0, 2, 1},    /* SEQ */
	{1, 2, 0},    /* SEQ + SMID */
	{0, 1, 2},    /* SEQ + ILEAVE */
	{-1, -1, -1}  /* SEQ + SMID + ILEAVE -> illegal combination */
};

/*
 * Array [language][message] with text string error messages that correspond
 * to return values from public functions in this library.
 */
#define NEMSG 9      /* number of error codes */
#define NEMSG_LANG 3 /* number of supported languages */
static const char *errmsg[NEMSG_LANG][NEMSG] = {
	/* English (JBG_EN) */
	{
		"Everything is ok",                             /* JBG_EOK */
		"Reached specified maximum size",               /* JBG_EOK_INTR */
		"Unexpected end of data",                       /* JBG_EAGAIN */
		"Not enough memory available",                  /* JBG_ENOMEM */
		"ABORT marker found",                           /* JBG_EABORT */
		"Unknown marker segment encountered",           /* JBG_EMARKER */
		"Incremental BIE does not fit to previous one", /* JBG_ENOCONT */
		"Invalid data encountered",                     /* JBG_EINVAL */
		"Unimplemented features used"                   /* JBG_EIMPL */
	},
	/* German (JBG_DE_8859_1) */
	{
		"Kein Problem aufgetreten",                          /* JBG_EOK */
		"Angegebene maximale Bildgr\366\337e erreicht",      /* JBG_EOK_INTR */
		"Unerwartetes Ende der Daten",                       /* JBG_EAGAIN */
		"Nicht gen\374gend Speicher vorhanden",              /* JBG_ENOMEM */
		"Es wurde eine Abbruch-Sequenz gefunden",            /* JBG_EABORT */
		"Eine unbekannte Markierungssequenz wurde gefunden", /* JBG_EMARKER */
		"Neue Daten passen nicht zu vorangegangenen Daten",  /* JBG_ENOCONT */
		"Es wurden ung\374ltige Daten gefunden",             /* JBG_EINVAL */
		"Noch nicht implementierte Optionen wurden benutzt"  /* JBG_EIMPL */
	},
	/* German (JBG_DE_UTF_8) */
	{
		"Kein Problem aufgetreten",                             /* JBG_EOK */
		"Angegebene maximale Bildgr\303\266\303\237e erreicht", /* JBG_EOK_INTR */
		"Unerwartetes Ende der Daten",                          /* JBG_EAGAIN */
		"Nicht gen\303\274gend Speicher vorhanden",             /* JBG_ENOMEM */
		"Es wurde eine Abbruch-Sequenz gefunden",               /* JBG_EABORT */
		"Eine unbekannte Markierungssequenz wurde gefunden",    /* JBG_EMARKER */
		"Neue Daten passen nicht zu vorangegangenen Daten",     /* JBG_ENOCONT */
		"Es wurden ung\303\274ltige Daten gefunden",            /* JBG_EINVAL */
		"Noch nicht implementierte Optionen wurden benutzt"     /* JBG_EIMPL */
	} };

/*
 * The following three functions are the only places in this code, were
 * C library memory management functions are called. The whole JBIG
 * library has been designed in order to allow multi-threaded
 * execution. No static or global variables are used, so all fuctions
 * are fully reentrant. However if you want to use this multi-thread
 * capability and your malloc, realloc and free are not reentrant,
 * then simply add the necessary semaphores or mutex primitives below.
 * In contrast to C's malloc() and realloc(), but like C's calloc(),
 * these functions take two parameters nmemb and size that are multiplied
 * before being passed on to the corresponding C function.
 * This we can catch all overflows during a size_t multiplication a
 * a single place.
 */

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1) /* largest value of size_t */
#endif

static void *checked_malloc(size_t nmemb, size_t size)
{
	void *p;

	/* Full manual exception handling is ugly here for performance
	 * reasons. If an adequate handling of lack of memory is required,
	 * then use C++ and throw a C++ exception instead of abort(). */

	 /* assert that nmemb * size <= SIZE_MAX */
	if (size > SIZE_MAX / nmemb)
		abort();

	p = malloc(nmemb * size);

	if (!p)
		abort();

#if 0
	fprintf(stderr, "%p = malloc(%lu * %lu)\n", p,
		(unsigned long)nmemb, (unsigned long)size);
#endif

	return p;
}

static void *checked_realloc(void *ptr, size_t nmemb, size_t size)
{
	void *p;

	/* Full manual exception handling is ugly here for performance
	 * reasons. If an adequate handling of lack of memory is required,
	 * then use C++ and throw a C++ exception here instead of abort(). */

	 /* assert that nmemb * size <= SIZE_MAX */
	if (size > SIZE_MAX / nmemb)
		abort();

	p = realloc(ptr, nmemb * size);

	if (!p)
		abort();

#if 0
	fprintf(stderr, "%p = realloc(%p, %lu * %lu)\n", p, ptr,
		(unsigned long)nmemb, (unsigned long)size);
#endif

	return p;
}

static void checked_free(void *ptr)
{
	free(ptr);

#if 0
	fprintf(stderr, "free(%p)\n", ptr);
#endif
}


ARITH void arith_decode_init(struct jbg_ardec_state *s, int reuse_st)
{
	int i;

	if (!reuse_st)
		for (i = 0; i < 4096; s->st[i++] = 0)
			;
	s->c = 0;
	s->a = 1;
	s->ct = 0;
	s->result = JBG_OK;
	s->startup = 1;
	return;
}

ARITH_INL int arith_decode(struct jbg_ardec_state *s, int cx)
{
	extern short jbg_lsz[];
	extern unsigned char jbg_nmps[], jbg_nlps[];
	register unsigned lsz, ss;
	register unsigned char *st;
	int pix;

	/* renormalization */
	while (s->a < 0x8000 || s->startup)
	{
		if (s->ct < 1 && s->result != JBG_READY)
		{
			/* first we have to move a new byte into s->c */
			if (s->pscd_ptr >= s->pscd_end)
			{
				s->result = JBG_MORE;
				return -1;
			}
			if (*s->pscd_ptr == 0xff)
				if (s->pscd_ptr + 1 >= s->pscd_end)
				{
					s->result = JBG_MARKER;
					return -1;
				}
				else
				{
					if (*(s->pscd_ptr + 1) == MARKER_STUFF)
					{
						s->c |= 0xffL << (8 - s->ct);
						s->ct += 8;
						s->pscd_ptr += 2;
						s->result = JBG_OK;
					}
					else
						s->result = JBG_READY;
				}
			else
			{
				s->c |= (long)*(s->pscd_ptr++) << (8 - s->ct);
				s->ct += 8;
				s->result = JBG_OK;
			}
		}
		s->c <<= 1;
		s->a <<= 1;
		--s->ct;
		if (s->a == 0x10000L)
			s->startup = 0;
	}

	st = s->st + cx;
	ss = *st & 0x7f;
	assert(ss < 113);
	lsz = jbg_lsz[ss];

#if 0
	fprintf(stderr, "cx = %d, mps = %d, st = %3d, lsz = 0x%04x, a = 0x%05lx, "
		"c = 0x%08lx, ct = %2d\n",
		cx, !!(s->st[cx] & 0x80), ss, lsz, s->a, s->c, s->ct);
#endif

	if ((s->c >> 16) < (s->a -= lsz))
		if (s->a & 0xffff8000L)
			return *st >> 7;
		else
		{
			/* MPS_EXCHANGE */
			if (s->a < lsz)
			{
				pix = 1 - (*st >> 7);
				/* Check whether MPS/LPS exchange is necessary
			 * and chose next probability estimator status */
				*st &= 0x80;
				*st ^= jbg_nlps[ss];
			}
			else
			{
				pix = *st >> 7;
				*st &= 0x80;
				*st |= jbg_nmps[ss];
			}
		}
	else
	{
		/* LPS_EXCHANGE */
		if (s->a < lsz)
		{
			s->c -= s->a << 16;
			s->a = lsz;
			pix = *st >> 7;
			*st &= 0x80;
			*st |= jbg_nmps[ss];
		}
		else
		{
			s->c -= s->a << 16;
			s->a = lsz;
			pix = 1 - (*st >> 7);
			/* Check whether MPS/LPS exchange is necessary
			 * and chose next probability estimator status */
			*st &= 0x80;
			*st ^= jbg_nlps[ss];
		}
	}

	return pix;
}

/*
 * Calculate y = ceil(x/2) applied n times, which is equivalent to
 * y = ceil(x/(2^n)). This function is used to
 * determine the number of pixels per row or column after n resolution
 * reductions. E.g. X[d-1] = jbg_ceil_half(X[d], 1) and X[0] =
 * jbg_ceil_half(X[d], d) as defined in clause 6.2.3 of T.82.
 */
unsigned long jbg_ceil_half(unsigned long x, int n)
{
	unsigned long mask;

	assert(n >= 0 && n < 32);
	mask = (1UL << n) - 1; /* the lowest n bits are 1 here */
	return (x >> n) + ((mask & x) != 0);
}


/*
 * Calculate the number of stripes, as defined in clause 6.2.3 of T.82.
 */
static unsigned long jbg_stripes(unsigned long l0, unsigned long yd,
	unsigned long d)
{
	unsigned long y0 = jbg_ceil_half(yd, d);

	return y0 / l0 + (y0 % l0 != 0);
}

/*
 * Convert the table which controls the deterministic prediction
 * process from the 1728 byte long DPTABLE format into the 6912 byte long
 * internal format.
 */
void jbg_dppriv2int(char *internal, const unsigned char *dptable)
{
	int i, j, k;
	int trans0[8] = { 1, 0, 3, 2, 7, 6, 5, 4 };
	int trans1[9] = { 1, 0, 3, 2, 8, 7, 6, 5, 4 };
	int trans2[11] = { 1, 0, 3, 2, 10, 9, 8, 7, 6, 5, 4 };
	int trans3[12] = { 1, 0, 3, 2, 11, 10, 9, 8, 7, 6, 5, 4 };

#define FILL_TABLE2(offset, len, trans)                           \
  for (i = 0; i < len; i++)                                       \
  {                                                               \
    k = 0;                                                        \
    for (j = 0; j < 8; j++)                                       \
      k |= ((i >> j) & 1) << trans[j];                            \
    internal[k + offset] =                                        \
        (dptable[(i + offset) >> 2] >> ((3 - (i & 3)) << 1)) & 3; \
  }

	FILL_TABLE2(0, 256, trans0);
	FILL_TABLE2(256, 512, trans1);
	FILL_TABLE2(768, 2048, trans2);
	FILL_TABLE2(2816, 4096, trans3);

	return;
}

/*
 * The constructor for a decoder
 */
void jbg_dec_init(struct jbg_dec_state *s)
{
	s->order = 0;
	s->d = -1;
	s->bie_len = 0;
	s->buf_len = 0;
	s->dppriv = NULL;
	s->xmax = 4294967295UL;
	s->ymax = 4294967295UL;
	s->dmax = 256;
	s->s = NULL;

	return;
}

/*
 * Specify a maximum image size for the decoder. If the JBIG file has
 * the order bit ILEAVE, but not the bit SEQ set, then the decoder
 * will abort to decode after the image has reached the maximal
 * resolution layer which is still not wider than xmax or higher than
 * ymax.
 */
void jbg_dec_maxsize(struct jbg_dec_state *s, unsigned long xmax,
	unsigned long ymax)
{
	if (xmax > 0)
		s->xmax = xmax;
	if (ymax > 0)
		s->ymax = ymax;

	return;
}

/*
 * Decode the new len PSDC bytes to which data points and add them to
 * the current stripe. Return the number of bytes which have actually
 * been read (this will be less than len if a marker segment was
 * part of the data or if the final byte was 0xff were this code
 * can not determine, whether we have a marker segment.
 */
static size_t decode_pscd(struct jbg_dec_state *s, unsigned char *data,
	size_t len)
{
	unsigned long stripe;
	unsigned int layer, plane;
	unsigned long hl, ll, y, hx, hy, lx, ly, hbpl, lbpl;
	unsigned char *hp, *lp1, *lp2, *p1, *q1;
	register unsigned long line_h1, line_h2, line_h3;
	register unsigned long line_l1, line_l2, line_l3;
	struct jbg_ardec_state *se;
	unsigned long x;
	long o;
	unsigned a;
	int n;
	int pix, cx = 0, slntp, tx;

	/* SDE loop variables */
	stripe = s->ii[iindex[s->order & 7][STRIPE]];
	layer = s->ii[iindex[s->order & 7][LAYER]];
	plane = s->ii[iindex[s->order & 7][PLANE]];

	/* forward data to arithmetic decoder */
	se = s->s[plane] + layer - s->dl;
	se->pscd_ptr = data;
	se->pscd_end = data + len;

	/* number of lines per stripe in highres image */
	hl = s->l0 << layer;
	/* number of lines per stripe in lowres image */
	ll = hl >> 1;
	/* current line number in highres image */
	y = stripe * hl + s->i;
	/* number of pixels in highres image */
	hx = jbg_ceil_half(s->xd, s->d - layer);
	hy = jbg_ceil_half(s->yd, s->d - layer);
	/* number of pixels in lowres image */
	lx = jbg_ceil_half(hx, 1);
	ly = jbg_ceil_half(hy, 1);
	/* bytes per line in highres and lowres image */
	hbpl = jbg_ceil_half(hx, 3);
	lbpl = jbg_ceil_half(lx, 3);
	/* pointer to highres and lowres image bytes */
	hp = s->lhp[layer & 1][plane] + (stripe * hl + s->i) * hbpl +
		(s->x >> 3);
	lp2 = s->lhp[(layer - 1) & 1][plane] + (stripe * ll + (s->i >> 1)) * lbpl +
		(s->x >> 4);
	lp1 = lp2 + lbpl;

	/* restore a few local variables */
	line_h1 = s->line_h1;
	line_h2 = s->line_h2;
	line_h3 = s->line_h3;
	line_l1 = s->line_l1;
	line_l2 = s->line_l2;
	line_l3 = s->line_l3;
	x = s->x;

	if (s->x == 0 && s->i == 0 &&
		(stripe == 0 || s->reset[plane][layer - s->dl]))
	{
		s->tx[plane][layer - s->dl] = s->ty[plane][layer - s->dl] = 0;
		if (s->pseudo)
			s->lntp[plane][layer - s->dl] = 1;
	}

#ifdef DEBUG
	if (s->x == 0 && s->i == 0 && s->pseudo)
		fprintf(stderr, "decode_pscd(%p, %p, %ld): s/d/p = %2lu/%2u/%2u\n",
		(void *)s, (void *)data, (long)len, stripe, layer, plane);
#endif

	if (layer == 0)
	{

		/*
		 *  Decode lowest resolution layer
		 */
		 /* dnhua：解码低层 20191128 */
		for (; s->i < hl && y < hy; s->i++, y++)
		{

			/* adaptive template changes */
			/* dnhua：自适应像素（AT）20191128 */
			if (x == 0)
				for (n = 0; n < s->at_moves; n++)
					if (s->at_line[n] == s->i)
					{
						s->tx[plane][layer - s->dl] = s->at_tx[n];
						s->ty[plane][layer - s->dl] = s->at_ty[n];
#ifdef DEBUG
						fprintf(stderr, "ATMOVE: line=%lu, tx=%d, ty=%d.\n", s->i,
							s->tx[plane][layer - s->dl], s->ty[plane][layer - s->dl]);
#endif
					}
			tx = s->tx[plane][layer - s->dl];
			assert(tx >= 0); /* i.e., tx can safely be cast to unsigned */

			/* typical prediction */
			/* dnhua：典型预测模块（TP）20191128 */
			if (s->options & JBG_TPBON && s->pseudo)
			{
				slntp = arith_decode(se, (s->options & JBG_LRLTWO) ? TPB2CX : TPB3CX);
				if (se->result == JBG_MORE || se->result == JBG_MARKER)
					goto leave;
				s->lntp[plane][layer - s->dl] =
					!(slntp ^ s->lntp[plane][layer - s->dl]);
				if (s->lntp[plane][layer - s->dl])
				{
					/* this line is 'not typical' and has to be coded completely */
					s->pseudo = 0;
				}
				else
				{
					/* this line is 'typical' (i.e. identical to the previous one) */
					p1 = hp;
					if (s->i == 0 && (stripe == 0 || s->reset[plane][layer - s->dl]))
						while (p1 < hp + hbpl)
							*p1++ = 0;
					else
					{
						q1 = hp - hbpl;
						while (q1 < hp)
							*p1++ = *q1++;
					}
					hp += hbpl;
					continue;
				}
			}

			/*
			 * Layout of the variables line_h1, line_h2, line_h3, which contain
			 * as bits the neighbour pixels of the currently decoded pixel X:
			 *
			 *                     76543210 76543210 76543210 76543210     line_h3
			 *                     76543210 76543210 76543210 76543210     line_h2
			 *   76543210 76543210 76543210 76543210 X                     line_h1
			 */

			if (x == 0)
			{
				line_h1 = line_h2 = line_h3 = 0;
				if (s->i > 0 || (y > 0 && !s->reset[plane][layer - s->dl]))
					line_h2 = (long)*(hp - hbpl) << 8;
				if (s->i > 1 || (y > 1 && !s->reset[plane][layer - s->dl]))
					line_h3 = (long)*(hp - hbpl - hbpl) << 8;
			}

			/*
			 * Another tiny JBIG standard bug:
			 *
			 * While implementing the line_h3 handling here, I discovered
			 * another problem with the ITU-T T.82(1993 E) specification.
			 * This might be a somewhat pathological case, however. The
			 * standard is unclear about how a decoder should behave in the
			 * following situation:
			 *
			 * Assume we are in layer 0 and all stripes are single lines
			 * (L0=1 allowed by table 9). We are now decoding the first (and
			 * only) line of the third stripe. Assume, the first stripe was
			 * terminated by SDRST and the second stripe was terminated by
			 * SDNORM. While decoding the only line of the third stripe with
			 * the three-line template, we need access to pixels from the
			 * previous two stripes. We know that the previous stripe
			 * terminated with SDNROM, so we access the pixel from the
			 * second stripe. But do we have to replace the pixels from the
			 * first stripe by background pixels, because this stripe ended
			 * with SDRST? The standard, especially clause 6.2.5 does never
			 * mention this case, so the behaviour is undefined here. My
			 * current implementation remembers only the marker used to
			 * terminate the previous stripe. In the above example, the
			 * pixels of the first stripe are accessed despite the fact that
			 * this stripe ended with SDRST. An alternative (only slightly
			 * more complicated) implementation would be to remember the end
			 * marker (SDNORM or SDRST) of the previous two stripes in a
			 * plane/layer and to act accordingly when accessing the two
			 * previous lines. What am I supposed to do here?
			 *
			 * As the standard is unclear about the correct behaviour in the
			 * situation of the above example, I strongly suggest to avoid
			 * the following situation while encoding data with JBIG:
			 *
			 *   LRLTWO = 0, L0=1 and both SDNORM and SDRST appear in layer 0.
			 *
			 * I guess that only a very few if any encoders will switch
			 * between SDNORM and SDRST, so let us hope that this ambiguity
			 * in the standard will never cause any interoperability
			 * problems.
			 *
			 * Markus Kuhn -- 1995-04-30
			 */

			 /* decode line */
			while (x < hx)
			{
				if ((x & 7) == 0)
				{
					if (x < hbpl * 8 - 8 &&
						(s->i > 0 || (y > 0 && !s->reset[plane][layer - s->dl])))
					{
						line_h2 |= *(hp - hbpl + 1);
						if (s->i > 1 || (y > 1 && !s->reset[plane][layer - s->dl]))
							line_h3 |= *(hp - hbpl - hbpl + 1);
					}
				}
				if (s->options & JBG_LRLTWO)
				{
					/* two line template */
					/* dnhua：两线模板 20191128 */
					do
					{
						if (tx)
						{
							if ((unsigned)tx > x)
								a = 0;
							else if (tx < 8)
								a = ((line_h1 >> (tx - 5)) & 0x010);
							else
							{
								o = (x - tx) - (x & ~7L);
								a = (hp[o >> 3] >> (7 - (o & 7))) & 1;
								a <<= 4;
							}
							assert(tx > 31 ||
								a == ((line_h1 >> (tx - 5)) & 0x010));
							pix = arith_decode(se, (((line_h2 >> 9) & 0x3e0) | a |
								(line_h1 & 0x00f)));
						}
						else
							pix = arith_decode(se, (((line_h2 >> 9) & 0x3f0) |
							(line_h1 & 0x00f)));
						if (se->result == JBG_MORE || se->result == JBG_MARKER)
							goto leave;
						line_h1 = (line_h1 << 1) | pix;
						line_h2 <<= 1;
					} while ((++x & 7) && x < hx);
				}
				else
				{
					/* three line template */
					/* dnhua：三线模板 20191128 */
					do
					{
						if (tx)
						{
							if ((unsigned)tx > x)
								a = 0;
							else if (tx < 8)
								a = ((line_h1 >> (tx - 3)) & 0x004);
							else
							{
								o = (x - tx) - (x & ~7L);
								a = (hp[o >> 3] >> (7 - (o & 7))) & 1;
								a <<= 2;
							}
							assert(tx > 31 ||
								a == ((line_h1 >> (tx - 3)) & 0x004));
							pix = arith_decode(se, (((line_h3 >> 7) & 0x380) |
								((line_h2 >> 11) & 0x078) | a |
								(line_h1 & 0x003)));
						}
						else
							pix = arith_decode(se, (((line_h3 >> 7) & 0x380) |
							((line_h2 >> 11) & 0x07c) |
								(line_h1 & 0x003)));
						if (se->result == JBG_MORE || se->result == JBG_MARKER)
							goto leave;

						line_h1 = (line_h1 << 1) | pix;
						line_h2 <<= 1;
						line_h3 <<= 1;
					} while ((++x & 7) && x < hx);
				} /* if (s->options & JBG_LRLTWO) */
				*hp++ = line_h1;
			} /* while */
			*(hp - 1) <<= hbpl * 8 - hx;
			x = 0;
			s->pseudo = 1;
		} /* for (i = ...) */
	}
	else
	{

		/*
		 *  Decode differential layer
		 */

		for (; s->i < hl && y < hy; s->i++, y++)
		{

			/* adaptive template changes */
			if (x == 0)
				for (n = 0; n < s->at_moves; n++)
					if (s->at_line[n] == s->i)
					{
						s->tx[plane][layer - s->dl] = s->at_tx[n];
						s->ty[plane][layer - s->dl] = s->at_ty[n];
#ifdef DEBUG
						fprintf(stderr, "ATMOVE: line=%lu, tx=%d, ty=%d.\n", s->i,
							s->tx[plane][layer - s->dl], s->ty[plane][layer - s->dl]);
#endif
					}
			tx = s->tx[plane][layer - s->dl];

			/* handle lower border of low-resolution image */
			if ((s->i >> 1) >= ll - 1 || (y >> 1) >= ly - 1)
				lp1 = lp2;

			/* typical prediction */
			if (s->options & JBG_TPDON && s->pseudo)
			{
				s->lntp[plane][layer - s->dl] = arith_decode(se, TPDCX);
				if (se->result == JBG_MORE || se->result == JBG_MARKER)
					goto leave;
				s->pseudo = 0;
			}

			/*
			 * Layout of the variables line_h1, line_h2, line_h3, which contain
			 * as bits the high resolution neighbour pixels of the currently
			 * decoded highres pixel X:
			 *
			 *                     76543210 76543210 76543210 76543210     line_h3
			 *                     76543210 76543210 76543210 76543210     line_h2
			 *   76543210 76543210 76543210 76543210 X                     line_h1
			 *
			 * Layout of the variables line_l1, line_l2, line_l3, which contain
			 * the low resolution pixels near the currently decoded pixel as bits.
			 * The lowres pixel in which the currently coded highres pixel is
			 * located is marked as Y:
			 *
			 *                     76543210 76543210 76543210 76543210     line_l3
			 *                     76543210 76543210 Y6543210 76543210     line_l2
			 *                     76543210 76543210 76543210 76543210     line_l1
			 */

			if (x == 0)
			{
				line_h1 = line_h2 = line_h3 = line_l1 = line_l2 = line_l3 = 0;
				if (s->i > 0 || (y > 0 && !s->reset[plane][layer - s->dl]))
				{
					line_h2 = (long)*(hp - hbpl) << 8;
					if (s->i > 1 || (y > 1 && !s->reset[plane][layer - s->dl]))
						line_h3 = (long)*(hp - hbpl - hbpl) << 8;
				}
				if (s->i > 1 || (y > 1 && !s->reset[plane][layer - s->dl]))
					line_l3 = (long)*(lp2 - lbpl) << 8;
				line_l2 = (long)*lp2 << 8;
				line_l1 = (long)*lp1 << 8;
			}

			/* decode line */
			while (x < hx)
			{
				if ((x & 15) == 0)
					if ((x >> 1) < lbpl * 8 - 8)
					{
						line_l1 |= *(lp1 + 1);
						line_l2 |= *(lp2 + 1);
						if (s->i > 1 ||
							(y > 1 && !s->reset[plane][layer - s->dl]))
							line_l3 |= *(lp2 - lbpl + 1);
					}
				do
				{

					assert(hp - (s->lhp[layer & 1][plane] + (stripe * hl + s->i) * hbpl) == (ptrdiff_t)x >> 3);
					assert(lp2 - (s->lhp[(layer - 1) & 1][plane] + (stripe * ll + (s->i >> 1)) * lbpl) == (ptrdiff_t)x >> 4);

					if ((x & 7) == 0)
						if (x < hbpl * 8 - 8)
						{
							if (s->i > 0 || (y > 0 && !s->reset[plane][layer - s->dl]))
							{
								line_h2 |= *(hp + 1 - hbpl);
								if (s->i > 1 || (y > 1 && !s->reset[plane][layer - s->dl]))
									line_h3 |= *(hp + 1 - hbpl - hbpl);
							}
						}
					do
					{
						if (!s->lntp[plane][layer - s->dl])
							cx = (((line_l3 >> 14) & 0x007) |
							((line_l2 >> 11) & 0x038) |
								((line_l1 >> 8) & 0x1c0));
						if (!s->lntp[plane][layer - s->dl] &&
							(cx == 0x000 || cx == 0x1ff))
						{
							/* pixels are typical and have not to be decoded */
							do
							{
								line_h1 = (line_h1 << 1) | (cx & 1);
							} while ((++x & 1) && x < hx);
							line_h2 <<= 2;
							line_h3 <<= 2;
						}
						else
							do
							{

								/* deterministic prediction */
								/* dnhua：确定性预测（DP）*/
								if (s->options & JBG_DPON)
									if ((y & 1) == 0)
										if ((x & 1) == 0)
											/* phase 0 */
											pix = s->dppriv[((line_l3 >> 15) & 0x003) |
											((line_l2 >> 13) & 0x00c) |
											((line_h1 << 4) & 0x010) |
											((line_h2 >> 9) & 0x0e0)];
										else
											/* phase 1 */
											pix = s->dppriv[(((line_l3 >> 15) & 0x003) |
											((line_l2 >> 13) & 0x00c) |
												((line_h1 << 4) & 0x030) |
												((line_h2 >> 9) & 0x1c0)) +
											256];
									else if ((x & 1) == 0)
										/* phase 2 */
										pix = s->dppriv[(((line_l3 >> 15) & 0x003) |
										((line_l2 >> 13) & 0x00c) |
											((line_h1 << 4) & 0x010) |
											((line_h2 >> 9) & 0x0e0) |
											((line_h3 >> 6) & 0x700)) +
										768];
									else
										/* phase 3 */
										pix = s->dppriv[(((line_l3 >> 15) & 0x003) |
										((line_l2 >> 13) & 0x00c) |
											((line_h1 << 4) & 0x030) |
											((line_h2 >> 9) & 0x1c0) |
											((line_h3 >> 6) & 0xe00)) +
										2816];
								else
									pix = 2;

								if (pix & 2)
								{
									if (tx)
										cx = ((line_h1 & 0x003) |
										(((line_h1 << 2) >> (tx - 3)) & 0x010) |
											((line_h2 >> 12) & 0x00c) |
											((line_h3 >> 10) & 0x020));
									else
										cx = ((line_h1 & 0x003) |
										((line_h2 >> 12) & 0x01c) |
											((line_h3 >> 10) & 0x020));
									if (x & 1)
										cx |= (((line_l2 >> 8) & 0x0c0) |
										((line_l1 >> 6) & 0x300)) |
											(1UL << 10);
									else
										cx |= (((line_l2 >> 9) & 0x0c0) |
										((line_l1 >> 7) & 0x300));
									cx |= (y & 1) << 11;

									pix = arith_decode(se, cx);
									if (se->result == JBG_MORE || se->result == JBG_MARKER)
										goto leave;
								}

								line_h1 = (line_h1 << 1) | pix;
								line_h2 <<= 1;
								line_h3 <<= 1;

							} while ((++x & 1) && x < hx);
							line_l1 <<= 1;
							line_l2 <<= 1;
							line_l3 <<= 1;
					} while ((x & 7) && x < hx);
					*hp++ = line_h1;
				} while ((x & 15) && x < hx);
				++lp1;
				++lp2;
			} /* while */
			x = 0;

			*(hp - 1) <<= hbpl * 8 - hx;
			if ((s->i & 1) == 0)
			{
				/* low resolution pixels are used twice */
				lp1 -= lbpl;
				lp2 -= lbpl;
			}
			else
				s->pseudo = 1;

		} /* for (i = ...) */
	}

leave:

	/* save a few local variables */
	s->line_h1 = line_h1;
	s->line_h2 = line_h2;
	s->line_h3 = line_h3;
	s->line_l1 = line_l1;
	s->line_l2 = line_l2;
	s->line_l3 = line_l3;
	s->x = x;

	return se->pscd_ptr - data;
}

/*
 * Provide a new BIE fragment to the decoder.
 *
 * If cnt is not NULL, then *cnt will contain after the call the
 * number of actually read bytes. If the data was not complete, then
 * the return value will be JBG_EAGAIN and *cnt == len. In case this
 * function has returned with JBG_EOK, then it has reached the end of
 * a BIE but it can be called again with data from the next BIE if
 * there exists one in order to get to a higher resolution layer. In
 * case the return value was JBG_EOK_INTR then this function can be
 * called again with the rest of the BIE, because parsing the BIE has
 * been interrupted by a jbg_dec_maxsize() specification. In both
 * cases the remaining len - *cnt bytes of the previous block will
 * have to passed to this function again (if len > *cnt). In case of
 * any other return value than JBG_EOK, JBG_EOK_INTR or JBG_EAGAIN, a
 * serious problem has occured and the only function you should call
 * is jbg_dec_free() in order to remove the mess (and probably
 * jbg_strerror() in order to find out what to tell the user).
 */
int jbg_dec_in(struct jbg_dec_state *s, unsigned char *data, size_t len,
	size_t *cnt)
{
	int i, j, required_length;
	unsigned long x, y;
	unsigned long is[3], ie[3];
	extern char jbg_dptable[];
	size_t dummy_cnt;

	if (!cnt)
		cnt = &dummy_cnt;
	*cnt = 0;
	if (len < 1)
		return JBG_EAGAIN;

	/* read in 20-byte BIH */
	if (s->bie_len < 20)
	{
		while (s->bie_len < 20 && *cnt < len)
			s->buffer[s->bie_len++] = data[(*cnt)++];
		if (s->bie_len < 20)
			return JBG_EAGAIN;
		if (s->buffer[1] < s->buffer[0])
			return JBG_EINVAL;
		/* test whether this looks like a valid JBIG header at all */
		if (s->buffer[3] != 0 || (s->buffer[18] & 0xf0) != 0 ||
			(s->buffer[19] & 0x80) != 0)
			return JBG_EINVAL;
		if (s->buffer[0] != s->d + 1)
			return JBG_ENOCONT;
		s->dl = s->buffer[0];
		s->d = s->buffer[1];
		if (s->dl == 0)
			s->planes = s->buffer[2];
		else if (s->planes != s->buffer[2])
			return JBG_ENOCONT;
		x = (((long)s->buffer[4] << 24) | ((long)s->buffer[5] << 16) |
			((long)s->buffer[6] << 8) | (long)s->buffer[7]);
		y = (((long)s->buffer[8] << 24) | ((long)s->buffer[9] << 16) |
			((long)s->buffer[10] << 8) | (long)s->buffer[11]);
		if (s->dl != 0 && ((s->xd << (s->d - s->dl + 1)) != x &&
			(s->yd << (s->d - s->dl + 1)) != y))
			return JBG_ENOCONT;
		s->xd = x;
		s->yd = y;
		s->l0 = (((long)s->buffer[12] << 24) | ((long)s->buffer[13] << 16) |
			((long)s->buffer[14] << 8) | (long)s->buffer[15]);
		/* ITU-T T.85 trick not directly supported by decoder; for full
		 * T.85 compatibility with respect to all NEWLEN marker scenarios,
		 * preprocess BIE with jbg_newlen() before passing it to the decoder. */
		if (s->yd == 0xffffffff)
			return JBG_EIMPL;
		if (!s->planes || !s->xd || !s->yd || !s->l0)
			return JBG_EINVAL;
		/* prevent uint32 overflow: s->l0 * 2 ^ s->d < 2 ^ 32 */
		if (s->d > 31 || (s->d != 0 && s->l0 >= (1UL << (32 - s->d))))
			return JBG_EIMPL;
		s->mx = s->buffer[16];
		if (s->mx > 127)
			return JBG_EINVAL;
		s->my = s->buffer[17];
#if 0
		if (s->my > 0)
			return JBG_EIMPL;
#endif
		s->order = s->buffer[18];
		if (iindex[s->order & 7][0] < 0)
			return JBG_EINVAL;
		/* HITOLO and SEQ currently not yet implemented */
		if (s->dl != s->d && (s->order & JBG_HITOLO || s->order & JBG_SEQ))
			return JBG_EIMPL;
		s->options = s->buffer[19];

		/* calculate number of stripes that will be required */
		/* dnhua：关键步骤1，计算条数目 20191128 */
		s->stripes = jbg_stripes(s->l0, s->yd, s->d);

		/* some initialization */
		/* dnhua：一些初始化工作，无性能瓶颈 20191128 */
		s->ii[iindex[s->order & 7][STRIPE]] = 0;
		s->ii[iindex[s->order & 7][LAYER]] = s->dl;
		s->ii[iindex[s->order & 7][PLANE]] = 0;
		if (s->dl == 0)
		{
			s->s = (struct jbg_ardec_state **)
				checked_malloc(s->planes, sizeof(struct jbg_ardec_state *));
			s->tx = (int **)checked_malloc(s->planes, sizeof(int *));
			s->ty = (int **)checked_malloc(s->planes, sizeof(int *));
			s->reset = (int **)checked_malloc(s->planes, sizeof(int *));
			s->lntp = (int **)checked_malloc(s->planes, sizeof(int *));
			s->lhp[0] = (unsigned char **)
				checked_malloc(s->planes, sizeof(unsigned char *));
			s->lhp[1] = (unsigned char **)
				checked_malloc(s->planes, sizeof(unsigned char *));
			for (i = 0; i < s->planes; i++)
			{
				s->s[i] = (struct jbg_ardec_state *)
					checked_malloc(s->d - s->dl + 1, sizeof(struct jbg_ardec_state));
				s->tx[i] = (int *)checked_malloc(s->d - s->dl + 1, sizeof(int));
				s->ty[i] = (int *)checked_malloc(s->d - s->dl + 1, sizeof(int));
				s->reset[i] = (int *)checked_malloc(s->d - s->dl + 1, sizeof(int));
				s->lntp[i] = (int *)checked_malloc(s->d - s->dl + 1, sizeof(int));
				s->lhp[s->d & 1][i] = (unsigned char *)
					checked_malloc(s->yd, jbg_ceil_half(s->xd, 3));
				s->lhp[(s->d - 1) & 1][i] = (unsigned char *)
					checked_malloc(jbg_ceil_half(s->yd, 1), jbg_ceil_half(s->xd, 1 + 3));
			}
		}
		else
		{
			for (i = 0; i < s->planes; i++)
			{
				s->s[i] = (struct jbg_ardec_state *)
					checked_realloc(s->s[i], s->d - s->dl + 1,
						sizeof(struct jbg_ardec_state));
				s->tx[i] = (int *)checked_realloc(s->tx[i],
					s->d - s->dl + 1, sizeof(int));
				s->ty[i] = (int *)checked_realloc(s->ty[i],
					s->d - s->dl + 1, sizeof(int));
				s->reset[i] = (int *)checked_realloc(s->reset[i],
					s->d - s->dl + 1, sizeof(int));
				s->lntp[i] = (int *)checked_realloc(s->lntp[i],
					s->d - s->dl + 1, sizeof(int));
				s->lhp[s->d & 1][i] = (unsigned char *)
					checked_realloc(s->lhp[s->d & 1][i],
						s->yd, jbg_ceil_half(s->xd, 3));
				s->lhp[(s->d - 1) & 1][i] = (unsigned char *)
					checked_realloc(s->lhp[(s->d - 1) & 1][i],
						jbg_ceil_half(s->yd, 1), jbg_ceil_half(s->xd, 1 + 3));
			}
		}
		for (i = 0; i < s->planes; i++)
			for (j = 0; j <= s->d - s->dl; j++)
				arith_decode_init(s->s[i] + j, 0);
		if (s->dl == 0 || (s->options & JBG_DPON && !(s->options & JBG_DPPRIV)))
			s->dppriv = jbg_dptable;
		s->comment_skip = 0;
		s->buf_len = 0;
		s->x = 0;
		s->i = 0;
		s->pseudo = 1;
		s->at_moves = 0;
	}

	/* read in DPTABLE */
	if (s->bie_len < 20 + 1728 &&
		(s->options & (JBG_DPON | JBG_DPPRIV | JBG_DPLAST)) ==
		(JBG_DPON | JBG_DPPRIV))
	{
		assert(s->bie_len >= 20);
		while (s->bie_len < 20 + 1728 && *cnt < len)
			s->buffer[s->bie_len++ - 20] = data[(*cnt)++];
		if (s->bie_len < 20 + 1728)
			return JBG_EAGAIN;
		if (!s->dppriv || s->dppriv == jbg_dptable)
			s->dppriv = (char *)checked_malloc(1728, sizeof(char)); /* dnhua：dppriv-optional private deterministic prediction table 20191128 */
		jbg_dppriv2int(s->dppriv, s->buffer);                     /* dnhua：buffer->dppriv  1728byte->6912byte  20191128 */
	}

	/*
	 * BID processing loop
	 */
	 /* dnhua：解码主程序 20191128 */
	while (*cnt < len)
	{

		/* process floating marker segments */

		/* skip COMMENT contents */
		if (s->comment_skip)
		{
			if (s->comment_skip <= len - *cnt)
			{
				*cnt += s->comment_skip;
				s->comment_skip = 0;
			}
			else
			{
				s->comment_skip -= len - *cnt;
				*cnt = len;
			}
			continue;
		}

		/* load complete marker segments into s->buffer for processing */
		if (s->buf_len > 0)
		{
			assert(s->buffer[0] == MARKER_ESC);
			while (s->buf_len < 2 && *cnt < len)
				s->buffer[s->buf_len++] = data[(*cnt)++];
			if (s->buf_len < 2)
				continue;
			switch (s->buffer[1])
			{
			case MARKER_COMMENT:
				required_length = 6;
				break;
			case MARKER_ATMOVE:
				required_length = 8;
				break;
			case MARKER_NEWLEN:
				required_length = 6;
				break;
			case MARKER_ABORT:
			case MARKER_SDNORM:
			case MARKER_SDRST:
				required_length = 2;
				break;
			case MARKER_STUFF:
				/* forward stuffed 0xff to arithmetic decoder */
				s->buf_len = 0;
				decode_pscd(s, s->buffer, 2);
				continue;
			default:
				return JBG_EMARKER;
			}
			while (s->buf_len < required_length && *cnt < len)
				s->buffer[s->buf_len++] = data[(*cnt)++];
			if (s->buf_len < required_length)
				continue;
			/* now the buffer is filled with exactly one marker segment */
			switch (s->buffer[1])
			{
			case MARKER_COMMENT:
				s->comment_skip =
					(((long)s->buffer[2] << 24) | ((long)s->buffer[3] << 16) |
					((long)s->buffer[4] << 8) | (long)s->buffer[5]);
				break;
			case MARKER_ATMOVE:
				if (s->at_moves < JBG_ATMOVES_MAX)
				{
					s->at_line[s->at_moves] =
						(((long)s->buffer[2] << 24) | ((long)s->buffer[3] << 16) |
						((long)s->buffer[4] << 8) | (long)s->buffer[5]);
					s->at_tx[s->at_moves] = (signed char)s->buffer[6];
					s->at_ty[s->at_moves] = s->buffer[7];
					if (s->at_tx[s->at_moves] < -(int)s->mx ||
						s->at_tx[s->at_moves] > (int)s->mx ||
						s->at_ty[s->at_moves] > (int)s->my ||
						(s->at_ty[s->at_moves] == 0 && s->at_tx[s->at_moves] < 0))
						return JBG_EINVAL;
					if (s->at_ty[s->at_moves] != 0)
						return JBG_EIMPL;
					s->at_moves++;
				}
				else
					return JBG_EIMPL;
				break;
			case MARKER_NEWLEN:
				y = (((long)s->buffer[2] << 24) | ((long)s->buffer[3] << 16) |
					((long)s->buffer[4] << 8) | (long)s->buffer[5]);
				if (y > s->yd || !(s->options & JBG_VLENGTH))
					return JBG_EINVAL;
				s->yd = y;
				/* calculate again number of stripes that will be required */
				s->stripes = jbg_stripes(s->l0, s->yd, s->d);
				break;
			case MARKER_ABORT:
				return JBG_EABORT;

			case MARKER_SDNORM:
			case MARKER_SDRST:
				/* decode final pixels based on trailing zero bytes */
				decode_pscd(s, s->buffer, 2);

				arith_decode_init(s->s[s->ii[iindex[s->order & 7][PLANE]]] +
					s->ii[iindex[s->order & 7][LAYER]] - s->dl,
					s->ii[iindex[s->order & 7][STRIPE]] != s->stripes - 1 && s->buffer[1] != MARKER_SDRST);

				s->reset[s->ii[iindex[s->order & 7][PLANE]]]
					[s->ii[iindex[s->order & 7][LAYER]] - s->dl] =
					(s->buffer[1] == MARKER_SDRST);

				/* prepare for next SDE */
				s->x = 0;
				s->i = 0;
				s->pseudo = 1;
				s->at_moves = 0;

				/* increment layer/stripe/plane loop variables */
				/* start and end value for each loop: */
				is[iindex[s->order & 7][STRIPE]] = 0;
				ie[iindex[s->order & 7][STRIPE]] = s->stripes - 1;
				is[iindex[s->order & 7][LAYER]] = s->dl;
				ie[iindex[s->order & 7][LAYER]] = s->d;
				is[iindex[s->order & 7][PLANE]] = 0;
				ie[iindex[s->order & 7][PLANE]] = s->planes - 1;
				i = 2; /* index to innermost loop */
				do
				{
					j = 0; /* carry flag */
					if (++s->ii[i] > ie[i])
					{
						/* handling overflow of loop variable */
						j = 1;
						if (i > 0)
							s->ii[i] = is[i];
					}
				} while (--i >= 0 && j);

				s->buf_len = 0;

				/* check whether this have been all SDEs */
				if (j)
				{
#ifdef DEBUG
					fprintf(stderr, "This was the final SDE in this BIE, "
						"%d bytes left.\n",
						len - *cnt);
#endif
					s->bie_len = 0;
					return JBG_EOK;
				}

				/* check whether we have to abort because of xmax/ymax */
				if (iindex[s->order & 7][LAYER] == 0 && i < 0)
				{
					/* LAYER is the outermost loop and we have just gone to next layer */
					if (jbg_ceil_half(s->xd, s->d - s->ii[0]) > s->xmax ||
						jbg_ceil_half(s->yd, s->d - s->ii[0]) > s->ymax)
					{
						s->xmax = 4294967295UL;
						s->ymax = 4294967295UL;
						return JBG_EOK_INTR;
					}
					if (s->ii[0] > (unsigned long)s->dmax)
					{
						s->dmax = 256;
						return JBG_EOK_INTR;
					}
				}

				break;
			}
			s->buf_len = 0;
		}
		else if (data[*cnt] == MARKER_ESC)
			s->buffer[s->buf_len++] = data[(*cnt)++];

		else
		{

			/* we have found PSCD bytes */
			*cnt += decode_pscd(s, data + *cnt, len - *cnt);
			if (*cnt < len && data[*cnt] != 0xff)
			{
#ifdef DEBUG
				fprintf(stderr, "PSCD was longer than expected, unread bytes "
					"%02x %02x %02x %02x ...\n",
					data[*cnt], data[*cnt + 1],
					data[*cnt + 2], data[*cnt + 3]);
#endif
				return JBG_EINVAL;
			}
		}
	} /* of BID processing loop 'while (*cnt < len) ...' */

	return JBG_EAGAIN;
}

/*
 * After jbg_dec_in() returned JBG_EOK or JBG_EOK_INTR, you can call this
 * function in order to find out the width of the image.
 */
long jbg_dec_getwidth(const struct jbg_dec_state *s)
{
	if (s->d < 0)
		return -1;
	if (iindex[s->order & 7][LAYER] == 0)
	{
		if (s->ii[0] < 1)
			return -1;
		else
			return jbg_ceil_half(s->xd, s->d - (s->ii[0] - 1));
	}

	return s->xd;
}

/*
 * After jbg_dec_in() returned JBG_EOK or JBG_EOK_INTR, you can call this
 * function in order to find out the height of the image.
 */
long jbg_dec_getheight(const struct jbg_dec_state *s)
{
	if (s->d < 0)
		return -1;
	if (iindex[s->order & 7][LAYER] == 0)
	{
		if (s->ii[0] < 1)
			return -1;
		else
			return jbg_ceil_half(s->yd, s->d - (s->ii[0] - 1));
	}

	return s->yd;
}

/*
 * After jbg_dec_in() returned JBG_EOK or JBG_EOK_INTR, you can call this
 * function in order to get a pointer to the image.
 */
unsigned char *jbg_dec_getimage(const struct jbg_dec_state *s, int plane)
{
	if (s->d < 0)
		return NULL;
	if (iindex[s->order & 7][LAYER] == 0)
	{
		if (s->ii[0] < 1)
			return NULL;
		else
			return s->lhp[(s->ii[0] - 1) & 1][plane];
	}

	return s->lhp[s->d & 1][plane];
}

/*
 * After jbg_dec_in() returned JBG_EOK or JBG_EOK_INTR, you can call
 * this function in order to find out the size in bytes of one
 * bitplane of the image.
 */
long jbg_dec_getsize(const struct jbg_dec_state *s)
{
	if (s->d < 0)
		return -1;
	if (iindex[s->order & 7][LAYER] == 0)
	{
		if (s->ii[0] < 1)
			return -1;
		else
			return jbg_ceil_half(s->xd, s->d - (s->ii[0] - 1) + 3) *
			jbg_ceil_half(s->yd, s->d - (s->ii[0] - 1));
	}

	return jbg_ceil_half(s->xd, 3) * s->yd;
}

/*
 * After jbg_dec_in() returned JBG_EOK or JBG_EOK_INTR, you can call
 * this function in order to find out the size of the image that you
 * can retrieve with jbg_merge_planes().
 */
long jbg_dec_getsize_merged(const struct jbg_dec_state *s)
{
	if (s->d < 0)
		return -1;
	if (iindex[s->order & 7][LAYER] == 0)
	{
		if (s->ii[0] < 1)
			return -1;
		else
			return jbg_ceil_half(s->xd, s->d - (s->ii[0] - 1)) *
			jbg_ceil_half(s->yd, s->d - (s->ii[0] - 1)) *
			((s->planes + 7) / 8);
	}

	return s->xd * s->yd * ((s->planes + 7) / 8);
}

/*
 * The destructor function which releases any resources obtained by the
 * other decoder functions.
 */
void jbg_dec_free(struct jbg_dec_state *s)
{
	int i;
	extern char jbg_dptable[];

	if (s->d < 0 || s->s == NULL)
		return;
	s->d = -2;

	for (i = 0; i < s->planes; i++)
	{
		checked_free(s->s[i]);
		checked_free(s->tx[i]);
		checked_free(s->ty[i]);
		checked_free(s->reset[i]);
		checked_free(s->lntp[i]);
		checked_free(s->lhp[0][i]);
		checked_free(s->lhp[1][i]);
	}

	checked_free(s->s);
	checked_free(s->tx);
	checked_free(s->ty);
	checked_free(s->reset);
	checked_free(s->lntp);
	checked_free(s->lhp[0]);
	checked_free(s->lhp[1]);
	if (s->dppriv && s->dppriv != jbg_dptable)
		checked_free(s->dppriv);

	s->s = NULL;

	return;
}