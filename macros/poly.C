#include <stdio.h>
#include <stdlib.h>
#include "sslib.h"

void qurt(Complex a, Complex b, Complex c, Complex x[])
{
  Complex d, d1, d2, w;

  if(a.r == 0. && a.i == 0.)
    {
      fprintf(stderr, "Warning : a = 0  in qurt()\n");
      if(b.r == 0. && b.i == 0.)
	{
	  fprintf(stderr, "Error : a = b = 0  in qurt()\n");
	  x[0] = x[1] = tocomplex(0., 0.);
	  return;
	}
      x[0] = cdiv(c, b);
      x[1] = tocomplex(0., 0.);
      return;
    }
  d = cmul1(a, -4.);
  d = cmul(d, c);
  w = cmul(b, b);
  d = cadd(d, w);
  if(d.r == 0. && d.i == 0.)
    {
      w = tocomplex(-b.r, -b.i);
      d = cmul1(a, 2.);
      x[0] = x[1] = cdiv(w, d);
      return;
    }
  d = csqrt(d);
  d1 = d2 = tocomplex(-b.r, -b.i);
  d1 = cadd(d1, d);
  d2 = csub(d2, d);
  if(cabslt(d2) > cabslt(d1))
    {
      w = d1;
      d1 = d2;
      d2 = w;
    }
  w = cmul1(a, 2.);
  x[0] = cdiv(d1, w);
  w = cmul1(c, 2.);
  x[1] = cdiv(w, d1);
  return;
}

void carda(double a[], Complex x[])
{
  int i;
  double alfa, beta, b1, b2, b3, d, p, q, theta, u3, v3, w;

  if(a[0] == 0.)
    {
      fprintf(stderr, "Error : a[0] = 0  in carda()\n");
      return;
    }
  b1 = *(a + 1) / *a / 3.;
  b2 = *(a + 2)/ *a;
  b3 = *(a + 3)/ *a;
  p = b2 / 3. - b1 * b1;
  q = -(b1 * (2. * b1 * b1 - b2) + b3);
  d = q * q + 4. * p * p * p;
  w = sqrt(fabs(d));
  if(d >= 0.)
    {
      u3 = (q + w) / 2.;
      v3 = (q - w) / 2.;
      alfa = cbrt(u3);
      beta = cbrt(v3);
      (*x).r = alfa + beta - b1;
      (*x).i = 0.;
      (*(x + 1)).r = - (alfa + beta) / 2. - b1;
      (*(x + 1)).i = 0.86602540378443864 * (alfa - beta);
      *(x + 2) = conj(*(x + 1));
      return;
    }
  theta = M_PI_2;
  if(q != 0.)theta = atan(w / q);
  w = 2. * sqrt(-p);
  (*x).r = w * cos(theta / 3.) - b1;
  (*(x + 1)).r = - (w * cos((M_PI - theta) / 3.)) - b1;
  (*(x + 2)).r = - (w * cos((M_PI + theta) / 3.)) - b1;
  (*x).i = (*(x + 1)).i = (*(x + 2)).i = 0.;
  return;
}


double newton(double a[], int n, double eps, int iter)
{
  int i, ic, m;
  double a1, d, dx, p, q, r, w, wx;

  a1 = fabs(*a);
  if(a1 == 0. || n < 2 || eps <= 0. || iter < 1)
    {
      fprintf(stderr, "Error : Illegal parameter in newton()\n");
      return 0;
    }
  for(i = 1, wx = *a; i <= n; i++)
    {
      w = *(a + i);
      wx += w;
      if (a1 < fabs(w))a1 = fabs(w);
    }
  wx = - wx / (double)(n + 1);
  for(ic = 1; ic <= iter; ic++)
    {
      p = *a;
      for(i = 1; i <= n; i++)p = p * wx + *(a + i);
      m = n;
      q = *a * (double)m;
      for(i = 1; i < n; i++)
	{
	  m--;
	  q = q * wx + (double)m * *(a + i);
	}
      if(q != 0.)dx = - p / q;
      else
	{
	  m = n;
	  r = *a * (double)(m * (n - 1));
	  for(i = 1; i < n - 1; i++)
	    {
	      m--;
	      r = r * wx + (double)(m * (m - 1)) * *(a + i);
	    }
	  d = p * (p - 2. * r);
	  if(d <= 0.)dx = - q / p;
	  else
	    {
	      w = sqrt(d);
	      if(q >= 0.)dx = 2. * p / (- q - w);
	      else dx = 2. * p / (- q + w);
	    }
	}
      wx += dx;
      if(fabs(dx) / a1 <= eps)return wx;
    }
  fprintf(stderr, "Error : No convergence in newton()\n");
  return wx;
}

void bairs(double a[], int n, double eps, int iter, Complex x[])
{
  int flag, i, ic, j, k, nm, nm1, nn;
  double amax, bk, b1, b2, ck, c1, c2, c3, d, p, q, w, xx;

  for(i = 0; i < n; i++)*(x + i) = tocomplex(0., 0.);
  if(*a == 0. || n < 2 || eps <= 0. || iter < 1)
    {
      fprintf(stderr, "Error : Illegal parameter  in bairs()\n");
      return;
    }
  nm = n;
  nm1 = n - 1;
  if((n % 2) != 0)
    {
      xx = newton(a, n, eps, iter);
      *(x + nm1) = tocomplex(xx, 0.);
      w = 0.;
      for(i = 0; i <= nm; i++)
	{
	  w = w * xx + *(a + i);
	  *(a + i) = w;
	}
      nm = nm1;
      nm1--;
    }
  amax = fabs(*a);
  for(i = 1; i <= nm; i++)
    if (amax < fabs(*(a + i)))amax = fabs(*(a + i));
  i = nm1;
  while(i >= 3)
    {
      p = q = 1.;
      flag = 0;
      for(ic = 0; ic < iter; ic++)
	{
	  b1 = b2 = c1 = c2 = 0.;
	  nn = i + 1;
	  for(k = 0; k <= nn; k++)
	    {
	      bk = *(a + k) - (p * b1 + q * b2);
	      ck = bk       - (p * c1 + q * c2);
	      if(k == (nn - 3))c3 = ck;
	      if(k != nn)
		{
		  b2 = b1;
		  b1 = bk;
		  c2 = c1;
		  c1 = ck;
		}
	    }
	  c1 = - (p * c2 + q * c3);
	  d = c2 * c2 - c1 * c3;
	  if(d != 0.)
	    {
	      p += (b1 * c2 - bk * c3) / d;
	      q += (bk * c2 - b1 * c1) / d;
	    }
	  else
	    {
	      p++;
	      q++;
	    }
	  if((fabs(b1) + fabs(bk + p * b1)) / amax <= eps)
	    {
	      flag = 1;
	      break;
	    }
	}
      if(!flag)
	{
	  fprintf(stderr, "Error : No convergence in bairs()\n");
	  return;
	}
      b1 = b2 = 0.;
      j = i - 1;
      for(k = 0; k <= j; k++)
	{
	  bk = *(a + k) - (p * b1 + q * b2);
	  *(a + k) = bk;
	  b2 = b1;
	  b1 = bk;
	}
      *(a + i) = p;
      *(a + i + 1) = q;
      i -= 2;
    }
  *(a + 1) = *(a + 1) / *a;
  *(a + 2) = *(a + 2) / *a;
  for(i = 1; i <= nm1; i += 2)
    {
      d = *(a + i) * *(a + i)- 4. * *(a + i + 1);
      w = sqrt(fabs(d));
      if(d >= 0.)
	{
	  *(x + i - 1) = tocomplex(- (*(a + i) - w) / 2., 0.);
	  *(x + i)     = tocomplex(- (*(a + i) + w) / 2., 0.);
	}
      else
	{
	  *(x + i - 1) = tocomplex(- *(a + i) / 2., w / 2.);
	  *(x + i)     = conj(*(x + i - 1));
	}
    }
  return;
}



int cnewton(Complex a[], Complex r[], int n, double eps, int iter)
{
  Complex *d, w, x, p, q;
  int root, nn, i, ic;
  double ww, amax, am;

  d = (Complex *)malloc((n + 1)*sizeof(Complex));
  if(d == NULL)
    {
      fprintf(stderr, "Error : Out of Memory in cnwtn()\n");
      return(997);
    }

  if(ceq(a[0], tocomplex(0.0, 0.0)) || n < 2 || eps <= 0.0 || iter < 1)
    {
      fprintf(stderr, "Error : Illigal parameter in cnwtn()\n");
      free((char *)d);
      return(999);
    }

  root = 0;
  nn = n;

  w = a[0];
  if(cne(w, tocomplex(1.0, 0.0)))
    {
      r[0] = tocomplex(1.0, 0.0);
      for(i = 1; i <= n; i++)r[i] = cdiv(a[i],w);
      /* normalizing coefficients */
    }
  else for(i = 0; i <= n; i++)r[i] = a[i];

  while(ceq(r[nn], tocomplex(0.0, 0.0)))
    {
      root++;
      nn--;
    }

  while(nn > 1)
    {
      amax = 1.0;
      for(i = 1; i <= nn; i++)
	{
	  ww = r[i].r * r[i].r + r[i].i * r[i].i;
	  if (amax < ww) amax = ww;
	}

      for(i = 0, am = nn; i < nn; i++, am--)
	d[i] = tocomplex(am * r[i].r, am * r[i].i);

      x = tocomplex(0.0, -1.0);

      for(ic = 0; ic < iter; ic++)
	{
	  p = cfunc(r, nn, x);
	  if(sqrt((p.r * p.r + p.i * p.i) / amax) <= eps)break;
	  q = cfunc(d, nn - 1, x); 
	  w = cdiv(p, q);
	  x = csub(x, w);
	}
      if(ic >= iter)
	{
	  free((char *)d);
	  if (nn == n)return(-1);/* no convergence , no root */
	  return(root);/* no convergence , any root */
	}
      r[nn] = x;
      nn--;
      root++;
      for(i = 1; i <= nn; i++)
	{
	  w = cmul(r[i - 1], x);
	  r[i] = cadd(r[i], w);
	}
    }
  w = cdiv(r[1], r[0]);
  r[1] = tocomplex(-w.r, -w.i);
  for(i = 1; i <= n; i++)r[i - 1] = r[i];
  free((char *)d);
  return 0;
}

int dka(Complex a[], Complex r[], int n, double eps, int iter)
{
  Complex b, *c, f1, f2, *p, *q, *ratio, ri, z, zmax;
  double cf, di, dn, p0, p1, rmax;
  int flag, i, j, k;

  if(n < 1 || eps <= 0. || iter < 1)
    {
      fprintf(stderr, "Error : Illegal parameter  in dka()\n");
      return 999;
    }
  ratio = (Complex *)malloc((n + 1) * sizeof(Complex));
  if(ratio == NULL)
    {
      fprintf(stderr, "Error : Out of Memory  in dka()!!\n");
      return 997;
    }
  c = (Complex *)malloc((n + 1) * sizeof(Complex));
  if(c == NULL)
    {
      fprintf(stderr, "Error : Out of Memory  in dka()!!\n");
      free((char *)ratio);
      return 997;
    }


  /* Normalize coefficients */

  z = *a;
  if(z.r != 1. || z.i != 0.) {
     for(i = 0, p = a, q = c; i <= n; i++)*q++ = cdiv(*p++, z);
  }
  else 
    for(i = 0, p = a, q = c; i <= n; i++)*q++ = *p++;

  /* Get initial roots */

  dn = (double)n;
  rmax = 0.;
  b.r = (*(c + 1)).r / dn;
  b.i = (*(c + 1)).i / dn;
  for(i = 2, di = 2., p = c + 2; i <= n; i++, di++, p++)
    {
      cf = pow(cabslt(*p), 1.0 / di);
      if(rmax < cf)rmax = cf;
    }

  p0 = 1.5 / dn;
  p1 = 2.0 * M_PI / dn;
  zmax = tocomplex(rmax, 0.);
  for(i = 1, di = 0.; i <= n; i++, di++)
    {
      z = cexp(tocomplex(0.0, p0 + di * p1));
      *(r + i) = csub(cmul(z, zmax), b);
    }

  /* Main loop */

  k = 0;
  while(k++ < iter)
    {
      for(i = 1; i <= n; i++)
	{
	  f1 = cfunc(c, n, (ri = *(r + i)));
	  f2 = tocomplex(1.0, 0.0);
	  for(j = 1; j <= n; j++)
	    if (j != i)f2 = cmul(f2, csub(ri, *(r + j)));
	  *(ratio + i) = cdiv(f1, f2);
	  *(r + i) = csub(ri, *(ratio + i));
	}

      /* test convergence */

      flag = 1;
      for(i = 1; i <= n; i++)
	{
	  if(cabslt(*(ratio + i)) > eps)
	    {
	      flag = 0;
	      break;
	    }
	}
      if(flag)
	{
	  free((char *)c);
	  free((char *)ratio);
	  return 0;
	}
    }
  fprintf(stderr, "Error :  No Convergence in dka()\n");
  free((char *)c);
  free((char *)ratio);
  return -1;
} 
