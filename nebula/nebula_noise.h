// Trigonometric approximation adapted from fastapprox.

#ifndef NEBULA_NOISE_H
#define NEBULA_NOISE_H

#include "../external/ciglet/ciglet.h"

static inline float
nebula_noise_sin (double x)
{
  static const float fouroverpi = 1.2732395447351627f;
  static const float fouroverpisq = 0.40528473456935109f;
  static const float q = 0.78444488374548933f;
  union { float f; uint32_t i; } p = { 0.20363937680730309f };
  union { float f; uint32_t i; } r = { 0.015124940802184233f };
  union { float f; uint32_t i; } s = { -0.0032225901625579573f };

  union { float f; uint32_t i; } vx = { x };
  uint32_t sign = vx.i & 0x80000000;
  vx.i = vx.i & 0x7FFFFFFF;

  float qpprox = fouroverpi * x - copysignf(fouroverpisq, x) * vx.f * vx.f;
  float qpproxsq = qpprox * qpprox;

  p.i |= sign;
  r.i |= sign;
  s.i ^= sign;

  return q * qpprox + qpproxsq * (p.f + qpproxsq * (r.f + qpproxsq * s.f));
}

static inline float
nebula_noise_sinfull (double x)
{
  static const float twopi = 6.2831853071795865f;
  static const float invtwopi = 0.15915494309189534f;

  int k = x * invtwopi;
  float half = (x < 0) ? -0.5f : 0.5f;
  return nebula_noise_sin ((half + k) * twopi - x);
}

static inline float
nebula_noise_cos (double x)
{
  static const float halfpi = 1.5707963267948966f;
  return nebula_noise_sinfull (x + halfpi);
}

static inline float nebula_normal(void)
{
    float u1 = (float)rand() / RAND_MAX;
    float u2 = (float)rand() / RAND_MAX;
    return sqrt(-2.0 * log_2(u1)) * nebula_noise_cos(2.0 * M_PI * u2);
}


#endif
