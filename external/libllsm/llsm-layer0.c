/*
  libllsm - Low Level Speech Model
  ===
  Copyright (c) 2015-2017 Kanru Hua.

  libllsm is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libllsm is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with libllsm. If not, see <http://www.gnu.org/licenses/>.
*/

//#define DEBUG_EXPORT_WAV

#include "llsm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "envelope.h"
#include "math-funcs.h"
#include "../libgvps/gvps.h"

#ifdef DEBUG_EXPORT_WAV
static void gain_write(FP_TYPE* x, int nx, int fs, char* path, FP_TYPE gain) {
  FP_TYPE* normx = calloc(nx, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++) normx[i] = x[i] * gain;
  wavwrite(normx, nx, fs, 16, path);
  free(normx);
}
#endif

llsm_parameters llsm_init(int nnosband) {
  llsm_parameters ret;
  ret.a_nhop = 256;
  ret.a_nhar = 80;
  ret.a_nhare = 4;
  ret.a_nnos = 54;
  ret.a_nosf = -1.0;
  ret.a_f0refine = 1;
  ret.a_suppress_subharmonics = 0;
  ret.a_unused1 = 1;
  ret.a_unused2 = 1;
  ret.a_wsize = 2.2;
  ret.a_tfft = 0.04;
  ret.a_mvf = 8000.0;
  ret.a_noswarp = 15000.0;
  ret.a_nnosband = nnosband;
  ret.a_nosbandf = calloc(nnosband - 1, sizeof(FP_TYPE));
  ret.a_nosbandf[0] = 2000;
  for(int i = 1; i < nnosband - 1; i ++) ret.a_nosbandf[i] = ret.a_nosbandf[i - 1] * 2.0;
  ret.s_fs = 0;
  ret.s_use_sinu_phase = 0;
  ret.s_noiseonly = 0;
  ret.s_n0 = 0;
  ret.s_n1 = 0;
  return ret;
}

void llsm_deinit(llsm_parameters dst) {
  if(dst.a_nosbandf != NULL)
    free(dst.a_nosbandf);
}

// compute magnitude, phase and delayed phase spectrogram on input signal
void spectrogram_analyze(FP_TYPE* x, int nx, int* t, int nfrm, int* winsize, int nfft, FP_TYPE fs,
  const char* wintype, FP_TYPE** spectrogram, FP_TYPE** phasegram, FP_TYPE** phasegram_d, FP_TYPE* dc) {
  
# pragma omp parallel for
  for(int n = 0; n < nfrm; n ++) {
    FP_TYPE fftbuff[65536];
    FP_TYPE inoutbuff[65536];
    FP_TYPE* xbuff = inoutbuff + 0;
    FP_TYPE* ybuffr = inoutbuff + nfft;
    FP_TYPE* ybuffi = inoutbuff + nfft * 2;
    int winlen = min(nfft, winsize[n]);
    int tn = t[n];
    
    FP_TYPE* w = NULL;
    if(! strcmp(wintype, "blackman_harris"))
      w = blackman_harris(winlen);
    else if(! strcmp(wintype, "hamming"))
      w = hamming(winlen);
    else if(! strcmp(wintype, "hanning"))
      w = hanning(winlen);
    else if(! strcmp(wintype, "mltsine"))
      w = mltsine(winlen);
    else if(! strcmp(wintype, "blackman"))
      w = blackman(winlen);
    FP_TYPE norm_factor = 2.0 / sumfp(w, winlen);

    FP_TYPE* xfrm, * xfrm_d;
    FP_TYPE* spec_magn, * spec_phse_, * spec_phse, * spec_phse_d;
    xfrm = xfrm_d = spec_magn = spec_phse_ = spec_phse = spec_phse_d = NULL;
    
    xfrm = fetch_frame(x, nx, tn, winlen);
    FP_TYPE mean_xfrm = sumfp(xfrm, winlen) / winlen;
    if(dc != NULL) dc[n] = minfp(xfrm, winlen);
    for(int i = 0; i < winlen; i ++)
      xfrm[i] -= mean_xfrm;
    if(phasegram_d) {
      xfrm_d = fetch_frame(x, nx, tn - 1, winlen);
      FP_TYPE mean_xfrm_d = sumfp(xfrm_d, winlen) / winlen;
      for(int i = 0; i < winlen; i ++)
        xfrm_d[i] -= mean_xfrm_d;
    }
    
    for(int i = 0; i < winlen; i ++) {
      xfrm[i] *= w[i];
      if(phasegram_d) xfrm_d[i] *= w[i];
    }
    
    // current frame
    memset(xbuff, 0, nfft * sizeof(FP_TYPE));
    for(int i = 0; i < winlen / 2; i ++) {
      xbuff[i] = xfrm[i + winlen / 2];
      xbuff[nfft - winlen / 2 + i] = xfrm[i];
    }
    fft(xbuff, NULL, ybuffr, ybuffi, nfft, fftbuff);
    
    spec_magn = abscplx(ybuffr, ybuffi, nfft / 2);
    if(phasegram) {
      spec_phse_ = argcplx(ybuffr, ybuffi, nfft / 2);
      spec_phse = unwrap(spec_phse_, nfft / 2); free(spec_phse_);
    }
    
    // delayed frame
    if(phasegram_d) {
      memset(xbuff, 0, nfft * sizeof(FP_TYPE));
      for(int i = 0; i < winlen / 2; i ++) {
        xbuff[i] = xfrm_d[i + winlen / 2];
        xbuff[nfft - winlen / 2 + i] = xfrm_d[i];
      }
      fft(xbuff, NULL, ybuffr, ybuffi, nfft, fftbuff);
    
      spec_phse_ = argcplx(ybuffr, ybuffi, nfft / 2);
      spec_phse_d = unwrap(spec_phse_, nfft / 2); free(spec_phse_);
    }
    
    for(int i = 0; i < nfft / 2; i ++) {
      spectrogram[n][i] = log_3(spec_magn[i] * norm_factor);
      if(isnan(spectrogram[n][i]) || isinf(spectrogram[n][i]))
        spectrogram[n][i] = -100.0;
      if(phasegram) phasegram[n][i] = spec_phse[i];
      if(phasegram_d) phasegram_d[n][i] = spec_phse_d[i];
    }

    free(spec_magn);
    if(spec_phse) free(spec_phse);
    if(spec_phse_d) free(spec_phse_d);
    free(xfrm);
    if(xfrm_d) free(xfrm_d);
    free(w);
  }
}

// a simple routine for refining F0 trajectory
// It is recommended to use a dedicated external routine for F0 refinement instead.
static FP_TYPE* refine_f0(llsm_parameters param, int nfft, int fs, FP_TYPE* f0, int nf0,
  FP_TYPE** spectrogram, FP_TYPE** phasegram, FP_TYPE** phasegram_d) {
  
  FP_TYPE* rf0 = calloc(nf0, sizeof(FP_TYPE));
  if(param.a_f0refine == 0) {
    for(int t = 0; t < nf0; t ++)
      rf0[t] = f0[t];
    return rf0;
  }

  for(int t = 0; t < nf0; t ++) {
    if(f0[t] <= 0) {
      rf0[t] = 0;
      continue;
    }
    int i_l = floor(f0[t] / fs * nfft * 0.7);
    int i_h = ceil(f0[t] / fs * nfft * 1.3);
    int i_max = round(f0[t] / fs * nfft);
    FP_TYPE peak = spectrogram[t][i_max];
    
    for(int i = i_l; i <= i_h; i ++)
      if(spectrogram[t][i] > peak) {
        peak = spectrogram[t][i];
        i_max = i;
      }
    
    FP_TYPE peak_freq;
    qifft(spectrogram[t], i_max, & peak_freq);
    peak_freq *= (FP_TYPE)fs / nfft;
    
    FP_TYPE p   = phasegram  [t][i_max];
    FP_TYPE p_d = phasegram_d[t][i_max];
    p   -= floor(p   / 2.0 / M_PI) * 2.0 * M_PI;
    p_d -= floor(p_d / 2.0 / M_PI) * 2.0 * M_PI;
    if(p < p_d)
      p += 2.0 * M_PI;
    
    rf0[t] = (p - p_d) / 2.0 / M_PI * fs;
    if(fabs(rf0[t] / fs * nfft - i_max) > 1)
      rf0[t] = peak_freq;
  }
  
  // take out spurious fluctuations
  for(int t = 1; t < nf0 - 1; t ++) {
    if(rf0[t] <= 0) continue;
    if(rf0[t + 1] > 0 && fabs(rf0[t + 1] - rf0[t]) > 10) {
      if(rf0[t - 1] == 0) rf0[t] = rf0[t + 1];
      else rf0[t] = (rf0[t - 1] + rf0[t + 1]) / 2.0;
    }
  }
  
  return rf0;
}

// estimate harmonic parameters by peak-picking
static int harmonic_peak_picking(llsm_parameters param, FP_TYPE* spectrum, FP_TYPE* phasetrum,
  int nfft, int fs, int nhar, FP_TYPE f0, FP_TYPE* freq, FP_TYPE* ampl, FP_TYPE* phse) {
  const FP_TYPE hardev = 0.3;

  if(f0 < 10) return 0;
  for(int i = 1; i <= nhar; i ++) {
    if(f0 * i > param.a_mvf || f0 * i > fs / 2)
      return i - 1;
    
    int l_idx = round(f0 * (i - hardev) / fs * nfft);
    int u_idx = round(f0 * (i + hardev) / fs * nfft);
    l_idx = l_idx < 1 ? 1 : l_idx;
    u_idx = u_idx > nfft / 2 - 1 ? nfft / 2 - 1 : u_idx;

    int peak_bin = cig_find_peak(spectrum, l_idx, u_idx, 1);
    
    FP_TYPE peak_freq, peak_ampl;
    peak_ampl = qifft(spectrum, peak_bin, & peak_freq);
    freq[i - 1] = peak_freq * fs / nfft;
    ampl[i - 1] = exp_3(peak_ampl);
    phse[i - 1] = linterp(phasetrum[(int)peak_freq], phasetrum[(int)peak_freq + 1], fmod(peak_freq, 1.0));
    
    if(isnan(ampl[i - 1])) // peak amplitude goes nan if one of the bins = -INF
      ampl[i - 1] = 0;
  }
  return nhar;
}

// generate a bunch of pitch-synchronous time marks over the entire signal
static int* gen_ps_instants(FP_TYPE* f0, int nhop, int nf0, int nx, int fs, int* ni) {
  FP_TYPE t = 0;
  *ni = 0;
  int ci = 100;
  int* instants = calloc(100, sizeof(int));
  while(t < (FP_TYPE)nx / fs) {
    int i = t * fs / nhop;
    instants[(*ni) ++] = t * fs;
    if((*ni) * 2 > ci) {
      ci *= 2;
      instants = realloc(instants, ci * sizeof(int));
    }
    t += 1.0 / max(50, f0[i >= nf0 ? nf0 - 1 : i]);
  }
  return instants;
}

// recolor noise by a time-varying warped PSD.
static FP_TYPE* filter_noise(llsm_parameters param, llsm_conf conf, FP_TYPE* x, int nx,
  FP_TYPE** warpped_spectrogram) {
  int nhop = conf.nhop;
  int nwin = 2 * nhop;
  if(nwin < 512)
    nwin = 512;
  int ns = nwin / 2 + 1;
  int step = nwin / 4;
  if(nhop <= step)
    step = nhop;

  int nfrm = conf.nfrm;
  int frm0 = (param.s_n0 - nwin) / step;
  if(frm0 < 0)
    frm0 = 0;
  int frm1 = nhop / step * nfrm;
  if(param.s_n1 != 0) {
    int hi = (param.s_n1 + nwin) / step;
    if(frm1 <= hi)
      hi = frm1;
    frm1 = hi;
  }

  FP_TYPE* y = calloc(nfrm * nhop + 4 * nwin, sizeof(FP_TYPE));
  FP_TYPE* w = hanning(nwin);
  FP_TYPE hop_norm = 0;
  for(int i = 0; i < nwin; i += step)
    hop_norm += w[i];
  for(int i = 0; i < nwin; i ++)
    w[i] = sqrt(w[i]);
  FP_TYPE weight_factor = 0.5 * sumfp(w, nwin);

  FP_TYPE* realbuff = calloc(nwin, sizeof(FP_TYPE));
  FP_TYPE* imagbuff = calloc(nwin, sizeof(FP_TYPE));
  FP_TYPE* yfrm = calloc(nwin, sizeof(FP_TYPE));
  FP_TYPE* fftbuff = calloc(nwin * 2, sizeof(FP_TYPE));
  FP_TYPE* freqwarp = llsm_warp_freq(0, conf.nosf, conf.nnos, conf.noswarp);

  for(int i = frm0; i < frm1; i ++) {
    int t = step * i;
    FP_TYPE* xwin = malloc(nwin * sizeof(FP_TYPE));

    for(int j = 0; j < nwin; j ++) {
      int isrc = t + j - nwin / 2;
      FP_TYPE v = (isrc >= 0 && isrc < nx) ? x[isrc] : 0;
      xwin[j] = v * w[j];
    }
    fft(xwin, NULL, realbuff, imagbuff, nwin, fftbuff);

    FP_TYPE* xA = abscplx(realbuff, imagbuff, ns);
    for(int j = 0; j < ns; j ++)
      xA[j] = xA[j] / weight_factor;
    FP_TYPE* xE = llsm_geometric_envelope(xA, nwin, param.s_fs, freqwarp, conf.nnos);
    int gidx = step * i / nhop;
    for(int j = 0; j < conf.nnos; j ++)
      xE[j] = warpped_spectrogram[gidx][j] - log(xE[j] + M_EPS);
    FP_TYPE* spec = llsm_spectrum_from_envelope(freqwarp, xE, conf.nnos, nwin, param.s_fs);

    for(int j = 0; j < ns; j ++) {
      FP_TYPE jgain = exp_2(spec[j]) + M_EPS;
      realbuff[j] *= jgain;
      imagbuff[j] *= jgain;
    }
    free(xE);
    free(xA);

    complete_symm (realbuff, nwin);
    complete_asymm(imagbuff, nwin);
    ifft(realbuff, imagbuff, yfrm, NULL, nwin, fftbuff);
    for(int j = 0; j < nwin; j ++) {
      int idx = t + j - nwin / 2;
      if(idx >= 0)
        y[idx] += w[j] * yfrm[j] / hop_norm;
    }
    free(spec);
    free(xwin);
  }

  free(w);
  free(yfrm);
  free(realbuff);
  free(imagbuff);
  free(fftbuff);
  free(freqwarp);
  return y;
}

// subtract lower envelope from a periodic signal
static void subtract_minimum_envelope(FP_TYPE* x, int nx, FP_TYPE* f0, int nhop, int nfrm, int fs) {
  int ninstant = 0;
  int* env_instants = gen_ps_instants(f0, nhop, nfrm, nx, fs, & ninstant);
  FP_TYPE* fenv_instants = calloc(ninstant, sizeof(FP_TYPE));
  int* env_winlen = calloc(ninstant, sizeof(int));
  for(int i = 0; i < ninstant; i ++) env_winlen[i] = fs / max(100, f0[i >= nfrm ? nfrm - 1 : i]) * 1.5;
  FP_TYPE* env_samples = llsm_nonuniform_envelope(x, nx, env_instants, env_winlen, ninstant, 0);
  FP_TYPE* iota = calloc(nx, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++) iota[i] = i;
  for(int i = 0; i < ninstant; i ++) fenv_instants[i] = env_instants[i];
  FP_TYPE* env = interp1(fenv_instants, env_samples, ninstant, iota, nx);
  for(int i = 0; i < nx; i ++) x[i] -= env[i];

  free(env_instants);
  free(fenv_instants);
  free(env_winlen);
  free(env_samples);
  free(iota);
  free(env);
}

// lengthen a noise signal by duplication
static FP_TYPE* stretch_static_noise(FP_TYPE* x, int nx, int ny, int novl) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  for(int i = 0; i < min(nx, ny); i ++)
    y[i] = x[i];
  if(ny <= nx) return y;
  
  int top = nx;
  while(1) {
    for(int i = 0; i < novl; i ++) {
      FP_TYPE r = (FP_TYPE)i / novl;
      y[top - novl + i] *= 1.0 - r;
      y[top - novl + i] += x[i] * r;
    }
    for(int i = 0; i < nx - novl; i ++) {
      if(top + i >= ny) return y;
      y[top + i] = x[i + novl];
    }
    top += nx - novl;
  }
  
  return y;
}

llsm_frame* llsm_create_frame(int nhar, int nhare, int nnos, int nchannel) {
  llsm_frame* ret = malloc(sizeof(llsm_frame));

  ret -> sinu = malloc(sizeof(llsm_sinframe));
  ret -> sinu -> ampl = calloc(nhar, sizeof(FP_TYPE));
  ret -> sinu -> phse = calloc(nhar, sizeof(FP_TYPE));
  ret -> sinu -> nhar = nhar;
  
  ret -> noise = malloc(sizeof(llsm_nosframe));
  ret -> noise -> eenv = calloc(nchannel, sizeof(llsm_sinframe*));
  ret -> noise -> emin = calloc(nchannel, sizeof(FP_TYPE));
  ret -> noise -> spec = calloc(nnos, sizeof(FP_TYPE));
  ret -> noise -> nchannel = nchannel;
  ret -> noise -> nnos = nnos;
  
  for(int i = 0; i < nchannel; i ++) {
    ret -> noise -> eenv[i] = malloc(sizeof(llsm_sinframe));
    llsm_sinframe* ieenv = ret -> noise -> eenv[i];
    ieenv -> ampl = calloc(nhare, sizeof(FP_TYPE));
    ieenv -> phse = calloc(nhare, sizeof(FP_TYPE));
    ieenv -> nhar = nhare;
  }

  ret -> wave = NULL;
  ret -> f0 = 0;
  return ret;
}

llsm_layer0* llsm_create_empty_layer0(llsm_conf conf) {
  llsm_layer0* ret = calloc(1, sizeof(llsm_layer0));
  ret -> conf = conf;
  ret -> conf.nosbandf = calloc(conf.nnosband - 1, sizeof(FP_TYPE));
  for(int i = 0; i < conf.nnosband - 1; i ++)
    ret -> conf.nosbandf[i] = conf.nosbandf[i];
  if(conf.nfrm > 0)
    ret -> frames = calloc(conf.nfrm, sizeof(llsm_frame*));
  else
    ret -> frames = NULL;
  return ret;
}

llsm_wave* llsm_create_wave(FP_TYPE* samples, int nsample, int fs) {
  llsm_wave* ret = calloc(1, sizeof(llsm_wave));
  ret -> samples = calloc(nsample, sizeof(FP_TYPE));
  ret -> nsample = nsample;
  ret -> fs = fs;
  if(samples != NULL)
    for(int i = 0; i < nsample; i ++)
      ret -> samples[i] = samples[i];
  return ret;
}

void llsm_copy_wave(llsm_wave* dst, llsm_wave* src) {
  if(dst -> nsample < src -> nsample)
    dst -> samples = realloc(dst -> samples, src -> nsample * sizeof(FP_TYPE));
  for(int i = 0; i < src -> nsample; i ++)
    dst -> samples[i] = src -> samples[i];
  dst -> nsample = src -> nsample;
  dst -> fs = src -> fs;
  dst -> t0 = src -> t0;
}

void llsm_copy_frame(llsm_frame* dst, llsm_frame* src) {
  llsm_copy_sinframe(dst -> sinu, src -> sinu);
  llsm_copy_nosframe(dst -> noise, src -> noise);
  if(src -> wave != NULL) {
    if(dst -> wave == NULL)
      dst -> wave = llsm_create_wave(src -> wave -> samples,
        src -> wave -> nsample, src -> wave -> fs);
    llsm_copy_wave(dst -> wave, src -> wave);
  }
  dst -> f0 = src -> f0;
}

void llsm_copy_sinframe(llsm_sinframe* dst, llsm_sinframe* src) {
  if(dst -> nhar < src -> nhar) {
    dst -> ampl = realloc(dst -> ampl, src -> nhar * sizeof(FP_TYPE));
    dst -> phse = realloc(dst -> phse, src -> nhar * sizeof(FP_TYPE));
  }
  dst -> nhar = src -> nhar;
  for(int i = 0; i < dst -> nhar; i ++) {
    dst -> ampl[i] = src -> ampl[i];
    dst -> phse[i] = src -> phse[i];
  }
}

void llsm_copy_nosframe(llsm_nosframe* dst, llsm_nosframe* src) {
  if(dst -> nnos < src -> nnos)
    dst -> spec = realloc(dst -> spec, src -> nnos * sizeof(FP_TYPE));
  dst -> nnos = src -> nnos;
  for(int i = 0; i < dst -> nnos; i ++)
    dst -> spec[i] = src -> spec[i];
  if(dst -> nchannel < src -> nchannel) {
    dst -> emin = realloc(dst -> emin, src -> nchannel * sizeof(FP_TYPE));
    dst -> eenv = realloc(dst -> eenv, src -> nchannel * sizeof(llsm_sinframe*));
    for(int i = dst -> nchannel; i < src -> nchannel; i ++) {
      dst -> eenv[i] = malloc(sizeof(llsm_sinframe));
      dst -> eenv[i] -> nhar = 0;
      dst -> eenv[i] -> ampl = NULL;
      dst -> eenv[i] -> phse = NULL;
    }
  } else if(dst -> nchannel > src -> nchannel)
    for(int i = src -> nchannel; i < dst -> nchannel; i ++) {
      free(dst -> eenv[i] -> ampl);
      free(dst -> eenv[i] -> phse);
      free(dst -> eenv[i]);
    }
  dst -> nchannel = src -> nchannel;
  for(int i = 0; i < dst -> nchannel; i ++) {
    dst -> emin[i] = src -> emin[i];
    llsm_copy_sinframe(dst -> eenv[i], src -> eenv[i]);
  }
}

void llsm_delete_frame(llsm_frame* dst) {
  if(dst == NULL) return;
  if(dst -> sinu != NULL) {
    free(dst -> sinu -> ampl);
    free(dst -> sinu -> phse);
    free(dst -> sinu);
  }
  if(dst -> noise != NULL) {
    for(int i = 0; i < dst -> noise -> nchannel; i ++) {
      llsm_sinframe* ieenv = dst -> noise -> eenv[i];
      if(ieenv == NULL) continue;
      free(ieenv -> ampl);
      free(ieenv -> phse);
      free(ieenv);
    }
    free(dst -> noise -> eenv);
    free(dst -> noise -> emin);
    free(dst -> noise -> spec);
    free(dst -> noise);
  }
  if(dst -> wave != NULL) {
    free(dst -> wave -> samples);
    free(dst -> wave);
  }
  free(dst);
}

void llsm_delete_layer0(llsm_layer0* dst) {
  if(dst == NULL) return;
  int nfrm = dst -> conf.nfrm;
  if(dst -> conf.nosbandf != NULL) free(dst -> conf.nosbandf);
  if(dst -> frames != NULL) {
    for(int i = 0; i < nfrm; i ++)
      llsm_delete_frame(dst -> frames[i]);
    free(dst -> frames);
  }
  free(dst);
}

void llsm_delete_output(llsm_output* dst) {
  if(dst == NULL) return;
  free(dst -> y);
  free(dst -> y_sin);
  free(dst -> y_nos);
  free(dst -> instants);
  free(dst);
}

// Viterbi neighborhood from context[1].
static int llsm_viterbi_ntran(void* context, int t) {
  (void)t;
  return ((int*)context)[1];
}

static FP_TYPE llsm_gaussian_window_sample(void* context, int ds, int t) {
  (void)t;
  FP_TYPE x = (FP_TYPE)ds / (FP_TYPE)((int*)context)[0];
  return (FP_TYPE)exp((double)(x * -x * 100.0f));
}

// mean abs wrapped residual of successive harmonic phases vs candidate offset
static FP_TYPE llsm_phase_incoherence(FP_TYPE* phse, int nhar, FP_TYPE cand) {
  FP_TYPE sum = 0.0f;
  for(int k = 1; k < nhar; k ++) {
    FP_TYPE d = phse[k] - phse[k - 1] - phse[0] + cand;
    FP_TYPE n = (FP_TYPE)round((double)d / (2.0 * M_PI));
    d = d - (n + n) * (FP_TYPE)M_PI;
    sum += (FP_TYPE)fabs((double)d);
  }
  return sum / (FP_TYPE)(nhar - 1);
}

// 32-state circular Viterbi smoothing of per-frame harmonic phase offset
static FP_TYPE* llsm_smooth_phase_offset(FP_TYPE** phse_har, llsm_layer0* model) {
  int nfrm = model -> conf.nfrm;
  int nstate = 32;
  FP_TYPE* axis = calloc(nstate, sizeof(FP_TYPE));
  for(int i = 0; i < nstate; i ++)
    axis[i] = (FP_TYPE)i * 0.03125f * (FP_TYPE)M_PI
            + (FP_TYPE)i * 0.03125f * (FP_TYPE)M_PI
            - (FP_TYPE)M_PI;

  gvps_obsrv* obs = gvps_obsrv_create_full(nfrm, nstate);
  for(int i = 0; i < nfrm; i ++) {
    int nhar = model -> frames[i] -> sinu -> nhar;
    for(int k = 0; k < nstate; k ++) {
      if(nhar > 3) {
        int n = nhar > 10 ? 10 : nhar;
        FP_TYPE e = -llsm_phase_incoherence(phse_har[i], n, axis[k]);
        obs -> slice[i] -> pair[k].p = (FP_TYPE)exp((double)e);
      } else {
        obs -> slice[i] -> pair[k].p = 1.0f;
      }
    }
  }

  int ctx[2] = {nstate, 4};
  int* dst = calloc(nfrm, sizeof(int));
  gvps_full_circular(dst, ctx, obs, llsm_gaussian_window_sample, llsm_viterbi_ntran);

  FP_TYPE* out = calloc(nfrm, sizeof(FP_TYPE));
  FP_TYPE* sbuf = calloc(nfrm, sizeof(FP_TYPE));
  FP_TYPE* cbuf = calloc(nfrm, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE voiced = model -> frames[i] -> f0 > 0.0f ? 1.0f : 0.0f;
    out[i] = axis[dst[i]];
    sbuf[i] = (FP_TYPE)sin((double)out[i]) * voiced;
    cbuf[i] = (FP_TYPE)cos((double)out[i]) * voiced;
  }
  FP_TYPE* savg = cig_moving_avg(sbuf, nfrm, 5.0f);
  FP_TYPE* cavg = cig_moving_avg(cbuf, nfrm, 5.0f);
  for(int i = 0; i < nfrm; i ++) {
    if(model -> frames[i] -> f0 > 0.0f)
      out[i] = (FP_TYPE)atan2((double)cavg[i], (double)savg[i]); // atan2(cos, sin)
  }

  gvps_obsrv_free(obs);
  free(axis);
  free(dst);
  free(sbuf);
  free(cbuf);
  free(savg);
  free(cavg);
  return out;
}

// map per-frame F0 and first-harmonic phase onto glottal-epoch sample indices
static int* llsm_phase_to_instants(
    FP_TYPE* ph0, FP_TYPE* ph_off, llsm_layer0* model, FP_TYPE fs, int* ninst) {
  int nfrm = model -> conf.nfrm;
  int nhop = model -> conf.nhop;
  int* inst = calloc((size_t)nhop * (size_t)nfrm, sizeof(int));
  FP_TYPE cursor = 0.0f;

  *ninst = 0;
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE f0 = model -> frames[i] -> f0;
    if(f0 <= 0.0f) {
      cursor = (FP_TYPE)(nhop * i);
      continue;
    }
    FP_TYPE t0 = fs / f0;
    FP_TYPE ph = ph0[i];
    ph = ph - (FP_TYPE)(2.0 * M_PI) * (FP_TYPE)round((double)ph / (2.0 * M_PI));
    FP_TYPE epoch = (FP_TYPE)(nhop * i)
          - ph * 0.5f / (FP_TYPE)M_PI * t0
          + 0.5f * ph_off[i] / (FP_TYPE)M_PI * t0;
    FP_TYPE delta = epoch - cursor;
    int n = (int)round((double)delta / (double)t0);
    if(n > 0) {
      FP_TYPE step = delta / (FP_TYPE)n;
      for(int k = 0; k < n; k ++) {
        cursor += step;
        inst[*ninst] = (int)round((double)cursor);
        *ninst += 1;
      }
      cursor = epoch;
    }
  }
  return inst;
}

static void llsm_layer0_analyze_noise(llsm_parameters param, llsm_layer0* model, FP_TYPE* x,
    FP_TYPE* residual, int nx, int fs, int* analysis_instant, FP_TYPE* rf0, int nf0) {
  int nfft = pow(2, ceil(log2(fs * param.a_tfft)));
  int ns = nfft / 2 + 1;
  FP_TYPE nosf = param.a_nosf > 0.0 ? param.a_nosf : fs / 2.0;
  int* winsize = calloc(nf0, sizeof(int));
  
  FP_TYPE* tmpfreq = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE* tmpampl = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE* tmpphse = calloc(param.a_nhar, sizeof(FP_TYPE));

  for(int i = 0; i < nf0; i ++) winsize[i] = param.a_nhop * 4;
  
  // estimate time-varying PSD

  FP_TYPE** noise_spectrogram = (FP_TYPE**)malloc2d(nf0, ns, sizeof(FP_TYPE));
  spectrogram_analyze(residual, nx, analysis_instant, nf0, winsize, nfft, fs, "hanning",
    noise_spectrogram, NULL, NULL, NULL);
  
  FP_TYPE* freqwarp = llsm_warp_freq(0, nosf, param.a_nnos, param.a_noswarp);
  for(int t = 0; t < nf0; t ++) {
    for(int j = 0; j < ns; j ++)
      noise_spectrogram[t][j] = exp(noise_spectrogram[t][j]);
    FP_TYPE* wrpspec = llsm_geometric_envelope(noise_spectrogram[t], nfft, fs, freqwarp, param.a_nnos);
    for(int k = 0; k < param.a_nnos; k ++)
      model -> frames[t] -> noise -> spec[k] = log(wrpspec[k]);
    model -> frames[t] -> noise -> nnos = param.a_nnos;
    free(wrpspec);
  }
  free(freqwarp);
  free2d(noise_spectrogram, nf0);

  // extract subband envelopes

  for(int i = 0; i < nf0; i ++) {
    if(rf0[i] == 0) winsize[i] = model -> conf.nhop * 2;
    else winsize[i] = round(fs / rf0[i] * 1.8) * 2;
  }
  // for each noise channel
  for(int b = 0; b < param.a_nnosband; b ++) {
    int filtord = LLSM_CHEBY_ORDER * 2;
    FP_TYPE* b_filtered = NULL;
    FP_TYPE favg = 0;
    if(b == 0) {
      b_filtered = chebyfilt(residual, nx, param.a_nosbandf[b] / fs * 2.0, 0, "lowpass");
      favg = param.a_nosbandf[b] / 2.0;
      filtord -= LLSM_CHEBY_ORDER;
    } else
    if(b == param.a_nnosband - 1) {
      b_filtered = chebyfilt((param.a_nosbandf[b - 1] > 6000.0 ? x : residual), nx,
        param.a_nosbandf[b - 1] / fs * 2.0, 0, "highpass");
      favg = (param.a_nosbandf[b - 1] + fs) / 2.0;
      filtord -= LLSM_CHEBY_ORDER;
    } else {
      b_filtered = chebyfilt((param.a_nosbandf[b - 1] > 6000.0 ? x : residual), nx,
        param.a_nosbandf[b - 1] / fs * 2.0, param.a_nosbandf[b] / fs * 2.0, "bandpass");
      favg = (param.a_nosbandf[b - 1] + param.a_nosbandf[b]) / 2.0;
    }

    FP_TYPE* b_env = moving_rms(b_filtered, nx, fs / favg * 5);
    free(b_filtered);
    
#   ifdef DEBUG_EXPORT_WAV
    char* strdup(const char*);
    char* writepath = strdup("./env .wav");
    writepath[5] = '0' + b;
    gain_write(b_env, nx, fs, writepath, 10);
    free(writepath);
#   endif

    FP_TYPE** b_spectrogram = (FP_TYPE**)malloc2d(nf0, nfft / 2 + 1, sizeof(FP_TYPE));
    FP_TYPE** b_phasegram   = (FP_TYPE**)malloc2d(nf0, nfft / 2 + 1, sizeof(FP_TYPE));
    FP_TYPE* emin = calloc(nf0, sizeof(FP_TYPE));
    spectrogram_analyze(b_env, nx, analysis_instant, nf0, winsize, nfft, fs, "blackman",
      b_spectrogram, b_phasegram, NULL, emin);
    for(int i = 0; i < nf0; i ++)
      model -> frames[i] -> noise -> emin[b] = emin[i];
    free(b_env);
    free(emin);
    
    for(int i = 0; i < nf0; i ++) {
      if(rf0[i] <= 0.0) continue;
      llsm_sinframe* ieenv = model -> frames[i] -> noise -> eenv[b];
      int inhare = harmonic_peak_picking(param, b_spectrogram[i], b_phasegram[i],
        nfft, fs, param.a_nhare, rf0[i], tmpfreq, tmpampl, tmpphse);
      for(int k = 0; k < inhare; k ++) {
        ieenv -> ampl[k] = tmpampl[k];
        ieenv -> phse[k] = tmpphse[k];
      }
    }
    free2d(b_spectrogram, nf0);
    free2d(b_phasegram, nf0);
  }
  
  free(winsize);
  free(tmpfreq);
  free(tmpampl);
  free(tmpphse);
}

llsm_layer0* llsm_layer0_analyze(llsm_parameters param, FP_TYPE* x, int nx, int fs,
    FP_TYPE* f0, int nf0, FP_TYPE** xap) {
  llsm_layer0* model = malloc(sizeof(llsm_layer0));
  int nfft = pow(2, ceil(log2(fs * param.a_tfft)));
  int ns = nfft / 2 + 1;
  FP_TYPE nosf = param.a_nosf > 0.0 ? param.a_nosf : fs / 2.0;

  model -> conf.nfrm = nf0;
  model -> conf.nhop = param.a_nhop;
  model -> conf.thop = (FP_TYPE)param.a_nhop / fs;
  model -> conf.nhar = param.a_nhar;
  model -> conf.nhare = param.a_nhare;
  model -> conf.nnos = param.a_nnos;
  model -> conf.nosf = nosf;
  model -> conf.mvf = param.a_mvf;
  model -> conf.noswarp = param.a_noswarp;
  model -> conf.nosbandf = calloc(param.a_nnosband - 1, sizeof(FP_TYPE));
  model -> conf.nnosband = param.a_nnosband;
  for(int i = 0; i < param.a_nnosband - 1; i ++)
    model -> conf.nosbandf[i] = param.a_nosbandf[i];

  // compute spectrogram

  FP_TYPE** spectrogram = (FP_TYPE**)malloc2d(nf0, ns, sizeof(FP_TYPE));
  FP_TYPE** phasegram   = (FP_TYPE**)malloc2d(nf0, ns, sizeof(FP_TYPE));
  FP_TYPE** phasegram_d = (FP_TYPE**)malloc2d(nf0, ns, sizeof(FP_TYPE));

  int* winsize = calloc(nf0, sizeof(int));
  int* analysis_instant = calloc(nf0, sizeof(int));
  for(int i = 0; i < nf0; i ++) {
    if(f0[i] == 0) winsize[i] = model -> conf.nhop * 2;
    else winsize[i] = round(fs / f0[i] * param.a_wsize) * 2;
    analysis_instant[i] = model -> conf.nhop * i;
  }
  spectrogram_analyze(x, nx, analysis_instant, nf0, winsize, nfft, fs, "blackman", spectrogram, phasegram, phasegram_d, NULL);

  FP_TYPE* rf0 = refine_f0(param, nfft, fs, f0, nf0, spectrogram, phasegram, phasegram_d);

  // harmonic analysis

  model -> frames = calloc(nf0, sizeof(llsm_frame*));
  FP_TYPE* tmpfreq = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE* tmpampl = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE* tmpphse = calloc(param.a_nhar, sizeof(FP_TYPE));
  FP_TYPE** matfreq = calloc(nf0, sizeof(FP_TYPE*));
  
  for(int i = 0; i < nf0; i ++) {
    if(rf0[i] <= 0.0) {
      model -> frames[i] = llsm_create_frame(0, param.a_nhare, param.a_nnos, param.a_nnosband);
      model -> frames[i] -> f0 = 0.0;
      continue;
    }
    
    int inhar = harmonic_peak_picking(param, spectrogram[i], phasegram[i],
      nfft, fs, param.a_nhar, rf0[i], tmpfreq, tmpampl, tmpphse);
    if(inhar == 0) {
      model -> frames[i] = llsm_create_frame(0, param.a_nhare, param.a_nnos, param.a_nnosband);
      model -> frames[i] -> f0 = 0.0;
      continue;
    }

    matfreq[i] = malloc(inhar * sizeof(FP_TYPE));
    model -> frames[i] = llsm_create_frame(inhar, param.a_nhare, param.a_nnos, param.a_nnosband);
    model -> frames[i] -> f0 = rf0[i];

    for(int k = 0; k < inhar; k ++) {
      model -> frames[i] -> sinu -> ampl[k] = tmpampl[k];
      model -> frames[i] -> sinu -> phse[k] = tmpphse[k];
      matfreq[i][k] = tmpfreq[k];
    }
  }

  free2d(spectrogram, nf0);
  free2d(phasegram, nf0);
  free2d(phasegram_d, nf0);
  free(tmpfreq);
  free(tmpampl);
  free(tmpphse);

  // resynthesis of harmonic component

  FP_TYPE* resynth = calloc(nx + param.a_nhop, sizeof(FP_TYPE));
  int nrwin = 2 * param.a_nhop;
  FP_TYPE* rwin = hanning(nrwin);
  for(int i = 0; i < nf0; i ++) {
    int tn = i * param.a_nhop;
    FP_TYPE* rfrm = gensins(
      matfreq[i], model -> frames[i] -> sinu -> ampl, model -> frames[i] -> sinu -> phse,
      model -> frames[i] -> sinu -> nhar, fs, nrwin);
    for(int j = 0; j < nrwin; j ++)
      if(tn + j - nrwin / 2 > 0)
        resynth[tn + j - nrwin / 2] += rfrm[j] * rwin[j];
    free(rfrm);
  }
  free(rwin);
  free2d(matfreq, nf0);
# ifdef DEBUG_EXPORT_WAV
  wavwrite(resynth, nx, fs, 16, "sin.wav");
# endif
  
  // extract noise/aperiodic component by subtraction

  for(int i = 0; i < nx; i ++)
    resynth[i] = x[i] - resynth[i];
  if(xap != NULL) {
    *xap = calloc(nx, sizeof(FP_TYPE));
    for(int i = 0; i < nx; i ++)
      xap[0][i] = resynth[i];
  }

# ifdef DEBUG_EXPORT_WAV
  wavwrite(resynth, nx, fs, 16, "res.wav");
# endif

  // analyze noise component
  llsm_layer0_analyze_noise(param, model, x, resynth, nx, fs, analysis_instant, rf0, nf0);

  free(resynth);
  free(rf0);
  free(winsize);
  free(analysis_instant);

  return model;
}

void llsm_layer0_phaseshift(llsm_layer0* dst, FP_TYPE* phaseshift) {
  int nfrm = dst -> conf.nfrm;
  for(int i = 0; i < nfrm; i ++) {
    llsm_sinframe* isinu = dst -> frames[i] -> sinu;
    if(isinu -> nhar <= 0) continue;
    
    for(int j = 0; j < isinu -> nhar; j ++)
      isinu -> phse[j] = wrap(isinu -> phse[j] + phaseshift[i] * (1.0 + j));
    for(int b = 0; b < dst -> conf.nnosband; b ++) {
      llsm_sinframe* ieenv = dst -> frames[i] -> noise -> eenv[b];
      for(int j = 0; j < ieenv -> nhar; j ++)
        ieenv -> phse[j] = wrap(ieenv -> phse[j] + phaseshift[i] * (1.0 + j));
    }
  }
}

// interpolate two adjacent layer-0 frames over one hop
static FP_TYPE* llsm_gensin_overlap_or_synth(
    FP_TYPE* amp0, FP_TYPE* ph0, int n0,
    FP_TYPE* amp1, FP_TYPE* ph1, int n1,
    FP_TYPE f0_0, FP_TYPE f0_1, int fs, int count)
{
  FP_TYPE* y = calloc(count, sizeof(FP_TYPE));
  int n = n0 > n1 ? n0 : n1;
  int i, t;
  FP_TYPE count_over_fs = (FP_TYPE)count / (FP_TYPE)fs;

  if(f0_0 <= 0.0) f0_0 = f0_1;
  if(f0_1 <= 0.0) f0_1 = f0_0;
  for(i = 0; i < n; i ++) {
    FP_TYPE k = (FP_TYPE)i + 1.0;
    FP_TYPE fmid = (k * f0_0 + (k * f0_1 - k * f0_0) * 0.5)
      * 2.0 * M_PI * count_over_fs;
    FP_TYPE p0, p1, p1u;
    FP_TYPE a0, a1;
    if(i >= n0) {
      a0 = 0.0;
      p0 = 0.0;
    } else {
      a0 = amp0[i];
      p0 = ph0[i];
    }
    if(i >= n1) {
      a1 = 0.0;
      p1 = 0.0;
    } else {
      a1 = amp1[i];
      p1 = ph1[i];
    }
    if(i >= n0) p0 = p1 - fmid;
    if(i >= n1) p1 = p0 + fmid;
    p1u = wrap(p1 - (fmid + p0)) + (fmid + p0);
    for(t = 0; t < count; t ++) {
      FP_TYPE u = (FP_TYPE)t / (FP_TYPE)count;
      FP_TYPE a = a0 + (a1 - a0) * u;
      FP_TYPE ph = p0 + u * (p1u - p0);
      y[t] += a * cos_2(ph);
    }
  }
  return y;
}

llsm_output* llsm_layer0_synthesize(llsm_parameters param, llsm_layer0* model) {
  int nfrm = model -> conf.nfrm;
  int nhop = model -> conf.nhop;
  int fs = param.s_fs;
  int nfft = nhop * 4;
  llsm_output* ret = malloc(sizeof(llsm_output));
  ret -> fs = fs;
  ret -> instants = NULL;
  ret -> ninst = 0;
  int* ny = & ret -> ny;
  int frm0 = max(0, (param.s_n0 - nfft) / nhop);
  int frm1 = param.s_n1 == 0 ? nfrm : min(nfrm, (param.s_n1 + nfft) / nhop);

  *ny = nfrm * nhop + nfft;
  FP_TYPE* y_sin = calloc(*ny, sizeof(FP_TYPE));
  FP_TYPE* f0 = calloc(nfrm, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++)
    f0[i] = model -> frames[i] -> f0;

  // shifted harmonic phases; kept for harmonic synth and noise-envelope shift from ph0 (moresampler only)
  FP_TYPE* ph0 = calloc(nfrm, sizeof(FP_TYPE));
  FP_TYPE** phse_har = calloc(nfrm, sizeof(FP_TYPE*));
  {
    FP_TYPE cumulative = 0.0f;
    FP_TYPE* ph_off;
    FP_TYPE* ph_inst;
    int i, j;
    for(i = 0; i < nfrm; i ++) {
      llsm_sinframe* sinu = model -> frames[i] -> sinu;
      FP_TYPE f_curr;
      FP_TYPE f_next;
      int nhar;
      if(sinu -> nhar <= 0)
        continue;
      nhar = sinu -> nhar;
      if(param.s_use_sinu_phase)
        ph0[i] = sinu -> phse[0];
      else
        ph0[i] = wrap(cumulative);
      f_curr = f0[i];
      f_next = f_curr;
      if(i != nfrm - 1 && f0[i + 1] != 0.0f)
        f_next = f0[i + 1];
      if(f_curr == 0.0f && i != nfrm - 1)
        f_curr = f0[i + 1];
      cumulative += (FP_TYPE)(2.0 * M_PI) * (0.5f * (f_curr + f_next))
                  * (FP_TYPE)nhop / (FP_TYPE)fs;
      phse_har[i] = calloc(nhar, sizeof(FP_TYPE));
      for(j = 0; j < nhar; j ++)
        phse_har[i][j] = (FP_TYPE)(j + 1) * (ph0[i] - sinu -> phse[0])
                       + sinu -> phse[j];
    }
    ph_off = llsm_smooth_phase_offset(phse_har, model);
    ph_inst = calloc(nfrm, sizeof(FP_TYPE));
    for(i = 0; i < nfrm; i ++) {
      if(phse_har[i] != NULL)
        ph_inst[i] = phse_har[i][0];
      else
        ph_inst[i] = 0.0f;
    }
    ret -> instants = llsm_phase_to_instants(
        ph_inst, ph_off, model, (FP_TYPE)fs, & ret -> ninst);
    free(ph_off);
    free(ph_inst);
  }

  // synthesize harmonic component
  if(! param.s_noiseonly) {
    for(int i = frm0; i < min(nfrm - 1, frm1); i ++) {
      int tn = i * nhop;
      llsm_sinframe* sinu0 = model -> frames[i] -> sinu;
      llsm_sinframe* sinu1 = model -> frames[i + 1] -> sinu;
      FP_TYPE* sfrm = llsm_gensin_overlap_or_synth(
        sinu0 -> ampl, phse_har[i], sinu0 -> nhar,
        sinu1 -> ampl, phse_har[i + 1], sinu1 -> nhar,
        f0[i], f0[i + 1], fs, nhop);
      for(int j = 0; j < nhop; j ++)
        y_sin[tn + j] = sfrm[j];
      free(sfrm);
    }
  }
  
  FP_TYPE* s = white_noise(1.0, fs); // one-second-long noise template
  FP_TYPE* noise_excitation = calloc(*ny, sizeof(FP_TYPE));
  FP_TYPE* ola_window = hanning(nhop * 2);

  // for each noise channel
  for(int b = 0; b < model -> conf.nnosband; b ++) {
    FP_TYPE* b_env = calloc(*ny, sizeof(FP_TYPE));
    FP_TYPE* b_env_mix = calloc(*ny, sizeof(FP_TYPE));
    FP_TYPE** phse_eenv = calloc(nfrm, sizeof(FP_TYPE*));

    for(int i = 0; i < nfrm; i ++) {
      llsm_sinframe* sinu = model -> frames[i] -> sinu;
      llsm_sinframe* ieenv = model -> frames[i] -> noise -> eenv[b];
      phse_eenv[i] = calloc(model -> conf.nhare, sizeof(FP_TYPE));
      if(sinu -> nhar > 0) {
        for(int k = 0; k < ieenv -> nhar; k ++)
          phse_eenv[i][k] = ((FP_TYPE)k + 1.0f)
              * (ph0[i] - sinu -> phse[0]) + ieenv -> phse[k];
      }
    }
    
    // generate noise envelope
    for(int i = frm0; i < frm1; i ++) {
      if(f0[i] <= 0.0) continue;
      llsm_sinframe* ieenv = model -> frames[i] -> noise -> eenv[b];
      int tn = i * nhop;
      FP_TYPE ifreq[32];
      for(int k = 0; k < ieenv -> nhar; k ++) ifreq[k] = f0[i] * (k + 1.0);
      FP_TYPE* b_env_frame = gensins(ifreq, ieenv -> ampl, phse_eenv[i], ieenv -> nhar, fs, nhop * 2);
      for(int j = 0; j < nhop * 2; j ++)
        if(tn + j - nhop > 0)
          b_env[tn + j - nhop] += b_env_frame[j] * ola_window[j];
      free(b_env_frame);
    }
    for(int i = 0; i < nfrm; i ++)
      free(phse_eenv[i]);
    free(phse_eenv);
    
    // generate band-limited noise
    int filtord = LLSM_CHEBY_ORDER * 2;
    FP_TYPE* b_template = NULL;
    FP_TYPE* b_filtered = NULL;
    if(b == 0) {
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b] / fs * 2.0, 0, "lowpass");
      filtord -= LLSM_CHEBY_ORDER;
    } else
    if(b == model -> conf.nnosband - 1) {
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b - 1] / fs * 2.0,
        model -> conf.nosf / fs * 2.0, "bandpass");
    } else {
      b_template = chebyfilt(s, fs, model -> conf.nosbandf[b - 1] / fs * 2.0,
        model -> conf.nosbandf[b] / fs * 2.0, "bandpass");
    }
    b_filtered = stretch_static_noise(b_template, fs, *ny, 256);
    free(b_template);
    subtract_minimum_envelope(b_env, *ny, f0, nhop, nfrm, fs);
    
    // normalize filtered noise and recover envelope at DC
    FP_TYPE* b_normalized = calloc(*ny, sizeof(FP_TYPE));
    for(int i = frm0; i < frm1; i ++) {
      FP_TYPE* hfrm = fetch_frame(b_filtered, *ny, i * nhop, nhop * 2);
      FP_TYPE* efrm = fetch_frame(b_env, *ny, i * nhop, nhop * 2);
      FP_TYPE havg = 0;
      for(int j = 0; j < nhop * 2; j ++)
        havg += hfrm[j] * hfrm[j];
      havg /= nhop * 2;
      for(int j = 0; j < nhop * 2; j ++)
        hfrm[j] *= sqrt(1.0 / (havg + 1e-8)); // in case havg is very small
      
      for(int j = 0; j < nhop * 2; j ++)
        efrm[j] += model -> frames[i] -> noise -> emin[b];
      for(int j = 0; j < nhop * 2; j ++) {
        efrm[j] = max(1e-7, efrm[j]);
      }
      
      for(int j = 0; j < nhop * 2; j ++)
        if(i * nhop + j - nhop > 0) {
          b_normalized[i * nhop + j - nhop] += hfrm[j] * ola_window[j];
          b_env_mix   [i * nhop + j - nhop] += efrm[j] * ola_window[j];
        }
      
      free(hfrm);
      free(efrm);
    }

    // modulate noise by envelope
    for(int i = 0; i < *ny; i ++)
      noise_excitation[i] += b_normalized[i] * b_env_mix[i];
    
    free(b_filtered);
    free(b_normalized);
    free(b_env_mix);
    free(b_env);
  }
  free(ph0);
  for(int i = 0; i < nfrm; i ++)
    free(phse_har[i]);
  free(phse_har);
  free(f0);
  free(s);

  // re-color excitation
  FP_TYPE** noisegram = calloc(nfrm, sizeof(FP_TYPE*));
  for(int i = 0; i < nfrm; i ++)
    noisegram[i] = model -> frames[i] -> noise -> spec;
  FP_TYPE* y_nos = filter_noise(param, model -> conf, noise_excitation, *ny, noisegram);
  free(noisegram);
# ifdef DEBUG_EXPORT_WAV
  wavwrite(noise_excitation, *ny, fs, 16, "synexc.wav");
  wavwrite(y_nos, *ny, fs, 16, "synres.wav");
  wavwrite(y_sin, *ny, fs, 16, "synsin.wav");
# endif
  free(noise_excitation);

  // mix optional per-frame PCM into y_nos (fade at overlap), then y = y_sin + y_nos
  {
    FP_TYPE* y_wave = calloc(*ny, sizeof(FP_TYPE));
    int i;
    for(i = frm0; i < frm1; i ++) {
      llsm_wave* wave = model -> frames[i] -> wave;
      FP_TYPE ratio;
      int out_n = 0;
      FP_TYPE* resampled;
      int fade_len;
      int start;
      int j;
      if(wave == NULL)
        continue;
      ratio = (FP_TYPE)(fs / wave -> fs);
      resampled = cig_rresample(wave -> samples, wave -> nsample, ratio, & out_n);
      start = i * nhop + (int)round((double)(wave -> t0 * (FP_TYPE)fs));
      fade_len = out_n > 0x7FF ? 0x400 : out_n / 2;
      for(j = 0; j < fade_len; j ++) {
        FP_TYPE g = (FP_TYPE)j / (FP_TYPE)fade_len;
        resampled[j] *= g;
        resampled[out_n - j - 1] *= g;
      }
      for(j = 0; j < out_n; j ++) {
        int pos = start + j;
        if(pos <= 0 || pos >= *ny)
          continue;
        if(j < fade_len)
          y_nos[pos] *= 1.0f - (FP_TYPE)j / (FP_TYPE)fade_len;
        else if(j > out_n - fade_len)
          y_nos[pos] *= (FP_TYPE)(j - out_n + fade_len) / (FP_TYPE)fade_len;
        else
          y_nos[pos] = 0.0f;
        y_wave[pos] += resampled[j];
      }
      free(resampled);
    }
    for(i = 0; i < *ny; i ++)
      y_nos[i] += y_wave[i];
    {
      FP_TYPE* y_mix = calloc(*ny, sizeof(FP_TYPE));
      for(i = 0; i < *ny; i ++)
        y_mix[i] = y_sin[i] + y_nos[i];
      ret -> y = y_mix;
    }
    ret -> y_sin = y_sin;
    ret -> y_nos = y_nos;
    free(y_wave);
  }

  free(ola_window);
  return ret;
}
