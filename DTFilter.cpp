/*
 * DTFilter — Recursive Domain Transform Filter for Nuke (NDK)
 *
 * Based on Gastal & Oliveira, SIGGRAPH 2011.
*/

// #if DD_IMAGE_VERSION_MAJOR >= 15
  // Since VS 2022 this macro is required to not crash Nuke
  // https://stackoverflow.com/questions/78598141/first-stdmutexlock-crashes-in-application-built-with-latest-visual-studio
#define _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR
// #endif

#include <mutex>
#include <algorithm>
#include <cmath>
#include <vector>

#include "DDImage/Iop.h"
#include "DDImage/Row.h"
#include "DDImage/Knobs.h"
#include "DDImage/Channel.h"
#include "DDImage/Tile.h"
#include "DDImage/Format.h"
#include "DDImage/Hash.h"



using namespace DD::Image;

namespace {

static const char* const pulldownKnobEntries[] = {"RF", "NC", 0};
// Small helper to clamp index access
static inline int clampi(int v, int lo, int hi) { return std::max(lo, std::min(v, hi)); }

// Layout helper to compute linear index into [H][W][C] contiguous storage
// TODO Change this to pre compute as much as possible
static inline size_t idx(int x, int y, int c, int W, int H, int C) {
  (void)H;
  return (static_cast<size_t>(y) * W * C) + (static_cast<size_t>(x) * C) + c;
}

// Compute dct derivatives (domain transform) along X and Y
static void compute_dct(const float* joint,
                        int W, int H, int C,
                        float sigma_s, float sigma_r,
                        std::vector<float>& dctx,
                        std::vector<float>& dcty)
{
  const float ratio = sigma_s / std::max(1e-8f, sigma_r);
  const int width = W - 1;
  dctx.assign(static_cast<size_t>(H) * std::max(0, width), 1.f);
  dcty.assign(std::max(0, H - 1) * static_cast<size_t>(W), 1.f);

  // horizontal diffs
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < width; ++x) {
      float accum = 0.f;
      for (int c = 0; c < C; ++c) {
        const float a = joint[idx(x + 1, y, c, W, H, C)];
        const float b = joint[idx(x, y, c, W, H, C)];
        accum += std::fabs(a - b);
      }
      dctx[y*width + x] = 1.f + ratio * accum;
    }
  }
  // vertical diffs
  for (int y = 0; y < H - 1; ++y) {
    for (int x = 0; x < W; ++x) {
      float accum = 0.f;
      for (int c = 0; c < C; ++c) {
        const float a = joint[idx(x, y + 1, c, W, H, C)];
        const float b = joint[idx(x, y, c, W, H, C)];
        accum += std::fabs(a - b);
      }
      dcty[y * W + x] = 1.f + ratio * accum;
    }
  }
}

static inline double sigmaH_for_iter(double sigma_s, int K, int i)
{
  // i in [0..K-1]
  const double denom = std::sqrt(std::pow(4.0, K) - 1.0);
  return sigma_s * std::sqrt(3.0) * std::pow(2.0, K - i - 1) / denom;
}

static void recursive_horizontal(float* out, const std::vector<float>& dctx,
                                 int W, int H, int C, double sigma_H)
{
  if (W <= 0 || H <= 0 || C <= 0) return;
  const double a = std::exp(-std::sqrt(2.0) / sigma_H);
  const int width = W - 1;

  // Precompute V = a^(dct)
  std::vector<double> V(static_cast<size_t>(H) * std::max(0, width));
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < width; ++x) {
      V[y * (width) + x] = std::pow(a, static_cast<double>(dctx[y*width + x]));
    }
  }

  // Forward & backward passes per row
  for (int y = 0; y < H; ++y) {
    // forward
    for (int x = 1; x < W; ++x) {
      const double p = V[y*width + (x-1)];
      for (int c = 0; c < C; ++c) {
        const size_t i0 = idx(x, y, c, W, H, C);
        const size_t i1 = idx(x - 1, y, c, W, H, C);
        const double val = static_cast<double>(out[i0]);
        const double prev = static_cast<double>(out[i1]);
        out[i0] = static_cast<float>(val + p * (prev - val));
      }
    }
    // backward
    for (int x = W - 2; x >= 0; --x) {
      const double p = V[y*width + x];
      for (int c = 0; c < C; ++c) {
        const size_t i0 = idx(x, y, c, W, H, C);
        const size_t i1 = idx(x + 1, y, c, W, H, C);
        const double val = static_cast<double>(out[i0]);
        const double next = static_cast<double>(out[i1]);
        out[i0] = static_cast<float>(val + p * (next - val));
      }
    }
  }
}

static void recursive_vertical(float* out, const std::vector<float>& dcty,
                               int W, int H, int C, double sigma_H)
{
  if (W <= 0 || H <= 0 || C <= 0) return;
  const double a = std::exp(-std::sqrt(2.0) / sigma_H);

  // Precompute V = a^(dct)
  std::vector<double> V(std::max(0, H - 1) * static_cast<size_t>(W));
  for (int y = 0; y < H - 1; ++y) {
    for (int x = 0; x < W; ++x) {
      V[y * W + x] = std::pow(a, static_cast<double>(dcty[y * W + x]));
    }
  }

  // Forward & backward passes per column
  for (int x = 0; x < W; ++x) {
    // forward
    for (int y = 1; y < H; ++y) {
      const double p = V[(y - 1) * W + x];
      for (int c = 0; c < C; ++c) {
        const size_t i0 = idx(x, y, c, W, H, C);
        const size_t i1 = idx(x, y - 1, c, W, H, C);
        const double val = static_cast<double>(out[i0]);
        const double prev = static_cast<double>(out[i1]);
        out[i0] = static_cast<float>(val + p * (prev - val));
      }
    }
    // backward
    for (int y = H - 2; y >= 0; --y) {
      const double p = V[y * W + x];
      for (int c = 0; c < C; ++c) {
        const size_t i0 = idx(x, y, c, W, H, C);
        const size_t i1 = idx(x, y + 1, c, W, H, C);
        const double val = static_cast<double>(out[i0]);
        const double next = static_cast<double>(out[i1]);
        out[i0] = static_cast<float>(val + p * (next - val));
      }
    }
  }
}

static void domain_transform_filter(const float* src, float* dst, const float* joint,
                                    int W, int H, int C,
                                    float sigma_s, float sigma_r, int K)
{
  std::vector<float> dctx, dcty;
  compute_dct(joint, W, H, C, sigma_s, sigma_r, dctx, dcty);

  const size_t N = static_cast<size_t>(W) * H * C;
  std::copy(src, src + N, dst);

  // Apply recursive filtering K iterations
  for (int i = 0; i < K; ++i) {
    const double sigmaH = sigmaH_for_iter(sigma_s, K, i);
    recursive_horizontal(dst, dctx, W, H, C, sigmaH);
    recursive_vertical(dst, dcty, W, H, C, sigmaH);
  }
}

}

class DTFilterOp : public Iop {
public:
  DTFilterOp(Node* node)
    : Iop(node),
    sigma_s_(30.0f), sigma_r_(0.4f), iterations_(3),
    has_rgb_(false), computed_(false)
  {}

  void knobs(Knob_Callback f) override {
    Float_knob(f, &sigma_s_, "Sigma_s");
    Tooltip(f, "Spatial standard deviation (pixels) controlling the extent of smoothing.");
    Float_knob(f, &sigma_r_, "Sigma_r");
    Tooltip(f, "Range standard deviation controlling edge sensitivity (larger = more smoothing across edges).");
    Int_knob(f, &iterations_, "Iterations");
    Tooltip(f, "Number of recursive filtering iterations (default to 3).");
  }
  
  const char* Class() const override { return "DTFilter"; }
  const char* node_help() const override {
    return "Recursive Domain Transform filter (edge-preserving).\n"
    "Implements Gastal & Oliveira 2011 using full-frame processing.";
  }
  
  
  void _validate(bool for_real) override {
    copy_info(); // same bbox/format/channels as input0
    out_channels_ = input0().channels();
  }
  
  
  void _request(int x, int y, int r, int t, ChannelMask channels, int count) override {
    input0().request(x, y, r, t, channels, count);
  }
  
  void _open() override {
    const Format format = input0().format();
    bx_ = format.x(); by_ = format.y(); br_ = format.r(); bt_ = format.t();
    W_ = br_ - bx_; H_ = bt_ - by_;
    
    ChannelSet inChans = input0().info().channels();
    has_rgb_ = inChans.contains(Chan_Red) && inChans.contains(Chan_Green) && inChans.contains(Chan_Blue);
    
    if (W_ <= 0 || H_ <= 0 || !has_rgb_) {
      std::lock_guard<std::mutex> lock(mutex_);
      src_.clear(); out_.clear(); joint_.clear();
      computed_ = true;
      return;
    }
    
    const int C = 3;
    src_.assign(static_cast<size_t>(W_) * H_ * C, 0.f);
    joint_.assign(static_cast<size_t>(W_) * H_ * C, 0.f);
    out_.assign(static_cast<size_t>(W_) * H_ * C, 0.f);
    
    ChannelSet readChans;
    readChans += Chan_Red; readChans += Chan_Green; readChans += Chan_Blue;
    Interest interest(input0(), bx_, by_, br_, bt_, readChans, true);
    interest.unlock();
    
    for (int y = by_; y < bt_; ++y) {
      Row row(bx_, br_);
      row.get(input0(), y, bx_, br_, readChans);
      if (aborted()) return;
      
      const float* R = row[Chan_Red]   + bx_;
      const float* G = row[Chan_Green] + bx_;
      const float* B = row[Chan_Blue]  + bx_;
      
      for (int x = bx_; x < br_; ++x) {
        const int lx = x - bx_;
        const int ly = y - by_;
        src_[idx(lx, ly, 0, W_, H_, C)] = R[x];
        src_[idx(lx, ly, 1, W_, H_, C)] = G[x];
        src_[idx(lx, ly, 2, W_, H_, C)] = B[x];
        
        // joint == src (non‑guided).
        // TODO implement join filtering as stated in paper.
        joint_[idx(lx, ly, 0, W_, H_, C)] = R[x];
        joint_[idx(lx, ly, 1, W_, H_, C)] = G[x];
        joint_[idx(lx, ly, 2, W_, H_, C)] = B[x];
      }
    }
    
    domain_transform_filter(src_.data(), out_.data(), joint_.data(), W_, H_, C,
    sigma_s_, sigma_r_, std::max(1, iterations_));
    
    std::lock_guard<std::mutex> lock(mutex_);
    computed_ = true;
  }
  
  void engine(int y, int l, int r, ChannelMask channels, Row& row) override {
    row.erase(channels);
    
    if (!has_rgb_ || W_ <= 0 || H_ <= 0) {
      row.get(input0(), y, l, r, channels);
      return;
    }
    
    if (!computed_) {
      row.get(input0(), y, l, r, channels);
      return;
    }
    
    foreach (z, channels) {
      float* outp = row.writable(z) + l;
      if (z == Chan_Red || z == Chan_Green || z == Chan_Blue) {
        const int c = (z == Chan_Red) ? 0 : (z == Chan_Green ? 1 : 2);
        for (int x = l; x < r; ++x) {
          const int lx = clampi(x - bx_, 0, W_ - 1);
          const int ly = clampi(y - by_, 0, H_ - 1);
          *outp++ = out_[idx(lx, ly, c, W_, H_, 3)];
        }
      } else {
        Row in(bx_, br_);
        in.get(input0(), y, l, r, ChannelSet(z));
        const float* srcp = in[z] + l;
        for (int x = l; x < r; ++x) *outp++ = *srcp++;
      } 
    }
  }
  
private:
  float sigma_s_;
  float sigma_r_;
  int   iterations_;

  int bx_ = 0, by_ = 0, br_ = 0, bt_ = 0;
  int W_ = 0, H_ = 0;
  bool has_rgb_ = false;
  
  std::vector<float> src_, joint_, out_;
  
  std::mutex mutex_;
  bool computed_ = false;
  Hash last_hash_;
};

static Iop* build_DTFilter(Node* node) { return new DTFilterOp(node); }
static const Iop::Description d_DTFilter("DTFilter", "Filter/DTFilter", build_DTFilter);
