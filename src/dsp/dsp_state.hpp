#pragma once
#include <numbers>
#include <complex>
#include "pluginshared/simd.hpp"
#include "global.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace warpcore {
template <simd::IsSimdFloat SimdT>
struct SvfLaneN {
    std::array<float, global::kMaxPoles> g{};
    std::array<float, global::kMaxPoles> d{};

    struct SvfState {
       SimdT s1_re_l;
       SimdT s1_re_r;
       SimdT s1_im_l;
       SimdT s1_im_r;
       SimdT s2_re_l;
       SimdT s2_re_r;
       SimdT s2_im_l;
       SimdT s2_im_r;
    };
    simd::Array256<SvfState, global::kMaxBands * global::kMaxPoles / simd::LaneSize<SimdT>> state{};

    void Reset() noexcept {
        state.fill(SvfState{});
    }

    void SetPole(int pole_idx, float w, float Q, float fmul) noexcept {
        float r2 = 1.0f / Q;
        float g_val = std::tan(w / 2.0f) * fmul;
        g[pole_idx] = g_val;
        d[pole_idx] = 1.0f / (1.0f + r2 * g_val + g_val * g_val);
    }
};

enum class FreqDistrbution {
    k0_n,
    k1_n,
    k0_2n,
    k1_2n,
};

struct ProcessorState {
    float fs{};
    int num_warps{};
    int poles{};
    float base_mix{};
    bool pitch_affect{};
    FreqDistrbution freq_distribution{};

    // smooth
    int smooth_samples{};
    int total_smooth_samples{};
    std::array<float, global::kMaxPoles> w_{};
    std::array<float, global::kMaxPoles> q_{};
    std::array<float, global::kMaxPoles> w_inc_{};
    std::array<float, global::kMaxPoles> q_inc_{};
    std::array<float, global::kMaxPoles> last_w_{};
    std::array<float, global::kMaxPoles> last_q_{};
    float pre_osc_phase_inc_inc_{};
    float post_osc_phase_inc_inc_{};
    float last_pre_osc_phase_inc_{};
    float last_post_osc_phase_inc_{};

    // complex sine generator
    float pre_osc_phase{};
    float pre_osc_phase_inc{};
    float post_osc_phase{};
    float post_osc_phase_inc{};

    // complex lowpass filter
    float analog_fmul{};
    SvfLaneN<simd::Float128> svf128;
    SvfLaneN<simd::Float256> svf256;

    std::complex<float> band0_s1[global::kMaxPoles * 2]{};
    std::complex<float> band0_s2[global::kMaxPoles * 2]{};

    
    void SetFreq(float w_lowpass, int num_filters) noexcept {
        constexpr float atten = 0.5f;
        constexpr float square_epsi = (1.0f - atten * atten) / atten;
        analog_fmul = 1.0f / std::pow(square_epsi, 0.25f / static_cast<float>(num_filters));

        int n = 2 * num_filters;
        for (int k = 1; k <= num_filters; ++k) {
            float phi =
                (2.0f * static_cast<float>(k) - 1.0f) * std::numbers::pi_v<float> / (2.0f * static_cast<float>(n));
            float Q = 1.0f / (2.0f * std::sin(phi));
            w_[k - 1] = w_lowpass;
            q_[k - 1] = Q;
        }
    }

    void StopSmooth() noexcept {
        smooth_samples = 0;
        last_q_ = q_;
        last_w_ = w_;
        q_inc_.fill(0.0f);
        w_inc_.fill(0.0f);

        pre_osc_phase_inc_inc_ = 0.0f;
        post_osc_phase_inc_inc_ = 0.0f;
        last_pre_osc_phase_inc_ = pre_osc_phase_inc;
        last_post_osc_phase_inc_ = post_osc_phase_inc;
    }

    void BeginSmooth() noexcept {
        smooth_samples = total_smooth_samples;
        float inv_samples = 1.0f / static_cast<float>(smooth_samples);
        pre_osc_phase_inc_inc_ = (pre_osc_phase_inc - last_pre_osc_phase_inc_) * inv_samples;
        post_osc_phase_inc_inc_ = (post_osc_phase_inc - last_post_osc_phase_inc_) * inv_samples;
        for (int i = 0; i < poles; ++i) {
            q_inc_[i] = (q_[i] - last_q_[i]) * inv_samples;
            w_inc_[i] = (w_[i] - last_w_[i]) * inv_samples;
        }
    }
};

struct Param {
    int bands;
    float f_low;
    float f_high;
    float filter_scale;
    int filter_order;
    float pitch_shift;
    float base_mix;
    bool pitch_affect;
    FreqDistrbution freq_distribution;
};

struct ProcessorDsp {
    void(*init)(ProcessorState& state, float fs) noexcept;
    void(*reset)(ProcessorState& state) noexcept;
    void(*update)(ProcessorState& state, const Param& p) noexcept;
    void(*process)(ProcessorState& state, float* left, float* right, int num_samples) noexcept;

    const char* name;

    bool IsValid() const noexcept {
        return init != nullptr;
    }
};

ProcessorDsp GetProcessorDsp() noexcept;
}

#pragma GCC diagnostic pop
