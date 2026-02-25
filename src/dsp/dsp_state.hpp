#pragma once
#include <numbers>
#include <complex>
#include "pluginshared/simd.hpp"
#include "global.hpp"

namespace warpcore {
template <simd::IsSimdFloat SimdT>
struct SvfLaneN {
    float g[global::kMaxPoles];
    float d[global::kMaxPoles];

    struct SvfState {
       SimdT s1_re;
       SimdT s1_im;
       SimdT s2_re;
       SimdT s2_im;
    };
    simd::Array256<SvfState, global::kMaxBands * global::kMaxPoles / simd::LaneSize<SimdT>> state;

    void Reset() noexcept {
        state.fill(SvfState{});
    }

    void SetPole(int pole_idx, float w, float Q) noexcept {
        float r2 = 1.0f / Q;
        float g_val = std::tan(w / 2.0f);
        g[pole_idx] = g_val;
        d[pole_idx] = 1.0f / (1.0f + r2 * g_val + g_val * g_val);
    }

    void SetFreq(float w_lowpass, int num_filters) noexcept {
        int n = 2 * num_filters;
        for (int k = 1; k <= num_filters; ++k) {
            float phi =
                (2.0f * static_cast<float>(k) - 1.0f) * std::numbers::pi_v<float> / (2.0f * static_cast<float>(n));
            float Q = 1.0f / (2.0f * std::sin(phi));
            SetPole(k - 1, w_lowpass, Q);
        }
    }
};

struct ProcessorState {
    float fs{};
    int num_warps{};
    int poles{};
    float base_mix{};

    // complex sine generator
    float osc_base_freq{};

    float pre_osc_phase{};
    float pre_osc_phase_inc{};

    float post_osc_phase{};
    float post_osc_phase_inc{};

    // complex lowpass filter
    union {
        SvfLaneN<simd::Float128> svf128;
        SvfLaneN<simd::Float256> svf256;
    };

    std::complex<float> band0_s1[global::kMaxPoles];
    std::complex<float> band0_s2[global::kMaxPoles];
};

struct Param {
    int bands;
    float f_low;
    float f_high;
    float filter_scale;
    int filter_order;
    float post_osc_freq_mul;
    float base_mix;
};

struct ProcessorDsp {
    void(*init)(ProcessorState& state, float fs) noexcept;
    void(*reset)(ProcessorState& state) noexcept;
    void(*update)(ProcessorState& state, const Param& p) noexcept;
    void(*process)(ProcessorState& state, float* left, float* right, int num_samples) noexcept;

    bool IsValid() const noexcept {
        return init != nullptr && reset != nullptr && update != nullptr && process != nullptr;
    }
};

ProcessorDsp GetProcessorDsp() noexcept;
}
