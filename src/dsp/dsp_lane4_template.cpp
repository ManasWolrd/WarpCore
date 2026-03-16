#include "dsp_state.hpp"

#include <cassert>
#include <complex>

namespace warpcore {
template <FreqDistrbution kFreqMode, int kPoles>
static void ProcessInternal(warpcore::ProcessorState& state, float* left, float* right, int num_samples) noexcept {
    constexpr simd::Array128<simd::Float128, 4> kBandGainLut{
        simd::Float128{1.0f, 1.0f, 1.0f, 1.0f},
        simd::Float128{1.0f, 0.0f, 0.0f, 0.0f},
        simd::Float128{1.0f, 1.0f, 0.0f, 0.0f},
        simd::Float128{1.0f, 1.0f, 1.0f, 0.0f},
    };
    simd::Float128 tail_gain = kBandGainLut[state.num_warps & 3];

    float band0_dry_mix = state.base_mix;
    float band0_wet_mix = 1.0f - band0_dry_mix;
    int simd_loop_count = (state.num_warps + 3) / 4;

    if (state.num_warps <= 4) {
        tail_gain[0] = band0_wet_mix;
    }

    for (int i = 0; i < num_samples; i++) {
        // -------------------- tick complex sine generators --------------------
        state.pre_osc_phase += state.pre_osc_phase_inc;
        state.pre_osc_phase -= std::floor(state.pre_osc_phase);

        state.post_osc_phase += state.post_osc_phase_inc;
        state.post_osc_phase -= std::floor(state.post_osc_phase);

        // e^jwt
        constexpr float twopi = 2.0f * std::numbers::pi_v<float>;
        std::complex<float> pre_osc_f32 = {std::cos(state.pre_osc_phase * twopi),
                                           std::sin(state.pre_osc_phase * twopi)};
        std::complex<float> post_osc_f32 = {std::cos(state.post_osc_phase * twopi),
                                            std::sin(state.post_osc_phase * twopi)};

        const auto pre_osc_f32_0 = std::complex{1.0f, 0.0f};
        const auto pre_osc_f32_1 = pre_osc_f32;
        const auto pre_osc_f32_2 = pre_osc_f32 * pre_osc_f32;
        const auto pre_osc_f32_3 = pre_osc_f32 * pre_osc_f32 * pre_osc_f32;
        const auto pre_osc_f32_4 = pre_osc_f32 * pre_osc_f32 * pre_osc_f32 * pre_osc_f32;
        const auto pre_osc_f32_5 = pre_osc_f32_1 * pre_osc_f32_4;
        const auto pre_osc_f32_6 = pre_osc_f32_2 * pre_osc_f32_4;
        const auto pre_osc_f32_7 = pre_osc_f32_3 * pre_osc_f32_4;
        const auto pre_osc_f32_8 = pre_osc_f32_4 * pre_osc_f32_4;

        const auto post_osc_f32_0 = std::complex{1.0f, 0.0f};
        const auto post_osc_f32_1 = post_osc_f32;
        const auto post_osc_f32_2 = post_osc_f32 * post_osc_f32;
        const auto post_osc_f32_3 = post_osc_f32 * post_osc_f32 * post_osc_f32;
        const auto post_osc_f32_4 = post_osc_f32 * post_osc_f32 * post_osc_f32 * post_osc_f32;
        const auto post_osc_f32_5 = post_osc_f32_1 * post_osc_f32_4;
        const auto post_osc_f32_6 = post_osc_f32_2 * post_osc_f32_4;
        const auto post_osc_f32_7 = post_osc_f32_3 * post_osc_f32_4;
        const auto post_osc_f32_8 = post_osc_f32_4 * post_osc_f32_4;

        simd::Complex128 pre_osc;
        simd::Complex128 post_osc;
        simd::Complex128 pre_osc_n_val;
        simd::Complex128 post_osc_n_val;

        if constexpr (kFreqMode == FreqDistrbution::k0_n) {
            pre_osc = simd::Complex128{
                .re = simd::BroadcastF128(pre_osc_f32_4.real()),
                .im = simd::BroadcastF128(pre_osc_f32_4.imag()),
            };
            post_osc = simd::Complex128{
                .re = simd::BroadcastF128(post_osc_f32_4.real()),
                .im = simd::BroadcastF128(post_osc_f32_4.imag()),
            };

            pre_osc_n_val = simd::Complex128{
                .re = {pre_osc_f32_0.real(), pre_osc_f32_1.real(), pre_osc_f32_2.real(), pre_osc_f32_3.real()},
                .im = {pre_osc_f32_0.imag(), pre_osc_f32_1.imag(), pre_osc_f32_2.imag(), pre_osc_f32_3.imag()},
            };
            post_osc_n_val = simd::Complex128{
                .re = {post_osc_f32_0.real(), post_osc_f32_1.real(), post_osc_f32_2.real(), post_osc_f32_3.real()},
                .im = {post_osc_f32_0.imag(), post_osc_f32_1.imag(), post_osc_f32_2.imag(), post_osc_f32_3.imag()},
            };
        }
        else if constexpr (kFreqMode == FreqDistrbution::k1_n) {
            pre_osc = simd::Complex128{
                .re = simd::BroadcastF128(pre_osc_f32_8.real()),
                .im = simd::BroadcastF128(pre_osc_f32_8.imag()),
            };
            post_osc = simd::Complex128{
                .re = simd::BroadcastF128(post_osc_f32_8.real()),
                .im = simd::BroadcastF128(post_osc_f32_8.imag()),
            };

            pre_osc_n_val = simd::Complex128{
                .re = {pre_osc_f32_1.real(), pre_osc_f32_2.real(), pre_osc_f32_3.real(), pre_osc_f32_4.real()},
                .im = {pre_osc_f32_1.imag(), pre_osc_f32_2.imag(), pre_osc_f32_3.imag(), pre_osc_f32_4.imag()},
            };
            post_osc_n_val = simd::Complex128{
                .re = {post_osc_f32_1.real(), post_osc_f32_2.real(), post_osc_f32_3.real(), post_osc_f32_4.real()},
                .im = {post_osc_f32_1.imag(), post_osc_f32_2.imag(), post_osc_f32_3.imag(), post_osc_f32_4.imag()},
            };
        }
        else if constexpr (kFreqMode == FreqDistrbution::k0_2n) {
            pre_osc = simd::Complex128{
                .re = simd::BroadcastF128(pre_osc_f32_8.real()),
                .im = simd::BroadcastF128(pre_osc_f32_8.imag()),
            };
            post_osc = simd::Complex128{
                .re = simd::BroadcastF128(post_osc_f32_8.real()),
                .im = simd::BroadcastF128(post_osc_f32_8.imag()),
            };

            pre_osc_n_val = simd::Complex128{
                .re = {pre_osc_f32_0.real(), pre_osc_f32_2.real(), pre_osc_f32_4.real(), pre_osc_f32_6.real()},
                .im = {pre_osc_f32_0.imag(), pre_osc_f32_2.imag(), pre_osc_f32_4.imag(), pre_osc_f32_6.imag()},
            };
            post_osc_n_val = simd::Complex128{
                .re = {post_osc_f32_0.real(), post_osc_f32_2.real(), post_osc_f32_4.real(), post_osc_f32_6.real()},
                .im = {post_osc_f32_0.imag(), post_osc_f32_2.imag(), post_osc_f32_4.imag(), post_osc_f32_6.imag()},
            };
        }
        else if constexpr (kFreqMode == FreqDistrbution::k1_2n) {
            pre_osc = simd::Complex128{
                .re = simd::BroadcastF128(pre_osc_f32_8.real()),
                .im = simd::BroadcastF128(pre_osc_f32_8.imag()),
            };
            post_osc = simd::Complex128{
                .re = simd::BroadcastF128(post_osc_f32_8.real()),
                .im = simd::BroadcastF128(post_osc_f32_8.imag()),
            };

            pre_osc_n_val = simd::Complex128{
                .re = {pre_osc_f32_1.real(), pre_osc_f32_3.real(), pre_osc_f32_5.real(), pre_osc_f32_7.real()},
                .im = {pre_osc_f32_1.imag(), pre_osc_f32_3.imag(), pre_osc_f32_5.imag(), pre_osc_f32_7.imag()},
            };
            post_osc_n_val = simd::Complex128{
                .re = {post_osc_f32_1.real(), post_osc_f32_3.real(), post_osc_f32_5.real(), post_osc_f32_7.real()},
                .im = {post_osc_f32_1.imag(), post_osc_f32_3.imag(), post_osc_f32_5.imag(), post_osc_f32_7.imag()},
            };
        }

        // -------------------- process first band --------------------
        float x_left = left[i];
        float x_right = right[i];

        std::complex<float> cpx_x_left = x_left * pre_osc_f32;
        std::complex<float> cpx_x_right = x_right * pre_osc_f32;
#pragma unroll
        for (int k = 0; k < kPoles; ++k) {
            const float gk = state.svf128.g[k];
            const float dk = state.svf128.d[k];
            auto& s1_l = state.band0_s1[2 * k];
            auto& s1_r = state.band0_s1[2 * k + 1];
            auto& s2_l = state.band0_s2[2 * k];
            auto& s2_r = state.band0_s2[2 * k + 1];

            auto bp_l = dk * (gk * (cpx_x_left - s2_l) + s1_l);
            auto bp_r = dk * (gk * (cpx_x_right - s2_r) + s1_r);
            auto v1_l = bp_l - s1_l;
            auto v1_r = bp_r - s1_r;
            auto v2_l = gk * bp_l;
            auto v2_r = gk * bp_r;
            auto lp_l = v2_l + s2_l;
            auto lp_r = v2_r + s2_r;

            s1_l = bp_l + v1_l;
            s1_r = bp_r + v1_r;
            s2_l = lp_l + v2_l;
            s2_r = lp_r + v2_r;
            cpx_x_left = lp_l;
            cpx_x_right = lp_r;
        }
        float y_l = (cpx_x_left * std::conj(post_osc_f32)).real() * band0_dry_mix;
        float y_r = (cpx_x_right * std::conj(post_osc_f32)).real() * band0_dry_mix;

        // -------------------- process bands --------------------
        auto* svf_state = state.svf128.state.data();
        simd::Float128 band_gain{band0_wet_mix, 1.0f, 1.0f, 1.0f};

        for (int j = 0; j < simd_loop_count - 1; ++j) {
            // std::complex<float> tmp = x * pre_osc_n_val;
            // pre_osc_n_val *= pre_osc;
            auto tmp_l = simd::BroadcastF128(x_left) * pre_osc_n_val;
            auto tmp_r = simd::BroadcastF128(x_right) * pre_osc_n_val;
            pre_osc_n_val *= pre_osc;

#pragma unroll
            for (int k = 0; k < kPoles; ++k) {
                const float gk = state.svf128.g[k];
                const float dk = state.svf128.d[k];

                auto s1_re_l = svf_state->s1_re_l;
                auto s1_im_l = svf_state->s1_im_l;
                auto s2_re_l = svf_state->s2_re_l;
                auto s2_im_l = svf_state->s2_im_l;

                auto s1_re_r = svf_state->s1_re_r;
                auto s1_im_r = svf_state->s1_im_r;
                auto s2_re_r = svf_state->s2_re_r;
                auto s2_im_r = svf_state->s2_im_r;

                auto bp_re_l = dk * (gk * (tmp_l.re - s2_re_l) + s1_re_l);
                auto bp_im_l = dk * (gk * (tmp_l.im - s2_im_l) + s1_im_l);
                auto v1_re_l = bp_re_l - s1_re_l;
                auto v1_im_l = bp_im_l - s1_im_l;
                auto v2_re_l = gk * bp_re_l;
                auto v2_im_l = gk * bp_im_l;
                auto lp_re_l = v2_re_l + s2_re_l;
                auto lp_im_l = v2_im_l + s2_im_l;

                auto bp_re_r = dk * (gk * (tmp_r.re - s2_re_r) + s1_re_r);
                auto bp_im_r = dk * (gk * (tmp_r.im - s2_im_r) + s1_im_r);
                auto v1_re_r = bp_re_r - s1_re_r;
                auto v1_im_r = bp_im_r - s1_im_r;
                auto v2_re_r = gk * bp_re_r;
                auto v2_im_r = gk * bp_im_r;
                auto lp_re_r = v2_re_r + s2_re_r;
                auto lp_im_r = v2_im_r + s2_im_r;

                svf_state->s1_re_l = bp_re_l + v1_re_l;
                svf_state->s1_im_l = bp_im_l + v1_im_l;
                svf_state->s2_re_l = lp_re_l + v2_re_l;
                svf_state->s2_im_l = lp_im_l + v2_im_l;
                svf_state->s1_re_r = bp_re_r + v1_re_r;
                svf_state->s1_im_r = bp_im_r + v1_im_r;
                svf_state->s2_re_r = lp_re_r + v2_re_r;
                svf_state->s2_im_r = lp_im_r + v2_im_r;

                ++svf_state;

                tmp_l.re = lp_re_l;
                tmp_l.im = lp_im_l;
                tmp_r.re = lp_re_r;
                tmp_r.im = lp_im_r;
            }

            // y += (tmp * post_osc_n_val).real();
            // post_osc_n_val *= post_osc;
            auto band_out_l = tmp_l * post_osc_n_val;
            auto band_out_r = tmp_r * post_osc_n_val;
            y_l += simd::ReduceAdd(band_out_l.re * band_gain);
            y_r += simd::ReduceAdd(band_out_r.re * band_gain);
            post_osc_n_val *= post_osc;

            band_gain = simd::BroadcastF128(1.0f);
        }

        // -------------------- here we have: 1/2/3/4 --------------------
        // std::complex<float> tmp = x * pre_osc_n_val;
        // pre_osc_n_val *= pre_osc;
        auto tmp_l = simd::BroadcastF128(x_left) * pre_osc_n_val;
        auto tmp_r = simd::BroadcastF128(x_right) * pre_osc_n_val;

#pragma unroll
        for (int k = 0; k < kPoles; ++k) {
            const float gk = state.svf128.g[k];
            const float dk = state.svf128.d[k];

            auto s1_re_l = svf_state->s1_re_l;
            auto s1_im_l = svf_state->s1_im_l;
            auto s2_re_l = svf_state->s2_re_l;
            auto s2_im_l = svf_state->s2_im_l;

            auto s1_re_r = svf_state->s1_re_r;
            auto s1_im_r = svf_state->s1_im_r;
            auto s2_re_r = svf_state->s2_re_r;
            auto s2_im_r = svf_state->s2_im_r;

            auto bp_re_l = dk * (gk * (tmp_l.re - s2_re_l) + s1_re_l);
            auto bp_im_l = dk * (gk * (tmp_l.im - s2_im_l) + s1_im_l);
            auto v1_re_l = bp_re_l - s1_re_l;
            auto v1_im_l = bp_im_l - s1_im_l;
            auto v2_re_l = gk * bp_re_l;
            auto v2_im_l = gk * bp_im_l;
            auto lp_re_l = v2_re_l + s2_re_l;
            auto lp_im_l = v2_im_l + s2_im_l;

            auto bp_re_r = dk * (gk * (tmp_r.re - s2_re_r) + s1_re_r);
            auto bp_im_r = dk * (gk * (tmp_r.im - s2_im_r) + s1_im_r);
            auto v1_re_r = bp_re_r - s1_re_r;
            auto v1_im_r = bp_im_r - s1_im_r;
            auto v2_re_r = gk * bp_re_r;
            auto v2_im_r = gk * bp_im_r;
            auto lp_re_r = v2_re_r + s2_re_r;
            auto lp_im_r = v2_im_r + s2_im_r;

            svf_state->s1_re_l = bp_re_l + v1_re_l;
            svf_state->s1_im_l = bp_im_l + v1_im_l;
            svf_state->s2_re_l = lp_re_l + v2_re_l;
            svf_state->s2_im_l = lp_im_l + v2_im_l;
            svf_state->s1_re_r = bp_re_r + v1_re_r;
            svf_state->s1_im_r = bp_im_r + v1_im_r;
            svf_state->s2_re_r = lp_re_r + v2_re_r;
            svf_state->s2_im_r = lp_im_r + v2_im_r;

            ++svf_state;

            tmp_l.re = lp_re_l;
            tmp_l.im = lp_im_l;
            tmp_r.re = lp_re_r;
            tmp_r.im = lp_im_r;
        }

        // y += (tmp * post_osc_n_val).real();
        // post_osc_n_val *= post_osc;
        auto band_out_l = tmp_l * post_osc_n_val;
        auto band_out_r = tmp_r * post_osc_n_val;
        y_l += simd::ReduceAdd(band_out_l.re * tail_gain);
        y_r += simd::ReduceAdd(band_out_r.re * tail_gain);

        left[i] = y_l;
        right[i] = y_r;
    }
}

template <FreqDistrbution kFreqMode>
static void ProcessPoles(warpcore::ProcessorState& state, float* left, float* right, int num_samples) noexcept {
    switch (state.poles) {
        case 1:
            ProcessInternal<kFreqMode, 1>(state, left, right, num_samples);
            break;
        case 2:
            ProcessInternal<kFreqMode, 2>(state, left, right, num_samples);
            break;
        case 3:
            ProcessInternal<kFreqMode, 3>(state, left, right, num_samples);
            break;
        case 4:
            ProcessInternal<kFreqMode, 4>(state, left, right, num_samples);
            break;
        case 5:
            ProcessInternal<kFreqMode, 5>(state, left, right, num_samples);
            break;
        case 6:
            ProcessInternal<kFreqMode, 6>(state, left, right, num_samples);
            break;
        case 7:
            ProcessInternal<kFreqMode, 7>(state, left, right, num_samples);
            break;
        case 8:
            ProcessInternal<kFreqMode, 8>(state, left, right, num_samples);
            break;
        default:
            assert(false);
            break;
    }
}

// ----------------------------------------
// dsp processor
// ----------------------------------------

static void Init(warpcore::ProcessorState& state, float fs) noexcept {
    state.fs = fs;
}

static void Reset(warpcore::ProcessorState& state) noexcept {
    state.svf128.Reset();
    state.pre_osc_phase = 0.0f;
    state.post_osc_phase = 0.0f;

    std::fill_n(state.band0_s1, global::kMaxPoles * 2, std::complex<float>{});
    std::fill_n(state.band0_s2, global::kMaxPoles * 2, std::complex<float>{});
}

static void Update(warpcore::ProcessorState& state, const warpcore::Param& p) noexcept {
    state.num_warps = p.bands;
    state.base_mix = p.base_mix;
    state.pitch_affect = p.pitch_affect;
    state.freq_distribution = p.freq_distribution;

        float fhigh = p.f_high;
    float fshit = p.pitch_affect ? -p.pitch_shift : p.pitch_shift;
    fshit = std::exp2(fshit / 12.0f);

    if (fhigh > 20000.0f) {
        fhigh = state.fs / 2;
    }
    fhigh = std::min(fhigh, state.fs / 2);

    if (state.pitch_affect && p.pitch_shift < 0.0f) {
        // formant anti-alasing
        fhigh /= fshit;
    }

    float f_first_band_stop = fhigh / static_cast<float>(p.bands);
    float f_first_band_center = f_first_band_stop / 2;

    if (!state.pitch_affect && p.pitch_shift > 0.0f) {
        // pitch anti-alasing
        f_first_band_stop *= fshit;
        state.num_warps = static_cast<int>(fhigh / f_first_band_stop);
    }

    if (state.freq_distribution == FreqDistrbution::k0_n || state.freq_distribution == FreqDistrbution::k1_n) {
        f_first_band_center *= 2;
    }

    state.pre_osc_phase_inc = f_first_band_center / state.fs;
    state.post_osc_phase_inc = state.pre_osc_phase_inc * fshit;

    if (p.pitch_affect) {
        std::swap(state.pre_osc_phase_inc, state.post_osc_phase_inc);
    }

    // butterworth lowpass
    float wbase = f_first_band_center * 2 * std::numbers::pi_v<float> / state.fs;
    if (state.freq_distribution == FreqDistrbution::k0_n || state.freq_distribution == FreqDistrbution::k1_n) {
        if (p.pitch_affect) {
            wbase *= fshit;
        }
    }
    else {
        if (!p.pitch_affect) {
            wbase *= fshit;
        }
    }
    float filter_w = wbase * p.filter_scale;
    filter_w = std::min(filter_w, std::numbers::pi_v<float> - 0.1f);

    if (state.poles != p.filter_order) {
        state.svf256.Reset();
    }
    state.poles = p.filter_order;
    state.svf256.SetFreq(filter_w, p.filter_order);
}

static void Process(warpcore::ProcessorState& state, float* left, float* right, int num_samples) noexcept {
    switch (state.freq_distribution) {
        case FreqDistrbution::k0_n:
            ProcessPoles<FreqDistrbution::k0_n>(state, left, right, num_samples);
            break;
        case FreqDistrbution::k1_n:
            ProcessPoles<FreqDistrbution::k1_n>(state, left, right, num_samples);
            break;
        case FreqDistrbution::k0_2n:
            ProcessPoles<FreqDistrbution::k0_2n>(state, left, right, num_samples);
            break;
        case FreqDistrbution::k1_2n:
            ProcessPoles<FreqDistrbution::k1_2n>(state, left, right, num_samples);
            break;
    }
}

// ----------------------------------------
// export
// ----------------------------------------

#ifndef DSP_EXPORT_NAME
#error "不应该编译这个文件,在其他cpp包含这个cpp并定义DSP_EXPORT_NAME=`dsp_dispatch.cpp里的变量`"
#endif

ProcessorDsp DSP_EXPORT_NAME{Init, Reset, Update, Process};
} // namespace warpcore
