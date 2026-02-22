#pragma once
#include <cmath>
#include <complex>
#include <numbers>

#include "warpcore.hpp"

namespace warpcore {
class SvfTPT {
public:
    void Reset() noexcept {
        s1_ = 0;
        s2_ = 0;
    }

    /**
     * @param w [0, pi]
     * @param r2 <0:不稳定 =0:无阻尼 0~1:复共轭极点 >1:分裂成两个单极点
     */
    void SetCoeffSVF(float w, float r2) noexcept {
        g_ = std::tan(w / 2);
        R2_ = r2;
        d_ = 1 / (1 + r2 * g_ + g_ * g_);
    }

    void SetCoeffQ(float w, float Q) noexcept {
        SetCoeffSVF(w, 1 / Q);
    }

    std::complex<float> TickLowpass(std::complex<float> x) noexcept {
        auto bp = d_ * (g_ * (x - s2_) + s1_);
        auto v1 = bp - s1_;
        auto v2 = g_ * bp;
        auto lp = v2 + s2_;
        s1_ = bp + v1;
        s2_ = lp + v2;
        return lp;
    }
private:
    float R2_{};
    float g_{};
    float d_{};
    std::complex<float> s1_{};
    std::complex<float> s2_{};
};

class VicSineOsc {
public:
    void Reset() noexcept {
        u_ = 1.0f;
        v_ = 0.0f;
    }

    void KeepAmp() noexcept {
        float g = 1.0f / std::sqrt(u_ * u_ + v_ * v_);
        u_ *= g;
        v_ *= g;
    }

    float Tick() noexcept {
        float w = u_ - k1_ * v_;
        v_ = v_ + k2_ * w;
        u_ = w - k1_ * v_;
        return v_;
    }

    /**
     * @param f hz
     * @param fs sample rate
     */
    void SetFreq(float f, float fs) noexcept {
        auto omega = f / fs * std::numbers::pi_v<float> * 2.0f;
        k1_ = std::tan(omega / 2.0f);
        k2_ = 2 * k1_ / (1 + k1_ * k1_);
    }

    /**
     * @param w [-pi, pi]
     */
    void SetFreq(float w) noexcept {
        k1_ = std::tan(w / 2.0f);
        k2_ = 2 * k1_ / (1 + k1_ * k1_);
    }

    std::complex<float> GetCpx() const noexcept {
        return {u_, v_};
    }
private:
    float k1_{};
    float k2_{};
    float u_{1.0f};
    float v_{};
};

class SingleWarp {
public:
    float Tick(std::complex<float> osc, float x) noexcept {
        auto x1 = x * osc;
        auto x2 = svf_.TickLowpass(x1);
        x2 = svf2_.TickLowpass(x2);
        // x2 = svf3_.TickLowpass(x2);
        // x2 = svf4_.TickLowpass(x2);
        auto x3 = x2 * osc;
        return x3.real();
    }

    void SetFreq(float fcenter) noexcept {
        svf_.SetCoeffQ(fcenter, 0.50979558f);
        svf2_.SetCoeffQ(fcenter, 0.60134489f);
        svf3_.SetCoeffQ(fcenter, 0.89997622f);
        svf4_.SetCoeffQ(fcenter, 2.5629154f);
    }

    void Reset() noexcept {
        svf_.Reset();
    }
private:
    SvfTPT svf_;
    SvfTPT svf2_;
    SvfTPT svf3_;
    SvfTPT svf4_;
};

class WarpCoreScalar final : public IWarpCore {
public:
    ~WarpCoreScalar() = default;

    void Init(float fs) override {
        fs_ = fs;
    }

    void Reset() noexcept override {
        osc_phase_ = 0;
        for (auto& svf : svfs_) {
            svf.Reset();
        }
    }

    void Process(float* left, float* right, int num_samples) noexcept override {
        for (int i = 0; i < num_samples; i++) {
            osc_phase_ += osc_phase_inc_;
            osc_phase_ -= std::floor(osc_phase_);

            float x = left[i];
            float y = 0.0f;
            for (int j = 0; j < num_warps_; ++j) {
                float f_mul = static_cast<float>(j) + 1.0f;
                float mod_phase = osc_phase_ * f_mul;
                float re = _PolyCos(mod_phase);
                float im = _PolySin(mod_phase);
                std::complex<float> osc = {re, im};

                auto tmp = x * osc;
                for (int k = 0; k < poles_; ++k) {
                    tmp = svfs_[k].TickLowpass(j, tmp);
                }
                tmp = tmp * osc;
                y += tmp.real();
            }

            left[i] = y;
        }
        std::copy_n(left, num_samples, right);
    }

    void Update(const IWarpCore::Param& p) noexcept override {
        num_warps_ = p.bands;

        float finc = p.f_high / static_cast<float>(p.bands);
        float fbase = finc / 2;
        osc_base_freq_ = fbase;
        osc_phase_inc_ = fbase / fs_;

        // butterworth lowpass
        float wbase = Freq2Omega(fbase);
        float filter_w = wbase * p.filter_scale;
        filter_w = std::min(filter_w, std::numbers::pi_v<float> / 2 - 1e-5f);

        poles_ = p.filter_order;
        int n = 2 * poles_;
        for (int k = 1; k <= poles_; ++k) {
            float phi =
                (2.0f * static_cast<float>(k) - 1.0f) * std::numbers::pi_v<float> / (2.0f * static_cast<float>(n));
            float Q = 1.0f / (2.0f * std::sin(phi));

            int index = k - 1;
            svfs_[index].Set(filter_w, Q);
        }
    }

    float Freq2Omega(float f) noexcept {
        return 2 * std::numbers::pi_v<float> * f / fs_;
    }
private:
    static float _PolyCos(float phase01) noexcept {
        float x = phase01 - std::floor(phase01);
        x = 2 * std::abs(x - 0.5f) - 0.5f;
        float const x2 = x * x;
        float u = -0.540347434104161f * x2 + 2.535656174488765f;
        u = u * x2 - 5.166512943349853f;
        u = u * x2 + 3.141592653589793f;
        return u * x;
    }
    static float _PolySin(float phase01) noexcept {
        return _PolyCos(phase01 - 0.25f);
    }

    float fs_{};
    int num_warps_{};
    int poles_{};

    // complex sine generator
    float osc_phase_{};
    float osc_phase_inc_{};
    float osc_base_freq_{};

    // shared filter
    struct Svf {
        float svf_r2_{};
        float svf_g_{};
        float svf_d_{};
        std::complex<float> svf_s1_[IWarpCore::kMaxBands];
        std::complex<float> svf_s2_[IWarpCore::kMaxBands];

        void Reset() noexcept {
            std::fill_n(svf_s1_, IWarpCore::kMaxBands, 0);
            std::fill_n(svf_s2_, IWarpCore::kMaxBands, 0);
        }

        void Set(float w, float Q) noexcept {
            svf_r2_ = 1 / Q;
            svf_g_ = std::tan(w / 2);
            svf_d_ = 1 / (1 + svf_r2_ * svf_g_ + svf_g_ * svf_g_);
        }

        std::complex<float> TickLowpass(int i, std::complex<float> x) noexcept {
            auto bp = svf_d_ * (svf_g_ * (x - svf_s2_[i]) + svf_s1_[i]);
            auto v1 = bp - svf_s1_[i];
            auto v2 = svf_g_ * bp;
            auto lp = v2 + svf_s2_[i];
            svf_s1_[i] = bp + v1;
            svf_s2_[i] = lp + v2;
            return lp;
        }
    };
    Svf svfs_[IWarpCore::kMaxPoles];
};
} // namespace warpcore
