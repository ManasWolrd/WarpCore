#include "plugin_ui.hpp"

#include "../PluginProcessor.h"

PluginUi::PluginUi(EmptyAudioProcessor& p)
    : preset_(*p.preset_manager_) {
    addAndMakeVisible(preset_);

    auto& apvt = *p.value_tree_;
    warp_.BindParam(apvt, "warp");
    warp_.slider.setTooltip("Split spectrum to {this value} bands");
    addAndMakeVisible(warp_);

    f_low_.BindParam(apvt, "scale");
    f_low_.slider.setTooltip("Set the filter's bandwidth of each band.\nlower sounds metallic, higher sounds clicky(downsampled)");
    addAndMakeVisible(f_low_);

    f_high_.BindParam(apvt, "f_high");
    f_high_.slider.setTooltip("Set the warp spectrum ceil frequency.\nHigher frequency components will be silenced");
    addAndMakeVisible(f_high_);

    poles_.BindParam(apvt, "poles");
    poles_.slider.setTooltip("Set the filter's poles");
    addAndMakeVisible(poles_);

    base_mix_.BindParam(apvt, "base_mix");
    base_mix_.slider.setTooltip("Set the band0's warp behavior.\nThis value is used to hear the pitch in low warp value.\n0 = no warp, 1 = warp");
    addAndMakeVisible(base_mix_);

    pitch_.BindParam(apvt, "pitch");
    pitch_.slider.setTooltip("Set the warped signal's pitch.\nThis is used to simulate PiWarp's pitch parameter");
    addAndMakeVisible(pitch_);

    setSize(480, 130);
}

void PluginUi::resized() {
    auto b = getLocalBounds();
    preset_.setBounds(b.removeFromTop(30));

    auto line = b.removeFromTop(100);
    warp_.setBounds(line.removeFromLeft(80));
    f_low_.setBounds(line.removeFromLeft(80));
    f_high_.setBounds(line.removeFromLeft(80));
    poles_.setBounds(line.removeFromLeft(80));
    pitch_.setBounds(line.removeFromLeft(80));
    base_mix_.setBounds(line.removeFromLeft(80));
}

void PluginUi::paint(juce::Graphics& g) {
    
}
