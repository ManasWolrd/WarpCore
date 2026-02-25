#include "plugin_ui.hpp"

#include "../PluginProcessor.h"

PluginUi::PluginUi(EmptyAudioProcessor& p)
    : preset_(*p.preset_manager_) {
    addAndMakeVisible(preset_);

    auto& apvt = *p.value_tree_;
    warp_.BindParam(apvt, "warp");
    addAndMakeVisible(warp_);

    f_low_.BindParam(apvt, "scale");
    addAndMakeVisible(f_low_);

    f_high_.BindParam(apvt, "f_high");
    addAndMakeVisible(f_high_);

    poles_.BindParam(apvt, "poles");
    addAndMakeVisible(poles_);

    base_mix_.BindParam(apvt, "base_mix");
    addAndMakeVisible(base_mix_);

    pitch_.BindParam(apvt, "pitch");
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
