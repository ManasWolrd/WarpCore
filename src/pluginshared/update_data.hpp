#pragma once
#include <juce_core/juce_core.h>
#include "global.hpp"

namespace pluginshared {
class UpdateData {
public:
    static constexpr int kNetworkTimeout = 500; // ms

    UpdateData() {
        is_thread_run_ = update_thread_.startThread();
    }

    ~UpdateData() {
        if (is_thread_run_) {
            update_thread_.Quit();
        }
    }

    void BeginCheck() {
        update_thread_.BeginCheck();
    }

    bool IsComplete() const {
        return update_thread_.is_complete_;
    }

    juce::String GetUpdateMessage() const {
        return update_thread_.GetUpdateMessage();
    }

    juce::String GetButtonLabel() const {
        return update_thread_.GetButtonLabel();
    }

    bool HaveNewVersion() const {
        return update_thread_.have_new_version_;
    }

    static juce::String GetPluginReleaseUrl() {
        return juce::String::formatted("https://github.com/%s/%s/releases/latest", global::kPluginRepoOwnerName,
                                       global::kPluginRepoName);
    }
private:
    class UpdateThread : public juce::Thread {
    public:
        UpdateThread()
            : juce::Thread("version check") {}

        void run() override {
            for (;;) {
                get_sem_.wait();
                if (threadShouldExit()) {
                    return;
                }

                juce::String api_url = juce::String::formatted("https://api.github.com/repos/%s/%s/releases/latest",
                                                               global::kPluginRepoOwnerName, global::kPluginRepoName);
                juce::URL url(api_url);

                auto op =
                    juce::URL::InputStreamOptions{juce::URL::ParameterHandling::inAddress}
                        .withConnectionTimeoutMs(kNetworkTimeout)
                        .withExtraHeaders(juce::String::formatted("User-Agent: %s_Updater", global::kPluginRepoName));
                auto stream = url.createInputStream(op);

                if (!stream) {
                    SetUpdateMessage("Network error or API rate limit");
                    SetButtonLabel("ok");
                    is_complete_ = true;
                    continue;
                }

                // 解析 JSON
                auto response = stream->readEntireStreamAsString();
                auto json = juce::JSON::fromString(response);

                if (json.isObject()) {
                    // 获取最新版本号，例如 "v1.2.0"
                    juce::String latest_version = json.getProperty("tag_name", "").toString();

                    // GitHub 标签通常带 'v'，对比前处理掉
                    juce::String current_version = JucePlugin_VersionString;

                    if (latest_version.startsWithIgnoreCase("v")) latest_version = latest_version.substring(1);

                    if (latest_version != current_version) {
                        juce::String msg;
                        msg << "New version available: " << latest_version << "\n"
                            << "Current version: " << current_version;
                        SetUpdateMessage(msg);
                        SetButtonLabel("ok");
                        have_new_version_ = true;
                    }
                    else {
                        SetUpdateMessage("You are up to date.");
                        SetButtonLabel("ok");
                        have_new_version_ = false;
                    }
                }
                else {
                    SetUpdateMessage("Payload error from GitHub");
                    SetButtonLabel("ok");
                }

                is_complete_ = true;
            }
        }

        void BeginCheck() {
            is_complete_ = false;
            have_new_version_ = false;
            get_sem_.signal();
            {
                juce::ScopedLock lock{data_lock_};
                update_message_.clear();
                button_text_.clear();
            }
        }

        void Quit() {
            get_sem_.signal();
            stopThread(-1);
        }

        juce::String GetUpdateMessage() const {
            juce::ScopedLock lock{data_lock_};
            return update_message_;
        }

        juce::String GetButtonLabel() const {
            juce::ScopedLock lock{data_lock_};
            return button_text_;
        }
    private:
        void SetUpdateMessage(const juce::String& message) {
            juce::ScopedLock lock{data_lock_};
            update_message_ = message;
        }

        void SetButtonLabel(const juce::String& label) {
            juce::ScopedLock lock{data_lock_};
            button_text_ = label;
        }

        friend class UpdateData;
        std::atomic<bool> is_complete_{false};
        std::atomic<bool> have_new_version_{false};
        juce::WaitableEvent get_sem_{0};
        juce::CriticalSection data_lock_;
        juce::String update_message_;
        juce::String button_text_;
    };

    UpdateThread update_thread_;
    bool is_thread_run_{};
};
} // namespace pluginshared
