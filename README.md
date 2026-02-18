# WarpCore
建议先了解一下 prosoniq piwarp 或者 zynaptiq wormhole 插件的效果样子  
WarpCore是一个在**没有反编译**的情况下[构思](doc/warp.md)出的`多段频谱反转`插件，类似于piwarp/wormhole的`时域局部频谱反转`效果  

> [!WARNING]
> 此插件正在开发中

## todolist
- [ ] 双声道处理
- [ ] 优化cpu
- [ ] 修改频率参数适配可变采样率
- [ ] 添加分频器来不对其他频率进行反转
- [ ] 添加warp tilt/pitch带来的共振峰移动效果

# 克隆到本地
```bash
git clone --recurse-submodules https://github.com/ManasWorld/plugin-template.git
```
或者
```bash
git clone https://github.com/ManasWorld/plugin-template.git
cd l-model-plugin-template
git submodule update --init --recursive
```

## MacOS授权

```bash
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.component
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.vst3
sudo xattr -dr com.apple.quarantine /path/to/your/plugins/plugin_name.lv2
```
## 构建

```bash
git clone --recurse-submodules https://github.com/ManasWorld/plugin-template.git

# windows
cmake -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Release -S . -B ./build

# linux
sudo apt update
sudo apt-get install libx11-dev libfreetype-dev libfontconfig1-dev libasound2-dev libxrandr-dev libxinerama-dev libxcursor-dev
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -S . -B ./build
cmake --build ./build --config Release

# macOS
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" -S . -B ./build
cmake --build ./build --config Release
```
