# WarpCore
尝试重现 prosoniq piwarp 或者 zynaptiq wormhole 插件的效果  
WarpCore是一个`多段频谱反转`插件，类似于piwarp/wormhole的`时域局部频谱反转`效果  

> [!WARNING]
> 此插件正在开发中

## featrue
- [x] 可配置极点数滤波器
    包括: 滤波器极点数量 滤波器截止频率缩放
- [ ] 内置带通分频器
    包括: 低分频和高分频(最后再做)
- [ ] 可调整带宽
    包括: 必须是递增的函数的带宽分配
- [x] 共振峰移动
    应该是加在输出的振荡器上面

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
