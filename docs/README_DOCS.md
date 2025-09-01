# s_mmpbsa 文档指南

## 文档概述

本文档是s_mmpbsa项目的官方文档，提供了项目的详细介绍、安装指南、使用方法和API参考等内容。文档使用Sphinx构建，可以生成多种格式（HTML、PDF、EPUB等）。

## 文档结构

本文档包含以下几个主要部分：

- **首页**：项目概述和主要功能介绍
- **简介**：MM/PB-SA方法的基本原理和s_mmpbsa的优势
- **安装**：系统要求和安装步骤
- **快速入门**：基本使用流程和示例
- **使用指南**：详细的功能说明和参数设置
- **API文档**：开发者参考文档
- **常见问题**：常见问题解答

## 构建文档

要构建本地文档，您需要先安装Python和Sphinx。以下是构建文档的步骤：

### 1. 安装依赖

首先，安装文档构建所需的Python依赖：

```bash
# 在项目根目录下执行
pip install -r docs/requirements.txt
```

### 2. 构建HTML文档

#### 在Linux/Mac系统上：

```bash
cd docs
make html
```

或者使用提供的简化脚本：

```bash
cd docs
chmod +x build_docs.sh
./build_docs.sh
```

#### 在Windows系统上：

```batch
cd docs
make.bat html
```

### 3. 查看文档

构建完成后，可以在`docs/_build/html`目录中找到生成的HTML文档。使用浏览器打开`index.html`文件即可查看完整文档。

## 构建其他格式的文档

除了HTML格式，您还可以构建其他格式的文档：

### PDF格式

```bash
# Linux/Mac
cd docs
make latexpdf

# Windows
cd docs
make.bat latexpdf
```

注意：构建PDF文档需要安装LaTeX发行版（如TeX Live、MiKTeX等）。

### EPUB格式

```bash
# Linux/Mac
cd docs
make epub

# Windows
cd docs
make.bat epub
```

## 在ReadTheDocs上查看文档

项目文档也托管在ReadTheDocs平台上，您可以直接在浏览器中访问：

- 稳定版：https://s-mmpbsa.readthedocs.io/zh_CN/stable/
- 开发版：https://s-mmpbsa.readthedocs.io/zh_CN/latest/

## 贡献文档

如果您发现文档中有错误或有改进建议，欢迎贡献您的力量：

1. 提出Issue描述问题或建议
2. 提交Pull Request修复问题或添加新内容

文档使用reStructuredText格式编写，主要内容文件位于`docs`目录下。

## 联系方式

如有任何问题，请联系项目维护者：

- 电子邮件：email@example.com
- QQ群：123456789