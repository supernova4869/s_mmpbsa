#!/bin/bash

# 帮助用户安装依赖并生成s_mmpbsa文档的脚本

# 设置脚本以在出错时立即退出
set -e

# 打印欢迎信息
echo "欢迎使用s_mmpbsa文档构建脚本！"
echo "此脚本将帮助您安装必要的依赖并生成文档。"
echo ""

# 检查是否安装了Python
if ! command -v python3 &> /dev/null
then
    echo "错误：未找到Python 3。请先安装Python 3。"
exit 1
fi

# 检查是否安装了pip
if ! command -v pip3 &> /dev/null
then
    echo "错误：未找到pip3。请先安装pip3。"
exit 1
fi

# 安装依赖
echo "正在安装文档构建依赖..."
pip3 install -r requirements.txt

# 检查是否安装了sphinx-build
if ! command -v sphinx-build &> /dev/null
then
    echo "错误：安装依赖后仍未找到sphinx-build。请检查pip安装。"
exit 1
fi

# 构建HTML文档
echo ""
echo "开始构建HTML文档..."
echo ""
make html

# 检查构建是否成功
if [ -f _build/html/index.html ]
then
    echo ""
echo "文档构建成功！"
echo "您可以在以下路径找到生成的HTML文档："
echo "$(pwd)/_build/html/index.html"
echo ""
echo "要在浏览器中打开文档，请运行："
echo "在Linux上: xdg-open _build/html/index.html"
echo "在Mac上: open _build/html/index.html"
echo ""
echo "构建其他格式的文档："
echo "- PDF: make latexpdf (需要安装LaTeX)"
echo "- EPUB: make epub"
echo "- 手册页: make man"
echo ""
echo "要查看所有可用的构建目标，请运行：make help"
echo ""
echo "祝您使用愉快！"
else
    echo "错误：文档构建失败。请查看上面的错误信息。"
exit 1
fi