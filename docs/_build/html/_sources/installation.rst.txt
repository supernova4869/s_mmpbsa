安装
====

本文档介绍s_mmpbsa的安装步骤、环境要求和配置方法，帮助用户快速搭建运行环境。

系统要求
--------

### 操作系统

s_mmpbsa支持以下操作系统：

- **Windows**：Windows 7/8/10/11
- **Linux**：Ubuntu、Debian、CentOS、Rocky等主流Linux发行版

### 硬件要求

- **处理器**：多核处理器（推荐4核及以上）
- **内存**：至少4GB RAM（大型系统推荐8GB及以上）
- **磁盘空间**：至少500MB可用空间

软件依赖
--------

### 基本依赖

s_mmpbsa的核心功能需要以下软件：

- **Gromacs**：用于处理分子动力学轨迹文件。s_mmpbsa支持多个版本的Gromacs，但建议使用较新版本以获得最佳兼容性。

可选依赖：

- **Matplotlib**：用于绘制结果图表。如果您需要使用s_mmpbsa的分析和绘图功能，则需要安装。
- **APBS**：用于计算泊松-玻尔兹曼表面面积。s_mmpbsa内置了APBS内核，但也支持使用外部APBS程序。

### 分子对接重打分功能的特殊依赖

- **PyMOL**：可选软件，用于绘制B因子着色的结构。
- **Gaussian**：可选软件，用于RESP原子电荷计算的DFT计算。
- **Multiwfn**：可选程序（已内置），用于拟合RESP原子电荷。

安装步骤
--------

### Windows系统安装

1. **下载s_mmpbsa**
   
   从GitHub发布页面下载最新版本的Windows可执行文件：
   
   .. code-block:: powershell
      
      # 从GitHub下载s_mmpbsa.exe
      # 访问: https://github.com/supernova4869/s_mmpbsa/releases/latest
   
2. **添加到系统路径**
   
   将s_mmpbsa.exe所在的文件夹添加到系统环境变量PATH中，以便在任何位置都能运行s_mmpbsa。
   
3. **安装Gromacs**
   
   从Gromacs官方网站下载并安装适合您系统的Gromacs版本。
   
4. **安装可选依赖（如需使用分析功能）**
   
   .. code-block:: powershell
      
      # 安装matplotlib
      pip install matplotlib

### Linux系统安装

#### Ubuntu/Debian系统

1. **下载s_mmpbsa**
   
   .. code-block:: bash
      
      # 从GitHub下载最新版本
      wget https://github.com/supernova4869/s_mmpbsa/releases/latest/download/s_mmpbsa
      
      # 添加执行权限
      chmod +x s_mmpbsa
   
2. **添加到系统路径**
   
   .. code-block:: bash
      
      # 将s_mmpbsa移动到/usr/local/bin或其他已在PATH中的目录
      sudo mv s_mmpbsa /usr/local/bin/
   
3. **安装必要依赖**
   
   .. code-block:: bash
      
      # 安装Gromacs（如果尚未安装）
      sudo apt-get update
      sudo apt-get install gromacs
      
      # 安装matplotlib和其他必要的Python包
      sudo apt -y install python3-matplotlib build-essential python-pip

#### CentOS/Rocky系统

1. **下载s_mmpbsa**
   
   .. code-block:: bash
      
      # 从GitHub下载最新版本
      wget https://github.com/supernova4869/s_mmpbsa/releases/latest/download/s_mmpbsa
      
      # 添加执行权限
      chmod +x s_mmpbsa
   
2. **添加到系统路径**
   
   .. code-block:: bash
      
      # 将s_mmpbsa移动到/usr/local/bin或其他已在PATH中的目录
      sudo mv s_mmpbsa /usr/local/bin/
   
3. **安装必要依赖**
   
   .. code-block:: bash
      
      # 安装Gromacs（可能需要从源码编译或使用EPEL仓库）
      # 这里假设您已经安装了Gromacs
      
      # 安装matplotlib和其他必要的Python包
      sudo dnf -y install python3-matplotlib python-pip

验证安装
--------

安装完成后，可以通过以下方式验证s_mmpbsa是否正确安装：

.. code-block:: bash
   
   # 在命令行中运行
   s_mmpbsa --version
   
   # 或者直接运行s_mmpbsa
   s_mmpbsa

如果安装成功，您将看到s_mmpbsa的欢迎信息和版本号。

配置s_mmpbsa
-----------

s_mmpbsa的配置文件为`settings.ini`，该文件包含了程序的各种设置参数。您可以根据需要修改这些参数以优化程序性能或调整计算设置。

### 配置文件位置

- Windows系统：配置文件通常位于s_mmpbsa可执行文件所在的目录
- Linux系统：配置文件通常位于`~/.s_mmpbsa/`目录

### 主要配置参数

配置文件中包含以下主要参数：

- **gmx_path**：Gromacs程序的路径
- **apbs_path**：APBS程序的路径（如果使用外部APBS）
- **nkernels**：并行计算使用的核心数
- **debug_mode**：是否启用调试模式
- **r_cutoff**：非键相互作用的截断距离
- **elec_screen**：电荷筛选方法设置

### 配置文件示例

```ini
[General]
last_opened = ""
debug_mode = false
nkernels = 4

[Program]
gmx_path = "built-in"
apbs_path = "built-in"
delphi_path = ""

[Parameters]
r_cutoff = 1.2
elec_screen = 0

[Display]
font_size = 12
```

常见安装问题
----------

### Gromacs未找到

如果s_mmpbsa无法找到Gromacs程序，请确保：

1. Gromacs已正确安装
2. Gromacs的可执行文件所在目录已添加到系统环境变量PATH中
3. 在settings.ini中正确设置了gmx_path参数

### APBS相关错误

如果使用内置APBS内核出现问题，可以尝试：

1. 安装外部APBS程序
2. 在settings.ini中设置apbs_path参数指向外部APBS程序

### Python/matplotlib相关错误

如果在使用分析功能时出现Python或matplotlib相关错误，请确保：

1. 已安装正确版本的Python（推荐Python 3.6及以上）
2. 已安装matplotlib包：`pip install matplotlib`

### 性能问题

如果计算速度较慢，可以尝试：

1. 在settings.ini中增加nkernels参数的值，利用更多CPU核心
2. 对于大型系统，考虑增加计算的时间间隔（即减少分析的帧数）

获取帮助
--------

如果您在安装过程中遇到任何问题，可以：

- 查看GitHub仓库中的问题页面：https://github.com/supernova4869/s_mmpbsa/issues
- 联系开发者：zhangjiaxing7137@tju.edu.cn
- 加入QQ群：864191465