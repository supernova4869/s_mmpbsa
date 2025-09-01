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

- **Gromacs**：用于处理分子动力学轨迹文件。s_mmpbsa内置了Gromacs，但也支持使用外部Gromacs程序。支持多个版本的Gromacs，但建议使用较新版本以获得最佳兼容性。

可选依赖：

- **Matplotlib**：用于绘制结果图表。如果您需要使用s_mmpbsa的分析和绘图功能，则需要安装。
- **APBS**：用于计算泊松-玻尔兹曼表面面积。s_mmpbsa内置了APBS内核，但也支持使用外部APBS程序。
- **PyMOL**：用于绘制B因子着色的结构。

部署方法
--------

### Windows系统

1. **下载s_mmpbsa**
   
   从GitHub发布页面下载最新版本的Windows可执行文件：
   
   .. code-block:: powershell
      
      # 从GitHub下载s_mmpbsa.exe
      # 访问: https://github.com/supernova4869/s_mmpbsa/releases/latest
   
2. **添加到系统路径**
   
   将s_mmpbsa.exe所在的文件夹添加到系统环境变量PATH中，以便在任何位置都能运行s_mmpbsa。
   
3. **安装可选依赖（如需使用分析功能）**
   
   .. code-block:: powershell
      
      # 安装matplotlib
      pip install matplotlib

### Linux系统

1. **下载s_mmpbsa**
   
   .. code-block:: bash
      
      # 从GitHub下载最新版本
      wget https://github.com/supernova4869/s_mmpbsa/releases/latest/download/s_mmpbsa
      
      # 添加执行权限
      chmod +x s_mmpbsa
   
2. **添加到系统路径**
   
   .. code-block:: bash
      
      # 将s_mmpbsa所在的文件夹添加到系统环境变量PATH中，以便在任何位置都能运行s_mmpbsa。
      export PATH=$PATH:/path/to/s_mmpbsa/
   
3. **安装必要依赖**
   
   .. code-block:: bash
      
      # 安装matplotlib和其他必要的Python包
      
      #### Ubuntu/Debian系统
      sudo apt -y install python3-matplotlib build-essential python-pip
      #### CentOS/Rocky系统
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

- 配置文件通常位于s_mmpbsa可执行文件所在的目录
- 程序启动时会自行检查`settings.ini`的文件位置，优先级：当前目录 > 程序所在目录。
- 若找不到`settings.ini`，程序将使用默认设置。

### 主要配置参数

配置文件中包含以下主要参数：

- **gmx_path**：Gromacs程序的路径。如果为"built-in"，程序将使用/programs/gmx/中的gmx程序。
- **apbs_path**：APBS程序的路径。如果为"built-in"，程序将使用/programs/apbs/中的apbs程序。
- **nkernels**：并行计算使用的核心数
- **debug_mode**：是否启用调试模式(y/n)。启用后，中间文件不会删除。
- **r_cutoff**：非键相互作用的截断距离。0为不截断。
- **elec_screen**：静电屏蔽方法设置。0为不使用静电屏蔽。1为使用德拜-休克尔屏蔽。

常见问题
----------

### Gromacs未找到

如果s_mmpbsa无法找到Gromacs程序，请确保：

1. Gromacs已正确安装
2. Gromacs的可执行文件所在目录已添加到系统环境变量PATH中
3. 在settings.ini中正确设置了gmx_path参数

### APBS相关错误

如果使用内置APBS内核出现问题，可以尝试：

1. 确保内置APBS程序具有可执行权限
2. 安装外部APBS程序
3. 在settings.ini中设置apbs_path参数指向外部APBS程序

### Python/matplotlib相关错误

如果在使用分析功能时出现Python或matplotlib相关错误，请确保：

1. 已安装正确版本的Python（推荐Python 3.6及以上）
2. 已安装matplotlib包

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