API文档
=======

本文档为开发者提供s_mmpbsa的内部结构和主要功能的详细信息，帮助您理解、扩展或修改s_mmpbsa的功能。

项目结构
--------

s_mmpbsa的项目结构如下：

.. code-block:: bash
   
   s_mmpbsa/
   ├── src/
   │   ├── main.rs           # 程序入口点
   │   ├── mmpbsa.rs         # MM-PBSA计算的主要实现
   │   ├── analyzation.rs    # 结果分析功能
   │   ├── parse_tpr.rs      # TPR文件解析
   │   ├── parse_ndx.rs      # NDX文件解析
   │   ├── parse_xtc.rs      # XTC文件解析
   │   ├── pdb.rs            # PDB文件处理
   │   ├── pbsa.rs           # PB/SA计算
   │   ├── mm.rs             # MM计算
   │   ├── utils.rs          # 工具函数
   │   └── settings.rs       # 设置管理
   ├── examples/             # 示例文件
   ├── docs/                 # 文档
   ├── Cargo.toml            # Rust依赖管理
   └── README.md             # 项目说明

主要模块
--------

下面详细介绍s_mmpbsa的主要模块及其功能。

main模块
--------

main模块是程序的入口点，负责处理命令行参数、初始化环境并协调其他模块的工作。

**主要功能**:
- 解析命令行参数
- 初始化程序环境
- 加载输入文件
- 协调MM-PBSA计算
- 提供交互式命令行界面

**核心函数**:

.. code-block:: rust
   
   // 程序入口点
   fn main() {}
   
   // 显示欢迎信息
   fn welcome() {}
   
   // 验证文件是否有效
   fn confirm_file_validity(path: &str) -> bool {}
   
   // 获取内置程序路径
   fn get_built_in_gmx() -> String {}

mmpbsa模块
----------

mmpbsa模块实现了MM-PBSA计算的核心功能，包括能量计算、丙氨酸扫描等。

**主要功能**:
- 执行MM-PBSA计算
- 实现丙氨酸扫描
- 管理临时文件
- 协调MM和PB/SA计算

**核心函数**:

.. code-block:: rust
   
   // 执行MM-PBSA计算
   pub fn fun_mmpbsa_calculations(tpr_path: &str, ...) -> Result<SMResult, Box<dyn Error>> {}
   
   // 实现丙氨酸突变
   pub fn ala_mutate(tpr_path: &str, ...) -> Result<(), Box<dyn Error>> {}
   
   // 设置进度条样式
   fn set_style() -> indicatif::ProgressStyle {}
   
   // 计算MM-PBSA能量
   fn calculate_mmpbsa(...) -> Result<SMResult, Box<dyn Error>> {}
   
   // 计算MM能量
   fn calc_mm(...) -> Result<(Array1<f64>, Array1<f64>), Box<dyn Error>> {}
   
   // 计算PB/SA能量
   fn calc_pbsa(...) -> Result<(Array1<f64>, Array1<f64>), Box<dyn Error>> {}

analyzation模块
---------------

analyzation模块实现了结果分析功能，包括结果的处理、可视化和导出。

**主要功能**:
- 处理MM-PBSA计算结果
- 提供结果可视化
- 导出结果数据
- 支持各种分析操作

**核心数据结构和函数**:

.. code-block:: rust
   
   // 存储MM-PBSA计算结果的数据结构
   pub struct SMResult {
       pub dh: Array1<f64>,          // 焓变
       pub mm: Array1<f64>,          // MM能量
       pub pb: Array1<f64>,          // PB能量
       pub sa: Array1<f64>,          // SA能量
       pub time: Array1<f64>,        // 时间点
       pub residues: Vec<String>,    // 残基名称
       pub res_energy: Array2<f64>,  // 残基能量
       // ... 其他字段
   }
   
   // 分析功能主控制器
   pub fn analyze_controller(result: &SMResult, ...) -> Result<(), Box<dyn Error>> {}
   
   // 获取时间范围
   pub fn get_time_range(result: &SMResult) -> (f64, f64) {}
   
   // 获取时间点对应的索引
   pub fn get_time_index(result: &SMResult, time: f64) -> usize {}

parse_tpr模块
-------------

parse_tpr模块负责解析Gromacs的TPR文件，提取系统的拓扑信息和原子参数。

**主要功能**:
- 解析TPR文件格式
- 提取原子类型、电荷、质量等信息
- 构建系统拓扑结构
- 提供对拓扑数据的访问接口

parse_ndx模块
-------------

parse_ndx模块负责解析Gromacs的NDX文件，提取系统的分组信息。

**主要功能**:
- 解析NDX文件格式
- 提取分组名称和原子索引
- 提供对分组数据的访问接口

parse_xtc模块
-------------

parse_xtc模块负责解析Gromacs的XTC文件，提取系统的坐标信息。

**主要功能**:
- 解析XTC文件格式
- 提取原子坐标数据
- 支持轨迹的随机访问
- 处理大型轨迹文件

pdb模块
-------

pdb模块负责处理PDB文件，包括读取、修改和写入PDB文件。

**主要功能**:
- 读取PDB文件
- 修改PDB文件中的原子坐标和属性
- 写入PDB文件
- 支持将能量信息编码到PDB文件中

pbsa模块
--------

pbsa模块实现了PB和SA能量计算的功能，包括调用外部程序（如APBS）进行计算。

**主要功能**:
- 准备PB计算的输入文件
- 调用APBS进行PB计算
- 计算SA能量
- 处理PB/SA计算的结果

mm模块
------

mm模块实现了MM能量计算的功能，包括范德华和静电相互作用的计算。

**主要功能**:
- 计算范德华相互作用能
- 计算静电相互作用能
- 实现距离截断优化
- 支持并行计算

utils模块
---------

utils模块提供了各种通用工具函数，供其他模块使用。

**主要功能**:
- 文件操作
- 字符串处理
- 数学计算
- 系统调用

settings模块
------------

settings模块负责管理程序的设置，包括读取、修改和保存设置。

**主要功能**:
- 读取settings.ini文件
- 提供设置的访问接口
- 保存设置更改
- 管理程序路径配置

关键数据结构
-------------

SMResult结构体
~~~~~~~~~~~~~

SMResult结构体是s_mmpbsa的核心数据结构，用于存储MM-PBSA计算的结果。

**主要字段**:
- **dh**: 焓变数组
- **mm**: MM能量数组
- **pb**: PB能量数组
- **sa**: SA能量数组
- **time**: 时间点数组
- **residues**: 残基名称列表
- **res_energy**: 残基能量矩阵
- **atom_energy**: 原子能量矩阵

**主要方法**:
- **new()**: 创建SMResult实例
- **to_bin()**: 将结果序列化到二进制文件
- **from_bin()**: 从二进制文件反序列化结果

使用s_mmpbsa作为库
-----------------

s_mmpbsa也可以作为Rust库使用，供其他Rust程序调用其功能。

**示例代码**:

.. code-block:: rust
   
   use s_mmpbsa::mmpbsa::fun_mmpbsa_calculations;
   use s_mmpbsa::analyzation::SMResult;
   
   fn main() -> Result<(), Box<dyn std::error::Error>> {
       // 设置计算参数
       let tpr_path = "path/to/md.tpr";
       let xtc_path = "path/to/md_xtc.xtc";
       let ndx_path = "path/to/index.ndx";
       let rec_group = 0;  // 受体组索引
       let lig_group = 1;  // 配体组索引
       let time_interval = 1.0;  // 时间间隔（ns）
       let temp = 298.15;  // 温度（K）
       let conc = 0.15;    // 盐浓度（mol/L）
       
       // 执行MM-PBSA计算
       let result = fun_mmpbsa_calculations(
           tpr_path,
           xtc_path,
           ndx_path,
           rec_group,
           lig_group,
           time_interval,
           temp,
           conc,
       )?;
       
       // 处理计算结果
       println!("平均结合能: {:.2} kJ/mol", result.dh.mean().unwrap());
       
       // 保存结果到文件
       result.to_bin("result.sm")?;
       
       Ok(())
   }

扩展s_mmpbsa
-----------

如果您想扩展s_mmpbsa的功能，可以考虑以下几个方面：

1. **添加新的能量计算方法**：可以在mm模块和pbsa模块中添加新的能量计算方法。

2. **支持新的输入文件格式**：可以在parse_tpr、parse_ndx和parse_xtc模块中添加对新文件格式的支持。

3. **增强分析功能**：可以在analyzation模块中添加新的分析方法和可视化功能。

4. **优化性能**：可以优化计算核心，提高计算速度和内存使用效率。

5. **添加新的溶剂化模型**：可以添加对其他溶剂化模型的支持，如GB模型、3D-RISM等。

贡献指南
--------

如果您想为s_mmpbsa项目做出贡献，请遵循以下步骤：

1. Fork GitHub仓库
2. 创建您的特性分支
3. 提交您的更改
4. 推送到您的分支
5. 创建新的Pull Request

在提交代码前，请确保您的代码符合项目的编码规范，并且通过了所有测试。

更多信息
--------

- :doc:`usage`：使用指南
- :doc:`installation`：安装说明
- :doc:`faq`：常见问题解答