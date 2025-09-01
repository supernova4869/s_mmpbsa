.. s_mmpbsa documentation master file, created by
   sphinx-quickstart on 2025-09-01.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

欢迎使用 s_mmpbsa 文档
======================

**s_mmpbsa** 是一个用于计算生物分子结合自由能的高效工具，专门用于分析Gromacs轨迹，采用分子力学泊松-玻尔兹曼表面积（MM/PB-SA）方法。

.. toctree::
   :maxdepth: 2
   :caption: 内容概述

   introduction
   installation
   quick_start
   usage
   features
   examples
   api_reference
   faq
   contribution
   authors

简介
----

**s_mmpbsa** 提供了一个便捷的界面（类似 `Multiwfn`）来计算Gromacs轨迹的结合自由能。相比于其他同类工具，它具有安装简便、运行高效、跨平台等优势，且支持多种高级功能，如电荷筛选效应和构象熵的计算。

主要功能
--------

- **MD模拟结合能计算**：从分子动力学模拟结果计算生物分子间的结合自由能
- **分子对接结果重打分**：为分子对接结果提供更准确的结合能预测
- **蛋白质-配体复合物丙氨酸扫描**：分析蛋白质中关键残基对结合的贡献

特点
----

- 开源免费，遵循LGPL许可证
- 环境依赖少，Linux系统仅需Gromacs程序，绘图功能需Python环境
- 使用Rust语言开发，性能优异
- 交互式操作，无需编写参数文件
- 考虑电荷筛选效应，如文献 [J. Chem. Inf. Model. 2021, 61, 2454] 所述
- 考虑构象熵，如文献 [J. Chem. Phys. 2017, 146, 124124] 所述
- 可存储分析结果，便于进一步的可重复分析

开始使用
--------

请参考 :doc:`installation` 章节安装s_mmpbsa，然后查看 :doc:`quick_start` 章节了解基本使用流程。

对于详细的使用说明，请参阅 :doc:`usage` 章节，其中包含了各种功能的使用方法和实例。

获取帮助
--------

如果您在使用过程中遇到任何问题，或者有任何改进建议，请联系开发者或加入QQ群：

- **开发者**：张嘉兴博士 (zhangjiaxing7137@tju.edu.cn, 天津大学)
- **QQ群**：864191465

引用
----

如果您在研究工作中使用了s_mmpbsa，请按照以下格式引用：

.. code-block:: text

   Jiaxing Zhang, s_mmpbsa, Version [your version], https://github.com/supernova4869/s_mmpbsa (accessed on yy-mm-dd)

.. note::
   当s_mmpbsa的详细论文发表后，请引用相应的论文而非此处的网页。

索引与表格
----------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`