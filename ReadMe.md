# 东南大学数据结构与算法实践 Project 4
## 项目结构
* header 中存放着所有的.hpp文件
* source 中存放着所有的.cpp文件
## 附加说明
* 所有的.hpp都是对应同名的.cpp的头文件 声明在.hpp中，定义在.cpp中
* Buffer.cpp 定义了 Buffer数据结构，包含判空、判满等函数
* FileOperator.cpp 封装了Buffer 与文件之间交互 以及 其它文件操作有关的函数
* Info.cpp 自定义败者树中存放的数据结构，包含value(值),key(下标),isValid(当前数据在树中是否有效（尽量延长Run 的Size 所需）)
* LoserTree.cpp 败者树数据结构 含基础的败者树相关操作，最终胜者（最小的值存放在 loserArray[0]（Info）中）
* Timer.cpp 封装的Chrono计时器，记录程序运行时间
* main.cpp 程序入口，输入BufferSize 和 LoserTree's Size 以及归并路数k以及Merge BufferSize开始测试
## 使用说明
* 首先运行 $DataProducer.py$ 随机生成指定数量的数据
* 然后运行 $main.cpp$运行，输入相关参数进行测试
## 补充
* 默认随机生成的数据初始存放在 $Input.txt$ 中
* 最后排序的结果在在 $Output.txt$ 中
* 中间生成 Run 的数据在$Run $文件夹里
* 使用测试样例的话，请将 $DataSize= .txt$ 中的数据覆盖到 $Input.txt$ 中
