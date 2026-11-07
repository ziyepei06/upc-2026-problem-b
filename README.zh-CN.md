# UPC 2026 B 题 —— 由两个早期球位预测点球的球门面落点

**230 队** · 2026 年国际物理竞赛（University Physics Competition）· B 题
[English](README.md)

题目给出踢球点和两个实测球位（$t=0.05$ s、$t=0.10$ s）。本仓库包含论文
（LaTeX 源码与编译好的 PDF）、产出论文中每一个数字的分析代码，以及机器可读的结果文件。

## 结论

六个坐标决定六个加速度分量。我们用一个六参数横向力模型——初速 $v_0$（3 个）、
阻力系数 $C_D$（1 个）、两个横向系数——在假设 5 mm 测量噪声下，用自写的阻尼
Gauss–Newton 拟合，并把后验传播到球门面：

| 量 | 值 |
|---|---|
| $x=11$ m 处的落点 | $y=2.785$ m，$z=1.945$ m，$t=0.399$ s |
| 相对门框 | 距右立柱内侧 0.87 m，距横梁 0.50 m |
| 90 % 区间 | $y\in[2.455,3.121]$ m，$z\in[1.617,2.276]$ m |
| 整球入框概率 | 0.971 |
| 模型形式包络（四个完整模型） | $y\in[2.125,3.121]$ m |

主要发现：

- **两个观测点无法确定 $C_D$**：$C_D=0.15$–$0.50$ 时残差几乎不变（12–18 mm），
  落点只在 $z$ 向移动约 5 cm，$y$ 完全不动。
- **自旋与横向初速退化**：在实测升力尺度 $C_L\approx0.7\,S_p$ 下，10 rad/s 的
  自旋让落点横移约 7 cm；数据只在 $\Delta\chi^2\approx1.4$（$p\approx0.24$）上
  轻微偏好负自旋，不显著。
- **模型形式误差主导**：统计区间约 3.5 cm，而模型包络约 1 m。诚实的答案是
  一个带宽，而不是一个点。
- 整条流程用合成真值做了注入–回收验证，论文还研究了"再增加一个观测帧"能带来多少收益。

## 目录结构

```
main.tex  preamble.tex  config.tex  refs.bib   LaTeX 源码；config.tex = 队号与标题
sections/                                       每节一个文件（00-frontmatter … 11-ai-report）
figures/                                        论文使用的矢量 PDF 图
code/                                           分析代码（只用 numpy + matplotlib）
results/                                        论文中每个数字对应的 JSON 输出
notes/                                          推导笔记与"数字溯源"索引
scripts/                                        编译/字数/合规检查脚本
230B.pdf                                        提交版（230 队，B 题）
```

## 编译论文

```bash
make        # latexmk + 提交前检查（页数、摘要长度、命名）
make pdf    # 同时生成 230B.pdf
```

需要带 `latexmk`、`xelatex`/`pdflatex` 与 `bibtex` 的 TeX Live。

## 复现分析

```bash
cd code
python3 fit6.py      # 阻尼 Gauss–Newton 拟合 + 区间   -> ../results/fit6.json
python3 an_valid.py  # 合成真值的注入–回收验证
python3 mkfig_a.py   # 不依赖外部数据的图
```

代码只依赖 `numpy` 与 `matplotlib`：拟合、重要性采样和绘图都没有用 SciPy。
论文里的每个数字在 `notes/numbers.md` 中都能查到来处。

少数绘图脚本还会读参考数据集与平行管线的结果，这些数据不随仓库分发。需要时：

```bash
export UPC_ROOT=$HOME/upc-2026-data
# <UPC_ROOT>/2026/UPC 2026 Problem B 资料库/08_数据集/...
# <UPC_ROOT>/GPT/prediction_revision_v5/results/...
```

## 数据

两个实测球位与场地几何来自题目本身；所有派生量在 `results/`（JSON）。用于标定
测量噪声模型的外部跟踪数据集（SoccerNet-v3D、SkillCorner 样例、StatsBomb 公开数据）
在 `notes/` 中引用，不再分发。

## AI 使用说明

`config.tex` 中 `\aiusedtrue`，论文的 *Report on the Use of AI* 一节记录了所用工具
（DeepSeek Harness，模型 `deepseek-flash`）及分工。第 1 页遵守 UPC 规则：只有队号、
题号、标题和摘要，不出现姓名与学校。

## 许可与引用

MIT，见 `LICENSE` 与 `CITATION.cff`。
