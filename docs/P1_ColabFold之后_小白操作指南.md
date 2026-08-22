# PocketState-AE：ColabFold 完成后的 P1 操作指南

这份指南只处理已经存在的 ColabFold 结果，不修改、不覆盖、也不重新运行原有的 ColabFold/Slurm 脚本。

## 1. 这一步要完成什么

现有 ColabFold 输出仍只是原始预测文件。P1 后处理要完成：

1. 检查每条 FASTA 是否真的生成完整结果，防止“作业退出正常但某条序列失败”。
2. 冻结由 5YWY–7UR 和 7D7M–P2E 共同定义的 EP4 6 Å 口袋。
3. 对每个模型计算全局、7TM 核心和 38 个口袋残基的 pLDDT/PAE。
4. 排除缺少口袋残基、口袋置信度过低或存在严重 Cα 碰撞的模型。
5. 仅把 WT 构象写入 PocketState-AE 主训练清单；突变体不会被混入状态学习。
6. 从固定 19 个代表残基提取 196 维几何特征。
7. 加入 5YWY 失活态和 7D7M 活化态作为参照，计算口袋 RMSD。

## 2. 运行前确认三个路径

在原服务器配置中确认：

- `OUTPUT_ROOT`：ColabFold 输出总目录。
- `RUN_LABEL`：本次运行标签，例如 `screen`。
- `BATCH_MANIFEST`：提交任务时使用的 `manifest.tsv`。

实际结果目录应类似：

```text
$OUTPUT_ROOT/$RUN_LABEL/
├─ chunk_0001/
├─ chunk_0002/
└─ chunk_0003/
```

每个 chunk 内应存在 `*.done.txt`、`*_scores_rank_001_*.json` 和 `*_rank_001_*.pdb` 等文件。

## 3. 安装后处理依赖

进入项目根目录：

```bash
cd /你的绝对路径/PocketState-AE
python -m pip install -r requirements-p1.txt
```

这里只增加 NumPy；不会更改 ColabFold 容器。

## 4. 最省事的一条命令

先加载原有服务器配置：

```bash
source server/config.private.env
```

然后运行：

```bash
python scripts/run_p1_after_colabfold.py \
  --input-root "$OUTPUT_ROOT/$RUN_LABEL" \
  --batch-manifest "$BATCH_MANIFEST" \
  --run-label "$RUN_LABEL" \
  --mode sampling
```

如果现有任务是 WT 与突变体比较，而不是同一 WT 的多构象采样，把最后一项改成：

```bash
--mode variant_screen
```

无论使用哪种模式，突变体都不会进入 PocketState-AE 主训练清单。

## 5. 在 Slurm 上作为独立 CPU 作业运行

也可以使用新增加、但不影响旧脚本的独立入口：

```bash
source server/config.private.env
export POSTPROCESS_MODE=sampling
sbatch --dependency=afterok:你的预测作业ID \
  server/slurm/colabfold_postprocess.sbatch
```

如果预测作业早已完成，可直接提交：

```bash
sbatch server/slurm/colabfold_postprocess.sbatch
```

该作业是 CPU 任务，不会启动新的 ColabFold GPU 预测。

## 6. 输出在哪里

原始结果不会改变。新增结果位于：

```text
data/processed/
├─ pocket_definition/
│  ├─ pocket_residues.csv
│  └─ pocket_definition.json
├─ conformers/<RUN_LABEL>/
│  ├─ colabfold_models.csv
│  ├─ query_summary.csv
│  ├─ conformer_manifest.csv
│  ├─ provenance.json
│  └─ coordinates/
└─ features/<RUN_LABEL>/
   ├─ pocket_features.csv
   ├─ feature_qc.csv
   └─ feature_schema.json
```

### `colabfold_models.csv`

每行是一个 ColabFold 模型。重点查看：

- `qc_status`：`PASS`、`WARN` 或 `FAIL`。
- `qc_fail_reasons`：不通过原因。
- `core_mean_plddt`：19–345 位受体核心平均置信度。
- `pocket_mean_plddt`：38 个实验口袋残基平均置信度。
- `pocket_fraction_ge70`：口袋中 pLDDT≥70 的比例。
- `pocket_mean_pae`：口袋内部平均 PAE。
- `ca_clash_count`：严重 Cα 空间碰撞数量。

### `conformer_manifest.csv`

这是进入 PocketState-AE 特征提取的冻结清单。它只包含：

- WT；
- 口袋残基完整；
- 自动 QC 为 `PASS`，或只有可追踪警告的 `WARN`；
- 已复制到派生目录的结构。

### `pocket_features.csv`

前四列是模型身份和文件来源，后面 196 列才是 AE 输入特征：

- 19 个固定残基两两 Cα 距离：171 维；
- 19 个侧链代表原子到口袋中心的距离：19 维；
- 口袋回转半径、接触比例等汇总：6 维。

pLDDT、PAE、模型来源和参考状态标签没有放入 AE 输入，防止模型只按“质量”或“来源”聚类。

## 7. 必须执行的人工检查

自动脚本不能决定结构是否具有真实生物学意义。请人工完成：

1. 打开 `query_summary.csv`，确认每个输入 ID 都找到模型。
2. 检查所有 `WARN` 和 `FAIL` 原因。
3. 用 ChimeraX 查看保留结构，检查低置信尾部、链穿插和口袋区域。
4. 对比 5YWY、7D7M 与预测结构的口袋形态。
5. 确认突变体没有出现在 `conformer_manifest.csv`。
6. 检查不同 seed/rank 是否提供了真正不同的 WT 口袋构象，而不是重复结构。

## 8. PocketState-AE 训练前的 Go/No-Go

在真实结果同步回来之前，代码只能验证流程，不能声称已经获得有效构象集合。

训练前至少确认：

- `conformer_manifest.csv` 中有足够多的 WT 构象；
- 目标方案建议最终质量合格构象不少于 200；
- 构象差异不是由低 pLDDT、缺失残基或原子碰撞驱动；
- 5YWY 和 7D7M 的特征行完整；
- `feature_schema.json` 显示 `feature_count` 为 196；
- 同一个原始结构没有因复制或重复文件被多次计入。

如果当前结果只有“1 个 WT + 6 个单点突变体，每条一个模型”，它可以用于测试后处理和突变敏感性，但不能作为 PocketState-AE 主训练集。此时应停止在 P1，不要用 7 个样本强行训练自编码器。

## 9. 常见报错

### 缺少 `.done.txt`

说明对应序列未完整结束或结果没有完整同步。先检查 ColabFold 日志和服务器原目录，不要手工创建 `.done.txt`。

### pLDDT 长度不是 488

说明结果与输入序列不对应、输出不完整或使用了不同构建体。不要把它加入当前全长 PTGER4 集合。

### 没有 WT 构象通过 QC

查看 `colabfold_models.csv` 中的失败原因。不要降低阈值来强行制造训练数据。

### 只有三个 feature rows

通常表示只有一个 WT 模型通过，再加上 5YWY 和 7D7M 两个参照。这只能证明脚本贯通，不能训练 AE。

### NumPy 缺失

运行：

```bash
python -m pip install -r requirements-p1.txt
```

## 10. 本阶段完成标准

- 原始 ColabFold 输出保持不变；
- 所有 FASTA ID 都通过逐项产物验证；
- 38 位点口袋定义和输入哈希已冻结；
- `colabfold_models.csv`、`conformer_manifest.csv` 已生成；
- 196 维特征矩阵和 schema 已生成；
- WT 与突变体严格分开；
- 人工检查和样本量门槛通过后，才进入下一阶段 AE/PCA/RMSD 基线训练。
