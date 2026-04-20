import json
from pathlib import Path

files = [
    ('markdown', '# 02 整合与注释\n\n基于 `01_qc_preprocess.ipynb` 生成的 QC-clean `h5ad` 对象，完成：\n- 合并全部样本并重建未整合嵌入\n- 在可用时比较 `Harmony` 效果\n- 基于 marker 参考完成大类细胞注释\n- 对上皮细胞进行第一版恶性初筛\n- 导出 `malignant_epithelial_screen`、`CAF`、`myeloid`、`T_NK`、`Schwann_neural` 五条后续分析主线\n\n> 当前版本以“analysis-ready object”为目标，不追求最终论文级终版命名。所有 PNI 趋势解释保持探索性表述。'),
    ('markdown', '## 1. 参数、路径与依赖检查\n\n运行前请确认：\n1. `01_qc_preprocess.ipynb` 已运行完成，并产生 `results/intermediate/qc_clean/*.h5ad`\n2. `results/tables/qc/qc_clean_object_index.tsv` 存在且与对象文件一致\n3. `references/marker_reference.tsv` 可正常读取\n\n本 notebook 的约定：\n- Markdown 说明使用中文\n- 代码中的 warning / exception / runtime message 使用英文\n- `Harmony` 为可选步骤，缺包时自动跳过而不阻断主流程'),
]
for idx in range(1, 14):
    if idx in (1, 2):
        files.append(('code', Path(f'tmp_code_cell_{idx}.py').read_text(encoding='utf-8-sig')))
files.append(('markdown', '## 2. 合并对象并重建未整合嵌入\n\n本部分使用 `01` 的 QC-clean 对象直接合并，不回到原始矩阵。目标是建立：\n- 可复跑的未整合 baseline\n- 大类注释与恶性初筛所需的程序分数\n- 后续 subset 导出的统一主对象'))
for idx in (3, 4):
    files.append(('code', Path(f'tmp_code_cell_{idx}.py').read_text(encoding='utf-8-sig')))
files.append(('markdown', '## 3. 可选 `Harmony` 比较\n\n`Harmony` 不是本 notebook 的硬前提。若环境中已安装 `harmonypy`，则运行并与未整合 UMAP 对比；若未安装，则明确跳过并继续后续注释。'))
files.append(('code', Path('tmp_code_cell_5.py').read_text(encoding='utf-8-sig')))
files.append(('markdown', '## 4. 大类注释与 cluster marker 汇总\n\n第一版大类注释采用“marker 分数 + cluster 层级摘要”的组合逻辑：\n- 先对细胞层面打分\n- 再汇总到 cluster\n- 用最高支持证据给出 `cell_class_major`\n\n当前输出偏向可复跑和可交接，不取代后续人工复核。'))
for idx in (6, 7):
    files.append(('code', Path(f'tmp_code_cell_{idx}.py').read_text(encoding='utf-8-sig')))
files.append(('markdown', '## 5. 上皮细胞恶性初筛（v1）\n\n当前只做启发式初筛，不用 CNV 证据替代：\n- `likely_malignant`\n- `epithelial_undetermined`\n- `likely_normal_epithelial`\n\n判读综合 marker、机制程序分数与 `tumor/adjacent` 富集，不把增殖或 stress 极端簇直接等同恶性。'))
for idx in (8, 9):
    files.append(('code', Path(f'tmp_code_cell_{idx}.py').read_text(encoding='utf-8-sig')))
files.append(('markdown', '## 6. 五条 lineage subset 导出\n\n每个 subset 的目标是形成第一版可继续深挖的对象，而不是在本 notebook 中完成全部终版分析。\n\n统一策略：\n- 从主对象按 `cell_class_major` 或恶性标签切出子集\n- 在 subset 内重跑轻量嵌入与聚类\n- 根据对应 marker 程序分数给出 `subset_subtype_v1`\n- 导出对象、marker 表和关键图'))
for idx in (10, 11):
    files.append(('code', Path(f'tmp_code_cell_{idx}.py').read_text(encoding='utf-8-sig')))
files.append(('markdown', '## 7. 结果汇总与 handoff checklist\n\n本部分导出统一 summary，作为后续 `03+` notebook 的直接入口。'))
for idx in (12, 13):
    files.append(('code', Path(f'tmp_code_cell_{idx}.py').read_text(encoding='utf-8-sig')))

cells = []
for cell_type, source in files:
    cells.append({
        'cell_type': cell_type,
        'metadata': {},
        'source': source.splitlines(True),
        **({'execution_count': None, 'outputs': []} if cell_type == 'code' else {}),
    })

nb = {
    'cells': cells,
    'metadata': {
        'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
        'language_info': {'name': 'python', 'version': '3.11'},
    },
    'nbformat': 4,
    'nbformat_minor': 5,
}

Path('notebooks/02_integration_annotation.ipynb').write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding='utf-8')
print('Wrote notebooks/02_integration_annotation.ipynb')
