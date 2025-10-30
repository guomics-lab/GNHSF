#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve,
    roc_auc_score,
    precision_recall_curve,
    auc,
)

# 显示名称与颜色映射
TASKS = [
    # (subdir, display_name, color)
    ("humanprotein", "Human protein", "#86bd42"),
    ("microprotein", "Microbial protein", "#d0919f"),
    ("genus", "Genus", "#2b7c80"),
    ("species", "Species", "#837ba3"),
    ("kegg", "KO", "#ee9634"),
    ("cog", "COG", "#76919f"),
]

SPLITS = {
    "internal": "all_internal_preds.csv",
    "external": "external_ensemble_preds.csv",
}


def load_preds(csv_path: str):
    if not os.path.exists(csv_path):
        warnings.warn(f"File not found, skip: {csv_path}")
        return None

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        warnings.warn(f"Failed to read {csv_path}: {e}")
        return None

    # 期望列：sample, y_true, y_pred
    req_cols = {"sample", "y_true", "y_pred"}
    if not req_cols.issubset(set(df.columns)):
        warnings.warn(f"Missing required columns in {csv_path}, need {req_cols}")
        return None

    # 清理无效行
    df = df.dropna(subset=["y_true", "y_pred"]).copy()
    # 二分类标签
    try:
        y_true = df["y_true"].astype(int).to_numpy()
        y_pred = df["y_pred"].astype(float).to_numpy()
    except Exception as e:
        warnings.warn(f"Type casting error in {csv_path}: {e}")
        return None

    # 需要至少包含两类
    if len(np.unique(y_true)) < 2:
        warnings.warn(f"Only one class present in {csv_path}, cannot compute metrics. Skipped.")
        return None

    return y_true, y_pred


def compute_curves(y_true, y_pred):
    # ROC
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    roc_auc = roc_auc_score(y_true, y_pred)

    # PR
    precision, recall, _ = precision_recall_curve(y_true, y_pred)
    pr_auc = auc(recall, precision)  # PR-AUC（与average_precision_score不同，这里为曲线几何面积）

    # 正例基线（用于PR图的水平线）
    pos_rate = np.mean(y_true)

    return {
        "fpr": fpr,
        "tpr": tpr,
        "roc_auc": roc_auc,
        "precision": precision,
        "recall": recall,
        "pr_auc": pr_auc,
        "pos_rate": pos_rate,
    }


def plot_roc(all_results, title, outpath):
    plt.figure(figsize=(7, 6))
    ax = plt.gca()

    for res in all_results:
        ax.plot(
            res["fpr"],
            res["tpr"],
            label=f'{res["display_name"]} (AUC={res["roc_auc"]:.3f})',
            color=res["color"],
            linewidth=2.0,
        )

    # 去除网格线与 chance 对角线
    ax.grid(False)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(title)
    ax.legend(loc="lower right", frameon=False)
    plt.tight_layout()
    plt.savefig(outpath, format="pdf")
    plt.close()


def plot_pr(all_results, title, outpath):
    plt.figure(figsize=(7, 6))
    ax = plt.gca()

    # 绘制PR曲线
    for res in all_results:
        ax.plot(
            res["recall"],
            res["precision"],
            label=f'{res["display_name"]} (PR-AUC={res["pr_auc"]:.3f})',
            color=res["color"],
            linewidth=2.0,
        )

    # 仍保留辅助基线（正例率中位数），如不需要可注释掉以下块
    if all_results:
        median_pos_rate = float(np.median([r["pos_rate"] for r in all_results]))
        ax.hlines(
            median_pos_rate,
            xmin=0,
            xmax=1,
            colors="#cccccc",
            linestyles="--",
            linewidth=1.0,
            label=f"Baseline ~ {median_pos_rate:.3f}",
        )

    # 去除网格线
    ax.grid(False)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(title)
    ax.legend(loc="lower left", frameon=False)
    plt.tight_layout()
    plt.savefig(outpath, format="pdf")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot ROC (AUC) and PR (PR-AUC) curves from internal/external prediction CSVs."
    )
    parser.add_argument("--root", type=str, default="./output/", help="Project root containing task subfolders.")
    parser.add_argument("--outdir", type=str, default="./figures", help="Output directory for PDF figures.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 使用更清晰的风格
    try:
        plt.style.use("seaborn-v0_8-paper")
    except Exception:
        pass

    # 先收集两个 split 的结果
    results_by_split = {}
    for split, filename in SPLITS.items():
        results = []
        for subdir, display_name, color in TASKS:
            csv_path = os.path.join(args.root, subdir, filename)
            data = load_preds(csv_path)
            if data is None:
                continue
            y_true, y_pred = data
            curves = compute_curves(y_true, y_pred)
            curves.update(
                {
                    "display_name": display_name,
                    "color": color,
                    "subdir": subdir,
                    "csv": csv_path,
                }
            )
            results.append(curves)
        results_by_split[split] = results

    # 基于 external 的 ROC AUC 从高到低确定统一图例顺序
    external_results = results_by_split.get("external", [])
    if external_results:
        # 显示名称 -> external ROC AUC
        ext_auc_map = {r["display_name"]: r["roc_auc"] for r in external_results}
        sorted_names = sorted(ext_auc_map.keys(), key=lambda k: ext_auc_map[k], reverse=True)
    else:
        # 若无 external，可退回到 TASKS 原顺序
        sorted_names = [t[1] for t in TASKS]

    # 对每个 split 的结果按上述顺序排序（没有 external 的条目排在末尾，保持 TASKS 顺序）
    default_order = [t[1] for t in TASKS]
    for split, results in results_by_split.items():
        if not results:
            warnings.warn(f"No valid results for split: {split}. Skipping plots.")
            continue

        # name 排序键：先按 external 排序表位置，再按默认顺序兜底
        def sort_key(r):
            name = r["display_name"]
            return (
                sorted_names.index(name) if name in sorted_names else len(sorted_names) + default_order.index(name)
            )

        results.sort(key=sort_key)

        # 绘制并保存
        roc_title = f"{split.capitalize()} test ROC"
        pr_title = f"{split.capitalize()} test PR"

        roc_out = os.path.join(args.outdir, f"{split}_roc.pdf")
        pr_out = os.path.join(args.outdir, f"{split}_pr.pdf")

        plot_roc(results, roc_title, roc_out)
        plot_pr(results, pr_title, pr_out)

        # 控制台打印汇总（保持 external AUC 排序）
        print(f"[{split}] saved:", roc_out, pr_out)
        for r in results:
            print(
                f'  - {r["display_name"]}: AUC={r["roc_auc"]:.4f}, PR-AUC={r["pr_auc"]:.4f} (pos_rate={r["pos_rate"]:.4f}) from {r["csv"]}'
            )


if __name__ == "__main__":
    main()