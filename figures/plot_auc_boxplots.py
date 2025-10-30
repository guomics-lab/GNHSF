import sys
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

# 映射：列名(小写) -> (展示名称, 颜色)
COL_MAP = {
    "humanprotein_auc": ("Human protein", "#86bd42"),
    "microprotein_auc": ("Microbial protein", "#d0919f"),
    "genus_auc": ("Genus", "#2b7c80"),
    "species_auc": ("Species", "#837ba3"),
    "kegg_auc": ("KO", "#ee9634"),
    "cog_auc": ("COG", "#76919f"),
}

# 为了固定展示顺序，定义顺序列表（使用展示名称）
ORDER = ["Human protein", "Microbial protein", "Genus", "Species", "KO", "COG"]

def load_and_melt(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"File not found: {csv_path}")

    df = pd.read_csv(csv_path)
    if df.shape[1] < 7:
        raise ValueError(f"Expect at least 7 columns (seed + 6 AUC columns) in {csv_path}")

    # 统一小写列名，去除空白
    df.columns = [str(c).strip().lower() for c in df.columns]

    # 确定 seed 列（优先使用名为 'seed' 的列，否则使用第一列）
    seed_col = "seed" if "seed" in df.columns else df.columns[0]
    if seed_col != "seed":
        # 重命名为 seed 以便后续 melt
        df = df.rename(columns={seed_col: "seed"})
        seed_col = "seed"

    # 检查并收集存在的 AUC 列
    auc_cols_present = [c for c in COL_MAP.keys() if c in df.columns]
    missing = [c for c in COL_MAP.keys() if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {csv_path}: {missing}\n"
            f"Found columns: {list(df.columns)}"
        )

    # 转成长表
    long_df = df.melt(
        id_vars=[seed_col],
        value_vars=auc_cols_present,
        var_name="metric",
        value_name="AUC",
    )

    # 显示名称与颜色
    long_df["display_name"] = long_df["metric"].map(lambda k: COL_MAP[k][0])
    long_df["color"] = long_df["metric"].map(lambda k: COL_MAP[k][1])

    # 转为数值型
    long_df["AUC"] = pd.to_numeric(long_df["AUC"], errors="coerce")

    # 丢弃空值
    long_df = long_df.dropna(subset=["AUC"])

    return long_df

def plot_boxplot(long_df: pd.DataFrame, title: str, out_pdf: str) -> None:
    # 使用固定顺序与调色板
    palette = {name: color for key, (name, color) in COL_MAP.items()}

    plt.figure(figsize=(8, 5), dpi=150)
    ax = sns.boxplot(
        data=long_df,
        x="display_name",
        y="AUC",
        order=ORDER,
        palette=palette,
        width=0.6,
        fliersize=3,
        linewidth=1,
    )
    sns.stripplot(
        data=long_df,
        x="display_name",
        y="AUC",
        order=ORDER,
        color="black",
        size=2.5,
        alpha=0.4,
        jitter=0.15,
        dodge=False,
    )

    ax.set_xlabel("")
    ax.set_ylabel("AUC")
    ax.set_title(title)

    # 纵轴设置：范围0.6~0.8；minor刻度每0.05，major刻度每0.1，仅major显示标签
    ax.set_ylim(0.6, 0.8)
    ax.yaxis.set_major_locator(MultipleLocator(0.1))   # 0.6, 0.7, 0.8
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))  # 0.65, 0.75
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.tick_params(axis='y', which='major', length=5)
    ax.tick_params(axis='y', which='minor', length=3, labelsize=0)

    # 横坐标标签斜着显示，右对齐并锚定，减少重叠
    for label in ax.get_xticklabels():
        label.set_rotation(45)              # 角度可改为 30/60 等
        label.set_horizontalalignment('right')
        label.set_rotation_mode('anchor')

    plt.tight_layout()
    plt.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close()

def main():
    # 默认文件名
    external_csv = "./output/seed_auc_summary_external.csv"
    internal_csv = "./output/seed_auc_summary_internal.csv"

    # 允许命令行参数覆盖
    # 用法：python plot_seed_auc_boxplots.py external.csv internal.csv
    if len(sys.argv) >= 2:
        external_csv = sys.argv[1]
    if len(sys.argv) >= 3:
        internal_csv = sys.argv[2]

    # External
    ext_long = load_and_melt(external_csv)
    plot_boxplot(
        ext_long,
        title="External test set",
        out_pdf="seed_auc_boxplot_external.pdf",
    )

    # Internal
    int_long = load_and_melt(internal_csv)
    plot_boxplot(
        int_long,
        title="Internal test set",
        out_pdf="seed_auc_boxplot_internal.pdf",
    )

    print("Saved:")
    print(" - seed_auc_boxplot_external.pdf")
    print(" - seed_auc_boxplot_internal.pdf")

if __name__ == "__main__":
    main()