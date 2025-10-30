import json
import numpy as np
import pandas as pd
import joblib
from pathlib import Path
from sklearn.metrics import (roc_auc_score, roc_curve, precision_recall_curve, auc)
import xgboost as xgb
import matplotlib.pyplot as plt
import re
# ----------------------- 全局配置 -----------------------
RANDOM_SEED = 42 # 与训练一致
N_FOLDS = 5 # 5折
TYPE = 'microprotein' # 测试类型，与训练一致
LABEL_COL = 'dm_cl' # 标签列
CLINICAL_FEATURES = ["BMI", "sex", "age", "waist", "TG", "HDL", "hyper_cl"] # 与训练一致
# 训练输出根目录（与训练代码一致）
output_root = Path(f"glm_xgb_gain_freq_{RANDOM_SEED}_output_{TYPE}")
feature_dir = output_root / "selected_features"
model_dir = output_root / "models"
plot_dir = output_root / "plots_external" # 新目录用于外部测试绘图
plot_dir.mkdir(exist_ok=True)
# 测试集路径
TEST_DATA_ROOT = "./ML_matrix_20250911/metaproteome/independent_test_104/"
TEST_SAMPLE_INFO = f"{TEST_DATA_ROOT}/GNHSF_com_sample_inform_normalized_ML.csv" # 假设样本信息文件路径，根据实际情况调整
TEST_FEATURE_FILE = f"{TEST_DATA_ROOT}/GNHSFc_diann_IGC_humanswiss_{TYPE}_sample.tsv" # 假设Original矩阵路径，根据实际情况调整
# 输出目录用于外部测试结果
external_test_dir = output_root / "external_test"
external_test_dir.mkdir(exist_ok=True)
# ----------------------- 数据预处理函数（类似训练代码） -----------------------
def _to_num(series):
    return pd.to_numeric(series, errors='coerce')
def preprocess_sample_data(filepath):
    """
    预处理样本信息表：
    - dm_cl转换为二元变量（yes=1, no=0）
    - hyper_cl转换为二元变量（yes=1, no=0）
    - sex转换为0/1（female=0, male=1）
    - 数值列：BMI, age, waist, TG, HDL
    - 剔除临床特征缺失的样本
    """
    sampledf = pd.read_csv(
        filepath,
        na_values=["", "NA", "NaN", "nan", "None", "none", "N/A", "n/a"]
    )
    # 转换dm_cl为二元变量
    print("原始dm_cl值分布：", sampledf["dm_cl"].value_counts(dropna=False))
    sampledf["dm_cl"] = sampledf["dm_cl"].map({"yes": 1, "no": 0, "Yes": 1, "No": 0})
    print("转换后dm_cl值分布：", sampledf["dm_cl"].value_counts(dropna=False))
    # sex 映射到0/1
    if "sex" in sampledf.columns:
        sex_map = {
            "female": 0, "f": 0, "0": 0, 0: 0,
            "male": 1, "m": 1, "1": 1, 1: 1
        }
        sampledf["sex"] = sampledf["sex"].apply(
            lambda v: sex_map.get(str(v).strip().lower(), np.nan)
        )
    # hyper_cl 二值化
    if "hyper_cl" in sampledf.columns:
        sampledf["hyper_cl"] = sampledf["hyper_cl"].apply(
            lambda v: 1 if str(v).strip().lower() == "yes" else (0 if str(v).strip().lower() == "no" else np.nan)
        )
    # 数值列
    for col in ["BMI", "age", "waist", "TG", "HDL"]:
        if col in sampledf.columns:
            sampledf[col] = _to_num(sampledf[col])
    # 保留必要列并剔除缺失
    required_cols = ["sample", LABEL_COL] + CLINICAL_FEATURES
    missing_cols = [c for c in required_cols if c not in sampledf.columns]
    if missing_cols:
        raise ValueError(f"样本信息表缺少以下必要列: {missing_cols}")
    sampledf = sampledf[required_cols].dropna(subset=CLINICAL_FEATURES + [LABEL_COL])
    print(f"\n测试样本信息表摘要：\n", sampledf.describe(include="all"))
    return sampledf
def load_feature_data(filename):
    """加载测试特征数据（Original矩阵）"""
    df = pd.read_csv(
        filename,
        sep="\t",
        index_col=0
    )
  
    # 转置并处理列名
    dft = df.T.reset_index()
    dft.columns = ["sample"] + list(dft.columns[1:])
  
    # 简化列名
    def simplify_colname(col):
        return col.split(" ")[0] if " " in col else col
    dft.columns = [simplify_colname(col) for col in dft.columns]
  
    # 清理特殊字符
    dft.columns = [re.sub(r"\|", ".", col) for col in dft.columns]
    dft.columns = [re.sub(r"\-", ".", col) for col in dft.columns]
    dft.columns = [re.sub(r"\.$", "", col) for col in dft.columns]
  
    print(f"加载测试Original矩阵 | 类型：{TYPE} | 样本数：{len(dft)} | 特征数：{len(dft.columns)-1}")
    return dft
# ----------------------- 评估函数（类似训练代码） -----------------------
def evaluate_clf(y_true, y_pred, split_name):
    """评估函数"""
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    roc_auc = auc(fpr, tpr)
    prec, rec, _ = precision_recall_curve(y_true, y_pred)
    pr_auc = auc(rec, prec)
  
    plot_path = plot_dir / f"external_{split_name}_curves.png"
    plot_metrics(
        y_true=y_true,
        y_score=y_pred,
        tag=f"外部测试 | {split_name}\nAUC：{roc_auc:.3f} | PR-AUC：{pr_auc:.3f}",
        save_path=plot_path
    )
  
    return {
        'split': 'external_test',
        'auc': roc_auc,
        'pr_auc': pr_auc,
        'sample_count': len(y_true),
        'positive_count': (y_true == 1).sum()
    }
def plot_metrics(y_true, y_score, tag, save_path):
    """绘图函数"""
    fpr, tpr, _ = roc_curve(y_true, y_score)
    prec, rec, _ = precision_recall_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    pr_auc = auc(rec, prec)
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(fpr, tpr, lw=2, color='#2E86AB', label=f'{tag}')
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Random Guess')
    plt.title('ROC Curve')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
  
    plt.subplot(1, 2, 2)
    plt.plot(rec, prec, lw=2, color='#A23B72', label=f'{tag}')
    plt.title('Precision-Recall Curve')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='lower left')
  
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
# ----------------------- 主流程：外部测试 -----------------------
def main():
    # 加载测试样本信息
    test_sampledf = preprocess_sample_data(TEST_SAMPLE_INFO)
    all_test_samples = test_sampledf["sample"].unique()
    print(f"\n测试集总样本数: {len(all_test_samples)}")
  
    # 加载测试特征矩阵
    dft_test = load_feature_data(TEST_FEATURE_FILE)
  
    # 收集5折高频特征的并集
    omics_sets = []
    fold_features_dict = {}
    for fold in range(N_FOLDS):
        fold_feature_dir = feature_dir / f"fold_{fold}"
        high_freq_file = fold_feature_dir / "xgb_high_freq_features.csv"
        model_features_file = fold_feature_dir / "model_features.csv"
        if not high_freq_file.exists() or not model_features_file.exists():
            print(f" 折 {fold} - 特征文件不存在，跳过")
            continue
        high_freq_df = pd.read_csv(high_freq_file)
        high_freq_features = high_freq_df["feature"].tolist()
        omics_sets.append(set(high_freq_features))
        model_features_df = pd.read_csv(model_features_file)
        fold_model_features = model_features_df["feature"].tolist()
        fold_features_dict[fold] = fold_model_features
  
    if not omics_sets:
        print("无有效折特征，无法进行测试")
        return
  
    union_omics = list(set.union(*omics_sets))
    union_features = union_omics + CLINICAL_FEATURES
    print(f"5折omics特征并集大小: {len(union_omics)} | 总特征: {len(union_features)}")
  
    # 准备测试数据（使用并集特征）
    matched_cols = [col for col in dft_test.columns if col in union_omics or col == "sample"]
    dft_test_filtered = dft_test[matched_cols]
    df_test = pd.merge(
        dft_test_filtered,
        test_sampledf,
        on="sample",
        how="inner"
    )
    missing_features = [feat for feat in union_features if feat not in df_test.columns]
    if missing_features:
        print(f"测试集缺少 {len(missing_features)} 个特征，已填充为 NaN")
        for feat in missing_features:
            df_test[feat] = np.nan
    X_test = df_test[union_features].copy()
    y_test = df_test[LABEL_COL].astype(int).values
    test_samples_list = df_test["sample"].tolist()
    if len(X_test) < 5 or (y_test == 1).sum() == 0 or (y_test == 0).sum() == 0:
        print("测试集样本量不足或标签分布异常，结束")
        return
  
    # 对每个折模型进行预测（使用其自身特征子集）
    fold_preds = []
    for fold in range(N_FOLDS):
        if fold not in fold_features_dict:
            continue
        fold_model_features = fold_features_dict[fold]
        model_file = model_dir / f"fold_{fold}" / "final_xgb_model_best.pkl"
        if not model_file.exists():
            print(f" 折 {fold} - 模型文件不存在，跳过")
            continue
        final_clf = joblib.load(model_file)
        # 选择模型所需特征子集
        X_fold = X_test[fold_model_features].copy()
        preds = final_clf.predict_proba(X_fold)[:, 1]
        fold_preds.append(preds)
        print(f" 折 {fold} - 预测完成")
  
    if not fold_preds:
        print("无有效模型，无法计算指标")
        return
  
    # 平均预测
    avg_preds = np.mean(fold_preds, axis=0)
  
    # 评估
    metrics = evaluate_clf(
        y_true=y_test,
        y_pred=avg_preds,
        split_name='ensemble'
    )
  
    # 保存评估结果
    pd.DataFrame([metrics]).to_csv(
        external_test_dir / "external_ensemble_metrics.csv", index=False
    )
  
    # 保存预测结果
    pred_df = pd.DataFrame({
        'sample': test_samples_list,
        'y_true': y_test,
        'y_pred': avg_preds
    })
    pred_df.to_csv(external_test_dir / "external_ensemble_preds.csv", index=False)
  
    # 打印关键指标
    print(f"\n外部测试集AUC：{metrics['auc']:.3f} | PR-AUC：{metrics['pr_auc']:.3f}")
  
    print(f"\n核心结果路径：")
    print(f" - 外部测试评估指标: {external_test_dir}")
if __name__ == "__main__":
    main()