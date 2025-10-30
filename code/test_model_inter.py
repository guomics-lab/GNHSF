# -*- coding: gbk -*-
import json
import numpy as np
import pandas as pd
import joblib
from pathlib import Path
from sklearn.metrics import (roc_auc_score, roc_curve, precision_recall_curve, auc)
from sklearn.model_selection import StratifiedKFold
import xgboost as xgb
import matplotlib.pyplot as plt
import re
# ----------------------- 全局配置 -----------------------
RANDOM_SEED = 42
N_FOLDS = 5
TYPE = 'microprotein'
LABEL_COL = 'dm_cl'
CLINICAL_FEATURES = ["BMI", "sex", "age", "waist", "TG", "HDL", "hyper_cl"]
# 训练输出根目录
output_root = Path(f"glm_xgb_gain_freq_{RANDOM_SEED}_output_{TYPE}")
feature_dir = output_root / "selected_features"
model_dir = output_root / "models"
plot_dir = output_root / "plots_internal"
plot_dir.mkdir(parents=True, exist_ok=True)
# 内部数据路径
INTERNAL_DATA_ROOT = "/sunyingying/test/ML_matrix_20250911/metaproteome/train_internalVal_1385"
INTERNAL_SAMPLE_INFO = f"{INTERNAL_DATA_ROOT}/GNHSF_sample_inform_normalized_1385_2_yes_all_metadata.csv"
# 输出目录用于内部测试结果
internal_test_dir = output_root / "internal_test"
internal_test_dir.mkdir(parents=True, exist_ok=True)
# ----------------------- 数据预处理函数 -----------------------
def _to_num(series):
    return pd.to_numeric(series, errors='coerce')
def preprocess_sample_data(filepath):
    sampledf = pd.read_csv(
        filepath,
        na_values=["", "NA", "NaN", "nan", "None", "none", "N/A", "n/a"]
    )
    print("原始dm_cl值分布：", sampledf["dm_cl"].value_counts(dropna=False))
    sampledf["dm_cl"] = sampledf["dm_cl"].map({"yes": 1, "no": 0, "Yes": 1, "No": 0})
    print("转换后dm_cl值分布：", sampledf["dm_cl"].value_counts(dropna=False))
    if "sex" in sampledf.columns:
        sex_map = {
            "female": 0, "f": 0, "0": 0, 0: 0,
            "male": 1, "m": 1, "1": 1, 1: 1
        }
        sampledf["sex"] = sampledf["sex"].apply(
            lambda v: sex_map.get(str(v).strip().lower(), np.nan)
        )
    if "hyper_cl" in sampledf.columns:
        sampledf["hyper_cl"] = sampledf["hyper_cl"].apply(
            lambda v: 1 if str(v).strip().lower() == "yes" else (0 if str(v).strip().lower() == "no" else np.nan)
        )
    for col in ["BMI", "age", "waist", "TG", "HDL"]:
        if col in sampledf.columns:
            sampledf[col] = _to_num(sampledf[col])
    required_cols = ["sample", LABEL_COL] + CLINICAL_FEATURES
    missing_cols = [c for c in required_cols if c not in sampledf.columns]
    if missing_cols:
        raise ValueError(f"样本信息表缺少以下必要列: {missing_cols}")
    sampledf = sampledf[required_cols].dropna(subset=CLINICAL_FEATURES + [LABEL_COL])
    print(f"\n样本信息表摘要：\n", sampledf.describe(include="all"))
    return sampledf
def load_feature_data(type_, sample_indices=None, is_original=True, base_path=INTERNAL_DATA_ROOT):
    if is_original:
        filename = f"{base_path}/original_expression_matrix/GNHSF_diann_IGC_humanswiss_{type_}_sample_1385_NA90.tsv"
    else:
        filename = f"{base_path}/INT_expression_matrix/GNHSF_diann_IGC_humanswiss_{type_}_sample_1385_NA90_INT.tsv"
    df = pd.read_csv(
        filename,
        sep="\t",
        index_col=0
    )
    dft = df.T.reset_index()
    dft.columns = ["sample"] + list(dft.columns[1:])
    def simplify_colname(col):
        return col.split(" ")[0] if " " in col else col
    dft.columns = [simplify_colname(col) for col in dft.columns]
    dft.columns = [re.sub(r"\|", ".", col) for col in dft.columns]
    dft.columns = [re.sub(r"\-", ".", col) for col in dft.columns]
    dft.columns = [re.sub(r"\.$", "", col) for col in dft.columns]
    if sample_indices is not None:
        dft = dft[dft["sample"].isin(sample_indices)]
    print(f"加载{'original' if is_original else 'INT'}矩阵 | 类型：{type_} | 样本数：{len(dft)} | 特征数：{len(dft.columns)-1}")
    return dft
# ----------------------- 评估函数 -----------------------
def evaluate_clf(y_true, y_pred, split_name):
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    roc_auc = auc(fpr, tpr)
    prec, rec, _ = precision_recall_curve(y_true, y_pred)
    pr_auc = auc(rec, prec)
    plot_path = plot_dir / f"internal_{split_name}_curves.png"
    plot_metrics(
        y_true=y_true,
        y_score=y_pred,
        tag=f"内部测试 | {split_name}\nAUC：{roc_auc:.3f} | PR-AUC：{pr_auc:.3f}",
        save_path=plot_path
    )
    return {
        'split': 'internal_test',
        'auc': roc_auc,
        'pr_auc': pr_auc,
        'sample_count': len(y_true),
        'positive_count': (y_true == 1).sum()
    }
def plot_metrics(y_true, y_score, tag, save_path):
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
# ----------------------- 主流程：内部测试 -----------------------
def main():
    sampledf = preprocess_sample_data(INTERNAL_SAMPLE_INFO)
    all_samples = sampledf["sample"].unique()
    print(f"\n总样本数: {len(all_samples)}")
    sample_to_label = dict(zip(sampledf["sample"], sampledf[LABEL_COL]))
    y_stratify = np.array([sample_to_label[s] for s in all_samples])
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_SEED)
    folds = list(skf.split(all_samples, y_stratify))
    all_preds = []
    for fold in range(N_FOLDS):
        print(f"\n----- 内部测试折 {fold + 1}/{N_FOLDS} -----")
        train_val_idx, test_idx = folds[fold]
        test_samples = all_samples[test_idx].tolist()
        fold_feature_dir = feature_dir / f"fold_{fold}"
        high_freq_file = fold_feature_dir / "xgb_high_freq_features.csv"
        model_features_file = fold_feature_dir / "model_features.csv"
        if not high_freq_file.exists() or not model_features_file.exists():
            print(f" 折 {fold} - 特征文件不存在，跳过")
            continue
        high_freq_df = pd.read_csv(high_freq_file)
        high_freq_features = high_freq_df["feature"].tolist()
        if not high_freq_features:
            print(f" 折 {fold} - 无高频特征，跳过")
            continue
        model_features_df = pd.read_csv(model_features_file)
        fold_model_features = model_features_df["feature"].tolist()
        dft_test = load_feature_data(TYPE, test_samples, is_original=True, base_path=INTERNAL_DATA_ROOT)
        matched_cols = [col for col in dft_test.columns if col in high_freq_features or col == "sample"]
        dft_test_filtered = dft_test[matched_cols]
        clinical_cols = ["sample"] + CLINICAL_FEATURES + [LABEL_COL]
        df_clinical = sampledf[clinical_cols]
        df_test = pd.merge(
            dft_test_filtered,
            df_clinical,
            on="sample",
            how="inner"
        )
        missing_features = [feat for feat in fold_model_features if feat not in df_test.columns]
        if missing_features:
            print(f" 折 {fold} - 测试集缺少 {len(missing_features)} 个特征，已填充为 NaN")
            for feat in missing_features:
                df_test[feat] = np.nan
        X_test = df_test[fold_model_features].copy()
        y_test = df_test[LABEL_COL].astype(int).values
        test_samples_list = df_test["sample"].tolist()
        if len(X_test) < 5 or (y_test == 1).sum() == 0 or (y_test == 0).sum() == 0:
            print(f" 折 {fold} - 测试集样本量不足或标签分布异常，跳过")
            continue
        model_file = model_dir / f"fold_{fold}" / "final_xgb_model_best.pkl"
        if not model_file.exists():
            print(f" 折 {fold} - 模型文件不存在，跳过")
            continue
        final_clf = joblib.load(model_file)
        preds = final_clf.predict_proba(X_test)[:, 1]
        pred_df = pd.DataFrame({
            'sample': test_samples_list,
            'y_true': y_test,
            'y_pred': preds,
            'fold': fold
        })
        all_preds.append(pred_df)
        pred_df.to_csv(internal_test_dir / f"fold_{fold}_internal_preds.csv", index=False)
        print(f" 折 {fold} - 预测完成")
    if all_preds:
        all_preds_df = pd.concat(all_preds, ignore_index=True)
        all_preds_df.to_csv(internal_test_dir / "all_internal_preds.csv", index=False)
        y_all = all_preds_df['y_true'].values
        pred_all = all_preds_df['y_pred'].values
        metrics = evaluate_clf(
            y_true=y_all,
            y_pred=pred_all,
            split_name='overall'
        )
        pd.DataFrame([metrics]).to_csv(internal_test_dir / "internal_overall_metrics.csv", index=False)
        print("\n===== 所有折内部测试完成 =====")
        print(f"内部测试集整体AUC：{metrics['auc']:.3f} | PR-AUC：{metrics['pr_auc']:.3f}")
    print(f"\n核心结果路径：")
    print(f" - 内部测试评估指标: {internal_test_dir}")
if __name__ == "__main__":
    main()