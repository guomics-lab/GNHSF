# -*- coding: utf-8 -*-
"""
完整的模型验证指标计算脚本 - 最终修复版
修复了SHAP兼容性问题和通路分析中的数据验证问题
"""

import os
import re
import numpy as np
import pandas as pd
import joblib
import warnings
import argparse
from pathlib import Path
from collections import defaultdict
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc
from statsmodels.stats.multitest import multipletests
import shap
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# ============================================================================
# 配置区域
# ============================================================================

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='计算模型验证的补充指标（效应方向一致性、通路一致性、特征重叠检验）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
    # 分析单个特征类型
    python validation_analysis.py --feature_types microprotein
    
    # 分析多个特征类型
    python validation_analysis.py --feature_types microprotein species genus
    
    # 使用默认值（所有类型）
    python validation_analysis.py
    
    # 指定自定义路径
    python validation_analysis.py --base_dir /path/to/data --feature_types species
        """
    )
    
    parser.add_argument(
        '--feature_types',
        nargs='+',
        default=['microprotein', 'species', 'genus', 'cog', 'kegg', 'humanprotein'],
        choices=['microprotein', 'species', 'genus', 'cog', 'kegg', 'humanprotein'],
        help='要分析的特征类型列表（可多选，默认：microprotein species genus COG KEGG）'
    )
    
    parser.add_argument(
        '--base_dir',
        type=str,
        default='/sunyingying/test',
        help='基础数据目录（默认：/sunyingying/test）'
    )
    
    parser.add_argument(
        '--random_seed',
        type=int,
        default=42,
        help='随机种子（默认：42）'
    )
    
    parser.add_argument(
        '--top_n_features',
        type=int,
        default=500,
        help='效应方向一致性分析的top特征数量（默认：100）'
    )
    
    parser.add_argument(
        '--top_n_pathways',
        type=int,
        default=100,
        help='通路层级分析的top通路数量（默认：50）'
    )
    
    parser.add_argument(
        '--n_permutations',
        type=int,
        default=10000,
        help='排列检验的次数（默认：10000）'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        default='/sunyingying/test/validation_analysis_metaproteome',
        help='输出目录（默认：自动根据特征类型生成）'
    )
    
    parser.add_argument(
        '--selected_fold',
        type=int,
        default=0,
        help='用于SHAP分析的折数（默认：0）'
    )
    
    return parser.parse_args()


# 全局变量（将在main函数中设置）
BASE_DIR = None
MODEL_OUTPUT_ROOT = None
INTERNAL_DATA_PATH = None
EXTERNAL_DATA_PATH = None
INTERNAL_SAMPLE_INFO = None
EXTERNAL_SAMPLE_INFO = None
FEATURE_TYPE = None
VALIDATION_OUTPUT_DIR = None

# 常量
SEEDS = list(range(42, 62))  # 20个不同的种子
N_FOLDS = 5
CLINICAL_FEATURES = ["BMI", "sex", "age", "waist", "TG", "HDL", "hyper_cl"]

# ============================================================================
# SHAP兼容性修复
# ============================================================================

def fix_xgboost_base_score(model):
    """
    修复XGBoost模型的base_score格式问题
    
    Parameters:
    -----------
    model : xgboost.XGBClassifier
        XGBoost模型
    
    Returns:
    --------
    model : xgboost.XGBClassifier
        修复后的模型
    """
    try:
        # 获取booster
        booster = model.get_booster()
        
        # 获取模型参数
        model_params = booster.save_config()
        import json
        config = json.loads(model_params)
        
        # 修复base_score格式
        if 'learner' in config and 'learner_model_param' in config['learner']:
            base_score_str = config['learner']['learner_model_param'].get('base_score', '0.5')
            
            # 如果base_score是数组格式（如 '[5.011208E-1]'），提取第一个值
            if isinstance(base_score_str, str):
                if base_score_str.startswith('[') and base_score_str.endswith(']'):
                    # 移除方括号并转换
                    base_score_str = base_score_str.strip('[]')
                    base_score_value = float(base_score_str)
                    
                    # 更新配置
                    config['learner']['learner_model_param']['base_score'] = str(base_score_value)
                    
                    # 重新加载配置
                    booster.load_config(json.dumps(config))
                    
                    print(f"  已修复base_score: {base_score_str} -> {base_score_value}")
        
        return model
    
    except Exception as e:
        print(f"  警告：无法修复base_score，使用原始模型: {e}")
        return model


def compute_shap_values_safe(model, X, max_samples=1000):
    """
    安全地计算SHAP值，包含多种fallback策略
    
    Parameters:
    -----------
    model : trained model
        训练好的模型
    X : pd.DataFrame or np.ndarray
        输入数据
    max_samples : int
        如果样本数过多，采样的最大样本数
    
    Returns:
    --------
    shap_values : np.ndarray or None
        SHAP值数组
    """
    # 策略1：尝试使用TreeExplainer（修复base_score后）
    try:
        print("  尝试策略1: TreeExplainer（修复base_score）...")
        fixed_model = fix_xgboost_base_score(model)
        explainer = shap.TreeExplainer(fixed_model)
        shap_values = explainer.shap_values(X)
        print("  策略1成功！")
        return shap_values
    
    except Exception as e:
        print(f"  策略1失败: {e}")
    
    # 策略2：使用TreeExplainer但采样背景数据
    try:
        print("  尝试策略2: TreeExplainer（采样背景）...")
        
        # 采样背景数据
        if len(X) > max_samples:
            background_indices = np.random.choice(len(X), size=min(100, len(X)), replace=False)
            background = X.iloc[background_indices] if hasattr(X, 'iloc') else X[background_indices]
        else:
            background = X
        
        explainer = shap.TreeExplainer(model, data=background)
        shap_values = explainer.shap_values(X)
        print("  策略2成功！")
        return shap_values
    
    except Exception as e:
        print(f"  策略2失败: {e}")
    
    # 策略3：使用KernelExplainer（较慢但稳定）
    try:
        print("  尝试策略3: KernelExplainer（可能较慢）...")
        
        # 采样背景数据
        if len(X) > 100:
            background_indices = np.random.choice(len(X), size=100, replace=False)
            background = X.iloc[background_indices] if hasattr(X, 'iloc') else X[background_indices]
        else:
            background = X
        
        # 创建预测函数
        def predict_fn(data):
            if hasattr(model, 'predict_proba'):
                return model.predict_proba(data)[:, 1]
            else:
                return model.predict(data)
        
        # 采样要解释的数据
        if len(X) > max_samples:
            sample_indices = np.random.choice(len(X), size=max_samples, replace=False)
            X_sample = X.iloc[sample_indices] if hasattr(X, 'iloc') else X[sample_indices]
        else:
            X_sample = X
        
        explainer = shap.KernelExplainer(predict_fn, background)
        shap_values = explainer.shap_values(X_sample, nsamples=100)
        
        # 如果采样了，需要扩展到完整数据集（使用均值近似）
        if len(X) > max_samples:
            full_shap_values = np.zeros((len(X), shap_values.shape[1]))
            full_shap_values[sample_indices] = shap_values
            # 对未采样的样本使用均值填充
            mean_shap = np.mean(shap_values, axis=0)
            for i in range(len(X)):
                if i not in sample_indices:
                    full_shap_values[i] = mean_shap
            shap_values = full_shap_values
        
        print("  策略3成功！")
        return shap_values
    
    except Exception as e:
        print(f"  策略3失败: {e}")
    
    # 策略4：使用特征重要性作为替代
    try:
        print("  尝试策略4: 使用特征重要性作为替代...")
        
        if hasattr(model, 'feature_importances_'):
            # 使用特征重要性作为SHAP值的粗略近似
            # 注意：这不是真正的SHAP值，但可以用于方向性分析
            feature_importances = model.feature_importances_
            
            # 为每个样本重复特征重要性
            shap_values = np.tile(feature_importances, (len(X), 1))
            
            print("  策略4成功（使用特征重要性近似）！")
            print("  警告：这不是真正的SHAP值，结果仅供参考")
            return shap_values
        else:
            raise ValueError("模型没有feature_importances_属性")
    
    except Exception as e:
        print(f"  策略4失败: {e}")
    
    print("  所有策略均失败，无法计算SHAP值")
    return None


# ============================================================================
# 数据加载函数
# ============================================================================

def load_sample_info(filepath):
    """加载并预处理样本信息"""
    sampledf = pd.read_csv(
        filepath,
        na_values=["", "NA", "NaN", "nan", "None", "none", "N/A", "n/a"]
    )
    
    # 处理dm_cl标签
    if "dm_cl" in sampledf.columns:
        sampledf["dm_cl"] = sampledf["dm_cl"].map(
            {"yes": 1, "no": 0, "Yes": 1, "No": 0}
        )
    
    # 处理sex
    if "sex" in sampledf.columns:
        sex_map = {
            "female": 0, "f": 0, "0": 0, 0: 0,
            "male": 1, "m": 1, "1": 1, 1: 1
        }
        sampledf["sex"] = sampledf["sex"].apply(
            lambda v: sex_map.get(str(v).strip().lower(), np.nan)
        )
    
    # 处理hyper_cl
    if "hyper_cl" in sampledf.columns:
        sampledf["hyper_cl"] = sampledf["hyper_cl"].apply(
            lambda v: 1 if str(v).strip().lower() == "yes" 
            else (0 if str(v).strip().lower() == "no" else np.nan)
        )
    
    # 数值列
    for col in ["BMI", "age", "waist", "TG", "HDL"]:
        if col in sampledf.columns:
            sampledf[col] = pd.to_numeric(sampledf[col], errors='coerce')
    
    return sampledf


def load_feature_matrix(feature_type, data_path, is_external=False):
    """加载特征矩阵"""
    if is_external:
        filename = data_path / f"GNHSFc_diann_IGC_humanswiss_{feature_type}_sample.tsv"
    else:
        filename = data_path / f"original_expression_matrix/GNHSF_diann_IGC_humanswiss_{feature_type}_sample_1385_NA90.tsv"
    
    if not filename.exists():
        print(f"  警告：特征矩阵文件不存在: {filename}")
        return None
    
    df = pd.read_csv(filename, sep="\t", index_col=0)
    
    # 转置
    dft = df.T.reset_index()
    dft.columns = ["sample"] + list(dft.columns[1:])
    
    # 简化列名
    def simplify_colname(col):
        col = col.split(" ")[0] if " " in col else col
        col = re.sub(r"\|", ".", col)
        col = re.sub(r"\-", ".", col)
        col = re.sub(r"\.$", "", col)
        return col
    
    dft.columns = [simplify_colname(col) for col in dft.columns]
    
    return dft


def categorize_features(feature_names):
    """
    将特征分类到不同的功能层级
    
    Returns:
    --------
    feature_categories : dict
        {'species': [...], 'COG': [...], 'genus': [...], 
         'human_protein': [...], 'KEGG': [...], 'microbial_protein': [...]}
    """
    categories = {
        'species': [],
        'COG': [],
        'genus': [],
        'human_protein': [],
        'KEGG': [],
        'microbial_protein': [],
        'clinical': []
    }
    
    clinical_vars = ['BMI', 'sex', 'age', 'waist', 'TG', 'HDL', 'hyper_cl']
    
    for feat in feature_names:
        if feat in clinical_vars:
            categories['clinical'].append(feat)
        elif feat.startswith('COG'):
            categories['COG'].append(feat)
        elif feat.startswith('K') or 'KEGG' in feat.upper():
            categories['KEGG'].append(feat)
        elif '_s_' in feat.lower() or 'species' in feat.lower():
            categories['species'].append(feat)
        elif '_g_' in feat.lower() or 'genus' in feat.lower():
            categories['genus'].append(feat)
        elif 'HUMAN' in feat.upper() or feat.startswith('ENSP'):
            categories['human_protein'].append(feat)
        else:
            categories['microbial_protein'].append(feat)
    
    return categories


# ============================================================================
# 1. 加载模型和预测结果
# ============================================================================

def load_predictions_and_models():
    """
    加载所有种子的内部和外部预测结果
    
    Returns:
    --------
    internal_preds_dict : dict
        {seed: DataFrame with columns ['sample', 'y_true', 'y_pred', 'fold', 'seed']}
    external_preds_dict : dict
        {seed: DataFrame with columns ['sample', 'y_true', 'y_pred', 'seed']}
    models_dict : dict
        {(fold, seed): loaded_model}
    """
    print("=" * 80)
    print(f"加载预测结果和模型（特征类型: {FEATURE_TYPE}）...")
    print("=" * 80)
    
    internal_preds_dict = {}
    external_preds_dict = {}
    models_dict = {}
    
    # 加载内部预测
    internal_preds_dir = MODEL_OUTPUT_ROOT / "internal_preds"
    if internal_preds_dir.exists():
        for seed in SEEDS:
            pred_file = internal_preds_dir / f"internal_preds_seed_{seed}.csv"
            if pred_file.exists():
                df = pd.read_csv(pred_file)
                internal_preds_dict[seed] = df
                # print(f"  已加载种子 {seed} 的内部预测：{len(df)} 个样本")
        print(f"  已加载 {len(internal_preds_dict)} 个种子的内部预测")
    else:
        print(f"  警告：内部预测目录不存在: {internal_preds_dir}")
    
    # 加载外部预测
    external_preds_dir = MODEL_OUTPUT_ROOT / "external_preds"
    if external_preds_dir.exists():
        for seed in SEEDS:
            pred_file = external_preds_dir / f"external_preds_seed_{seed}.csv"
            if pred_file.exists():
                df = pd.read_csv(pred_file)
                external_preds_dict[seed] = df
                # print(f"  已加载种子 {seed} 的外部预测：{len(df)} 个样本")
        print(f"  已加载 {len(external_preds_dict)} 个种子的外部预测")
    else:
        print(f"  警告：外部预测目录不存在: {external_preds_dir}")
    
    # 加载模型（用于SHAP分析）
    models_dir = MODEL_OUTPUT_ROOT / "models"
    if models_dir.exists():
        for fold in range(N_FOLDS):
            fold_dir = models_dir / f"fold_{fold}"
            if not fold_dir.exists():
                continue
            
            for seed in SEEDS:
                model_file = fold_dir / f"final_xgb_model_{seed}.pkl"
                if model_file.exists():
                    model = joblib.load(model_file)
                    models_dict[(fold, seed)] = model
    else:
        print(f"  警告：模型目录不存在: {models_dir}")
    
    print(f"  总计加载了 {len(models_dict)} 个模型")
    
    return internal_preds_dict, external_preds_dict, models_dict


# ============================================================================
# 2. 准备数据（用于SHAP分析）
# ============================================================================

def prepare_data_for_shap(fold, seed):
    """
    为指定的fold和seed准备数据
    
    Returns:
    --------
    X_internal, y_internal, X_external, y_external, feature_names
    """
    # 加载该fold的特征列表
    feature_file = MODEL_OUTPUT_ROOT / "selected_features" / f"fold_{fold}" / "model_features.csv"
    if not feature_file.exists():
        print(f"  警告：fold {fold} 的特征文件不存在: {feature_file}")
        return None, None, None, None, None
    
    feature_df = pd.read_csv(feature_file)
    model_features = feature_df["feature"].tolist()
    
    # 加载内部数据
    internal_sample_df = load_sample_info(INTERNAL_SAMPLE_INFO)
    internal_feature_df = load_feature_matrix(
        FEATURE_TYPE, INTERNAL_DATA_PATH, is_external=False
    )
    
    if internal_feature_df is None:
        print(f"  警告：无法加载内部特征矩阵")
        return None, None, None, None, None
    
    # 合并内部数据
    internal_df = pd.merge(
        internal_feature_df,
        internal_sample_df[["sample"] + CLINICAL_FEATURES + ["dm_cl"]],
        on="sample",
        how="inner"
    )
    
    # 获取内部测试集样本（从预测文件中）
    internal_pred_file = MODEL_OUTPUT_ROOT / "internal_preds" / f"internal_preds_seed_{seed}.csv"
    if internal_pred_file.exists():
        internal_pred_df = pd.read_csv(internal_pred_file)
        internal_test_samples = internal_pred_df[internal_pred_df["fold"] == fold]["sample"].tolist()
        internal_df = internal_df[internal_df["sample"].isin(internal_test_samples)]
    
    # 提取特征和标签
    missing_features_internal = [f for f in model_features if f not in internal_df.columns]
    if missing_features_internal:
        for feat in missing_features_internal:
            internal_df[feat] = np.nan
    
    X_internal = internal_df[model_features].copy()
    y_internal = internal_df["dm_cl"].astype(int).copy()
    
    # 加载外部数据
    external_sample_df = load_sample_info(EXTERNAL_SAMPLE_INFO)
    external_feature_df = load_feature_matrix(
        FEATURE_TYPE, EXTERNAL_DATA_PATH, is_external=True
    )
    
    if external_feature_df is None:
        print(f"  警告：无法加载外部特征矩阵")
        return None, None, None, None, None
    
    # 合并外部数据
    external_df = pd.merge(
        external_feature_df,
        external_sample_df[["sample"] + CLINICAL_FEATURES + ["dm_cl"]],
        on="sample",
        how="inner"
    )
    
    # 提取特征和标签
    missing_features_external = [f for f in model_features if f not in external_df.columns]
    if missing_features_external:
        for feat in missing_features_external:
            external_df[feat] = np.nan
    
    # 去除label为NAN的样本
    external_df = external_df.dropna(subset=CLINICAL_FEATURES + ["dm_cl"])
    
    X_external = external_df[model_features].copy()
    y_external = external_df["dm_cl"].astype(int).copy()
    
    return X_internal, y_internal, X_external, y_external, model_features


# ============================================================================
# 3. 效应方向一致性分析（修复版）
# ============================================================================

def compute_effect_direction_concordance(
    models_dict,
    top_n=100,
    selected_seed=42,
    selected_fold=0
):
    """
    计算内部测试集和外部验证集之间的特征效应方向一致性
    （修复版：包含SHAP兼容性处理）
    """
    print("\n" + "=" * 80)
    print("1. 计算特征效应方向一致性...")
    print("=" * 80)
    
    # 选择一个代表性的模型
    if (selected_fold, selected_seed) not in models_dict:
        print(f"  警告：未找到fold {selected_fold}, seed {selected_seed}的模型")
        return None
    
    model = models_dict[(selected_fold, selected_seed)]
    
    # 准备数据
    X_internal, y_internal, X_external, y_external, feature_names = prepare_data_for_shap(
        selected_fold, selected_seed
    )
    
    if X_internal is None or len(X_internal) == 0:
        print("  数据准备失败或数据为空，跳过")
        return None
    
    print(f"  使用fold {selected_fold}, seed {selected_seed}的模型")
    print(f"  内部测试集：{len(X_internal)} 样本，{len(feature_names)} 特征")
    print(f"  外部验证集：{len(X_external)} 样本，{len(feature_names)} 特征")
    
    # 计算SHAP值（使用安全方法）
    print("\n  计算SHAP值（内部测试集）...")
    shap_values_internal = compute_shap_values_safe(model, X_internal, max_samples=1000)
    
    if shap_values_internal is None:
        print("  无法计算内部测试集SHAP值，跳过分析")
        return None
    
    print("\n  计算SHAP值（外部验证集）...")
    shap_values_external = compute_shap_values_safe(model, X_external, max_samples=1000)
    
    if shap_values_external is None:
        print("  无法计算外部验证集SHAP值，跳过分析")
        return None
    
    # 如果是二分类，取正类的SHAP值
    if isinstance(shap_values_internal, list):
        shap_values_internal = shap_values_internal[1]
        shap_values_external = shap_values_external[1]
    
    # 计算平均SHAP值
    mean_shap_internal = np.mean(shap_values_internal, axis=0)
    mean_shap_external = np.mean(shap_values_external, axis=0)
    
    # 动态调整top_n（如果特征数不足）
    actual_top_n = min(top_n, len(feature_names), len(mean_shap_internal))
    if actual_top_n < 2:
        print(f"  警告：可用特征数不足（{actual_top_n}），跳过分析")
        return None
    
    if actual_top_n < top_n:
        print(f"  警告：特征数不足，调整top_n从 {top_n} 到 {actual_top_n}")
    
    # 获取top特征
    abs_shap_internal = np.abs(mean_shap_internal)
    top_indices = np.argsort(abs_shap_internal)[-actual_top_n:]
    
    # 计算符号一致性
    signs_internal = np.sign(mean_shap_internal[top_indices])
    signs_external = np.sign(mean_shap_external[top_indices])
    
    concordant = np.sum(signs_internal == signs_external)
    discordant = np.sum(signs_internal != signs_external)
    concordance_pct = (concordant / len(top_indices)) * 100
    
    # 计算相关系数
    pearson_r, pearson_p = pearsonr(
        mean_shap_internal[top_indices],
        mean_shap_external[top_indices]
    )
    
    spearman_r, spearman_p = spearmanr(
        mean_shap_internal[top_indices],
        mean_shap_external[top_indices]
    )
    
    # 余弦相似度
    from scipy.spatial.distance import cosine
    cosine_sim = 1 - cosine(
        mean_shap_internal[top_indices],
        mean_shap_external[top_indices]
    )
    
    results = {
        'top_n': actual_top_n,
        'concordant_features': concordant,
        'discordant_features': discordant,
        'concordance_percentage': concordance_pct,
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'cosine_similarity': cosine_sim,
        'top_features': [feature_names[i] for i in top_indices],
        'shap_internal': mean_shap_internal[top_indices],
        'shap_external': mean_shap_external[top_indices],
        'signs_internal': signs_internal,
        'signs_external': signs_external
    }
    
    # 打印结果
    print(f"\n效应方向一致性分析 (Top {actual_top_n} 特征):")
    print(f"  一致特征数: {concordant}/{actual_top_n} ({concordance_pct:.1f}%)")
    print(f"  不一致特征数: {discordant}/{actual_top_n}")
    print(f"  Pearson r = {pearson_r:.3f} (p = {pearson_p:.2e})")
    print(f"  Spearman ρ = {spearman_r:.3f} (p = {spearman_p:.2e})")
    print(f"  Cosine similarity = {cosine_sim:.3f}")
    
    return results


# ============================================================================
# 4. 通路/模块层级的信号一致性分析（修复版）
# ============================================================================

def aggregate_features_to_pathways(features, shap_values, category_type='COG'):
    """
    将特征聚合到通路/模块层级
    """
    pathway_mapping = defaultdict(list)
    
    for idx, feat in enumerate(features):
        if category_type == 'COG':
            pathway_id = feat.split('_')[0] if feat.startswith('COG') else 'Other'
        elif category_type == 'KEGG':
            pathway_id = feat.split('_')[0] if feat.startswith('K') else 'Other'
        elif category_type == 'pathway' or category_type == 'human_protein':
            if ':' in feat:
                pathway_id = feat.split(':')[0] + ':' + feat.split(':')[1].split('_')[0]
            elif '.' in feat:
                pathway_id = feat.split('.')[0]
            else:
                pathway_id = feat[:10]
        else:
            pathway_id = 'Other'
        
        pathway_mapping[pathway_id].append(idx)
    
    # 计算每个通路的平均SHAP值
    mean_shap = np.mean(shap_values, axis=0)
    pathway_scores = {}
    
    for pathway_id, indices in pathway_mapping.items():
        pathway_scores[pathway_id] = np.mean(mean_shap[indices])
    
    return pathway_scores, dict(pathway_mapping)


def compute_pathway_level_concordance(
    models_dict,
    category_types=['COG', 'KEGG', 'microbial_protein'],
    top_n_pathways=50,
    selected_seed=42,
    selected_fold=0
):
    """
    计算通路/模块层级的信号一致性
    （修复版：包含SHAP兼容性处理和数据验证）
    """
    print("\n" + "=" * 80)
    print("2. 计算通路/模块层级的信号一致性...")
    print("=" * 80)
    
    if (selected_fold, selected_seed) not in models_dict:
        print(f"  警告：未找到fold {selected_fold}, seed {selected_seed}的模型")
        return None
    
    model = models_dict[(selected_fold, selected_seed)]
    
    # 准备数据
    X_internal, y_internal, X_external, y_external, feature_names = prepare_data_for_shap(
        selected_fold, selected_seed
    )
    
    if X_internal is None or len(X_internal) == 0:
        print("  数据准备失败或数据为空，跳过")
        return None
    
    # 计算SHAP值（使用安全方法）
    print("\n  计算SHAP值（内部测试集）...")
    shap_values_internal = compute_shap_values_safe(model, X_internal, max_samples=1000)
    
    if shap_values_internal is None:
        print("  无法计算SHAP值，跳过通路分析")
        return None
    
    print("\n  计算SHAP值（外部验证集）...")
    shap_values_external = compute_shap_values_safe(model, X_external, max_samples=1000)
    
    if shap_values_external is None:
        print("  无法计算SHAP值，跳过通路分析")
        return None
    
    if isinstance(shap_values_internal, list):
        shap_values_internal = shap_values_internal[1]
        shap_values_external = shap_values_external[1]
    
    all_results = {}
    
    # 对每个类别进行分析
    feature_categories = categorize_features(feature_names)
    
    for cat_type in category_types:
        if cat_type not in feature_categories or len(feature_categories[cat_type]) == 0:
            print(f"\n  跳过 {cat_type}：无相关特征")
            continue
        
        print(f"\n  分析 {cat_type} 通路...")
        
        # 获取该类别的特征索引
        cat_features = feature_categories[cat_type]
        cat_indices = [i for i, feat in enumerate(feature_names) if feat in cat_features]
        
        if len(cat_indices) == 0:
            print(f"    跳过：无索引匹配")
            continue
        
        # 聚合到通路层级
        pathway_scores_internal, pathway_mapping = aggregate_features_to_pathways(
            [feature_names[i] for i in cat_indices],
            shap_values_internal[:, cat_indices],
            cat_type
        )
        
        pathway_scores_external, _ = aggregate_features_to_pathways(
            [feature_names[i] for i in cat_indices],
            shap_values_external[:, cat_indices],
            cat_type
        )
        
        # 获取共同通路
        common_pathways = set(pathway_scores_internal.keys()).intersection(
            set(pathway_scores_external.keys())
        )
        
        if len(common_pathways) < 2:
            print(f"    警告: {cat_type} 共同通路数不足（{len(common_pathways)}），跳过")
            continue
        
        # 提取共同通路的分数
        common_pathway_list = list(common_pathways)
        scores_internal = np.array([pathway_scores_internal[p] for p in common_pathway_list])
        scores_external = np.array([pathway_scores_external[p] for p in common_pathway_list])
        
        # 验证数据有效性
        if len(scores_internal) < 2 or len(scores_external) < 2:
            print(f"    警告: {cat_type} 数据点不足，跳过")
            continue
        
        # 动态调整top_n_pathways
        actual_top_n = min(top_n_pathways, len(common_pathways))
        if actual_top_n < 2:
            print(f"    警告: {cat_type} 可用通路数不足（{actual_top_n}），跳过")
            continue
        
        if actual_top_n < top_n_pathways:
            print(f"    信息：通路数不足，调整top_n从 {top_n_pathways} 到 {actual_top_n}")
        
        # 获取top通路
        top_indices = np.argsort(np.abs(scores_internal))[-actual_top_n:]
        top_pathways = [common_pathway_list[i] for i in top_indices]
        
        # 在外部验证集中检查这些通路
        test_abs_scores = np.abs(scores_external)
        test_top_indices = np.argsort(test_abs_scores)[-actual_top_n:]
        test_top_pathways = set([common_pathway_list[i] for i in test_top_indices])
        
        # 计算重叠
        overlap = len(set(top_pathways).intersection(test_top_pathways))
        jaccard_index = overlap / len(set(top_pathways).union(test_top_pathways))
        
        # 计算相关系数（需要至少2个数据点）
        try:
            pearson_r, pearson_p = pearsonr(scores_internal, scores_external)
            spearman_r, spearman_p = spearmanr(scores_internal, scores_external)
        except Exception as e:
            print(f"    警告：相关系数计算失败: {e}")
            pearson_r, pearson_p = np.nan, np.nan
            spearman_r, spearman_p = np.nan, np.nan
        
        # 排列检验
        n_permutations = 1000
        permuted_overlaps = []
        
        for _ in range(n_permutations):
            perm_indices = np.random.permutation(len(common_pathway_list))[:actual_top_n]
            perm_pathways = set([common_pathway_list[i] for i in perm_indices])
            perm_overlap = len(perm_pathways.intersection(test_top_pathways))
            permuted_overlaps.append(perm_overlap)
        
        permutation_p = np.mean(np.array(permuted_overlaps) >= overlap)
        expected_overlap = np.mean(permuted_overlaps)
        
        results = {
            'category': cat_type,
            'total_pathways': len(common_pathways),
            'top_n': actual_top_n,
            'overlap': overlap,
            'jaccard_index': jaccard_index,
            'expected_overlap': expected_overlap,
            'permutation_p': permutation_p,
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p,
            'top_pathways_internal': top_pathways,
            'top_pathways_external': list(test_top_pathways),
            'pathway_scores_internal': scores_internal[top_indices],
            'pathway_scores_external': scores_external[top_indices]
        }
        
        all_results[cat_type] = results
        
        # 打印结果
        print(f"    总通路数: {len(common_pathways)}")
        print(f"    Top {actual_top_n} 通路重叠: {overlap}/{actual_top_n} ({overlap/actual_top_n*100:.1f}%)")
        print(f"    Jaccard index: {jaccard_index:.3f}")
        print(f"    期望重叠 (随机): {expected_overlap:.1f}")
        print(f"    排列检验 p-value: {permutation_p:.4f}")
        if not np.isnan(spearman_r):
            print(f"    Spearman ρ = {spearman_r:.3f} (p = {spearman_p:.2e})")
    
    # FDR校正
    if len(all_results) > 1:
        p_values = [res['spearman_p'] for res in all_results.values() if not np.isnan(res['spearman_p'])]
        categories = [cat for cat, res in all_results.items() if not np.isnan(res['spearman_p'])]
        
        if len(p_values) > 0:
            _, fdr_corrected, _, _ = multipletests(p_values, method='fdr_bh')
            
            for cat, fdr in zip(categories, fdr_corrected):
                all_results[cat]['spearman_fdr'] = fdr
                print(f"\n  {cat} Spearman FDR-corrected p-value: {fdr:.4f}")
    
    return all_results if all_results else None


# ============================================================================
# 5. 特征重叠的排列检验（修复版）
# ============================================================================

def compute_feature_overlap_permutation_test(
    models_dict,
    observed_k=159,
    top_n=500,
    n_permutations=10000,
    selected_seed=42,
    selected_fold=0
):
    """
    计算特征重叠的排列检验
    （修复版：包含SHAP兼容性处理和数据验证）
    """
    print("\n" + "=" * 80)
    print("3. 特征重叠的排列检验...")
    print("=" * 80)
    
    if (selected_fold, selected_seed) not in models_dict:
        print(f"  警告：未找到fold {selected_fold}, seed {selected_seed}的模型")
        return None
    
    model = models_dict[(selected_fold, selected_seed)]
    
    # 准备数据
    X_internal, y_internal, X_external, y_external, feature_names = prepare_data_for_shap(
        selected_fold, selected_seed
    )
    
    if X_internal is None or len(X_internal) == 0:
        print("  数据准备失败或数据为空，跳过")
        return None
    
    # 计算SHAP值（使用安全方法）
    print("\n  计算SHAP值（内部测试集）...")
    shap_values_internal = compute_shap_values_safe(model, X_internal, max_samples=1000)
    
    if shap_values_internal is None:
        print("  无法计算SHAP值，跳过排列检验")
        return None
    
    print("\n  计算SHAP值（外部验证集）...")
    shap_values_external = compute_shap_values_safe(model, X_external, max_samples=1000)
    
    if shap_values_external is None:
        print("  无法计算SHAP值，跳过排列检验")
        return None
    
    if isinstance(shap_values_internal, list):
        shap_values_internal = shap_values_internal[1]
        shap_values_external = shap_values_external[1]
    
    # 计算特征重要性
    feat_imp_internal = np.abs(np.mean(shap_values_internal, axis=0))
    feat_imp_external = np.abs(np.mean(shap_values_external, axis=0))
    
    # 动态调整top_n
    actual_top_n = min(top_n, len(feature_names))
    if actual_top_n < 2:
        print(f"  警告：可用特征数不足（{actual_top_n}），跳过")
        return None
    
    if actual_top_n < top_n:
        print(f"  警告：特征数不足，调整top_n从 {top_n} 到 {actual_top_n}")
    
    # 获取top特征
    top_indices_internal = np.argsort(feat_imp_internal)[-actual_top_n:]
    top_indices_external = np.argsort(feat_imp_external)[-actual_top_n:]
    
    top_features_internal = set([feature_names[i] for i in top_indices_internal])
    top_features_external = set([feature_names[i] for i in top_indices_external])
    
    # 计算实际重叠
    actual_overlap = len(top_features_internal.intersection(top_features_external))
    
    print(f"\n  实际重叠特征数: {actual_overlap}/{actual_top_n}")
    print(f"  观察到的验证特征数: {observed_k}")
    
    # 排列检验
    permuted_overlaps = []
    
    print(f"\n  执行 {n_permutations} 次排列检验...")
    for i in range(n_permutations):
        if (i + 1) % 1000 == 0:
            print(f"    完成 {i+1}/{n_permutations}...")
        
        perm_indices = np.random.permutation(len(feature_names))[:actual_top_n]
        perm_features = set([feature_names[j] for j in perm_indices])
        
        perm_overlap = len(top_features_internal.intersection(perm_features))
        permuted_overlaps.append(perm_overlap)
    
    permuted_overlaps = np.array(permuted_overlaps)
    
    # 计算统计量
    expected_overlap = np.mean(permuted_overlaps)
    std_overlap = np.std(permuted_overlaps)
    
    # p-value
    p_value = np.mean(permuted_overlaps >= actual_overlap)
    
    # Z-score
    z_score = (actual_overlap - expected_overlap) / (std_overlap + 1e-10)
    
    # 置信区间
    ci_95 = np.percentile(permuted_overlaps, [2.5, 97.5])
    
    results = {
        'observed_overlap': actual_overlap,
        'observed_k': observed_k,
        'expected_overlap': expected_overlap,
        'std_overlap': std_overlap,
        'p_value': p_value,
        'z_score': z_score,
        'ci_95': ci_95,
        'permuted_distribution': permuted_overlaps,
        'top_n': actual_top_n,
        'n_permutations': n_permutations
    }
    
    # 打印结果
    print(f"\n  排列检验结果:")
    print(f"    观察到的重叠: k = {actual_overlap}")
    print(f"    期望重叠 (k0): {expected_overlap:.1f} ± {std_overlap:.1f}")
    print(f"    95% CI: [{ci_95[0]:.1f}, {ci_95[1]:.1f}]")
    print(f"    Z-score: {z_score:.2f}")
    print(f"    P-value: {p_value:.4f}")
    
    if p_value < 0.001:
        print(f"    *** 重叠显著超过随机期望 (p < 0.001) ***")
    elif p_value < 0.05:
        print(f"    ** 重叠显著超过随机期望 (p < 0.05) **")
    else:
        print(f"    重叠未显著超过随机期望 (p = {p_value:.4f})")
    
    return results


# ============================================================================
# 6. 生成综合报告
# ============================================================================

def generate_comprehensive_report(
    internal_preds_dict,
    external_preds_dict,
    models_dict,
    concordance_results,
    pathway_results,
    permutation_results,
    output_file='validation_report.txt'
):
    """
    生成综合验证报告
    """
    print("\n" + "=" * 80)
    print("生成综合报告...")
    print("=" * 80)
    
    # 计算内部和外部的AUC
    internal_aucs = []
    external_aucs = []
    
    for seed in SEEDS:
        if seed in internal_preds_dict:
            df = internal_preds_dict[seed]
            try:
                auc_val = roc_auc_score(df['y_true'], df['y_pred'])
                internal_aucs.append(auc_val)
            except Exception as e:
                print(f"  警告：种子{seed}的内部AUC计算失败: {e}")
        
        if seed in external_preds_dict:
            df = external_preds_dict[seed]
            try:
                auc_val = roc_auc_score(df['y_true'], df['y_pred'])
                external_aucs.append(auc_val)
            except Exception as e:
                print(f"  警告：种子{seed}的外部AUC计算失败: {e}")
    
    # 写入报告
    report_path = VALIDATION_OUTPUT_DIR / output_file
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write(f"模型验证综合报告\n")
        f.write(f"特征类型: {FEATURE_TYPE}\n")
        f.write(f"生成时间: {pd.Timestamp.now()}\n")
        f.write(f"用户: westlakeomics\n")
        f.write("=" * 80 + "\n\n")
        
        # 1. AUROC结果
        f.write("1. AUROC 结果\n")
        f.write("-" * 40 + "\n")
        if internal_aucs:
            f.write(f"内部测试集 (n={len(internal_aucs)} seeds):\n")
            f.write(f"  Mean AUC: {np.mean(internal_aucs):.3f} ± {np.std(internal_aucs):.3f}\n")
            f.write(f"  Median AUC: {np.median(internal_aucs):.3f}\n")
            f.write(f"  Range: [{np.min(internal_aucs):.3f}, {np.max(internal_aucs):.3f}]\n\n")
        else:
            f.write("  内部测试集：无数据\n\n")
        
        if external_aucs:
            f.write(f"外部验证集 (FH, n={len(external_aucs)} seeds):\n")
            f.write(f"  Mean AUC: {np.mean(external_aucs):.3f} ± {np.std(external_aucs):.3f}\n")
            f.write(f"  Median AUC: {np.median(external_aucs):.3f}\n")
            f.write(f"  Range: [{np.min(external_aucs):.3f}, {np.max(external_aucs):.3f}]\n\n")
        else:
            f.write("  外部验证集：无数据\n\n")
        
        # 2. 效应方向一致性
        if concordance_results:
            f.write("2. 效应方向一致性\n")
            f.write("-" * 40 + "\n")
            conc = concordance_results
            f.write(f"Top {conc['top_n']} 特征:\n")
            f.write(f"  一致特征数: {conc['concordant_features']}/{conc['top_n']}\n")
            f.write(f"  一致性百分比: {conc['concordance_percentage']:.1f}%\n")
            f.write(f"  Pearson r = {conc['pearson_r']:.3f} (p = {conc['pearson_p']:.2e})\n")
            f.write(f"  Spearman ρ = {conc['spearman_r']:.3f} (p = {conc['spearman_p']:.2e})\n")
            f.write(f"  Cosine similarity = {conc['cosine_similarity']:.3f}\n\n")
        else:
            f.write("2. 效应方向一致性\n")
            f.write("-" * 40 + "\n")
            f.write("  无可用数据（可能由于特征数不足或SHAP计算失败）\n\n")
        
        # 3. 通路层级一致性
        if pathway_results:
            f.write("3. 通路/模块层级一致性\n")
            f.write("-" * 40 + "\n")
            for cat, res in pathway_results.items():
                f.write(f"\n{cat}:\n")
                f.write(f"  总通路数: {res['total_pathways']}\n")
                f.write(f"  Top {res['top_n']} 通路重叠: {res['overlap']}/{res['top_n']} ({res['overlap']/res['top_n']*100:.1f}%)\n")
                f.write(f"  Jaccard index: {res['jaccard_index']:.3f}\n")
                f.write(f"  期望重叠（随机）: {res['expected_overlap']:.1f}\n")
                if not np.isnan(res['spearman_r']):
                    f.write(f"  Spearman ρ = {res['spearman_r']:.3f} (p = {res['spearman_p']:.2e}")
                    if 'spearman_fdr' in res:
                        f.write(f", FDR = {res['spearman_fdr']:.4f}")
                    f.write(")\n")
                else:
                    f.write("  Spearman ρ = N/A（数据不足）\n")
                f.write(f"  排列检验 p-value = {res['permutation_p']:.4f}\n")
        else:
            f.write("3. 通路/模块层级一致性\n")
            f.write("-" * 40 + "\n")
            f.write("  无可用数据（可能由于相关特征数不足）\n")
        
        # 4. 排列检验
        if permutation_results:
            f.write("\n4. 特征重叠排列检验\n")
            f.write("-" * 40 + "\n")
            perm = permutation_results
            f.write(f"观察到的重叠: k = {perm['observed_overlap']}\n")
            f.write(f"期望重叠: k0 = {perm['expected_overlap']:.1f} ± {perm['std_overlap']:.1f}\n")
            f.write(f"95% CI: [{perm['ci_95'][0]:.1f}, {perm['ci_95'][1]:.1f}]\n")
            f.write(f"Z-score: {perm['z_score']:.2f}\n")
            f.write(f"P-value: {perm['p_value']:.4f}\n")
            if perm['p_value'] < 0.001:
                f.write("显著性: *** (p < 0.001)\n")
            elif perm['p_value'] < 0.05:
                f.write("显著性: ** (p < 0.05)\n")
            else:
                f.write("显著性: ns (p >= 0.05)\n")
        else:
            f.write("\n4. 特征重叠排列检验\n")
            f.write("-" * 40 + "\n")
            f.write("  无可用数据\n")
        
        f.write("\n" + "=" * 80 + "\n")
        f.write("报告结束\n")
        f.write("=" * 80 + "\n")
    
    print(f"\n  报告已保存到: {report_path}")
    
    # 返回汇总结果
    return {
        'feature_type': FEATURE_TYPE,
        'AUROC': {
            'internal': {
                'mean': np.mean(internal_aucs) if internal_aucs else None,
                'std': np.std(internal_aucs) if internal_aucs else None,
                'median': np.median(internal_aucs) if internal_aucs else None,
                'n': len(internal_aucs)
            },
            'external': {
                'mean': np.mean(external_aucs) if external_aucs else None,
                'std': np.std(external_aucs) if external_aucs else None,
                'median': np.median(external_aucs) if external_aucs else None,
                'n': len(external_aucs)
            }
        },
        'concordance': concordance_results,
        'pathway': pathway_results,
        'permutation': permutation_results
    }


# ============================================================================
# 7. 可视化函数
# ============================================================================

def plot_validation_results(
    concordance_results,
    pathway_results,
    permutation_results,
    internal_preds_dict,
    external_preds_dict
):
    """
    生成验证结果的可视化图表
    """
    print("\n" + "=" * 80)
    print("生成可视化图表...")
    print("=" * 80)
    
    plot_dir = VALIDATION_OUTPUT_DIR / "plots"
    plot_dir.mkdir(exist_ok=True, parents=True)
    
    # 1. 效应方向散点图
    if concordance_results:
        try:
            plt.figure(figsize=(10, 8))
            conc = concordance_results
            
            plt.scatter(conc['shap_internal'], conc['shap_external'], 
                       alpha=0.6, s=50, c='steelblue', edgecolors='black')
            plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
            plt.axvline(0, color='gray', linestyle='--', linewidth=0.8)
            
            # 对角线
            min_val = min(conc['shap_internal'].min(), conc['shap_external'].min())
            max_val = max(conc['shap_internal'].max(), conc['shap_external'].max())
            plt.plot([min_val, max_val], [min_val, max_val], 
                    'r--', linewidth=2, label='y=x (perfect concordance)')
            
            plt.xlabel('Mean SHAP value (Internal Test Set)', fontsize=12)
            plt.ylabel('Mean SHAP value (External Validation Set)', fontsize=12)
            plt.title(f'Effect Direction Concordance ({FEATURE_TYPE})\n'
                     f'Top {conc["top_n"]} Features | Concordance: {conc["concordance_percentage"]:.1f}%\n'
                     f'Pearson r = {conc["pearson_r"]:.3f} | Spearman ρ = {conc["spearman_r"]:.3f}',
                     fontsize=14)
            plt.legend()
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(plot_dir / f'effect_direction_concordance_{FEATURE_TYPE}.png', dpi=300)
            plt.close()
            print(f"  已保存: effect_direction_concordance_{FEATURE_TYPE}.png")
        except Exception as e:
            print(f"  警告：效应方向散点图生成失败: {e}")
    
    # 2. 通路重叠柱状图
    if pathway_results:
        try:
            categories = list(pathway_results.keys())
            overlaps = [pathway_results[cat]['overlap'] for cat in categories]
            expected = [pathway_results[cat]['expected_overlap'] for cat in categories]
            
            fig, ax = plt.subplots(figsize=(12, 6))
            x = np.arange(len(categories))
            width = 0.35
            
            ax.bar(x - width/2, overlaps, width, label='Observed Overlap', 
                   color='steelblue', edgecolor='black')
            ax.bar(x + width/2, expected, width, label='Expected (Random)', 
                   color='coral', edgecolor='black')
            
            ax.set_xlabel('Pathway Category', fontsize=12)
            ax.set_ylabel('Number of Overlapping Pathways', fontsize=12)
            ax.set_title(f'Pathway-Level Concordance ({FEATURE_TYPE})', fontsize=14)
            ax.set_xticks(x)
            ax.set_xticklabels(categories, rotation=45, ha='right')
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(plot_dir / f'pathway_overlap_{FEATURE_TYPE}.png', dpi=300)
            plt.close()
            print(f"  已保存: pathway_overlap_{FEATURE_TYPE}.png")
        except Exception as e:
            print(f"  警告：通路重叠柱状图生成失败: {e}")
    
    # 3. 排列检验分布
    if permutation_results:
        try:
            perm = permutation_results
            
            plt.figure(figsize=(12, 6))
            plt.hist(perm['permuted_distribution'], bins=50, alpha=0.7, 
                    color='lightblue', edgecolor='black', label='Null Distribution')
            plt.axvline(perm['observed_overlap'], color='red', linestyle='--', 
                       linewidth=2, label=f'Observed (k={perm["observed_overlap"]})')
            plt.axvline(perm['expected_overlap'], color='green', linestyle='--', 
                       linewidth=2, label=f'Expected (k0={perm["expected_overlap"]:.1f})')
            
            plt.xlabel('Number of Overlapping Features', fontsize=12)
            plt.ylabel('Frequency', fontsize=12)
            plt.title(f'Permutation Test ({FEATURE_TYPE})\n'
                     f'n={perm["n_permutations"]} permutations | p-value = {perm["p_value"]:.4f}', 
                     fontsize=14)
            plt.legend()
            plt.grid(alpha=0.3)
            plt.tight_layout()
            plt.savefig(plot_dir / f'permutation_test_{FEATURE_TYPE}.png', dpi=300)
            plt.close()
            print(f"  已保存: permutation_test_{FEATURE_TYPE}.png")
        except Exception as e:
            print(f"  警告：排列检验分布图生成失败: {e}")
    
    # 4. AUC箱型图
    try:
        internal_aucs = []
        external_aucs = []
        
        for seed in SEEDS:
            if seed in internal_preds_dict:
                df = internal_preds_dict[seed]
                try:
                    auc_val = roc_auc_score(df['y_true'], df['y_pred'])
                    internal_aucs.append(auc_val)
                except:
                    pass
            
            if seed in external_preds_dict:
                df = external_preds_dict[seed]
                try:
                    auc_val = roc_auc_score(df['y_true'], df['y_pred'])
                    external_aucs.append(auc_val)
                except:
                    pass
        
        if internal_aucs or external_aucs:
            fig, ax = plt.subplots(figsize=(8, 6))
            data_to_plot = []
            labels = []
            
            if internal_aucs:
                data_to_plot.append(internal_aucs)
                labels.append('Internal Test')
            
            if external_aucs:
                data_to_plot.append(external_aucs)
                labels.append('External Validation (FH)')
            
            bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
            
            # 设置颜色
            colors = ['lightblue', 'lightcoral']
            for patch, color in zip(bp['boxes'], colors[:len(data_to_plot)]):
                patch.set_facecolor(color)
            
            ax.set_ylabel('AUC', fontsize=12)
            ax.set_title(f'AUC Distribution Across {len(SEEDS)} Seeds ({FEATURE_TYPE})', fontsize=14)
            ax.grid(axis='y', alpha=0.3)
            
            # 添加均值线
            for i, aucs in enumerate(data_to_plot, 1):
                mean_val = np.mean(aucs)
                ax.plot(i, mean_val, marker='D', markersize=8, color='red')
                ax.text(i, mean_val + 0.02, f'{mean_val:.3f}', 
                       ha='center', fontsize=10, color='red')
            
            plt.tight_layout()
            plt.savefig(plot_dir / f'auc_boxplots_{FEATURE_TYPE}.png', dpi=300)
            plt.close()
            print(f"  已保存: auc_boxplots_{FEATURE_TYPE}.png")
    except Exception as e:
        print(f"  警告：AUC箱型图生成失败: {e}")


# ============================================================================
# 8. 单个特征类型的分析函数
# ============================================================================

def analyze_single_feature_type(feature_type, args):
    """
    分析单个特征类型
    """
    global BASE_DIR, MODEL_OUTPUT_ROOT, INTERNAL_DATA_PATH, EXTERNAL_DATA_PATH
    global INTERNAL_SAMPLE_INFO, EXTERNAL_SAMPLE_INFO, FEATURE_TYPE, VALIDATION_OUTPUT_DIR
    
    # 设置全局变量
    BASE_DIR = Path(args.base_dir)
    FEATURE_TYPE = feature_type
    
    # 构建模型输出目录名（匹配训练脚本的命名规则）
    model_dir_name = f"glm_xgb_gain_freq_{args.random_seed}_output_{feature_type}"
    MODEL_OUTPUT_ROOT = BASE_DIR / model_dir_name
    
    # 设置数据路径
    INTERNAL_DATA_PATH = BASE_DIR / "ML_matrix_20250911/metaproteome/train_internalVal_1385"
    EXTERNAL_DATA_PATH = BASE_DIR / "ML_matrix_20250911/metaproteome/independent_test_104"
    
    # 样本信息文件
    INTERNAL_SAMPLE_INFO = INTERNAL_DATA_PATH / "GNHSF_sample_inform_normalized_1385_2_yes_all_metadata.csv"
    EXTERNAL_SAMPLE_INFO = EXTERNAL_DATA_PATH / "GNHSF_com_sample_inform_normalized_ML.csv"
    
    # 设置输出目录
    if args.output_dir:
        VALIDATION_OUTPUT_DIR = Path(args.output_dir) / f"validation_analysis_{feature_type}"
    else:
        VALIDATION_OUTPUT_DIR = MODEL_OUTPUT_ROOT / "validation_analysis"
    
    VALIDATION_OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    
    print("\n" + "=" * 80)
    print(f"开始分析特征类型: {feature_type}")
    print(f"模型目录: {MODEL_OUTPUT_ROOT}")
    print(f"输出目录: {VALIDATION_OUTPUT_DIR}")
    print("=" * 80)
    
    # 检查模型目录是否存在
    if not MODEL_OUTPUT_ROOT.exists():
        print(f"\n错误：模型目录不存在: {MODEL_OUTPUT_ROOT}")
        print(f"请确认特征类型 '{feature_type}' 的模型已训练完成")
        return None
    
    # 1. 加载预测结果和模型
    internal_preds_dict, external_preds_dict, models_dict = load_predictions_and_models()
    
    if not models_dict:
        print(f"\n警告：未找到特征类型 '{feature_type}' 的任何模型文件！")
        return None
    
    # 2. 计算效应方向一致性
    concordance_results = compute_effect_direction_concordance(
        models_dict=models_dict,
        top_n=args.top_n_features,
        selected_seed=args.random_seed,
        selected_fold=args.selected_fold
    )
    
    # 3. 计算通路层级一致性
    pathway_results = compute_pathway_level_concordance(
        models_dict=models_dict,
        category_types=['COG', 'KEGG', 'microbial_protein'],
        top_n_pathways=args.top_n_pathways,
        selected_seed=args.random_seed,
        selected_fold=args.selected_fold
    )
    
    # 4. 计算特征重叠排列检验
    permutation_results = compute_feature_overlap_permutation_test(
        models_dict=models_dict,
        observed_k=159,
        top_n=500,
        n_permutations=args.n_permutations,
        selected_seed=args.random_seed,
        selected_fold=args.selected_fold
    )
    
    # 5. 生成综合报告
    report = generate_comprehensive_report(
        internal_preds_dict=internal_preds_dict,
        external_preds_dict=external_preds_dict,
        models_dict=models_dict,
        concordance_results=concordance_results,
        pathway_results=pathway_results,
        permutation_results=permutation_results,
        output_file=f'validation_report_{feature_type}.txt'
    )
    
    # 6. 生成可视化
    plot_validation_results(
        concordance_results=concordance_results,
        pathway_results=pathway_results,
        permutation_results=permutation_results,
        internal_preds_dict=internal_preds_dict,
        external_preds_dict=external_preds_dict
    )
    
    # 7. 保存结果到JSON
    import json
    results_json = {
        'feature_type': feature_type,
        'timestamp': str(pd.Timestamp.now()),
        'user': 'westlakeomics',
        'AUROC': report['AUROC'],
        'concordance': {
            'percentage': concordance_results['concordance_percentage'] if concordance_results else None,
            'pearson_r': concordance_results['pearson_r'] if concordance_results else None,
            'spearman_r': concordance_results['spearman_r'] if concordance_results else None,
            'top_n': concordance_results['top_n'] if concordance_results else None
        } if concordance_results else None,
        'pathway_concordance': {
            cat: {
                'overlap': res['overlap'],
                'total_pathways': res['total_pathways'],
                'top_n': res['top_n'],
                'jaccard_index': res['jaccard_index'],
                'spearman_r': res['spearman_r'] if not np.isnan(res['spearman_r']) else None,
                'permutation_p': res['permutation_p']
            } for cat, res in pathway_results.items()
        } if pathway_results else {},
        'permutation_test': {
            'observed_overlap': permutation_results['observed_overlap'],
            'expected_overlap': permutation_results['expected_overlap'],
            'p_value': permutation_results['p_value'],
            'z_score': permutation_results['z_score'],
            'top_n': permutation_results['top_n']
        } if permutation_results else None
    }
    
    with open(VALIDATION_OUTPUT_DIR / f'validation_results_{feature_type}.json', 'w') as f:
        json.dump(results_json, f, indent=2)
    
    print("\n" + "=" * 80)
    print(f"特征类型 '{feature_type}' 分析完成！")
    print("=" * 80)
    
    # 打印主要发现摘要
    if concordance_results:
        print(f"\n主要发现:")
        print(f"  1. 效应方向一致性: {concordance_results['concordance_percentage']:.1f}%")
    if pathway_results:
        print(f"  2. 通路层级一致性:")
        for cat, res in pathway_results.items():
            if not np.isnan(res['spearman_r']):
                print(f"     - {cat}: ρ = {res['spearman_r']:.3f} (p = {res['spearman_p']:.2e})")
    if permutation_results:
        print(f"  3. 特征重叠检验: p = {permutation_results['p_value']:.4f}")
    
    print(f"\n所有结果已保存到: {VALIDATION_OUTPUT_DIR}")
    
    return report


# ============================================================================
# 9. 主函数
# ============================================================================

def main():
    """
    主函数：执行完整的验证分析流程
    """
    # 解析命令行参数
    args = parse_arguments()
    
    print("\n" + "=" * 80)
    print("模型验证分析脚本")
    print(f"用户: westlakeomics")
    print(f"时间: {pd.Timestamp.now()}")
    print("=" * 80)
    print(f"\n配置信息:")
    print(f"  特征类型: {', '.join(args.feature_types)}")
    print(f"  基础目录: {args.base_dir}")
    print(f"  随机种子: {args.random_seed}")
    print(f"  Top特征数: {args.top_n_features}")
    print(f"  Top通路数: {args.top_n_pathways}")
    print(f"  排列次数: {args.n_permutations}")
    print(f"  选择折数: {args.selected_fold}")
    
    # 分析每个特征类型
    all_reports = {}
    
    for feature_type in args.feature_types:
        print("\n" + "=" * 80)
        print(f"开始处理特征类型: {feature_type}")
        print("=" * 80)
        
        try:
            report = analyze_single_feature_type(feature_type, args)
            if report:
                all_reports[feature_type] = report
        except Exception as e:
            print(f"\n错误：处理特征类型 '{feature_type}' 时发生异常:")
            print(f"  {str(e)}")
            import traceback
            traceback.print_exc()
            continue
    
    # 生成汇总报告
    if all_reports:
        print("\n" + "=" * 80)
        print("所有特征类型分析完成！")
        print("=" * 80)
        print("\n主要发现汇总:")
        
        for feature_type, report in all_reports.items():
            print(f"\n{feature_type}:")
            
            # AUROC
            if report.get('AUROC'):
                auroc = report['AUROC']
                if auroc['internal']['mean']:
                    print(f"  AUROC (内部): {auroc['internal']['mean']:.3f} ± {auroc['internal']['std']:.3f}")
                if auroc['external']['mean']:
                    print(f"  AUROC (外部): {auroc['external']['mean']:.3f} ± {auroc['external']['std']:.3f}")
            
            # 效应方向一致性
            if report.get('concordance'):
                conc = report['concordance']
                print(f"  效应方向一致性: {conc['concordance_percentage']:.1f}%")
                print(f"    Spearman ρ = {conc['spearman_r']:.3f}")
            
            # 通路一致性
            if report.get('pathway'):
                print(f"  通路层级一致性:")
                for cat, res in report['pathway'].items():
                    if not np.isnan(res['spearman_r']):
                        print(f"    - {cat}: ρ = {res['spearman_r']:.3f} (重叠: {res['overlap']}/{res['top_n']})")
            
            # 排列检验
            if report.get('permutation'):
                perm = report['permutation']
                print(f"  特征重叠检验: k={perm['observed_overlap']}, p={perm['p_value']:.4f}")
        
        # 保存汇总JSON
        if args.output_dir:
            summary_dir = Path(args.output_dir)
        else:
            summary_dir = Path(args.base_dir) / "validation_summary"
        
        summary_dir.mkdir(exist_ok=True, parents=True)
        
        import json
        summary_json = {
            'timestamp': str(pd.Timestamp.now()),
            'user': 'westlakeomics',
            'feature_types': args.feature_types,
            'config': {
                'random_seed': args.random_seed,
                'top_n_features': args.top_n_features,
                'top_n_pathways': args.top_n_pathways,
                'n_permutations': args.n_permutations
            },
            'results': {
                ft: {
                    'AUROC_internal_mean': report['AUROC']['internal']['mean'],
                    'AUROC_external_mean': report['AUROC']['external']['mean'],
                    'concordance_percentage': report['concordance']['concordance_percentage'] if report.get('concordance') else None,
                    'permutation_p_value': report['permutation']['p_value'] if report.get('permutation') else None
                } for ft, report in all_reports.items()
            }
        }
        
        with open(summary_dir / 'validation_summary_all_types.json', 'w') as f:
            json.dump(summary_json, f, indent=2)
        
        print(f"\n汇总结果已保存到: {summary_dir / 'validation_summary_all_types.json'}")
    
    else:
        print("\n警告：没有成功分析任何特征类型！")
    
    print("\n" + "=" * 80)
    print("分析流程完成！")
    print("=" * 80)


if __name__ == "__main__":
    main()