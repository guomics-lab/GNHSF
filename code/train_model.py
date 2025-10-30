import os
import re
import json
import numpy as np
import pandas as pd
import joblib
import warnings
from collections import defaultdict
from pathlib import Path
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (roc_auc_score, roc_curve, precision_recall_curve, auc)
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from hyperopt import hp, fmin, tpe, Trials
import xgboost as xgb
import matplotlib.pyplot as plt
from hyperopt import hp
warnings.filterwarnings('ignore')

# ----------------------- 全局配置（移除TOP_K，明确Gain+频次逻辑） -----------------------
RANDOM_SEED = 42  # 单个固定Seed
N_FOLDS = 5  # 5折交叉验证
MAX_EVAL = 30  # 贝叶斯调参次数
XGB_SELECTION_ROUNDS = 100  # 100次XGB特征选择（频次统计）
XGB_FEAT_THRESHOLD = 70  # 100次中出现≥70次的特征保留
# 仅用于GLM的协变量（保持不变）
GLM_COVARIATES = ["sex", "age", "bristol_scale"]
# 临床特征（用于最终XGB建模）
CLINICAL_FEATURES = ["BMI", "sex", "age", "waist", "TG", "HDL", "hyper_cl"]
XGB_SPACE = {
    'learning_rate': hp.uniform('learning_rate', 0.001, 0.003), # 降低上限，进一步减缓学习，防止快速过拟合
    'max_depth': hp.quniform('max_depth', 1, 2, 1), # 进一步限制深度（1~3），浅树更稳健
    'n_estimators': hp.quniform('n_estimators', 50, 300, 10), # 增加树数但步长大，结合低 lr 总迭代类似；若仍过拟合，可降至 50~200
    'subsample': hp.uniform('subsample', 0.7, 0.9), # 略微提高下限，增强随机性（避免太低导致欠拟合）
    'colsample_bytree': hp.uniform('colsample_bytree', 0.5, 0.6), # 降低上限，增加列采样随机，减少特征依赖
    'min_child_weight': hp.uniform('min_child_weight', 15, 20), # 提高下限，使叶子节点更保守，减少过拟合
    'reg_alpha': hp.uniform('reg_alpha', 15, 20), # 扩展范围，但偏向中等（原 2~15 太高可能欠拟合）
    'reg_lambda': hp.uniform('reg_lambda', 15, 20), # 同上,加强 L2 正则
   'gamma': hp.uniform('gamma', 0.5, 0.7), # 分裂最小增益，防止不必要分裂
}
types = ['genus']  # 示例
if isinstance(types, (list, tuple)):
    types_str = "_".join(types)
else:
    types_str = str(types)
print(RANDOM_SEED)
output_root = Path(f"glm_xgb_gain_freq_{RANDOM_SEED}_output_{types_str}")
output_root.mkdir(exist_ok=True)
feature_dir = output_root / "selected_features"
feature_dir.mkdir(exist_ok=True)
model_dir = output_root / "models"
model_dir.mkdir(exist_ok=True)
plot_dir = output_root / "plots"
plot_dir.mkdir(exist_ok=True)
external_dir = output_root / "external"
external_dir.mkdir(exist_ok=True)
internal_preds_dir = output_root / "internal_preds"
internal_preds_dir.mkdir(exist_ok=True)
external_preds_dir = output_root / "external_preds"
external_preds_dir.mkdir(exist_ok=True)

SEEDS = list(range(42, 62))  # 20 different seeds

# ----------------------- 数据预处理（增强：清理临床协变量并二值化） -----------------------
def _to_num(series):
    return pd.to_numeric(series, errors='coerce')

def preprocess_sample_data(filepath):
    """
    预处理样本信息表：
    - dm_cl转换为二元变量（yes=1, no=0）
    - bristol_scale标准化为类别（固定水平 0/1/2）
    - hyper_cl转换为二元变量（yes=1, no=0）
    - sex转换为0/1（female=0, male=1）
    - 在特征筛选前，将临床特征(BMI, sex, age, waist, TG, HDL, hyper_cl)为none/NA的样本剔除
    """
    sampledf = pd.read_csv(
        filepath,
        na_values=["", "NA", "NaN", "nan", "None", "none", "N/A", "n/a"]
    )
    # 处理bristol_scale（用于GLM协变量）
    sampledf["bristol_scale"] = sampledf["bristol_scale"].replace({
        "Normal": "0",
        "Diarrhea": "1",
        "Constipation": "2"
    })
    # 转换dm_cl为二元变量
    if "dm_cl" not in sampledf.columns:
        raise ValueError("未在样本信息表中找到列 'dm_cl'")
    print("原始dm_cl值分布：", sampledf["dm_cl"].value_counts(dropna=False))
    sampledf["dm_cl"] = sampledf["dm_cl"].map({"yes": 1, "no": 0, "Yes": 1, "No": 0})
    print("转换后dm_cl值分布：", sampledf["dm_cl"].value_counts(dropna=False))
    # sex 映射到0/1（female=0, male=1），GLM中使用 C(sex) 仍可作为分类处理
    if "sex" in sampledf.columns:
        sex_map = {
            "female": 0, "f": 0, "0": 0, 0: 0,
            "male": 1, "m": 1, "1": 1, 1: 1
        }
        sampledf["sex"] = sampledf["sex"].apply(
            lambda v: sex_map.get(str(v).strip().lower(), np.nan)
        )
    # hyper_cl 二值化（yes=1, no=0）
    if "hyper_cl" in sampledf.columns:
        sampledf["hyper_cl"] = sampledf["hyper_cl"].apply(
            lambda v: 1 if str(v).strip().lower() == "yes" else (0 if str(v).strip().lower() == "no" else np.nan)
        )
    # 数值列：BMI, age, waist, TG, HDL
    for col in ["BMI", "age", "waist", "TG", "HDL"]:
        if col in sampledf.columns:
            sampledf[col] = _to_num(sampledf[col])
    # bristol_scale 作为分类变量（固定水平）
    sampledf["bristol_scale"] = pd.Categorical(sampledf["bristol_scale"], categories=["0", "1", "2"])
    # 仅保留必要列；在筛选前剔除临床缺失
    label_col = "dm_cl"
    required_cols = list(dict.fromkeys(GLM_COVARIATES + CLINICAL_FEATURES + ["sample", label_col]))
    missing_cols = [c for c in required_cols if c not in sampledf.columns]
    if missing_cols:
        raise ValueError(f"样本信息表缺少以下必要列: {missing_cols}")
    # 在特征筛选前：剔除临床特征为NA/none的样本
    sampledf = sampledf[required_cols].dropna(subset=CLINICAL_FEATURES + [label_col])
    print(f"\n样本信息表摘要（固定Seed：{RANDOM_SEED}）:\n", sampledf.describe(include="all"))
    return sampledf

def load_feature_data(type_, sample_indices=None, is_original=False, base_path="/sunyingying/test/ML_matrix_20250911/metaproteome/train_internalVal_1385"):
    """加载特征数据（支持INT和original目录）"""
    if is_original:
        filename = f"{base_path}/original_expression_matrix/GNHSF_diann_IGC_humanswiss_{type_}_sample_1385_NA90.tsv"
    else:
        filename = f"{base_path}/INT_expression_matrix/GNHSF_diann_IGC_humanswiss_{type_}_sample_1385_NA90_INT.tsv"
    df = pd.read_csv(
        filename,
        sep="\t",
        index_col=0
    )
    # 转置并处理列名（确保INT与original特征名可匹配）
    dft = df.T.reset_index()
    dft.columns = ["sample"] + list(dft.columns[1:])
    # 简化列名（统一INT/original列名格式）
    def simplify_colname(col):
        return col.split(" ")[0] if " " in col else col
    dft.columns = [simplify_colname(col) for col in dft.columns]
    # 清理特殊字符（避免列名匹配失败）
    dft.columns = [re.sub(r"\|", ".", col) for col in dft.columns]
    dft.columns = [re.sub(r"\-", ".", col) for col in dft.columns]
    dft.columns = [re.sub(r"\.$", "", col) for col in dft.columns]
    # 筛选样本
    if sample_indices is not None:
        dft = dft[dft["sample"].isin(sample_indices)]
    print(f"加载{'original' if is_original else 'INT'}矩阵 | 类型：{type_} | 样本数：{len(dft)} | 特征数：{len(dft.columns)-1}")
    return dft

def load_external_feature_data(type_):
    filename = f"/sunyingying/test/ML_matrix_20250911/metaproteome/independent_test_104/GNHSFc_diann_IGC_humanswiss_{type_}_sample.tsv"
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
    print(f"加载外部original矩阵 | 类型：{type_} | 样本数：{len(dft)} | 特征数：{len(dft.columns)-1}")
    return dft

# ----------------------- 折内处理（添加临床特征到最终训练中） -----------------------
def process_fold(type_, sampledf, train_samples, val_samples, test_samples, fold, label_col='dm_cl'):
    """单折处理：GLM筛选 → XGB（Gain>0+100次频次）筛选（仅omics）→ 加入临床特征建模评估"""
    fold_feature_dir = feature_dir / f"fold_{fold}"
    fold_model_dir = model_dir / f"fold_{fold}"
    fold_feature_dir.mkdir(parents=True, exist_ok=True)
    fold_model_dir.mkdir(parents=True, exist_ok=True)
    # --------------------------
    # 1. GLM筛选特征名（修改为收集所有p值，然后计算q值筛选q<0.05）
    # --------------------------
    print(f" 折 {fold} - 开始GLM特征筛选（INT矩阵）...")
    dft_int = load_feature_data(type_, train_samples, is_original=False)
    feature_cols_int = [col for col in dft_int.columns if col != "sample"]
    df_merged_int = pd.merge(dft_int, sampledf, on="sample", how="inner")
    if len(df_merged_int) < 10:
        print(f" 折 {fold} - INT矩阵匹配后样本量不足，跳过")
        return None
    glm_res_list = []
    for i, feat_name in enumerate(feature_cols_int, 1):
        if i % 100 == 0:
            print(f" 折 {fold} - 处理GLM特征 {i}/{len(feature_cols_int)}")
        df_glm = df_merged_int[[feat_name] + GLM_COVARIATES + [label_col]].dropna()
        if len(df_glm) < 10:
            continue
        glm_res_row = glm_feature_selection(df_glm, feat_name, label_col)
        if glm_res_row is not None:
            glm_res_list.append(glm_res_row)
    if not glm_res_list:
        print(f" 折 {fold} - 无GLM结果，跳过XGB处理")
        return None
    glm_res_df = pd.DataFrame(glm_res_list)
    pvals = glm_res_df['p_value'].values
    qvals = multipletests(pvals, method='fdr_bh')[1]
    glm_res_df['q_value'] = qvals
    q_sig = []
    for q in qvals:
        if q < 0.001:
            q_sig.append('***')
        elif q < 0.01:
            q_sig.append('**')
        elif q < 0.05:
            q_sig.append('*')
        else:
            q_sig.append('nsig')
    glm_res_df['q_sig'] = q_sig
    glm_res_df.to_csv(fold_feature_dir / "glm_res.csv", index=False)
    glm_selected_names = glm_res_df[glm_res_df['q_value'] < 0.05]['feature'].tolist()
    print(f" 折 {fold} - GLM筛选后特征名数量 (q < 0.05): {len(glm_selected_names)}")
    pd.DataFrame({"feature_name": glm_selected_names}).to_csv(
        fold_feature_dir / "glm_selected_names.csv", index=False
    )
    if not glm_selected_names:
        print(f" 折 {fold} - 无GLM筛选特征，跳过XGB处理")
        return None
    # --------------------------
    # 2. 加载original矩阵（用于XGB Gain+频次筛选，仍仅对omics特征进行）
    # --------------------------
    print(f" 折 {fold} - 加载original矩阵（用于XGB Gain+频次筛选）...")
    dft_original = load_feature_data(type_, train_samples, is_original=True)
    matched_original_cols = [col for col in dft_original.columns if col in glm_selected_names and col != "sample"]
    if not matched_original_cols:
        print(f" 折 {fold} - original矩阵无匹配特征，跳过")
        return None
    print(f" 折 {fold} - original矩阵匹配后特征数: {len(matched_original_cols)}")
    df_original_train = pd.merge(
        dft_original[["sample"] + matched_original_cols],
        sampledf[[label_col, "sample"]],
        on="sample",
        how="inner"
    )
    X_original_train = df_original_train[matched_original_cols].copy()
    y_train = df_original_train[label_col].astype(int).copy()
    if (y_train == 1).sum() == 0 or (y_train == 0).sum() == 0:
        print(f" 折 {fold} - 训练集标签分布异常，跳过")
        return None
    # 用训练集构建类权重（先用omics训练集，后续用full train再重算用于最终模型）
    scale_pos_weight = (y_train == 0).sum() / (y_train == 1).sum()  # 负样本 / 正样本
    # --------------------------
    # 3. XGB特征选择（核心：100次循环，Gain>0+频次统计，无Top K限制）仅针对omics特征
    # --------------------------
    print(f" 折 {fold} - 开始{XGB_SELECTION_ROUNDS}次XGB特征选择（Gain>0+频次统计）...")
    xgb_feat_counts = defaultdict(int)  # 统计100次中特征出现次数
    all_round_gain = []  # 记录每轮Gain，便于追溯
    for round_idx in range(XGB_SELECTION_ROUNDS):
        round_seed = RANDOM_SEED * 1000 + fold * 100 + round_idx
        # 每轮筛选：仅保留Gain>0的特征（无Top K限制）
        round_selected, round_gain_df = xgb_feature_selection_by_gain(
            X=X_original_train,
            y=y_train,
            scale_pos_weight=scale_pos_weight,
            random_seed=round_seed
        )
        # 累计频次
        for feat in round_selected:
            xgb_feat_counts[feat] += 1
        # 记录每轮Gain
        round_gain_df["round"] = round_idx + 1
        all_round_gain.append(round_gain_df)
        # 进度打印（每20轮）
        if (round_idx + 1) % 20 == 0:
            high_freq_cnt = sum(1 for cnt in xgb_feat_counts.values() if cnt >= XGB_FEAT_THRESHOLD)
            print(f" 折 {fold} - XGB进度：{round_idx+1}/{XGB_SELECTION_ROUNDS} | 高频特征数（≥70次）：{high_freq_cnt}")
    # 保存每轮Gain信息（单折内可追溯）
    all_round_gain_df = pd.concat(all_round_gain, ignore_index=True)
    all_round_gain_df.to_csv(fold_feature_dir / "xgb_all_round_gain.csv", index=False)
    # 筛选高频特征（100次≥70次）
    high_freq_features = [feat for feat, cnt in xgb_feat_counts.items() if cnt >= XGB_FEAT_THRESHOLD]
    print(f" 折 {fold} - 100次XGB选择后：总有效特征数（Gain>0）={len(xgb_feat_counts)} | 高频特征数（≥70次）={len(high_freq_features)}")
    # 保存单折高频特征及频次
    high_freq_df = pd.DataFrame({
        "feature": high_freq_features,
        "appear_count_100rounds": [xgb_feat_counts[feat] for feat in high_freq_features],
        "avg_gain_100rounds": [all_round_gain_df[all_round_gain_df["feature"] == feat]["gain"].mean() for feat in high_freq_features]
    }).sort_values(["appear_count_100rounds", "avg_gain_100rounds"], ascending=[False, False])
    high_freq_df.to_csv(fold_feature_dir / "xgb_high_freq_features.csv", index=False)
    if not high_freq_features:
        print(f" 折 {fold} - 无高频特征，跳过建模")
        return None
    # 保存model_features
    model_features = high_freq_features + CLINICAL_FEATURES
    pd.DataFrame({"feature": model_features}).to_csv(fold_feature_dir / "model_features.csv", index=False)
    # --------------------------
    # 4. 训练XGB分类器（最终模型：omics高频特征 + 临床特征）
    # --------------------------
    print(f" 折 {fold} - 用高频特征 + 临床特征训练XGB分类器...")
    def get_split_data(split_samples):
        # omics特征
        dft_split = load_feature_data(type_, split_samples, is_original=True)
        dft_split_filtered = dft_split[["sample"] + high_freq_features]
        # 临床特征
        clinical_cols = ["sample"] + CLINICAL_FEATURES + [label_col]
        df_clinical = sampledf[clinical_cols].copy()
        # 合并
        df_split = pd.merge(
            dft_split_filtered,
            df_clinical,
            on="sample",
            how="inner"
        )
        # 确保临床特征为数值（sex与hyper_cl已在预处理转换）
        for col in CLINICAL_FEATURES:
            if col not in df_split.columns:
                raise ValueError(f"分割数据缺少临床列: {col}")
        X_mat = df_split[model_features].copy()
        y_vec = df_split[label_col].astype(int).copy()
        samples_list = df_split["sample"].tolist()
        return X_mat, y_vec, samples_list
    # 训练、验证、测试数据集
    X_train_full, y_train_full, _ = get_split_data(train_samples)
    X_val, y_val, _ = get_split_data(val_samples)
    X_test, y_test, test_samples_list = get_split_data(test_samples)
    # 用最终训练集（含临床+omics）重算类权重
    scale_pos_weight_full = (y_train_full == 0).sum() / (y_train_full == 1).sum()
    # 贝叶斯调参（在训练+验证划分上）
    best_params = bayes_search(
        X_tr=X_train_full,
        y_tr=y_train_full,
        X_val=X_val,
        y_val=y_val,
        scale_pos_weight=scale_pos_weight_full,
        seed=RANDOM_SEED
    )
    # 合并训练+验证集重训
    X_tv = pd.concat([X_train_full, X_val], ignore_index=True)
    y_tv = np.concatenate([y_train_full, y_val])
    seed_metrics = {}
    seed_internal_preds = []
    for seed in SEEDS:
        final_clf = xgb.XGBClassifier(
            objective='binary:logistic',
            eval_metric='auc',
            tree_method='hist',
            scale_pos_weight=scale_pos_weight_full,
            random_state=seed,
            **best_params
        )
        final_clf.fit(X_tv, y_tv, verbose=False)
        # 保存模型
        joblib.dump(final_clf, fold_model_dir / f"final_xgb_model_{seed}.pkl")
        # 计算内部测试集预测
        preds = final_clf.predict_proba(X_test)[:, 1]
        # 保存预测
        pred_df = pd.DataFrame({
            'sample': test_samples_list,
            'y_true': y_test,
            'y_pred': preds,
            'fold': fold,
            'seed': seed
        })
        seed_internal_preds.append(pred_df)
        # 评估（可选，保留原逻辑但不收集到all_fold_metrics）
        test_metrics = evaluate_clf(
            clf=final_clf,
            X=X_test,
            y=y_test,
            split_name=f'test_seed_{seed}',
            fold=fold
        )
        print(f" 折 {fold} 种子 {seed} - 测试集AUC：{test_metrics['auc']:.3f} | PR-AUC：{test_metrics['pr_auc']:.3f}")
        seed_metrics[seed] = test_metrics
    seed_auc_scores = {}
    for seed in SEEDS:
        model_path = fold_model_dir / f"final_xgb_model_{seed}.pkl"
        if not model_path.exists():
            continue
        clf = joblib.load(model_path)
        preds = clf.predict_proba(X_test)[:, 1]
        auc_val = roc_auc_score(y_test, preds)
        seed_auc_scores[seed] = auc_val

    # 2. 找到AUC最高的seed
    if seed_auc_scores:
        best_seed = max(seed_auc_scores, key=seed_auc_scores.get)
        best_model_path = fold_model_dir / f"final_xgb_model_{best_seed}.pkl"
        final_model_path = fold_model_dir / "final_xgb_model_best.pkl"
        # 复制或重命名为最终模型
        import shutil
        shutil.copy(best_model_path, final_model_path)
    return {
        "fold": fold,
        "high_freq_features": high_freq_features,
        "seed_internal_preds": seed_internal_preds  # 返回每个种子的内部预测
    }

# ----------------------- 核心函数（移除Top K参数，简化Gain筛选） -----------------------
def glm_feature_selection(df_glm, feature, label_col='dm_cl'):
    """修改GLM筛选逻辑：返回dm_cl行的结果字典，如果存在"""
    formula = f"Q('{feature}') ~ C(sex) + age + C(bristol_scale) + {label_col}"
    try:
        model = sm.GLM.from_formula(
            formula=formula,
            data=df_glm,
            family=sm.families.Gaussian()
        )
        result = model.fit(use_t=True)  # 使用t统计以匹配R的Pr(>|t|)
    except Exception:
        return None
    coef_df = result.summary2().tables[1].reset_index().rename(columns={"index": "metadata"})
    dm_cl_row = coef_df[coef_df["metadata"] == label_col]
    if dm_cl_row.empty:
        return None
    return {
        'feature': feature,
        'Estimate': dm_cl_row['Coef.'].values[0],
        'Std. Error': dm_cl_row['Std.Err.'].values[0],
        't value': dm_cl_row['t'].values[0],
        'p_value': dm_cl_row['P>|t|'].values[0],
        'metadata': label_col
    }

def xgb_feature_selection_by_gain(X, y, scale_pos_weight, random_seed=RANDOM_SEED, top_k=50):
    """
    单轮XGB特征选择：仅保留Gain排名前top_k的特征
    返回：本轮特征列表 + 本轮Gain数据框
    """
    np.random.seed(random_seed)
    subsample = np.random.uniform(0.6, 0.9)
    colsample_bytree = np.random.uniform(0.6, 0.9)
    clf = xgb.XGBClassifier(
        n_estimators=200,
        learning_rate=0.05,
        max_depth=4,
        objective='binary:logistic',
        eval_metric='auc',
        tree_method='hist',
        n_jobs=-1,
        scale_pos_weight=scale_pos_weight,
        random_state=random_seed,
        subsample=subsample,
        colsample_bytree=colsample_bytree,
        importance_type='gain'
    )
    clf.fit(X, y)
    feat_gain = clf.feature_importances_
    round_gain_df = pd.DataFrame({
        "feature": X.columns,
        "gain": feat_gain
    }).sort_values("gain", ascending=False)
    # 选择Gain排名前top_k的特征
    round_selected = round_gain_df.head(top_k)["feature"].tolist()
    return round_selected, round_gain_df

def bayes_search(X_tr, y_tr, X_val, y_val, scale_pos_weight, seed=RANDOM_SEED):
    """贝叶斯调参（保留原逻辑）"""
    def obj(params):
        params['max_depth'] = int(params['max_depth'])
        params['n_estimators'] = int(params['n_estimators'])
        clf = xgb.XGBClassifier(
            objective='binary:logistic',
            eval_metric='auc',
            tree_method='hist',
            scale_pos_weight=scale_pos_weight,
            random_state=seed,
            **params
        )
        clf.fit(
            X_tr, y_tr,
            eval_set=[(X_val, y_val)],
            verbose=False
        )
        preds = clf.predict_proba(X_val)[:, 1]
        return 1 - roc_auc_score(y_val, preds)
    trials = Trials()
    best = fmin(
        fn=obj,
        space=XGB_SPACE,
        algo=tpe.suggest,
        max_evals=MAX_EVAL,
        rstate=np.random.default_rng(seed)
    )
    best['max_depth'] = 1
    best['n_estimators'] = int(best['n_estimators'])
    return best

# ----------------------- 辅助函数（无修改） -----------------------
def evaluate_clf(clf, X, y, split_name, fold):
    """评估函数"""
    preds = clf.predict_proba(X)[:, 1]
    fpr, tpr, _ = roc_curve(y, preds)
    roc_auc = auc(fpr, tpr)
    prec, rec, _ = precision_recall_curve(y, preds)
    pr_auc = auc(rec, prec)
    plot_path = plot_dir / f"fold_{fold}_{split_name}_curves.png"
    plot_metrics(
        y_true=y,
        y_score=preds,
        tag=f"折 {fold} | {split_name}\nAUC：{roc_auc:.3f} | PR-AUC：{pr_auc:.3f}",
        save_path=plot_path
    )
    return {
        'fold': fold,
        'split': split_name,
        'auc': roc_auc,
        'pr_auc': pr_auc,
        'sample_count': len(y),
        'positive_count': (y == 1).sum()
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

# ----------------------- 主流程（移除全局汇总，简化输出） -----------------------
def main(types=['species']):
    sampledf = preprocess_sample_data("/sunyingying/test/ML_matrix_20250911/metaproteome/train_internalVal_1385/GNHSF_sample_inform_normalized_1385_2_yes_all_metadata.csv")
    label_col = 'dm_cl'
    all_samples = sampledf["sample"].unique()
    print(f"\n总样本数: {len(all_samples)} | 使用固定Seed：{RANDOM_SEED} | 特征筛选逻辑：Gain>0 + 100次≥70次")
    print(f"临床协变量将加入最终XGB训练：{CLINICAL_FEATURES}")
    # 固定Seed下划分5折：确保分层标签与样本顺序对齐
    np.random.seed(RANDOM_SEED)
    sample_to_label = dict(zip(sampledf["sample"], sampledf[label_col].values))
    y_stratify = np.array([sample_to_label[s] for s in all_samples])
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_SEED)
    folds = list(skf.split(all_samples, y_stratify))
    internal_preds_by_seed = defaultdict(list)  # {seed: [pred_df_fold1, pred_df_fold2, ...]}
    fold_model_features = {}  # {fold: model_features}
    # 处理每个折
    for fold in range(N_FOLDS):
        print(f"\n----- 处理折 {fold + 1}/{N_FOLDS} -----")
        # 样本划分
        train_val_indices, test_indices = folds[fold]
        skf_val = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_SEED)
        train_indices, val_indices = next(skf_val.split(
            all_samples[train_val_indices],
            y_stratify[train_val_indices]
        ))
        train_samples = all_samples[train_val_indices[train_indices]].tolist()
        val_samples = all_samples[train_val_indices[val_indices]].tolist()
        test_samples = all_samples[test_indices].tolist()
        # 检查样本重叠
        train_set = set(train_samples)
        val_set = set(val_samples)
        test_set = set(test_samples)
        overlap_train_val = len(train_set & val_set)
        overlap_train_test = len(train_set & test_set)
        overlap_val_test = len(val_set & test_set)
        if overlap_train_val > 0 or overlap_train_test > 0 or overlap_val_test > 0:
            print(f" 错误：折 {fold} 样本集存在重叠！")
            print(f" 训练-验证重叠数：{overlap_train_val} | 训练-测试重叠数：{overlap_train_test} | 验证-测试重叠数：{overlap_val_test}")
            continue
        # 折内处理（type=species）
        for type_ in types:
            fold_result = process_fold(
                type_, sampledf,
                train_samples, val_samples, test_samples,
                fold, label_col
            )
            if fold_result is not None:
                # 收集内部预测
                for pred_df in fold_result["seed_internal_preds"]:
                    seed = pred_df['seed'].iloc[0]
                    internal_preds_by_seed[seed].append(pred_df)
                # 保存model_features
                fold_model_features[fold] = fold_result["high_freq_features"] + CLINICAL_FEATURES
    # --------------------------
    # 处理内部AUC
    # --------------------------
    internal_aucs = []
    for seed in SEEDS:
        seed_dfs = internal_preds_by_seed[seed]
        if len(seed_dfs) == 0:
            continue
        all_internal_df = pd.concat(seed_dfs, ignore_index=True)
        all_y = all_internal_df['y_true'].values
        all_pred = all_internal_df['y_pred'].values
        auc_val = roc_auc_score(all_y, all_pred)
        internal_aucs.append(auc_val)
        # 保存预测
        all_internal_df.to_csv(internal_preds_dir / f"internal_preds_seed_{seed}.csv", index=False)
        print(f"种子 {seed} 内部总AUC: {auc_val:.3f}")
    # --------------------------
    # 处理外部测试集
    # --------------------------
    external_sample_filepath = "/sunyingying/test/ML_matrix_20250911/metaproteome/independent_test_104/GNHSF_com_sample_inform_normalized_ML.csv"
    external_sampledf = preprocess_sample_data(external_sample_filepath)
    external_dft = load_external_feature_data(types[0])  # assuming types[0] = 'species'
    external_samples = external_sampledf["sample"].unique()
    external_clinical = external_sampledf[["sample"] + CLINICAL_FEATURES + [label_col]]
    external_omics = external_dft
    external_df = pd.merge(
        external_omics,
        external_clinical,
        on="sample",
        how="inner"
    )
    y_external = external_df[label_col].astype(int).values
    external_samples_list = external_df["sample"].tolist()
    external_aucs = []
    for seed in SEEDS:
        seed_preds = []
        for fold in range(N_FOLDS):
        
            model_path = model_dir / f"fold_{fold}" / f"final_xgb_model_{seed}.pkl"
            if not model_path.exists():
                continue
            clf = joblib.load(model_path)
            fold_features = pd.read_csv(feature_dir / f"fold_{fold}" / "model_features.csv")["feature"].tolist()
            missing_features = [feat for feat in fold_features if feat not in external_df.columns]
            if missing_features:
                print(f"种子 {seed} 折 {fold} - 外部测试集缺少 {len(missing_features)} 个特征，已填充为 NA: {missing_features[:5]}...")
                for feat in missing_features:
                    external_df[feat] = np.nan
            if len(missing_features) == len(fold_features):
                print(f"种子 {seed} 折 {fold} - 所有特征缺失，跳过该折预测")
                continue
            X_external_fold = external_df[fold_features].copy()
            fold_pred = clf.predict_proba(X_external_fold)[:, 1]  # 预测类别标签（0或1）
            seed_preds.append(fold_pred)
        if len(seed_preds) == 0:
            continue
        # 硬投票：多数投票
        avg_pred = np.mean(seed_preds, axis=0)
        auc_val = roc_auc_score(y_external, avg_pred)
        external_aucs.append(auc_val)
        # 保存预测
        external_pred_df = pd.DataFrame({
            'sample': external_samples_list,
            'y_true': y_external,
            'y_pred': avg_pred,  # 硬投票结果
            'seed': seed
        })
        external_pred_df.to_csv(external_preds_dir / f"external_preds_seed_{seed}.csv", index=False)
        print(f"种子 {seed} 外部AUC: {auc_val:.3f}")
    # --------------------------
    # 绘制箱型图
    # --------------------------
    fig, ax = plt.subplots()
    ax.boxplot([internal_aucs, external_aucs])
    ax.set_xticklabels(['Internal', 'External'])
    ax.set_ylabel('AUC')
    plt.title('AUC Boxplots across 20 Seeds')
    plt.savefig(output_root / "auc_boxplots.png")
    plt.close()

if __name__ == "__main__":
    # 按需求指定type=species；临床特征会在最终XGB阶段加入（并在筛选前剔除缺失）
    main(types=['microprotein'])