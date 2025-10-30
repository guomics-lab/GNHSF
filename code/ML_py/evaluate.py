# -*- coding: gbk -*-
import os
import pandas as pd
from sklearn.metrics import roc_auc_score

# 设置你的预测结果文件夹路径
input_folder = "/sunyingying/test/glm_xgb_gain_freq_42_output_microprotein/internal_preds/"
output_file = os.path.join(input_folder, "seed_auc_summary.csv")

# 存放结果
results = []

# 遍历文件夹下所有csv文件
for fname in os.listdir(input_folder):
    if fname.endswith(".csv") and "seed" in fname:
        file_path = os.path.join(input_folder, fname)
        df = pd.read_csv(file_path)
        # 检查文件内是否有seed列（如果有多seed混合，需分组）
        if "seed" in df.columns and df["seed"].nunique() > 1:
            for seed, subdf in df.groupby("seed"):
                auc = roc_auc_score(subdf["y_true"], subdf["y_pred"])
                results.append({"seed": seed, "auc": auc})
        else:
            # seed作为文件名的一部分
            seed = df["seed"].iloc[0] if "seed" in df.columns else fname.split("seed_")[-1].split(".")[0]
            auc = roc_auc_score(df["y_true"], df["y_pred"])
            results.append({"seed": seed, "auc": auc})

# 保存到csv
result_df = pd.DataFrame(results)
result_df.to_csv(output_file, index=False)
print(f"AUC结果已保存到: {output_file}")