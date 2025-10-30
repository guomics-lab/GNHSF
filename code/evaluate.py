# -*- coding: gbk -*-
import os
import pandas as pd
from sklearn.metrics import roc_auc_score

# �������Ԥ�����ļ���·��
input_folder = "/sunyingying/test/glm_xgb_gain_freq_42_output_microprotein/internal_preds/"
output_file = os.path.join(input_folder, "seed_auc_summary.csv")

# ��Ž��
results = []

# �����ļ���������csv�ļ�
for fname in os.listdir(input_folder):
    if fname.endswith(".csv") and "seed" in fname:
        file_path = os.path.join(input_folder, fname)
        df = pd.read_csv(file_path)
        # ����ļ����Ƿ���seed�У�����ж�seed��ϣ�����飩
        if "seed" in df.columns and df["seed"].nunique() > 1:
            for seed, subdf in df.groupby("seed"):
                auc = roc_auc_score(subdf["y_true"], subdf["y_pred"])
                results.append({"seed": seed, "auc": auc})
        else:
            # seed��Ϊ�ļ�����һ����
            seed = df["seed"].iloc[0] if "seed" in df.columns else fname.split("seed_")[-1].split(".")[0]
            auc = roc_auc_score(df["y_true"], df["y_pred"])
            results.append({"seed": seed, "auc": auc})

# ���浽csv
result_df = pd.DataFrame(results)
result_df.to_csv(output_file, index=False)
print(f"AUC����ѱ��浽: {output_file}")