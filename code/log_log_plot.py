import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ===== ファイル名を指定 =====
file_path = "/home/shishida/SMAC/.d/err_convection_v2.d"  # 実際のファイル名に合わせて変更

# ===== データ読み込み =====
data = pd.read_csv(file_path, delim_whitespace=True, header=None)
data.columns = ["N", "RMSE", "extra"][:len(data.columns)]

# ===== 対数を取る =====
x = np.log10(data["N"])
y = np.log10(data["RMSE"])

# ===== 最小二乗法による直線フィッティング =====
coeff = np.polyfit(x, y, 1)  # y = a*x + b
slope, intercept = coeff
fit_y = np.polyval(coeff, x)

# ===== グラフ描画 =====
plt.figure(figsize=(6, 4))
plt.loglog(data["N"], data["RMSE"], 'o', label="RMSE data")
plt.loglog(data["N"], 10**fit_y, '--', label=f"Fit: slope = {slope:.2f}")
plt.xlabel("Grid points")
plt.ylabel("RMSE")
plt.title("convection error (log-log scale)")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()

# ===== 傾きを出力 =====
print(f"近似直線の傾き（収束次数）: {slope:.3f}")
