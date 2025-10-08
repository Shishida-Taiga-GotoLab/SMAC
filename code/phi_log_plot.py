import numpy as np
import matplotlib.pyplot as plt

# 格子点数 Nx（与えられたデータ）
Nx = np.array([  8,   16,   32,   64,  128 ])

# 誤差（L2ノルム）
error = np.array([
    0.92250397E-01,
    0.27122280E-01,
    0.73149612E-02,
    0.18987094E-02,
    0.48372348E-03
])

# 格子幅 h ≈ 1/(Nx-1)
h = 1.0 / (Nx - 1)

# 収束次数 p を最小二乗で推定（log-log直線フィット）
logh = np.log(h)
loge = np.log(error)
p, logC = np.polyfit(logh, loge, 1)  # 1次フィット
C = np.exp(logC)

print(f"推定された収束次数 p ≈ {p:.4f}")
print(f"近似式： error ≈ {C:.4e} * h^{p:.2f}")

# --- プロット ---
plt.figure(figsize=(6,5))
plt.loglog(h, error, 'o-', label='数値データ', basex=10, basey=10)

# フィット直線
h_fit = np.linspace(h.min(), h.max(), 100)
error_fit = C * h_fit**p
plt.loglog(h_fit, error_fit, '--', label=f'フィット: p ≈ {p:.2f}')

# 見た目の調整
plt.gca().invert_xaxis()  # h が小さい方を右にする
plt.xlabel('h (格子幅)')
plt.ylabel('L2誤差')
plt.title('誤差収束の log–log プロット')
plt.grid(True, which='both', ls=':')
plt.legend()
plt.tight_layout()
plt.show()
