import numpy as np
import matplotlib.pyplot as plt
import japanize_matplotlib

# ファイルの読み込み
data = np.loadtxt('//home/shishida/SMAC/.d/phi(x = 0.5).d')

# 列を分離
y_data = data[:, 0]
phi_data = data[:, 1]

# 理論式 φ = 0.25 * y * (1 - y)
y_theory = np.linspace(0, 1, 200)
phi_theory = 0.25 * y_theory * (1 - y_theory)

# グラフ描画
plt.figure(figsize=(7,5))
plt.plot(y_data, phi_data, 'o', label='データ点(x=0.5)', markersize=5)
plt.plot(y_theory, phi_theory, '-', label=r'$\phi = 0.25\,y(1 - y)$', linewidth=2)

plt.xlabel('y')
plt.ylabel('φ')
plt.title('データ点と理論式の比較')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
