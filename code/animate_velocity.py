import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import re

# === 設定 ===
# 速度ファイルのパターン（例: uv_00010.d, uv_00020.d, ...）
file_pattern = "/home/shishida/SMAC/.d/uv_*.d"
# ファイル一覧をソート（ステップ番号順）
files = sorted(glob.glob(file_pattern),
               key=lambda x: int(re.findall(r'\d+', x)[-1]))

# 10ステップごとに抽出
files = files[::1]  # 例: 10ステップごとにする場合は [::10]

# === データ読み込み関数 ===
def load_velocity(filename):
    data = np.loadtxt(filename)
    i = data[:, 0].astype(int)
    j = data[:, 1].astype(int)
    u = data[:, 2]
    v = data[:, 3]

    Nx = i.max()
    Ny = j.max()

    # 格子状に並び替え
    U = np.zeros((Ny, Nx))
    V = np.zeros((Ny, Nx))
    for k in range(len(i)):
        U[j[k]-1, i[k]-1] = u[k]
        V[j[k]-1, i[k]-1] = v[k]

    return U, V

# === 最初のファイルから格子を取得 ===
U0, V0 = load_velocity(files[0])
Ny, Nx = U0.shape
x = np.arange(1, Nx + 1)
y = np.arange(1, Ny + 1)
X, Y = np.meshgrid(x, y)

# === 図の設定 ===
fig, ax = plt.subplots(figsize=(6, 6))
qv = ax.quiver(X, Y, U0, V0, scale=10)
ax.set_title(f"Velocity field: {files[0]}")
ax.set_xlabel("i")
ax.set_ylabel("j")

# === アニメーション更新関数 ===
def update(frame):
    U, V = load_velocity(files[frame])
    qv.set_UVC(U, V)
    ax.set_title(f"Velocity field: {files[frame]}")
    return qv,

ani = animation.FuncAnimation(fig, update, frames=len(files), interval=30, blit=False)

plt.tight_layout()
plt.show()
