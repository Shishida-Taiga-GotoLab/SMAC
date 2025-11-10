import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# === 設定 ===
# 出力ファイルの保存場所（あなたの環境に合わせて変更）
data_dir = "/home/shishida/SMAC/.d"
step = 100   # 表示したいステップ番号
u_file = os.path.join(data_dir, f"u_{step:05d}.d")
v_file = os.path.join(data_dir, f"v_{step:05d}.d")

# === データ読み込み関数 ===
def read_field(filename):
    """
    Fortran出力ファイル (i, j, x, y, value) を読み込み、
    x, y, value の2次元配列を返す。
    """
    data = np.loadtxt(filename)
    i = data[:, 0].astype(int)
    j = data[:, 1].astype(int)
    x = data[:, 2]
    y = data[:, 3]
    val = data[:, 4]

    Nx = i.max()
    Ny = j.max()
    # 格子点ごとの配列に整形
    X = np.zeros((Nx + 1, Ny + 1))
    Y = np.zeros((Nx + 1, Ny + 1))
    F = np.zeros((Nx + 1, Ny + 1))

    for ii, jj, xx, yy, vv in zip(i, j, x, y, val):
        X[ii, jj] = xx
        Y[ii, jj] = yy
        F[ii, jj] = vv

    return X, Y, F

# === 読み込み ===
print(f"Loading: {u_file}")
X_u, Y_u, U = read_field(u_file)
print(f"Loading: {v_file}")
X_v, Y_v, V = read_field(v_file)

# 格子の中心をそろえる（平均位置でベクトルを描く）
X = 0.5 * (X_u[:-1, :-1] + X_v[1:, 1:])
Y = 0.5 * (Y_u[:-1, :-1] + Y_v[1:, 1:])
U_c = 0.5 * (U[:-1, :-1] + U[1:, 1:])
V_c = 0.5 * (V[:-1, :-1] + V[1:, 1:])

# === 可視化 ===
plt.figure(figsize=(6, 3))
plt.quiver(X, Y, U_c, V_c, scale=2.0, pivot='mid')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f"Velocity field (step = {step})")
plt.axis('equal')
plt.tight_layout()
plt.show()
