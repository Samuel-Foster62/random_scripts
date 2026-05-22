import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# 1. Load CSV
# ---------------------------------------------------------
df = pd.read_csv("/home/iuh23949/Foster_S/scripts/figure_making/Oxy_vs_deoxy_Mb.csv", header=0)
df = df.apply(pd.to_numeric, errors="coerce")
df = df.dropna()
df.columns = ["x", "group2", "group3", "group4", "group5"]

wl = df["x"].values

# Masks
m_low  = wl < 470
m_mid  = (wl >= 470) & (wl <= 700)
m_high = wl > 700

# ---------------------------------------------------------
# 2. Load CIE 1931 colour matching functions
# ---------------------------------------------------------
cmf = pd.read_csv("/home/iuh23949/Foster_S/scripts/figure_making/CIE_xyz_1931_2deg.csv", header=0)  # wavelength, x_bar, y_bar, z_bar
cmf.columns = ["wavelength", "x_bar", "y_bar", "z_bar"]
cmf = cmf.apply(pd.to_numeric, errors="coerce").dropna()

cmf_w = cmf["wavelength"].values
cmf_x = cmf["x_bar"].values
cmf_y = cmf["y_bar"].values
cmf_z = cmf["z_bar"].values

# Clip wavelengths to CMF range
wl_clipped = np.clip(wl, cmf_w.min(), cmf_w.max())

# Interpolate
xbar = np.interp(wl_clipped, cmf_w, cmf_x)
ybar = np.interp(wl_clipped, cmf_w, cmf_y)
zbar = np.interp(wl_clipped, cmf_w, cmf_z)

# ---------------------------------------------------------
# 3. Convert absorbance → transmission
# ---------------------------------------------------------
def absorbance_to_transmission(A, k=10.0):
    return 10 ** (-k * A)

path_factor = 1000.0

T_oxy = absorbance_to_transmission(df["group2"].values, k=path_factor)
T_deoxy = absorbance_to_transmission(df["group3"].values, k=path_factor)

# ---------------------------------------------------------
# 4. Convert spectrum → XYZ
# ---------------------------------------------------------
def spectrum_to_xyz(T, xbar, ybar, zbar):
    X = np.sum(T * xbar)
    Y = np.sum(T * ybar)
    Z = np.sum(T * zbar)
    return np.array([X, Y, Z])

XYZ_oxy = spectrum_to_xyz(T_oxy, xbar, ybar, zbar)
XYZ_deoxy = spectrum_to_xyz(T_deoxy, xbar, ybar, zbar)

XYZ_stack = np.vstack([XYZ_oxy, XYZ_deoxy])
XYZ_stack /= np.max(XYZ_stack)

XYZ_oxy, XYZ_deoxy = XYZ_stack

# ---------------------------------------------------------
# 5. Convert XYZ → sRGB
# ---------------------------------------------------------
def xyz_to_srgb(XYZ):
    M = np.array([
        [ 3.2406, -1.5372, -0.4986],
        [-0.9689,  1.8758,  0.0415],
        [ 0.0557, -0.2040,  1.0570]
    ])

    rgb = np.dot(M, XYZ)
    rgb = np.nan_to_num(rgb)
    rgb = np.clip(rgb, 0, None)
    rgb = np.where(rgb <= 0.0031308,
                   12.92 * rgb,
                   1.055 * (rgb ** (1/2.4)) - 0.055)
    return np.clip(rgb, 0, 1)

oxy_color = tuple(xyz_to_srgb(XYZ_oxy))
deoxy_color = tuple(xyz_to_srgb(XYZ_deoxy)
)

print("Computed oxy-Hb colour:", oxy_color)
print("Computed deoxy-Hb colour:", deoxy_color)

# ---------------------------------------------------------
# 6. Plot
# ---------------------------------------------------------

sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))

# Oxy: group 2 outside 470–700, group 4 inside
sns.lineplot(x=wl[m_low],  y=df["group2"][m_low],  color=oxy_color,   label="Oxy-hemoglobin")
sns.lineplot(x=wl[m_mid],  y=df["group4"][m_mid],  color=oxy_color,   legend=False)
sns.lineplot(x=wl[m_high], y=df["group2"][m_high], color=oxy_color,   legend=False)

# Deoxy: group 3 outside 470–700, group 5 inside
sns.lineplot(x=wl[m_low],  y=df["group3"][m_low],  color=deoxy_color, label="Deoxy-hemoglobin")
sns.lineplot(x=wl[m_mid],  y=df["group5"][m_mid],  color=deoxy_color, legend=False)
sns.lineplot(x=wl[m_high], y=df["group3"][m_high], color=deoxy_color, legend=False)

plt.xlabel("Wavelength (nm)")
plt.ylabel("Absorbance")
plt.title("Oxy/Deoxy Hemoglobin with Region-Specific Scaling")
plt.legend()
plt.tight_layout()
plt.show()
