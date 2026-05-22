import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t
from matplotlib.table import Table
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import sys

dir_path = os.getcwd()
csv_file_name = sys.argv[1]
file_path = os.path.join(dir_path, csv_file_name)

df = pd.read_csv(file_path)
df_filtered = df[~df.apply(lambda row: row.astype(str).str.contains("dials.show").any(), axis=1)]

df_filtered['Mean_Intensity'] = pd.to_numeric(df_filtered['Mean_Intensity'], errors='coerce')
df_filtered['Dose'] = pd.to_numeric(df_filtered['Dose'], errors='coerce')
print("\nDataframe info after converting to numeric:")
df_filtered.info()

df_filtered = df_filtered.dropna(subset=['Dose', 'Mean_Intensity', 'Group-ID'])

def normalize_safely(x):
    if (x.max() - x.min()) > 0:
        return (x - x.min()) / (x.max() - x.min())
    else:
        return np.nan
    
df_filtered['normalised_intensity'] = df_filtered.groupby('Group-ID')['Mean_Intensity'].transform(normalize_safely) #lambda x: (x - x.min()) / (x.max() - x.min())

print("\nDataFrame info post normalization")
df_filtered.info
print(df_filtered.describe())
print(df_filtered.head(10))

df_filtered = df_filtered.dropna(subset=['normalised_intensity'])
def exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

#plotting func
groups = df_filtered.groupby('Group-ID')

cmap = cm.get_cmap('tab20', len(groups))
plt.figure(figsize=(12, 12))

table_data = [["Group", "a", "b", "c", "R²"]]

i=0
for name, group in groups:
    print(f"\nprocessing group: {name}")

    try:
        x = group['Dose'].values.astype(float)
        y = group['normalised_intensity'].values.astype(float)

        if len(x) < 3:
            print(f"Skipping group {group} due to insufficient valid data.")
            continue

        p0 = (1.0, 0.1, 0.0)    
        popt, pcov = curve_fit(exp_decay, x, y, p0=p0, maxfev=10000)
        x_fit = np.linspace(x.min(), x.max(), 100)
        y_fit = exp_decay(x_fit, *popt)

        #confidence intervals
        alpha = 0.05
        n = len(y)
        p = len(popt)
        dof = max(0, n - p)
        tval = t.ppf(1.0 - alpha / 2., dof)
        
        plt.scatter(x, y, label=name, color=cmap(i))
        plt.plot(x_fit, y_fit, color=cmap(i))

        #calc r²
        residuals = y - exp_decay(x, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        if ss_tot > 0:
            r_squared = 1 - (ss_res / ss_tot)
        else:
            r_squared = 1.0

        table_data.append([name] + [f"{param:.3f}" for param in popt] + [f"{r_squared:.3f}"])
        print(f"successfully fitted group: {name}")
    except Exception as e:
        print(f"Could not fit group {name}: {e}")
    #plt.scatter(group['Dose'], group['normalised_intensity'].astype(float), label=name, color=cmap(i))
    i += 1

plt.xlabel('Dose (MGy)')
plt.ylabel('Mean Intensity (normalised - Profile fitted)')
plt.title('Intensity decay per unit dose')

handles, labels = plt.gca().get_legend_handles_labels()
if handles:
    plt.legend(title='Group ID')
else:
    print("\nNo data to plot, legend will not be shown")

if len(table_data) > 1:
    the_table = plt.table(cellText=table_data,
                          colLabels=None,
                          cellLoc='center',
                          loc='bottom',
                          bbox=[0, -0.7, 1, 0.5])
    the_table.auto_set_font_size=False
    the_table.set_fontsize(10)
    the_table.scale(1, 1.5)
else:
    print("\nNo data to create table")

plt.subplots_adjust(left=0.2, bottom=0.2)
plt.grid(True)
plt.tight_layout()
plt.show()

### need to average on crystal size, plot error bands