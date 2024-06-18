import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

vu = np.loadtxt("../los_full/decompose/vu_3_for_landcover.xyz", usecols=2)
land = np.loadtxt("/Users/qi/leeds/datasets/rasters/landcover_PROBAV_LC100_global_v3.0.1_2019_classification_0.005.xyz", usecols=2)
df = pd.DataFrame({'vu': vu, 'land': land})
df.dropna(inplace=True)
df = df[df["land"] > 0]
df = df[df["land"] < 100]

# Create an array with the colors you want to use
colors = ["#ffbb22",
"#ffff4c",
"#f096ff",
"#fa0000",
"#b4b4b4",
"#f0f0f0",
"#0032c8",
"#0096a0"]

# Set your custom color palette
sns.set_palette(sns.color_palette(colors))

# extract mean and std by landcover
grouped_data = df.groupby("land")["vu"]
means = grouped_data.mean()
std_devs = grouped_data.std()

# fig, ax = plt.subplots(figsize=(4, 6)) # 9
fig, ax = plt.subplots(figsize=(3, 6)) # 9

sns.violinplot(data=df, x="vu", y="land", ax=ax, orient="h")
for i, id in enumerate(np.arange(20,100,10)):
    plt.text(-6.8, i+0.1, "{:.1f}%".format(
        len(df[df["land"] == id]) / len(df) * 100),
        horizontalalignment='left',
        weight='semibold')
    plt.text(6.9, i+0.1, "{:.1f} ± {:.1f}".format(
        means.iloc[i], std_devs.iloc[i]),
        horizontalalignment='right',
        weight='semibold')

plt.xlim(-7, 7)
plt.xlabel('Vu, mm/yr')
ax.tick_params(labelleft=False,labelright=True)
plt.ylabel('')
y_ticks_labels = ['Shrubs', 'Herb', 'Crops', 'Urban', 'Bare',
                  'Snow/Ice', 'Water', 'Wetland']
plt.yticks(ticks=np.arange(9), labels=y_ticks_labels) #, rotation=45
plt.vlines(0, -0.5, 7.5, linestyles=':')
plt.ylim(-0.5, 7.5)
plt.tight_layout()
plt.show()
fig.savefig('../los_full/decompose/land_cover_vu_correlation.png', format='PNG', dpi=300,
            bbox_inches='tight')


# for i in np.arange(100,127,1):
# # for i in np.arange(20,112,10):
#     perc = len(nonnan_land[nonnan_land == i]) / len(nonnan_land) * 100
#     if perc > 0.01:
#         print(i,"{:.3f}%".format(perc))

# 20 2.368%
# 30 62.706%
# 40 8.515%
# 50 0.740%
# 60 21.035%
# 70 1.288%
# 80 0.576%
# 90 0.716%
# 111 0.409%
# 113 0.131%
# 116 0.965%
# 121 0.022%
# 126 0.524%


