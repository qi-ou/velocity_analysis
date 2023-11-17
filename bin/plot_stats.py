# importing libraries
import matplotlib.pyplot as plt
import pandas as pd

# file = "~/leeds/projects/insar/tianshan_processing/frame_epoch_stats.csv"
file = "../frame_epoch_stats.csv"
df = pd.read_csv(file)
df_stats = df.describe().round(1)

fig, ax = plt.subplots(1,4, sharey='all', figsize=(10, 3))
ax[1].hist(df['initial_epoch'], bins=range(0, 300, 20), label="Start: {:.0f} ± {:.0f}".format(df_stats["initial_epoch"].loc["mean"], df_stats["initial_epoch"].loc["std"]), alpha=.5)
ax[1].hist(df['final_epoch'], bins=range(0, 300, 20), label='Final: {:.0f} ± {:.0f}'.format(df_stats["final_epoch"].loc["mean"], df_stats["final_epoch"].loc["std"]), alpha=.5) #, color='red'
ax[1].set_xlabel("Number of Epochs")
ax[1].set_title("Epoch")
ax[1].legend()

ax[0].hist(df['initial_ifg'], bins=range(0, 3000, 200), label="Start: {:.0f} ± {:.0f}".format(df_stats["initial_ifg"].loc["mean"], df_stats["initial_ifg"].loc["std"]), alpha=.5)
ax[0].hist(df['final_ifg'], bins=range(0, 3000, 200), label='Final: {:.0f} ± {:.0f}'.format(df_stats["final_ifg"].loc["mean"], df_stats["final_ifg"].loc["std"]), alpha=.5) #, color='red'
ax[0].set_xlabel("Number of IFGs")
ax[0].set_title("IFG")
ax[0].legend()
ax[0].set_ylabel("Number of Frames")


ax[2].hist(df['Initial_years'].astype(float), label="Start: {:.1f} ± {:.1f}".format(df_stats["Initial_years"].loc["mean"], df_stats["Initial_years"].loc["std"]), alpha=.5)
ax[2].hist(df['Final_years'].astype(float), label='Final: {:.1f} ± {:.1f}'.format(df_stats["Final_years"].loc["mean"], df_stats["Final_years"].loc["std"]), alpha=.5) #, color='red'
ax[2].set_xlabel("Years")
ax[2].set_title("Duration")
ax[2].legend()

ax[3].hist(df['initial_mean_bt'].astype(float), bins=range(50, 500, 24), label="Start: {:.0f} ± {:.0f}".format(df_stats["initial_mean_bt"].loc["mean"], df_stats["initial_mean_bt"].loc["std"]), alpha=.5)
ax[3].hist(df['final_mean_bt'].astype(float), bins=range(50, 500, 24), label='Final: {:.0f} ± {:.0f}'.format(df_stats["final_mean_bt"].loc["mean"], df_stats["final_mean_bt"].loc["std"]), alpha=.5) #, color='red'
ax[3].set_xlabel("Days")
ax[3].set_xlim([50, 500])
ax[3].set_title("Avg. Temp. Baseline")
ax[3].legend()

# ax[4].hist(df['initial_mean_avg_coh'].astype(float), label="Start: {:.1f} ± {:.1f}".format(df_stats["initial_mean_avg_coh"].loc["mean"], df_stats["initial_mean_avg_coh"].loc["std"]),  alpha=.5)
# ax[4].hist(df['final_mean_avg_coh'].astype(float), label='Final: {:.1f} ± {:.1f}'.format(df_stats["final_mean_avg_coh"].loc["mean"], df_stats["final_mean_avg_coh"].loc["std"]), alpha=.5) #, color='red'
# ax[4].set_xlabel("Coherence")
# ax[4].set_title("Avg_Coh")
# ax[4].legend()

# Showing the plot using plt.show()
plt.tight_layout()
fig.savefig('../frame_stats.png', format='PNG', dpi=300, bbox_inches='tight')
plt.show()
