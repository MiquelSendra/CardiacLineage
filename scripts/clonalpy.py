import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import patches
from matplotlib.colors import to_rgb
import os

def create_colormap(color):
    """
    Create a custom colormap where the first 30% is a solid color 
    that is 30% between white and the target color,
    and then transition smoothly from this 30% color to the full target color over the remaining 70%.
    """
    target_rgb = np.array(to_rgb(color))
    white = np.array([1, 1, 1])
    mid_rgb = white + 0.3 * (target_rgb - white)
    
    cdict = {
        'red': [(0.0, mid_rgb[0], mid_rgb[0]), (0.3, mid_rgb[0], mid_rgb[0]), (1.0, target_rgb[0], target_rgb[0])],
        'green': [(0.0, mid_rgb[1], mid_rgb[1]), (0.3, mid_rgb[1], mid_rgb[1]), (1.0, target_rgb[1], target_rgb[1])],
        'blue': [(0.0, mid_rgb[2], mid_rgb[2]), (0.3, mid_rgb[2], mid_rgb[2]), (1.0, target_rgb[2], target_rgb[2])]
    }
    
    return LinearSegmentedColormap(f'custom_colormap_{color}', cdict)

def contribution_heatmap(clusters_summary_norm, subtypes_df, subtypes_columns, subtype_colors, df, colors, figsize=(10, 10), out = '../figures'):
    """
    Plot heatmaps for the data.
    """

    colormaps = [create_colormap(color) for color in colors]
    fig, axes = plt.subplots(1, len(clusters_summary_norm.columns) + 2, figsize=figsize,
                             gridspec_kw={'width_ratios': [1] * len(clusters_summary_norm.columns) + [1/3, 1/3], 'wspace': 0.0})

    mask = clusters_summary_norm.isna()

    for idx, col in enumerate(clusters_summary_norm.columns):
        ax = axes[idx]
        sns.heatmap(clusters_summary_norm[[col]], fmt='g', linewidth=0.8, cmap=colormaps[idx % len(colormaps)], 
                    mask=mask[[col]], cbar=False, ax=ax)
        
        ax.set_xticks([]); ax.set_xticklabels([]); ax.set_xlabel('')
        ax.set_yticks([]); ax.set_ylabel(''); ax.set_title('')
        
        if idx == 0:
            ax.set_yticks(np.arange(0, len(clusters_summary_norm)) + 0.5)
            ax.set_yticklabels(clusters_summary_norm.index)
            ax.tick_params(axis='y', labelrotation=0)

    EmM_ax = axes[clusters_summary_norm.columns.get_loc('EmM')]
    for y, (row_index, subtypes) in enumerate(subtypes_df.iterrows()):
        for subtype, is_present in subtypes.items():
            if is_present > 0:
                rect = patches.Rectangle(
                    (subtypes_columns.index(subtype) / len(subtypes_columns), y), 1 / len(subtypes_columns), 1,
                    linewidth=1.5, edgecolor=subtype_colors[subtype], facecolor='none'
                )
                EmM_ax.add_patch(rect)

    bipotent_cmap = LinearSegmentedColormap.from_list('bipotent_cmap', ['white', 'lightblue'], N=2)
    sns.heatmap(df[['bipotent']], fmt='g', linewidth=0.8, cmap=bipotent_cmap, cbar=False, ax=axes[-2])
    axes[-2].set_xticks([]); axes[-2].set_yticks([]); axes[-2].set_xlabel(''); axes[-2].set_ylabel(''); axes[-2].set_title('')

    clusters_bilateral = df['bilateral'].replace({'YES': 1, 'NO': 0, '?': 0}).astype(int).to_frame()
    sns.heatmap(clusters_bilateral, fmt='g', linewidth=0.8, cmap=['white', 'darkorange'], cbar=False, ax=axes[-1])
    axes[-1].set_xticks([]); axes[-1].set_yticks([]); axes[-1].set_xlabel(''); axes[-1].set_ylabel(''); axes[-1].set_title('')

    for y, total_cells in enumerate(df['total']):
        axes[-1].text(1.35, y + 0.5, total_cells, verticalalignment='center', horizontalalignment='center', fontsize=12)

    for y in range(len(clusters_summary_norm)):
        axes[0].axhline(y + 1, color='dimgray', linewidth=0.8)

    plt.tight_layout()
    plt.savefig(os.path.join(out, "Tamoxifen_clone_contribution_heatmap.png"))
    plt.savefig(os.path.join(out, "Tamoxifen_clone_contribution_heatmap.svg"))



if __name__ == "__main__":
    # To test the functions if needed
    pass
