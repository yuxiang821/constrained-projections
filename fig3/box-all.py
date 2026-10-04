import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.ticker import AutoMinorLocator

# 全局设置字体
plt.rcParams['font.family'] = 'Helvetica'  # 如果Times New Roman可用，可改为'Times New Roman'


def box():
    # 读取数据
    f1 = xr.open_dataset('E:/north-future/ec/new-output/box/raw/all/hot-all.nc')
    f2 = xr.open_dataset('E:/north-future/ec/new-output/box/c/all/hot-all.nc')
    f3 = xr.open_dataset('E:/north-future/ec/new-output/box/raw/all/cold-all.nc')
    f4 = xr.open_dataset('E:/north-future/ec/new-output/box/c/all/cold-all.nc')
    f5 = xr.open_dataset('E:/north-future/ec/new-output/box/raw/all/pr-all.nc')
    f6 = xr.open_dataset('E:/north-future/ec/new-output/box/c/all/pr-all.nc')

    # 提取第二行数据（索引为1）
    t1 = f1['raw'].values[0, :]  # 提取热极端第二行
    t2 = f2['EC'].values[0, :]  # 提取约束后热极端第二行
    t3 = f3['raw'].values[0, :]  # 提取冷极端第二行
    t4 = f4['EC'].values[0, :]  # 提取约束后冷极端第二行
    t5 = f5['raw'].values[1, :]  # 提取降水第二行
    t6 = f6['EC'].values[1, :]  # 提取约束后降水第二行
    print(t5)
    labels = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))  # 增大图形尺寸

    width = 0.4
    spacing = 0.8
    positions = np.arange(1, 9, 2)

    # 定义箱型图属性
    boxprops1 = dict(facecolor='#ef8a62', color='#ef8a62')
    boxprops2 = dict(facecolor='#92c5de', color='#92c5de')
    medianprops = dict(color='black', linewidth=2)

    # 子图1: 热极端
    box1 = axes[0].boxplot(t1.T, patch_artist=True, labels=labels, positions=positions, widths=width, vert=True,
                           boxprops=boxprops1, medianprops=medianprops, whiskerprops=dict(linestyle='--'),
                           capprops=dict(linewidth=1.5), flierprops=dict(marker='o', markersize=5, alpha=0.5))
    box2 = axes[0].boxplot(t2.T, patch_artist=True, positions=positions + spacing, widths=width, vert=True,
                           boxprops=boxprops2, medianprops=medianprops, whiskerprops=dict(linestyle='--'),
                           capprops=dict(linewidth=1.5), flierprops=dict(marker='o', markersize=5, alpha=0.5))
    axes[0].legend([box1["boxes"][0], box2["boxes"][0]], ['Raw Projection', 'Constrained Projection'],
                   loc='upper center', bbox_to_anchor=(0.5, 1.15), fontsize=12, frameon=False)
    axes[0].set_title('(a) TXx', fontsize=14, loc='left', pad=10)
    axes[0].yaxis.set_minor_locator(AutoMinorLocator())

    # 子图2: 冷极端
    box1 = axes[1].boxplot(t3.T, patch_artist=True, labels=labels, positions=positions, widths=width, vert=True,
                           boxprops=boxprops1, medianprops=medianprops, whiskerprops=dict(linestyle='--'),
                           capprops=dict(linewidth=1.5), flierprops=dict(marker='o', markersize=5, alpha=0.5))
    box2 = axes[1].boxplot(t4.T, patch_artist=True, positions=positions + spacing, widths=width, vert=True,
                           boxprops=boxprops2, medianprops=medianprops, whiskerprops=dict(linestyle='--'),
                           capprops=dict(linewidth=1.5), flierprops=dict(marker='o', markersize=5, alpha=0.5))
    axes[1].set_title('(b) TNn', fontsize=14, loc='left', pad=10)
    axes[1].yaxis.set_minor_locator(AutoMinorLocator())

    # 子图3: 降水
    box1 = axes[2].boxplot(t5.T, patch_artist=True, labels=labels, positions=positions, widths=width, vert=True,
                           boxprops=boxprops1, medianprops=medianprops, whiskerprops=dict(linestyle='--'),
                           capprops=dict(linewidth=1.5), flierprops=dict(marker='o', markersize=5, alpha=0.5))
    box2 = axes[2].boxplot(t6.T, patch_artist=True, positions=positions + spacing, widths=width, vert=True,
                           boxprops=boxprops2, medianprops=medianprops, whiskerprops=dict(linestyle='--'),
                           capprops=dict(linewidth=1.5), flierprops=dict(marker='o', markersize=5, alpha=0.5))
    axes[2].set_title('(c) RX5day', fontsize=14, loc='left', pad=10)
    axes[2].yaxis.set_minor_locator(AutoMinorLocator())

     # ==================== 设置刻度 ====================
    for ax in axes:
        # X 轴刻度标签
        ax.set_xticks(positions + spacing / 2)
        ax.set_xticklabels(labels, fontsize=14)

        # Y 轴刻度标签
        ax.tick_params(axis='y', labelsize=14, length=7, width=0.8)

        ax.tick_params(axis='both', which='major',
                       direction='out', length=7, width=0.8,
                       top=True, right=True)

        ax.tick_params(axis='both', which='minor',
                       direction='out', length=3.5, width=0.6,
                       top=True, right=True)


    # for label in ax.get_yticklabels():
    #     # label.set_fontweight('bold')
    #     label.set_fontsize(14)

    plt.tight_layout()
    # plt.savefig('E:/north-future/new-p/original/fig5/box-all.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.show()


box()