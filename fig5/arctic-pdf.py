import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import AutoMinorLocator
import matplotlib  # Import matplotlib explicitly


def calculate_weighted_pdf(data, lat, inter_model_std=None):
    """计算面积加权的PDF，并评估显著增加的面积百分比"""
    # 确保数据是2D的
    if len(data.shape) != 2:
        raise ValueError("数据必须是2D数组 (纬度, 经度)")

    # 计算纬度权重（面积加权）
    weights = np.cos(np.deg2rad(lat))

    # 扩展权重到与数据相同的形状
    weights_expanded = np.outer(weights, np.ones(data.shape[1]))

    # 创建有效数据掩膜（非缺测值）
    mask = ~np.isnan(data)

    # 提取有效数据和权重
    valid_data = data[mask]
    valid_weights = weights_expanded[mask]

    # Handle cases where valid_data might be empty
    if valid_data.size == 0:
        # Return empty arrays and 0% if no valid data
        return np.array([]), np.array([]), 0.0

    # 标准化权重，使其总和为1
    valid_weights = valid_weights / np.sum(valid_weights)

    # 计算显著增加的面积百分比
    significant_area_percent = 0
    if inter_model_std is not None:
        # 修改为使用2倍标准差作为显著性阈值
        significant_mask = valid_data > inter_model_std
        # 计算相对于总面积的比例
        significant_area_percent = np.sum(valid_weights[significant_mask]) * 100

    # Use KDE to calculate PDF
    # Ensure there's enough data for KDE, otherwise it might fail
    if valid_data.size < 2:  # gaussian_kde needs at least 2 points
        # If not enough data, return empty arrays
        return np.array([]), np.array([]), significant_area_percent

    kde = stats.gaussian_kde(valid_data, weights=valid_weights, bw_method='silverman')

    # Generate evaluation points
    # Add a small buffer to the range to ensure KDE covers the tails
    data_min = np.min(valid_data)
    data_max = np.max(valid_data)
    x_range = data_max - data_min
    # Ensure x_range is not zero for single-point distributions
    if x_range == 0:
        x_range = 1.0  # Use a default range if all data points are identical

    x = np.linspace(data_min - 0.1 * x_range, data_max + 0.1 * x_range, 1000)

    # Ensure x is not empty after linspace, though unlikely with 1000 points
    if x.size == 0:
        return np.array([]), np.array([]), significant_area_percent

    y = kde(x)

    return x, y, significant_area_percent


def plot_separated_indices(data_before_list, data_after_list, lat,
                           std_list, indices_names, units_list,
                           output_file=None):

    # 直接设置绘图风格，而不是调用单独的函数
    plt.style.use('default')  # 重置风格
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 8,

        # 标题和标签
        'axes.labelsize': 8,  # 轴标签字体大小
        'axes.titlesize': 9,  # 标题字体大小
        'xtick.labelsize': 8,  # x轴刻度标签字体大小
        'ytick.labelsize': 8,  # y轴刻度标签字体大小
        'legend.fontsize': 8,  # 图例字体大小

        # 线条设置
        'lines.linewidth': 1.0,  # 线条宽度
        'axes.linewidth': 0.5,  # 轴线宽度
        'xtick.major.width': 0.5,  # x轴主刻度宽度
        'ytick.major.width': 0.5,  # y轴主刻度宽度
        'xtick.minor.width': 0.5,  # x轴次刻度宽度
        'ytick.minor.width': 0.5,  # y轴次刻度宽度

        # 刻度长度
        'xtick.major.size': 4,  # x轴主刻度长度
        'ytick.major.size': 4,  # y轴主刻度长度
        'xtick.minor.size': 2,  # x轴次刻度长度
        'ytick.minor.size': 2,  # y轴次刻度长度
    })

    # 创建图形 - 使用Nature双栏宽度（转换为英寸：183mm约等于7.2英寸）
    fig_width = 15  # 英寸
    fig_height = 4  # 比例保持约3:1
    fig, axes = plt.subplots(1, 3, figsize=(fig_width, fig_height))

    # 颜色方案 - 调整为约束前用红色，约束后用蓝色
    red_color = '#b2182b'  # 深红色（约束前）
    blue_color = '#2166ac'  # 深蓝色（约束后）
    red_fill = '#ef8a62'  # 浅红色填充
    blue_fill = '#92c5de'  # 浅蓝色填充

    # 绘制每个子图
    for i, (data_before, data_after, std, name, units) in enumerate(
            zip(data_before_list, data_after_list, std_list, indices_names, units_list)):

        # 计算约束前的PDF和显著增加区域
        x_before, y_before, sig_area_before = calculate_weighted_pdf(
            data_before, lat, inter_model_std=std
        )

        # 计算约束后的PDF和显著增加区域
        x_after, y_after, sig_area_after = calculate_weighted_pdf(
            data_after, lat, inter_model_std=std
        )

        # 设置当前子图
        ax = axes[i]

        # 处理没有有效数据的情况
        if x_before.size == 0 and x_after.size == 0:
            ax.set_title(f"{name} (No valid data)")
            ax.text(0.5, 0.5, "No valid data", transform=ax.transAxes,
                    ha='center', va='center', fontsize=8, color='gray')
            ax.set_xlabel(f"Change in {name} ({units})")
            ax.set_ylabel("Land fraction (%)")
            continue  # 跳过当前子图的后续绘制

        # 设置y轴范围为固定值
        y_max_values = [0.9,0.4, 0.25]  # 为三个子图分别设定最大值
        ax.set_ylim(0, y_max_values[i])

        # 按照正确的顺序绘制图形以确保适当的重叠
        # 首先绘制两条线
        if x_before.size > 0:
            ax.plot(x_before, y_before, color=red_color, linewidth=1.0,
                    label=f'Unconstrained: {sig_area_before:.1f}%')
        if x_after.size > 0:
            ax.plot(x_after, y_after, color=blue_color, linewidth=1.0,
                    label=f'Constrained: {sig_area_after:.1f}%')

        # 填充超过阈值的区域（未约束数据）
        if x_before.size > 0:
            sig_mask_before = x_before > std
            if np.any(sig_mask_before):  # 确保有点被选中
                ax.fill_between(x_before[sig_mask_before], 0, y_before[sig_mask_before],
                                color=red_fill, alpha=0.5)

        # 填充超过阈值的区域（约束数据）
        if x_after.size > 0:
            sig_mask_after = x_after > std
            if np.any(sig_mask_after):  # 确保有点被选中
                ax.fill_between(x_after[sig_mask_after], 0, y_after[sig_mask_after],
                                color=blue_fill, alpha=0.5)

        # 设置子图标题和标签
        panel_labels = ['a', 'b', 'c']
        ax.set_title(f"{panel_labels[i]} {name}", loc='left', fontweight='bold')
        ax.set_xlabel(f"Change in {name} ({units})")
        ax.set_ylabel("Land fraction (%)")

        # 设置刻度
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.grid(False)

        # 添加图例（简洁风格）
        ax.legend(loc='upper right', frameon=False, handlelength=1.5,
                  handletextpad=0.5, borderaxespad=0.3)

    # 调整子图之间的间距 - 减小间距使图形更紧凑
    plt.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.20, wspace=0.15)  # 减小wspace

    return fig, axes


def main():
    # 加载数据
    f1 = xr.open_dataset('E:/north-future/pdf/input/long/raw/raw-ssp370.nc')
    f2 = xr.open_dataset('E:/north-future/pdf/input/long/c/c-ssp370.nc')
    f3 = xr.open_dataset('E:/north-future/pdf/input/long/std/3index-ssp370-std.nc')

    # 提取原始极端指数
    txx_before = f1['index'].values[0, :, :]
    tnn_before = f1['index'].values[1, :, :]
    rx5day_before = f1['index'].values[2, :, :]

    # 提取约束后的极端指数
    txx_after = f2['index'].values[0, :, :]
    tnn_after = f2['index'].values[1, :, :]
    rx5day_after = f2['index'].values[2, :, :]

    # 提取模型间标准差
    txx_std_map = f3['index'].values[0, :, :]
    tnn_std_map = f3['index'].values[1, :, :]
    rx5day_std_map = f3['index'].values[2, :, :]

    # 获取纬度信息
    lat = f1['lat'].values

    # 计算模型间标准差的空间加权平均值
    lat_weights = np.cos(np.deg2rad(lat))
    weights_2d = np.outer(lat_weights, np.ones(txx_std_map.shape[1]))

    # 计算TXx的面积加权平均标准差
    mask_txx_std = ~np.isnan(txx_std_map)
    arctic_txx_std = np.average(txx_std_map[mask_txx_std], weights=weights_2d[mask_txx_std])

    # 计算TNn的面积加权平均标准差
    mask_tnn_std = ~np.isnan(tnn_std_map)
    arctic_tnn_std = np.average(tnn_std_map[mask_tnn_std], weights=weights_2d[mask_tnn_std])

    # 计算RX5day的面积加权平均标准差
    mask_rx5day_std = ~np.isnan(rx5day_std_map)
    arctic_rx5day_std = np.average(rx5day_std_map[mask_rx5day_std], weights=weights_2d[mask_rx5day_std])

    print(f"Arctic Hot Extreme (TXx) Area-Weighted Inter-Model Std: {arctic_txx_std:.2f}")
    print(f"Arctic Cold Extreme (TNn) Area-Weighted Inter-Model Std: {arctic_tnn_std:.2f}")
    print(f"Arctic Precipitation Extreme (RX5day) Area-Weighted Inter-Model Std: {arctic_rx5day_std:.2f}")

    # 创建组合图（顺序要与标签匹配）
    fig, axes = plot_separated_indices(
        [txx_before, tnn_before, rx5day_before],
        [txx_after, tnn_after, rx5day_after],
        lat,
        [arctic_txx_std, arctic_tnn_std, arctic_rx5day_std],
        ["TXx", "TNn", "RX5day"],
        ["°C", "°C", "mm"],
        output_file=None  # 先不保存
    )

    # 显示预览
    plt.savefig('E:/north-future/new-p/original/fig4/ssp370/long-ssp370.pdf', format='pdf', bbox_inches='tight')
    plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
