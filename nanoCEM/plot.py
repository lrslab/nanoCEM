import pandas as pd
import plotnine as p9
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
# plt.rcParams['font.sans-serif'] = ['Arial']
SIG_PCTL_RANGE = (2.5, 97.5)


# def draw_boxplot(df,results_path,pos,base_list,title):
#     item_list = ['Mean', 'STD', 'Median', 'Dwell_time']
#     plot_list=[]
#     for item in item_list:
#
#         sig_min, sig_max = np.percentile(df[item], SIG_PCTL_RANGE)
#         sig_diff = sig_max - sig_min
#         ylim_tuple = (sig_min - sig_diff * 0.1, sig_max + sig_diff * 0.1)
#
#         plot = p9.ggplot(df, p9.aes(x='position', y=item, fill='type')) \
#                + p9.theme_bw() \
#                +p9.geom_boxplot( outlier_shape='',position=p9.position_dodge(0.9),size=0.25) \
#                + p9.scale_fill_manual(values={"Sample": "#ff6f91", "Control": "#7389af"}) \
#                + p9.scale_x_discrete(labels=list(base_list)) \
#                + p9.theme(
#                     figure_size=(6, 3),
#                     panel_grid_minor=p9.element_blank(),
#                     axis_text=p9.element_text(size=13),
#                     axis_title=p9.element_text(size=13),
#                     title=p9.element_text(size=13),
#                     legend_position='none')\
#                + p9.ylim(ylim_tuple)\
#                + p9.labs(title=title, x=str(pos + 1), y=item)
#
#         plot.save(filename=results_path + "/" + item + "_boxplot.pdf", dpi=300)
#
# def draw_violin(df,results_path,pos,base_list,title):
#
#     item_list = ['Mean', 'STD', 'Median', 'Dwell_time']
#     for item in item_list:
#
#         sig_min, sig_max = np.percentile(df[item], SIG_PCTL_RANGE)
#         sig_diff = sig_max - sig_min
#         ylim_tuple = (sig_min - sig_diff * 0.1, sig_max + sig_diff * 0.1)
#
#         plot = p9.ggplot(df, p9.aes(x='position', y=item, fill='type')) \
#                + p9.geom_violin(style='left-right',position=p9.position_dodge(0),color='none',width=1.5) \
#                + p9.theme_bw() \
#                + p9.scale_fill_manual(values={"Sample": "#ff6f91", "Control": "#7389af"}) \
#                + p9.scale_x_discrete(labels=list(base_list)) \
#                + p9.theme(
#             figure_size=(6, 3),
#             panel_grid_minor=p9.element_blank(),
#             axis_text=p9.element_text(size=13),
#             axis_title=p9.element_text(size=13),
#             title=p9.element_text(size=13),
#             legend_position='none'
#         )\
#         + p9.ylim(ylim_tuple)
#
#         plot = plot + p9.labs(title=title, x=str(pos + 1), y=item)
#         # plot.render_matplotlib()
#         plot.save(filename=results_path + "/" + item + "_violin.pdf", dpi=300)
# print(plot)

def signal_plot(df, results_path, pos, base_list, title, plot_type,filter=False):
    item_list = ['Mean', 'STD', 'Median', 'Dwell time']


    if plot_type != 'merged':

        for item in item_list:
            sig_min, sig_max = np.percentile(df[item], SIG_PCTL_RANGE)
            sig_diff = sig_max - sig_min
            ylim_tuple = (sig_min - sig_diff * 0.1, sig_max + sig_diff * 0.1)

            plot = p9.ggplot(df, p9.aes(x='position', y=item, fill='type')) \
                   + p9.theme_bw() \
                   + p9.scale_fill_manual(values={"Sample": "#ff6f91", "Control": "#7389af"}) \
                   + p9.scale_x_discrete(labels=list(base_list)) \
                   + p9.theme(
                figure_size=(6, 3),
                panel_grid_minor=p9.element_blank(),
                axis_text=p9.element_text(size=13),
                axis_title=p9.element_text(size=13),
                title=p9.element_text(size=13),
                legend_position='none'
            ) \
                   + p9.ylim(ylim_tuple)
            plot = plot + p9.labs(title=title, x=str(pos + 1), y=item)
            if plot_type == 'boxplot':
                plot = plot + p9.geom_boxplot(outlier_shape='', position=p9.position_dodge(0.9), size=0.25)
            elif plot_type == 'violin_plot':
                plot = plot + p9.geom_violin(style='left-right', position=p9.position_dodge(0), color='none', width=1.5)
            else:
                raise Exception("Unsupported figure type!")
            # plot.render_matplotlib()
            if item == 'Dwell time':
                item = 'Dwell_time'
            plot.save(filename=results_path + "/" + item + "_" + plot_type + ".pdf", dpi=300)
    else:
        new_df = None
        for item in item_list:
            # collect data
            temp = df[[item, 'position', 'type']].copy()
            temp.columns = ['value', 'position', 'Group']
            temp.loc[:, 'stats'] = item
            if filter:
                sig_min, sig_max = np.percentile(temp['value'], SIG_PCTL_RANGE)
                sig_diff = sig_max - sig_min
                ylim_tuple = [sig_min - sig_diff * 0.1, sig_max + sig_diff * 0.1]
                temp = temp[(temp['value'] >= ylim_tuple[0]) & (temp['value'] <= ylim_tuple[1])]
            if new_df is None:
                new_df = temp
            else:
                new_df = pd.concat([new_df, temp], axis=0)

        plot = p9.ggplot(new_df, p9.aes(x='position', y="value", fill='Group')) \
               + p9.theme_bw() \
               + p9.scale_fill_manual(values={"Sample": "#ff6f91", "Control": "#7389af", "Single": "#a3abbd"}) \
               + p9.scale_x_discrete(labels=list(base_list)) \
               + p9.theme(
            figure_size=(8, 8),
            panel_grid_minor=p9.element_blank(),
            axis_text=p9.element_text(size=13),
            axis_title=p9.element_text(size=13),
            title=p9.element_text(size=13),
            legend_position='bottom',
            legend_title=p9.element_blank(),
            strip_text=p9.element_text(size=13),
            strip_background=p9.element_rect(alpha=0),

        ) \
               + p9.facet_grid('stats ~', scales='free_y')
        plot = plot + p9.labs(title=title, x=str(pos + 1), y='')

        if new_df['Group'].drop_duplicates().shape[0] == 1:
            plot2 = plot + p9.geom_violin(color='none', position=p9.position_dodge(0.9), width=1)
            plot2 = plot2 + p9.geom_boxplot(outlier_shape='', position=p9.position_dodge(0.9), size=0.2, width=0.1)
            plot2.save(filename=results_path + "/Merged_single.pdf", dpi=300)
        else:
            plot1 = plot + p9.geom_boxplot(outlier_shape='', position=p9.position_dodge(0.9), size=0.2, width=0.75,alpha = 0.8)
            plot1.save(filename=results_path + "/Merged_boxplot.pdf", dpi=300)
            plot2 = plot + p9.geom_violin(style='left-right', position=p9.position_dodge(0), color='none', width=1.5,alpha = 0.8)
            plot2.save(filename=results_path + "/Merged_violin.pdf", dpi=300)


def draw_signal(df, start, base,start_index,end_index):
    df = pd.DataFrame(df)
    df.columns = ['raw']
    start=start[start_index:end_index]
    df=df[start[0]:start[-1]]
    df = df.reset_index()
    plot = p9.ggplot(df, p9.aes(x='index', y="raw")) \
           + p9.theme_bw() \
           + p9.geom_line() \
           + p9.theme(
        figure_size=(6, 3),
        panel_grid_minor=p9.element_blank(),
        axis_text=p9.element_text(size=13),
        axis_title=p9.element_text(size=13),
        title=p9.element_text(size=13),
        strip_text=p9.element_text(size=13),
        legend_position='none',
        strip_background=p9.element_rect(alpha=0)
    )
    # plot.save(filename="/home/zhguo/Dropbox/@labfiles/projects/GUO/ONT_showcase_tool/plot/signal.pdf", dpi=300)
    for item in start:
        plot = plot + p9.geom_vline(xintercept=item, linetype='dashed', color='red')
    print(plot)
    # plot.save(filename="/home/zhguo/Dropbox/@labfiles/projects/GUO/ONT_showcase_tool/plot/norm_signal_aligned.pdf", dpi=300)

    # print(1)
