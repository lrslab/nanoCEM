import plotnine as p9
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.sans-serif'] = ['Arial']
SIG_PCTL_RANGE = (2.5, 97.5)
def draw_boxplot(df,results_path,pos,leng,chrom,base_list,aligned_num_wt,aligned_num_ivt,strand="+"):
    title = chrom + ':' + str(pos - leng + 1) + '-' + str(pos + leng + 2) + ':' + strand
    title = title + '   Sample:' + str(aligned_num_wt) + ' Control:' + str(aligned_num_ivt)
    item_list = ['Mean', 'STD', 'Median', 'Dwell_time']
    plot_list=[]
    for item in item_list:

        sig_min, sig_max = np.percentile(df[item], SIG_PCTL_RANGE)
        sig_diff = sig_max - sig_min
        ylim_tuple = (sig_min - sig_diff * 0.1, sig_max + sig_diff * 0.1)

        plot = p9.ggplot(df, p9.aes(x='position', y=item, fill='type')) \
               + p9.theme_bw() \
               +p9.geom_boxplot( outlier_shape='',position=p9.position_dodge(0.9),size=0.25) \
               + p9.scale_fill_manual(values={"Sample": "#ff6f91", "Control": "#7389af"}) \
               + p9.scale_x_discrete(labels=list(base_list)) \
               + p9.theme(
                    figure_size=(6, 3),
                    panel_grid_minor=p9.element_blank(),
                    axis_text=p9.element_text(size=13),
                    axis_title=p9.element_text(size=13),
                    title=p9.element_text(size=13),
                    legend_position='none')\
               + p9.ylim(ylim_tuple)\
               + p9.labs(title=title, x=str(pos + 1), y=item)

        plot.save(filename=results_path + "/" + item + "_boxplot.pdf", dpi=300)

def draw_volin(df,results_path,pos,leng,chrom,base_list,aligned_num_wt,aligned_num_ivt,strand="+"):
    title = chrom + ':' + str(pos - leng + 1) + '-' + str(pos + leng + 2) + ':' + strand
    title = title + '   Sample:' + str(aligned_num_wt) + '  Control:' + str(aligned_num_ivt)
    item_list = ['Mean', 'STD', 'Median', 'Dwell_time']
    for item in item_list:

        sig_min, sig_max = np.percentile(df[item], SIG_PCTL_RANGE)
        sig_diff = sig_max - sig_min
        ylim_tuple = (sig_min - sig_diff * 0.1, sig_max + sig_diff * 0.1)

        plot = p9.ggplot(df, p9.aes(x='position', y=item, fill='type')) \
               + p9.geom_violin(style='left-right',position=p9.position_dodge(0),color='none',width=1.5) \
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
        )\
        + p9.ylim(ylim_tuple)

        plot = plot + p9.labs(title=title, x=str(pos + 1), y=item)
        # plot.render_matplotlib()
        plot.save(filename=results_path + "/" + item + "_violin.pdf", dpi=300)
        # print(plot)