# this script defines a function that returns and plots a matrix in the form of a heatmap. The returned matrix can be normalised by
# either over the rows or over the columns, by dividing by either the sum or the max element
# normalisation over columns entirely recommended


import seaborn as sns
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np


# this heatmap generator works well for values between 0 and 1 (inclusive)
def heatmaps_seq(matrix,title=''):
    cmaps = OrderedDict()
    cmaps['Sequential (2)'] = ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
        'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
        'hot', 'afmhot', 'gist_heat', 'copper']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax1 = sns.heatmap(matrix, cmap=cmaps['Sequential (2)'][-8], yticklabels=False, xticklabels=False)
    # add diagonal
    # for i in range(len(matrix)):
    #     ax.add_patch(Rectangle((i, i), 1, 1, fill=False, edgecolor='blue', lw=3))
    ax1.set_title(title)
    # ax.set_ylabel(ylabel='PSMC output')
    fig.show()
    # fig.savefig(heatmaps_path + sys.argv[4])
    return None

# this heatmap works well for values between -1 and 1
def heatmaps_div(matrix,title=''):
    cmaps = OrderedDict()
    cmaps['Diverging'] = [
        'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
        'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax1 = sns.heatmap(matrix, cmap=cmaps['Diverging'][11], center=0,yticklabels=False, xticklabels=False)
    ax1.set_title(title)
    # for i in range(len(matrix)):
    #     ax.add_patch(Rectangle((i, i), 1, 1, fill=False, edgecolor='blue', lw=3))
    # ax.set_title(plot_title)
    # ax.set_ylabel(ylabel='PSMC output')
    fig.show()
    # fig.savefig(heatmaps_path + sys.argv[4])
    return None