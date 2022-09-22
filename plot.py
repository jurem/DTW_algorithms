import sys
import pandas as pd
import matplotlib.pyplot as plt

# folder where the results are stored
scenario = sys.argv[1]

# plot kind
kind = sys.argv[2]

# algorithms to be visualized
algs = sys.argv[3:]


def load_data(name):
    fn = scenario + '/' + name + '.csv'
    return (name, pd.read_csv(fn, sep=' ', header=1,
        names=['n', 'm', 'v', 'realtime', 'cputime']))


def plot_data(title, data):
    plt.figure(figsize=(16, 9), dpi=80)
    # plt.tight_layout()

    for lbl, df in data:
        vals = df['realtime']
        if kind == 'acc': vals = base / vals
        plt.plot(df['n'], vals, label=lbl)
    plt.legend(loc='best')

    plt.show()

if kind == 'acc':
    base = load_data('rect_fw')[1]['realtime']

dfs = []
for name in algs:
    df = load_data(name)
    dfs.append(df)

plot_data('DTW', dfs)