import sys
import pandas as pd
import matplotlib.pyplot as plt

# folder where the results are stored
scenario = sys.argv[1]

# algorithms to be visualized
algs = sys.argv[2:]


def load_data(name):
    fn = scenario + '/' + name + '.csv'
    return (name, pd.read_csv(fn, sep=' ', header=1,
        names=['n', 'm', 'v', 'realtime', 'cputime']))


def plot_data(title, data):
    plt.figure(figsize=(16, 9), dpi=80)
    # plt.tight_layout()

    for lbl, df in data:
        plt.plot(df['n'], df['realtime'], label=lbl)
    plt.legend(loc='best')

    plt.show()

dfs = []
for name in algs:
    df = load_data(name)
    dfs.append(df)

plot_data('DTW', dfs)