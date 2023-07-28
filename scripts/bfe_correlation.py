import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import kendalltau
from sklearn.metrics import mean_absolute_error, mean_squared_error

def plt_bfe_vs_spr(data,x_col, y_col, y_err, FigureName):
    """
        data: dataframe
       FigureName:  String; Name of the generated plot
    """

    lines = {'linestyle': 'None'}
    plt.rc('lines', **lines)
    fig, ax = plt.subplots(figsize=(8, 5))
    lim_min=int(min(data[x_col].tolist() + data[y_col].tolist()))-2
    lim_max=int(max(data[x_col].tolist() + data[y_col].tolist())) +2
    x = np.linspace(lim_min, lim_max, 1000)
    ax.plot(x, x + 0, linestyle='solid', color='black')
    ax.plot(x, x + 1, linestyle='dashed', color='grey')
    ax.plot(x, x - 1, linestyle='dashed', color='grey')
    ax.fill_between(x, x + 2, x - 2, facecolor='gold',alpha=0.2)
    ax.fill_between(x, x + 1, x - 1, facecolor='orange', alpha=0.5)
    ax.set_facecolor('ivory')
    ax.set_xlabel(f'{x_col} (kcal/mol)')
    ax.set_ylabel(f'{y_col} (Kcal/mol)')
    ax.set_title(f'{FigureName}')
    plt.ylim(lim_min,lim_max)
    plt.xlim(lim_min,lim_max)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid()
    plt.tight_layout()
    for i, row in data.iterrows():
        ax.scatter(row[x_col], row[y_col],marker= 'o', edgecolor="black", s=30)
        if bool(y_err):
            ax.errorbar(row[x_col], row[y_col], yerr=row[y_err], color='black')

##### Calculate statistics
    p = stats.pearsonr(data[y_col].values.tolist(), data[x_col].values.tolist())
    rmse = np.sqrt(mean_squared_error(data[y_col].values.tolist(), data[x_col].values.tolist()))
    mue = mean_absolute_error(data[y_col], data[x_col])
    taocoef= kendalltau(data[y_col], data[x_col])[0]
    ##### Print statiscs on plot
#tao
    ax.annotate(r'$τ=$' + str(round(taocoef, 2)), xy=(0, 0),
        xycoords='axes fraction', xytext=(0.05, 0.90),
        textcoords='axes fraction', size=12)
#rmse
    ax.annotate(r'$RMSE=$' + str(round(rmse, 2)), xy=(0, 0),
        xycoords='axes fraction', xytext=(0.05, 0.85),
        textcoords='axes fraction', size=12)
#mue
    ax.annotate(r'$MUE=$' + str(round(mue, 2)), xy=(0, 0),
        xycoords='axes fraction', xytext=(0.05, 0.8),
        textcoords='axes fraction', size=12)

#R^2
    ax.annotate(r'$R^2=$' + str(round(p[0]**2, 2)), xy=(0, 0),
        xycoords='axes fraction', xytext=(0.05, 0.95),
        textcoords='axes fraction', size=12)
    plt.savefig(f'{FigureName}_sim_vs_exp.png')
    #plt.show()


def main():
    """
    plot the x, y for the input_csv file
    """

    parser = argparse.ArgumentParser(description="Collect BFE throughput.")
    parser.add_argument(
        "--input_csv",
        dest="input_csv",
        action="store",
        help="the path to input csv file.",
    )

    parser.add_argument(
        "--x",
        dest="x",
        action="store",
        help="x column name",
    )

    parser.add_argument(
        "--y",
        dest="y",
        action="store",
        help="y column name",
    )

    parser.add_argument(
        "--y_error",
        dest="y_error",
        action="store",
        help="y error column name",
    )

    parser.add_argument(
        "--title",
        dest="title",
        action="store",
        help="plot title name",
    )

    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)



    plt_bfe_vs_spr(df,args.x,args.y,args.y_error, args.title)

if __name__ == "__main__":
    main()
