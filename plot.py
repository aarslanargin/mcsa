import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker
import argparse

def main(params):
    ion = params.ion

    df = pd.read_csv(f'./ions/{ion}-ec-qm.dat', sep=",", header=None)
    df_1 = pd.read_csv(f'./bind_ener/{ion}-ec-before-fit.txt', sep=",", header=None)
    df_2 = pd.read_csv(f'./bind_ener/{ion}-ec-after-fit.txt', sep=",", header=None)

    df.columns = ["x", "qm"]
    df_1.columns = ["x", "bfit"]
    df_2.columns = ["x", "afit"]
    fig, ax = plt.subplots()

    ax.plot(df['x'],df['qm'], label='Quantum')
    ax.plot(df_1['x'],df_1['bfit'], label='MM/Before Fit')
    ax.plot(df_2['x'],df_2['afit'], label='MM/After Fit')

    ax.set_xlim(-0.2,4.0)
    ax.set_ylim(df['qm'].min() - 5 ,0)
    ax.set_title(f'{ion}-EC Binding Energies')

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
    plt.xlabel('Ion Displacement (Å)')
    plt.ylabel('Binding Energies (kcal/mol)')
    plt.legend(loc = "lower right")
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plotting after/before fit binding energies')

    parser.add_argument('--ion', help='ion name')
    args = parser.parse_args()

    main(args)