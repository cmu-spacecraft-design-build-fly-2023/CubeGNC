import numpy as np
import pickle
import pandas as pd
import plotly.express as px

def pkl_loader(path_to_config_pkl):
    with open(path_to_config_pkl, 'rb') as f:
        data = pickle.load(f)
    return data

def mc_plot_momentum_magnitude_vs_time(mc_results, max_samples=500):
    J = np.array([[ 4.6e-03 , 2.0e-05 ,-3.0e-06],
                  [ 2.0e-05,  4.6e-03 ,-2.0e-05],
                  [-3.0e-06, -2.0e-05,  4.6e-03]])

    Ntrials = mc_results["T"].shape[0]
    h_average = np.zeros((max_samples + 1,))
    t_plot_average = np.zeros_like(h_average)
    long_data = []

    for mc_step in range(Ntrials):
        xhist = mc_results["X"][mc_step, :, :]
        thist = mc_results["T"][mc_step, 0, :]
        print(thist.shape)

        downsample = np.linspace(0, len(thist) - 1, max_samples + 1).astype(int)
        #omega = xhist[10:13, downsample]
        omega = xhist[10:13]
        h = J @ omega
        h_mag = np.linalg.norm(h, axis=0)
        #t_plot = thist[downsample] / (60 * 60)

        for i in range(len(h_mag)):
            if not (thist[i] == 0):
                long_data.append({"Trial": mc_step + 1, "time": thist[i], "variable": h_mag[i]})

    df = pd.DataFrame(long_data)
    print(df.size)

    fig = px.line(df, x="time", y="variable", color="Trial")
    fig.update_layout(title="B-cross control - MC analysis", xaxis_title="Time (seconds)", yaxis_title = "Magnitude of Momentum (Nms)")
    fig.show()

if __name__ == "__main__":
    x = pkl_loader('my_data_test_10_runs.pickle')
    print(x["X"].shape)
    mc_plot_momentum_magnitude_vs_time(x)
