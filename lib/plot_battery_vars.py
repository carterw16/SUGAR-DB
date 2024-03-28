import matplotlib.pyplot as plt
import numpy as np
import pickle

def plot_output_period(output, epoch_to_plot=0):
    num_outputs = len(output)
    num_periods = len(output['B'][0])
    x_values = np.linspace(1,num_periods, num_periods)
    breakpoint()

    fig, axs = plt.subplots(6, 2, figsize=(12, 24))  # Create a grid of 6 rows and 2 columns
    epochs_to_plot = [0, 1, 2, 3, 4, 5, 6]
    colors = plt.cm.viridis(np.linspace(0, 1, len(epochs_to_plot)))  # Generate colors for each epoch
    handles, labels = [], []  # Collect handles and labels for the legend  
    
    row = 0
    col = 0
    for i, (output_name, output_data) in enumerate(output.items()):
        # Plot for each epoch
        for j, epoch_to_plot in enumerate(epochs_to_plot):
            data_to_plot = [output_data[epoch_to_plot][k][0][0] for k in range(num_periods)]
            line = axs[row, col].plot(data_to_plot, color=colors[j], alpha=0.5)
            if i == 0:  # Collect handles and labels only once
                handles.append(line[0])
                labels.append(f'Epoch {epoch_to_plot}')
            
            axs[row, col].set_title(f'{output_name} - Epochs {repr(epochs_to_plot)}')
            axs[row, col].set_xlabel('Period')
            axs[row, col].set_ylabel('Value')
            
            # Set x-axis ticks to integers only
            axs[row, col].xaxis.set_major_locator(plt.MaxNLocator(integer=True))
            
        col += 1
        if col == 2:
            col = 0
            row += 1
            
    fig.legend(handles, labels, loc='upper right')    
    plt.tight_layout()
    plt.show()

def load_from_pickle(file_path):
    with open(file_path, 'rb') as pickle_file:
        output = pickle.load(pickle_file)
    
    return output

output = load_from_pickle('battery_outputs.pkl')
plot_output_period(output)
