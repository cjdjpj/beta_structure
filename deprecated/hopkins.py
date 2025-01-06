import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance

def create_clustered_data(num_samples_per_cluster):
    """
    Generate a synthetic clustered dataset.
    """
    cluster_centers = np.array([
        [1, 1],
        [5, 5],
        [1, 5],
        [5, 1]
    ])
    data = []
    for center in cluster_centers:
        cluster_data = center + np.random.randn(num_samples_per_cluster, 2) * 0.3
        data.append(cluster_data)
    data = np.vstack(data)
    return data

def hopkins_statistic(data, n=None):
    """
    Calculate the Hopkins statistic for a 2D numpy array.
    """
    
    num_samples = data.shape[0]
    
    if n is None:
        n = num_samples

    min_values = np.min(data, axis=0)
    max_values = np.max(data, axis=0)
    random_points = np.random.uniform(min_values, max_values, (n, data.shape[1]))

    dist_data = distance.cdist(data, data, 'euclidean')
    nearest_data_distances = np.min(dist_data + np.eye(num_samples) * np.max(dist_data), axis=1)

    dist_random = distance.cdist(random_points, data, 'euclidean')
    nearest_random_distances = np.min(dist_random, axis=1)

    H = np.sum(nearest_random_distances) / (np.sum(nearest_random_distances) + np.sum(nearest_data_distances))
    
    return H

def main():
    clustered_data = np.random.rand(1000, 2) #RANDOM DATA
    uniform_data = create_clustered_data(250)
    H_clustered = hopkins_statistic(clustered_data)
    H_uniform = hopkins_statistic(uniform_data)

    plt.figure(figsize=(16, 8))

    plt.subplot(1, 2, 1)
    plt.scatter(clustered_data[:, 0], clustered_data[:, 1], color='blue', alpha=0.6)
    plt.title(f'Clustered Data\nHopkins Statistic: {H_clustered:.2f}')
    plt.xlabel('Feature 1')
    plt.ylabel('Feature 2')

    plt.subplot(1, 2, 2)
    plt.scatter(uniform_data[:, 0], uniform_data[:, 1], color='orange', alpha=0.6)
    plt.title(f'Uniform Data\nHopkins Statistic: {H_uniform:.2f}')
    plt.xlabel('Feature 1')
    plt.ylabel('Feature 2')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

