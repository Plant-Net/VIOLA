####################################################################################################################################################################################################
##### IMPORTS
####################################################################################################################################################################################################

import os
import time
import random
from math import sqrt

import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN

import argparse
import warnings

warnings.filterwarnings('ignore')

print(f"sklearn: {__import__('sklearn').__version__}")
print(f"Tensorflow/Keras: {keras.__version__}")

####################################################################################################################################################################################################
## VAE FUNCTION
#####################################################################################################################################################################################################

tfk = tf.keras
tfkl = tfk.layers


def seed_everything(seed=779):
    """Set all seeds to ensure reproducibility."""
    tf.random.set_seed(seed)
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)


def make_features(data):
    """
    Select and normalize features for the model.
    
    :param data: Input DataFrame.
    :return: Normalized DataFrame with selected features.
    """
    
    features_data = data[
        ['SIFTval', 'PolyPhenVal', 'priPhCons', 'mamPhCons', 'verPhCons',
         'priPhyloP', 'mamPhyloP', 'verPhyloP', 'EncodeH3K27ac.max', 'EncodeH3K27ac.sum',
         'EncodeH3K4me1.max', 'EncodeH3K4me1.sum', 'EncodeH3K4me3.sum', 'EncodeH3K4me3.max',
         'EncodeDNase.sum', 'EncodeDNase.max', 'MMSp_acceptorIntron', 'MMSp_acceptor',
         'MMSp_exon', 'MMSp_donor', 'MMSp_donorIntron', 'Freq100bp', 'Rare100bp', 'Sngl100bp',
         'Freq1000bp', 'Rare1000bp', 'Sngl1000bp', 'Freq10000bp', 'Rare10000bp', 'Sngl10000bp',
         'RawScore', 'PHRED']
        ]

    features_data_treated = features_data.fillna(features_data.median())
    scaler = MinMaxScaler()
    print(scaler.fit(features_data_treated))

    scaled_data = scaler.transform(features_data_treated)
    return pd.DataFrame(scaled_data, columns=features_data_treated.columns)


def Vae_function(file, batch_size=16, epochs=50, kernel="lecun_normal",
                 kl_loss_weight=0.5, edl1=25, latent_size=16, ddl1=25):
    """
    Variational Autoencoder (VAE).
    
    :param file: Pre-processed DataFrame.
    :return: Latent space as a DataFrame.
    """

    seed_everything()

    def sampling(args):
        z_mean, z_log_var = args
        batch = tfk.backend.shape(z_mean)[0]
        dim = tfk.backend.int_shape(z_mean)[1]
        epsilon = tfk.backend.random_normal(shape=(batch, dim), seed=779)
        return z_mean + tf.keras.backend.exp(0.5 * z_log_var) * epsilon

    data = file
    n_cols = data.shape[1]
    start_time = time.time()

    # Encoder
    inputs = tfk.Input(shape=(n_cols,), name='encoder_input')
    hidden_encoder_1 = tfkl.Dense(edl1, kernel_initializer=kernel, name="hidden_encoder_1")(inputs)
    encoder_norm_1 = tfkl.BatchNormalization(name="encoder_norm_1")(hidden_encoder_1)
    hidden_encoder_1_activation = tfkl.ELU()(encoder_norm_1)

    z_mean = tfkl.Dense(latent_size, name='z_mean')(hidden_encoder_1_activation)
    z_log_var = tfkl.Dense(latent_size, name='z_log_var')(hidden_encoder_1_activation)
    z = tfkl.Lambda(sampling, output_shape=(latent_size,), name='z')([z_mean, z_log_var])

    encoder = tfk.Model(inputs, [z_mean, z_log_var, z], name='encoder')

    # Decoder
    latent_inputs = tfk.Input(shape=(latent_size,), name='z_sampling')
    hidden_decoder_5 = tfkl.Dense(ddl1, kernel_initializer=kernel, name="hidden_decoder_4")(latent_inputs)
    decoder_norm_5 = tfkl.BatchNormalization(name="decoder_norm_4")(hidden_decoder_5)
    hidden_encoder_5_activation = tfkl.ELU()(decoder_norm_5)
    outputs = tfkl.Dense(n_cols)(hidden_encoder_5_activation)

    decoder = tfk.Model(latent_inputs, outputs, name="decoder")

    # VAE
    outputs = decoder(encoder(inputs)[2])
    vae = tfk.Model(inputs, outputs=[outputs, z], name="vae")

    msle = tfk.losses.MSLE(inputs, outputs[0]) * n_cols
    kl_loss = 1 + z_log_var - tfk.backend.square(z_mean) - tfk.backend.exp(z_log_var)
    kl_loss = tfk.backend.sum(kl_loss, axis=-1) * -kl_loss_weight

    vae.add_loss(msle + kl_loss)
    vae.add_metric(msle, name='msle', aggregation='mean')
    vae.add_loss(kl_loss)
    vae.add_metric(kl_loss, name='kl_loss', aggregation='mean')
    vae.compile(optimizer="adam")
    vae.summary()

    # Training
    history = vae.fit(data, epochs=epochs, batch_size=batch_size)
    end_time = time.time() - start_time
    print(f"----- {end_time // 3600} hours {(end_time % 3600) // 60} minutes {round((end_time % 3600) % 60, 2)} seconds -----")

    pred, latent_space = vae.predict(data)
    latent_space = pd.DataFrame(latent_space)

    pred = pd.DataFrame(pred, index=data.index, columns=data.columns)

    # Generate random samples (not used later but kept for consistency)
    latent_samples = np.random.normal(size=(len(data), latent_size))
    generated_samples = decoder.predict(latent_samples)

    return latent_space

####################################################################################################################################################################################################
## DBSCAN
####################################################################################################################################################################################################

def aff(data):
    """Return only the first 16 columns (latent space)."""
    return data.iloc[:, :16]


def find_eps(d, n=2500):
    """
    Estimate epsilon value for DBSCAN.
    
    :param d: Latent space data.
    :param n: Number of points considered for averaging.
    :return: Estimated epsilon.
    """
    k = int(sqrt(d.shape[0]))
    neigh = NearestNeighbors(n_neighbors=k).fit(d)
    distances, _ = neigh.kneighbors(d)
    avg_distances = np.mean(distances, axis=1)
    y_values = np.sort(avg_distances)
    x_values = np.arange(len(y_values))
    y_values_subset = y_values[(x_values >= 0) & (x_values <= n)]
    return sum(y_values_subset) / len(y_values_subset)


def dbs(d, e=10, m=32):
    """
    Apply DBSCAN clustering.
    
    :param d: Data containing latent space.
    :param e: Epsilon value.
    :param m: min_samples (recommended = 2 * latent_size).
    """
    d_db = aff(d)
    d_std = StandardScaler().fit_transform(d_db)
    model = DBSCAN(eps=e, min_samples=m).fit(d_std)
    d["Cluster"] = model.fit_predict(d_std)
    print(d["Cluster"].value_counts())
    print(e)
    return d

####################################################################################################################################################################################################
## MAIN
####################################################################################################################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file_path", type=str, required=True, help="Path to the file to process")
    parser.add_argument("-o", "--output_folder_path", type=str, required=True, help="Path to the output folder")
    args = parser.parse_args()

    filepath = args.file_path
    output = args.output_folder_path

    if filepath.endswith(".csv"):
        filename = os.path.basename(filepath)
        name = filename.split("_")[0]
        print(name)

        data = pd.read_csv(filepath)

        # VAE
        ls = Vae_function(make_features(data))

        res = pd.concat([ls, data], axis=1, join="inner")
        res.to_csv(f"{output}{name}_rare_variant_latent_space.csv", index=False)

        # DBSCAN
        e = find_eps(aff(res), 1500)
        res_db = dbs(res, e)
        res_db = res_db[res_db["Cluster"] == -1]
        res_db.to_csv(f"{output}{name}_res_dbscan.csv", index=False)
    else:
        print("File needs to be in CSV format")
