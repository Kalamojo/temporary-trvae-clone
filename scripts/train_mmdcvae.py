import sys

import numpy as np
import scanpy as sc

import reptrvae

data_name = sys.argv[1]

if data_name == "haber":
    conditions = ["Control", "Hpoly.Day10"]
    target_conditions = ["Hpoly.Day10"]
    source_condition = "Control"
    target_condition = "Hpoly.Day10"
    labelencoder = {"Control": 0, "Hpoly.Day10": 1}
    cell_type_key = "cell_label"
    condition_key = "condition"
    if len(sys.argv) == 3:
        specific_celltype = sys.argv[2]
    else:
        specific_celltype = "Tuft"

elif data_name == "kang":
    conditions = ["control", "stimulated"]
    target_conditions = ["stimulated"]
    source_condition = "control"
    target_condition = "stimulated"
    labelencoder = {"control": 0, "stimulated": 1}
    cell_type_key = "cell_type"
    condition_key = "condition"
    if len(sys.argv) == 3:
        specific_celltype = sys.argv[2]
    else:
        specific_celltype = "NK"
else:
    raise Exception("InValid data name")

adata = sc.read(f"./data/{data_name}_normalized.h5ad")
adata = adata[adata.obs[condition_key].isin(conditions)]

if adata.shape[1] > 2000:
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var['highly_variable']]

train_adata, valid_adata = reptrvae.utils.train_test_split(adata, 0.80)

if specific_celltype == 'all':
    for specific_celltype in adata.obs[cell_type_key].unique().tolist():
        net_train_adata = train_adata[
            ~((train_adata.obs[cell_type_key] == specific_celltype) & (train_adata.obs[condition_key].isin(target_conditions)))]
        net_valid_adata = valid_adata[
            ~((valid_adata.obs[cell_type_key] == specific_celltype) & (valid_adata.obs[condition_key].isin(target_conditions)))]

        network = reptrvae.models.MMDCVAE(x_dimension=net_train_adata.shape[1],
                                          z_dimension=40,
                                          n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                          alpha=1e-5,
                                          beta=100,
                                          eta=100,
                                          clip_value=1e6,
                                          lambda_l1=0.0,
                                          lambda_l2=0.0,
                                          learning_rate=0.001,
                                          model_path=f"./models/MMDCVAE/best/{data_name}-{specific_celltype}/",
                                          dropout_rate=0.2,
                                          output_activation='relu')

        network.train(net_train_adata,
                      net_valid_adata,
                      labelencoder,
                      condition_key,
                      n_epochs=500,
                      batch_size=512,
                      verbose=2,
                      early_stop_limit=50,
                      lr_reducer=35,
                      shuffle=True,
                      save=True,
                      )

        encoder_labels, _ = reptrvae.tl.label_encoder(net_train_adata, labelencoder, condition_key)
        decoder_labels = encoder_labels
        mmd_adata = network.to_mmd_layer(net_train_adata, encoder_labels, decoder_labels)

        cell_type_adata = adata[adata.obs[cell_type_key] == specific_celltype]
        source_adata = cell_type_adata[cell_type_adata.obs[condition_key] == source_condition]
        source_labels = np.zeros(source_adata.shape[0]) + labelencoder[source_condition]
        target_labels = np.zeros(source_adata.shape[0]) + labelencoder[target_condition]

        pred_adata = network.predict(source_adata,
                                     encoder_labels=source_labels,
                                     decoder_labels=target_labels,
                                     )

        pred_adata.obs[condition_key] = [f"{specific_celltype}_pred_{target_condition}"] * pred_adata.shape[0]

        pred_adata.write_h5ad(f"./data/reconstructed/MMDCVAE_{data_name}/{specific_celltype}.h5ad")

#         reptrvae.pl.plot_umap(mmd_adata,
#                               condition_key, cell_type_key,
#                               frameon=False, path_to_save=f"./results/{data_name}/", model_name="MMDCVAE_MMD", ext="png")

else:
    net_train_adata = train_adata[
            ~((train_adata.obs[cell_type_key] == specific_celltype) & (train_adata.obs[condition_key].isin(target_conditions)))]
    net_valid_adata = valid_adata[
            ~((valid_adata.obs[cell_type_key] == specific_celltype) & (valid_adata.obs[condition_key].isin(target_conditions)))]

    network = reptrvae.models.MMDCVAE(x_dimension=net_train_adata.shape[1],
                                      z_dimension=40,
                                      n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                      alpha=1e-5,
                                      beta=100,
                                      eta=100,
                                      clip_value=1e6,
                                      lambda_l1=0.0,
                                      lambda_l2=0.0,
                                      learning_rate=0.001,
                                      model_path=f"./models/MMDCVAE/best/{data_name}-{specific_celltype}/",
                                      dropout_rate=0.2,
                                      output_activation='relu')

    network.train(net_train_adata,
                  net_valid_adata,
                  labelencoder,
                  condition_key,
                  n_epochs=500,
                  batch_size=512,
                  verbose=2,
                  early_stop_limit=50,
                  lr_reducer=35,
                  shuffle=True,
                  save=True,
                  )

    encoder_labels, _ = reptrvae.tl.label_encoder(net_train_adata, labelencoder, condition_key)
    decoder_labels = encoder_labels
    mmd_adata = network.to_mmd_layer(net_train_adata, encoder_labels, decoder_labels)

    cell_type_adata = adata[adata.obs[cell_type_key] == specific_celltype]
    source_adata = cell_type_adata[cell_type_adata.obs[condition_key] == source_condition]
    source_labels = np.zeros(source_adata.shape[0]) + labelencoder[source_condition]
    target_labels = np.zeros(source_adata.shape[0]) + labelencoder[target_condition]

    pred_adata = network.predict(source_adata,
                                 encoder_labels=source_labels,
                                 decoder_labels=target_labels,
                                 )

    pred_adata.obs[condition_key] = [f"{specific_celltype}_pred_{target_condition}"] * pred_adata.shape[0]

    pred_adata.write_h5ad(f"./data/reconstructed/MMDCVAE_{data_name}/{specific_celltype}.h5ad")

#     reptrvae.pl.plot_umap(mmd_adata,
#                           condition_key, cell_type_key,
#                           frameon=False, path_to_save=f"./results/{data_name}/", model_name="MMDCVAE_MMD", ext="png")