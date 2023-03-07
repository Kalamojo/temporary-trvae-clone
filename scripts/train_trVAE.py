import sys
import numpy as np
import scanpy as sc

import reptrvae

data_name = sys.argv[1]
# specific_cell_type = sys.argv[2]

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
    dname = data_name

elif data_name == "alzPro" or data_name == "alzPho":
    conditions = ["WT", "HET"]
    target_conditions = ["HET"]
    source_condition = "WT"
    target_condition = "HET"
    labelencoder = {"WT": 0, "HET": 1}
    cell_type_key = "Timepoint"
    condition_key = "Group"
    if len(sys.argv) == 3:
        specific_celltype = sys.argv[2]
    else:
        specific_celltype = "3m"
    dname = data_name

elif data_name == "alzProTrain" or data_name == "alzPhoTrain":
    conditions = ["WT", "HET"]
    target_conditions = ["HET"]
    source_condition = "WT"
    target_condition = "HET"
    labelencoder = {"WT": 0, "HET": 1}
    cell_type_key = "Validation"
    condition_key = "Group"
    if len(sys.argv) == 3:
        specific_celltype = sys.argv[2]
    else:
        specific_celltype = "Test"
    dname = data_name[:-5]

elif data_name == "alzProTime" or data_name == "alzPhoTime":
    conditions = ["3m", "6m", "9m"]
    target_conditions = ["9m"]
    source_condition = "9m"
    target_condition = "9m"
    labelencoder = {"3m": 0, "6m": 1, "9m": 2}
    cell_type_key = "Group"
    condition_key = "Timepoint"
    if len(sys.argv) == 3:
        specific_celltype = sys.argv[2]
    else:
        specific_celltype = "HET"
    dname = data_name[:-4]

elif data_name == "alzProTT" or data_name == "alzPhoTT":
    conditions = ["3m", "6m", "9m"]
    target_conditions = ["3m"]
    source_condition = "3m"
    target_condition = "3m"
    labelencoder = {"3m": 0, "6m": 1, "9m": 2}
    cell_type_key = "Validation"
    condition_key = "Timepoint"
    if len(sys.argv) == 3:
        specific_celltype = sys.argv[2]
    else:
        specific_celltype = "Test"
    dname = data_name[:-2]

elif data_name == "alzProSex":
    conditions = ["M", "F"]
    target_conditions = ["F"]
    source_condition = "M"
    target_condition = "F"
    labelencoder = {"M": 0, "F": 1}
    cell_type_key = "Timepoint"
    condition_key = "sex"
    if len(sys.argv) == 3:
        specific_celltype = sys.argv[2]
    else:
        specific_celltype = "3m"
    dname = data_name[:-3]

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
    dname = data_name
else:
    raise Exception("InValid data name")

adata = sc.read(f"./data/{dname}_normalized.h5ad")
#adata = sc.read(f"./data/{data_name}_count.h5ad")
adata = adata[adata.obs[condition_key].isin(conditions)]

#if adata.shape[1] > 2000:
#    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
#    adata = adata[:, adata.var['highly_variable']]


print(specific_celltype, sys.argv)
if specific_celltype == 'all':
    train_adata, valid_adata = reptrvae.utils.train_test_split(adata, 0.80)
    for specific_celltype in adata.obs[cell_type_key].unique().tolist():
        print("""
                VVV
                VVV
                VVV
                VVV
                VVV
                VVV""")
        print(specific_celltype)
        net_train_adata = train_adata[
            ~((train_adata.obs[cell_type_key] == specific_celltype) & (train_adata.obs[condition_key].isin(target_conditions)))]
        net_valid_adata = valid_adata[
            ~((valid_adata.obs[cell_type_key] == specific_celltype) & (valid_adata.obs[condition_key].isin(target_conditions)))]

        network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                        z_dimension=40,
                                        n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                        alpha=5e-5,
                                        beta=500,
                                        eta=100,
                                        clip_value=1e6,
                                        lambda_l1=0.0,
                                        lambda_l2=0.0,
                                        learning_rate=0.00005,
                                        model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                        dropout_rate=0.2,
                                        output_activation='relu')

        network.train(net_train_adata,
                        net_valid_adata,
                        labelencoder,
                        condition_key,
                        n_epochs=2000,
                        batch_size=4,
                        verbose=2,
                        early_stop_limit=500,
                        lr_reducer=250,
                        shuffle=True,
                        save=True,
                        retrain=True,
                        )

        train_labels, _ = reptrvae.tl.label_encoder(net_train_adata, labelencoder, condition_key)
        mmd_adata = network.to_mmd_layer(net_train_adata, train_labels, feed_fake=-1)

        cell_type_adata = adata[adata.obs[cell_type_key] == specific_celltype]
        source_adata = cell_type_adata[cell_type_adata.obs[condition_key] == source_condition]
        target_adata = cell_type_adata[cell_type_adata.obs[condition_key] == target_condition]
        source_labels = np.zeros(source_adata.shape[0]) + labelencoder[source_condition]
        target_labels = np.zeros(source_adata.shape[0]) + labelencoder[target_condition]

        pred_adata = network.predict(source_adata,
                                     encoder_labels=source_labels,
                                     decoder_labels=target_labels,
                                     )

        pred_adata.obs[condition_key] = [f"{source_condition}_to_{target_condition}"] * pred_adata.shape[0]
        pred_adata.obs[cell_type_key] = specific_celltype

        adata_to_write = pred_adata.concatenate(target_adata)
        adata_to_write.write_h5ad(f"./data/reconstructed/trVAE_{data_name}/{specific_celltype}.h5ad")
        # reptrvae.pl.plot_umap(mmd_adata,
        #                       condition_key, cell_type_key,
        #                       frameon=False, path_to_save=f"./results/{data_name}/", model_name="trVAE_MMD",
        #                       ext="png")
elif specific_celltype == '3m1':
    specific_celltype = specific_celltype[:-1]
    indices = np.arange(adata.shape[0])
    train_idx = np.concatenate([indices[:0], indices[16:]])
    test_idx = indices[0:16]
    train_adata = adata[train_idx, :]
    valid_adata = adata[test_idx, :]
    print("""
            VVV
            VVV
            VVV
            VVV
            VVV
            VVV""")
    print(specific_celltype)
    net_train_adata = train_adata[
        ~((train_adata.obs[cell_type_key] == specific_celltype) & (train_adata.obs[condition_key].isin(target_conditions)))]
    net_valid_adata = valid_adata[
        ~((valid_adata.obs[cell_type_key] == specific_celltype) & (valid_adata.obs[condition_key].isin(target_conditions)))]

    network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                    z_dimension=40,
                                    n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                    alpha=5e-5,
                                    beta=500,
                                    eta=100,
                                    clip_value=1e6,
                                    lambda_l1=0.0,
                                    lambda_l2=0.0,
                                    learning_rate=0.00005,
                                    model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                    dropout_rate=0.2,
                                    output_activation='relu')

    network.train(net_train_adata,
                    net_valid_adata,
                    labelencoder,
                    condition_key,
                    n_epochs=2000,
                    batch_size=4,
                    verbose=2,
                    early_stop_limit=500,
                    lr_reducer=250,
                    shuffle=True,
                    save=True,
                    retrain=True,
                    )

    train_labels, _ = reptrvae.tl.label_encoder(net_train_adata, labelencoder, condition_key)
    mmd_adata = network.to_mmd_layer(net_train_adata, train_labels, feed_fake=-1)

    cell_type_adata = adata[adata.obs[cell_type_key] == specific_celltype]
    source_adata = cell_type_adata[cell_type_adata.obs[condition_key] == source_condition]
    target_adata = cell_type_adata[cell_type_adata.obs[condition_key] == target_condition]
    source_labels = np.zeros(source_adata.shape[0]) + labelencoder[source_condition]
    target_labels = np.zeros(source_adata.shape[0]) + labelencoder[target_condition]

    pred_adata = network.predict(source_adata,
                                    encoder_labels=source_labels,
                                    decoder_labels=target_labels,
                                    )
    
    network.get_corrected(source_adata, labels=source_labels, return_z=True)

    pred_adata.obs[condition_key] = [f"{source_condition}_to_{target_condition}"] * pred_adata.shape[0]
    pred_adata.obs[cell_type_key] = specific_celltype

    adata_to_write = pred_adata.concatenate(target_adata)
    adata_to_write.write_h5ad(f"./data/reconstructed/trVAE_{data_name}/{specific_celltype}.h5ad")

    source_adata.write_h5ad(f"./data/reconstructed/trVAE_{data_name}/{specific_celltype}_latent.h5ad")

    reptrvae.pl.plot_umap(mmd_adata,
                          condition_key, cell_type_key,
                          frameon=False, path_to_save=f"./results/{data_name}/", model_name="trVAE_MMD",
                          ext="png")
elif specific_celltype == '6m1':
    specific_celltype = specific_celltype[:-1]
    indices = np.arange(adata.shape[0])
    train_idx = np.concatenate([indices[:16], indices[32:]])
    test_idx = indices[16:32]
    train_adata = adata[train_idx, :]
    valid_adata = adata[test_idx, :]
    print("""
            VVV
            VVV
            VVV
            VVV
            VVV
            VVV""")
    print(specific_celltype)
    net_train_adata = train_adata[
        ~((train_adata.obs[cell_type_key] == specific_celltype) & (train_adata.obs[condition_key].isin(target_conditions)))]
    net_valid_adata = valid_adata[
        ~((valid_adata.obs[cell_type_key] == specific_celltype) & (valid_adata.obs[condition_key].isin(target_conditions)))]

    network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                    z_dimension=40,
                                    n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                    alpha=5e-5,
                                    beta=500,
                                    eta=100,
                                    clip_value=1e6,
                                    lambda_l1=0.0,
                                    lambda_l2=0.0,
                                    learning_rate=0.001,
                                    model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                    dropout_rate=0.2,
                                    output_activation='relu')

    network.train(net_train_adata,
                    net_valid_adata,
                    labelencoder,
                    condition_key,
                    n_epochs=100,
                    batch_size=4,
                    verbose=2,
                    early_stop_limit=0,
                    lr_reducer=40,
                    shuffle=True,
                    save=True,
                    retrain=True,
                    )

    train_labels, _ = reptrvae.tl.label_encoder(net_train_adata, labelencoder, condition_key)
    mmd_adata = network.to_mmd_layer(net_train_adata, train_labels, feed_fake=-1)

    cell_type_adata = adata[adata.obs[cell_type_key] == specific_celltype]
    source_adata = cell_type_adata[cell_type_adata.obs[condition_key] == source_condition]
    target_adata = cell_type_adata[cell_type_adata.obs[condition_key] == target_condition]
    source_labels = np.zeros(source_adata.shape[0]) + labelencoder[source_condition]
    target_labels = np.zeros(source_adata.shape[0]) + labelencoder[target_condition]

    pred_adata = network.predict(source_adata,
                                    encoder_labels=source_labels,
                                    decoder_labels=target_labels,
                                    )

    pred_adata.obs[condition_key] = [f"{source_condition}_to_{target_condition}"] * pred_adata.shape[0]
    pred_adata.obs[cell_type_key] = specific_celltype

    adata_to_write = pred_adata.concatenate(target_adata)
    adata_to_write.write_h5ad(f"./data/reconstructed/trVAE_{data_name}/{specific_celltype}.h5ad")
    # reptrvae.pl.plot_umap(mmd_adata,
    #                       condition_key, cell_type_key,
    #                       frameon=False, path_to_save=f"./results/{data_name}/", model_name="trVAE_MMD",
    #                       ext="png")
elif specific_celltype == '9m1':
    specific_celltype = specific_celltype[:-1]
    indices = np.arange(adata.shape[0])
    train_idx = np.concatenate([indices[:32], indices[48:]])
    test_idx = indices[32:48]
    train_adata = adata[train_idx, :]
    valid_adata = adata[test_idx, :]
    print("""
            VVV
            VVV
            VVV
            VVV
            VVV
            VVV""")
    print(specific_celltype)
    net_train_adata = train_adata[
        ~((train_adata.obs[cell_type_key] == specific_celltype) & (train_adata.obs[condition_key].isin(target_conditions)))]
    net_valid_adata = valid_adata[
        ~((valid_adata.obs[cell_type_key] == specific_celltype) & (valid_adata.obs[condition_key].isin(target_conditions)))]

    network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                    z_dimension=40,
                                    n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                    alpha=5e-5,
                                    beta=500,
                                    eta=100,
                                    clip_value=1e6,
                                    lambda_l1=0.0,
                                    lambda_l2=0.0,
                                    learning_rate=0.001,
                                    model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                    dropout_rate=0.2,
                                    output_activation='relu')

    network.train(net_train_adata,
                    net_valid_adata,
                    labelencoder,
                    condition_key,
                    n_epochs=100,
                    batch_size=4,
                    verbose=2,
                    early_stop_limit=0,
                    lr_reducer=40,
                    shuffle=True,
                    save=True,
                    retrain=True,
                    )

    train_labels, _ = reptrvae.tl.label_encoder(net_train_adata, labelencoder, condition_key)
    mmd_adata = network.to_mmd_layer(net_train_adata, train_labels, feed_fake=-1)

    cell_type_adata = adata[adata.obs[cell_type_key] == specific_celltype]
    source_adata = cell_type_adata[cell_type_adata.obs[condition_key] == source_condition]
    target_adata = cell_type_adata[cell_type_adata.obs[condition_key] == target_condition]
    source_labels = np.zeros(source_adata.shape[0]) + labelencoder[source_condition]
    target_labels = np.zeros(source_adata.shape[0]) + labelencoder[target_condition]

    pred_adata = network.predict(source_adata,
                                    encoder_labels=source_labels,
                                    decoder_labels=target_labels,
                                    )

    pred_adata.obs[condition_key] = [f"{source_condition}_to_{target_condition}"] * pred_adata.shape[0]
    pred_adata.obs[cell_type_key] = specific_celltype

    adata_to_write = pred_adata.concatenate(target_adata)
    adata_to_write.write_h5ad(f"./data/reconstructed/trVAE_{data_name}/{specific_celltype}.h5ad")
    # reptrvae.pl.plot_umap(mmd_adata,
    #                       condition_key, cell_type_key,
    #                       frameon=False, path_to_save=f"./results/{data_name}/", model_name="trVAE_MMD",
    #                       ext="png")
else:
    train_adata, valid_adata = reptrvae.utils.train_test_split(adata, 0.80)
    net_train_adata = train_adata[
        ~((train_adata.obs[cell_type_key] == specific_celltype) & (train_adata.obs[condition_key].isin(target_conditions)))]
    net_valid_adata = valid_adata[
        ~((valid_adata.obs[cell_type_key] == specific_celltype) & (valid_adata.obs[condition_key].isin(target_conditions)))]

    network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                    z_dimension=40,
                                    n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                    alpha=5e-5,
                                    beta=500,
                                    eta=100,
                                    clip_value=1e6,
                                    lambda_l1=0.0,
                                    lambda_l2=0.0,
                                    learning_rate=0.00005,
                                    model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                    dropout_rate=0.2,
                                    output_activation='relu')
    #"""
    network.train(net_train_adata,
                  net_valid_adata,
                  labelencoder,
                  condition_key,
                  n_epochs=2000,
                  batch_size=4,
                  verbose=2,
                  early_stop_limit=500,
                  lr_reducer=250,
                  shuffle=True,
                  save=True,
                  retrain=True,
                  )
    #"""

    train_labels, _ = reptrvae.tl.label_encoder(net_train_adata, labelencoder, condition_key)
    mmd_adata = network.to_mmd_layer(net_train_adata, train_labels, feed_fake=-1)

    cell_type_adata = adata[adata.obs[cell_type_key] == specific_celltype]
    source_adata = cell_type_adata[cell_type_adata.obs[condition_key] == source_condition]
    target_adata = cell_type_adata[cell_type_adata.obs[condition_key] == target_condition]
    source_labels = np.zeros(source_adata.shape[0]) + labelencoder[source_condition]
    target_labels = np.zeros(source_adata.shape[0]) + labelencoder[target_condition]

    pred_adata = network.predict(source_adata,
                                 encoder_labels=source_labels,
                                 decoder_labels=target_labels,
                                 )

    pred_adata.obs[condition_key] = [f"{source_condition}_to_{target_condition}"] * pred_adata.shape[0]
    pred_adata.obs[cell_type_key] = specific_celltype

    adata_to_write = pred_adata.concatenate(target_adata)
    adata_to_write.write_h5ad(f"./data/reconstructed/trVAE_{data_name}/{specific_celltype}.h5ad")
    
    reptrvae.pl.plot_umap(mmd_adata,
                           condition_key, cell_type_key,
                           frameon=False, path_to_save=f"./results/{data_name}/", model_name="trVAE_MMD",
                           ext="png")
print("All done")