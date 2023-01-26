# Model Training Parameters

## AlzPro

### For Normalized Data

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
                early_stop_limit=300,
                lr_reducer=250,
                shuffle=True,
                save=True,
                retrain=True,
                )

### Reduced Normalized Data

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
                n_epochs=1000,
                batch_size=16,
                verbose=2,
                early_stop_limit=300,
                lr_reducer=150,
                shuffle=True,
                save=True,
                retrain=True,
                )

### For Count Data

#### Normal 1

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
                n_epochs=10000,
                batch_size=4,
                verbose=2,
                early_stop_limit=500,
                lr_reducer=40,
                shuffle=True,
                save=True,
                retrain=True,
                )

#### Normal 2

network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                z_dimension=40,
                                n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                alpha=5e-5,
                                beta=500,
                                eta=100,
                                clip_value=1e6,
                                lambda_l1=0.0,
                                lambda_l2=0.0,
                                learning_rate=0.00001,
                                model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                dropout_rate=0.2,
                                output_activation='relu')

network.train(net_train_adata,
                net_valid_adata,
                labelencoder,
                condition_key,
                n_epochs=10000,
                batch_size=8,
                verbose=2,
                early_stop_limit=300,
                lr_reducer=150,
                shuffle=True,
                save=True,
                retrain=True,
                )

#### Normal 2B

network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                z_dimension=40,
                                n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                alpha=5e-5,
                                beta=500,
                                eta=100,
                                clip_value=1e6,
                                lambda_l1=0.0,
                                lambda_l2=0.0,
                                learning_rate=0.00001,
                                model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                dropout_rate=0.2,
                                output_activation='relu')

network.train(net_train_adata,
                net_valid_adata,
                labelencoder,
                condition_key,
                n_epochs=10000,
                batch_size=8,
                verbose=2,
                early_stop_limit=600,
                lr_reducer=300,
                shuffle=True,
                save=True,
                retrain=True,
                )

#### Janky 1

network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                z_dimension=40,
                                n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                alpha=5e-5,
                                beta=500,
                                eta=100,
                                clip_value=1e6,
                                lambda_l1=0.0,
                                lambda_l2=0.0,
                                learning_rate=0.1,
                                model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                dropout_rate=0.2,
                                output_activation='relu')

network.train(net_train_adata,
                net_valid_adata,
                labelencoder,
                condition_key,
                n_epochs=1000,
                batch_size=2,
                verbose=2,
                early_stop_limit=300,
                lr_reducer=100,
                shuffle=True,
                save=True,
                retrain=True,
                )

#### Janky 2

network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                z_dimension=40,
                                n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                alpha=5e-5,
                                beta=500,
                                eta=100,
                                clip_value=1e6,
                                lambda_l1=0.0,
                                lambda_l2=0.0,
                                learning_rate=0.00001,
                                model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                dropout_rate=0.2,
                                output_activation='relu')

network.train(net_train_adata,
                net_valid_adata,
                labelencoder,
                condition_key,
                n_epochs=1000,
                batch_size=2,
                verbose=2,
                early_stop_limit=30,
                lr_reducer=25,
                shuffle=True,
                save=True,
                retrain=True,
                )

### Janky 3

network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                z_dimension=40,
                                n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                alpha=5e-5,
                                beta=500,
                                eta=100,
                                clip_value=1e6,
                                lambda_l1=0.0,
                                lambda_l2=0.0,
                                learning_rate=0.0001,
                                model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                dropout_rate=0.2,
                                output_activation='relu')

network.train(net_train_adata,
                net_valid_adata,
                labelencoder,
                condition_key,
                n_epochs=1000,
                batch_size=2,
                verbose=2,
                early_stop_limit=30,
                lr_reducer=25,
                shuffle=True,
                save=True,
                retrain=True,
                )

2272712
361786

#### 16o 1

network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                z_dimension=40,
                                n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                alpha=5e-5,
                                beta=500,
                                eta=100,
                                clip_value=1e6,
                                lambda_l1=0.0,
                                lambda_l2=0.0,
                                learning_rate=0.00001,
                                model_path=f"./models/trVAE/best/{data_name}-{specific_celltype}/",
                                dropout_rate=0.2,
                                output_activation='relu')

network.train(net_train_adata,
                net_valid_adata,
                labelencoder,
                condition_key,
                n_epochs=10000,
                batch_size=16,
                verbose=2,
                early_stop_limit=300,
                lr_reducer=150,
                shuffle=True,
                save=True,
                retrain=True,
                )

#### 32o 1

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
                n_epochs=1000,
                batch_size=32,
                verbose=2,
                early_stop_limit=400,
                lr_reducer=100,
                shuffle=True,
                save=True,
                retrain=True,
                )
