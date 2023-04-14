# NN Parameters

## Official

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
                lr_reducer=250,
                shuffle=True,
                save=True,
                retrain=True,
                )

## Low Fitting

network = reptrvae.models.trVAE(x_dimension=net_train_adata.shape[1],
                                z_dimension=10,
                                n_conditions=len(net_train_adata.obs[condition_key].unique()),
                                alpha=1e-5,
                                beta=50,
                                eta=10,
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
                n_epochs=10,
                batch_size=8,
                verbose=2,
                early_stop_limit=0,
                lr_reducer=0,
                shuffle=True,
                save=True,
                retrain=True,
                )
