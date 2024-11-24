from MultimodalSDK_loader.process_dataset import get_feature_matrix

dataset_name = "cmumosi" # eg. "cmumosi"
seq_len = 20 # eg. 20
X_train, y_train, X_val, y_val, X_test, y_test = get_feature_matrix(dataset_name, seq_len) 

print(X_train)