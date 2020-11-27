#DL models


# Import these libraries
from sklearn.neural_network import MLPRegressor
import numpy as np
import pandas as pd
from keras.models import Sequential, load_model
from keras.layers import Dense, Activation, Dropout
from keras.layers import Flatten, Conv1D, MaxPooling1D, LSTM
from keras.activations import relu, elu, linear, softmax
from keras.callbacks import EarlyStopping, Callback
from keras.optimizers import adam, Nadam, sgd
from keras.losses import mean_squared_error, categorical_crossentropy, logcosh
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from keras.wrappers.scikit_learn import KerasRegressor
from tensorflow.keras.callbacks import EarlyStopping
import time



# Data processing
def data_preprocess(dataset,target): 
    # Scaled data
    min_max_scaler = MinMaxScaler(feature_range = (0,1))
    np_scaled = min_max_scaler.fit_transform(dataset)
    X = pd.DataFrame(np_scaled)
    
    target_edit = pd.Series(target).values
    target_edit = target_edit.reshape(-1,1)
    np_scaled = min_max_scaler.fit_transform(target_edit)
    Y = pd.DataFrame(np_scaled)
    # Split data
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X,Y, test_size=0.2)
      
    scaler = StandardScaler()
    # Fit only to the training data
    scaler.fit(X_train)
        
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)
    
    X_train = pd.DataFrame(X_train)
    X_test = pd.DataFrame(X_test)
    return X_train, y_train, X_test, y_test

def cal_correlation(pred, y_test, target):
    pred = pred.reshape((-1,1))
    target_edit = pd.Series(target).values
    target_edit = target_edit.reshape(-1,1)
    min_max_scaler = MinMaxScaler(feature_range = (0,1))
    np_scaled = min_max_scaler.fit_transform(target_edit)
    target_pred = min_max_scaler.inverse_transform(pred)
    target_orig = min_max_scaler.inverse_transform(y_test)
    target_orig = target_orig[:,0]
    target_orig = pd.Series(target_orig)
    target_pred = target_pred[:,0]
    target_pred = pd.Series(target_pred)
    cor1 = target_orig.corr(target_pred, method='pearson')
    return cor1

# Multilayer perceptron
def multi_layer_perceptron(X_train, y_train, X_test, y_test):
    mlp = MLPRegressor(max_iter=200, early_stopping = True)   #---------
    #parameter spce for mlp
    parameter_space = {
        'hidden_layer_sizes': [(19,19,19), (19,38,19), (19, 38, 38, 19), (20, 20, 40, 40, 20), (38,38,38,19),  (50, 50, 38), (90,90,90), (120,90,90)],
        'activation': ['tanh', 'relu', 'linear', 'identity', 'logistic'],
        'solver': ['sgd', 'adam','lbfgs'],
        'alpha': [0.001, 0.05, 0.4], # regularization
        'learning_rate': ['constant','adaptive']
    }
    
    # Grid search with 5 fold cross validation
    clf = GridSearchCV(mlp, parameter_space, n_jobs=-1, cv=5)
    clf.fit(X_train, y_train.values.ravel())
    
    # Best parameter set
    print('Best parameters found:\n', clf.best_params_)
        
    # All results
    means = clf.cv_results_['mean_test_score']
    stds = clf.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, params))
    
    print(clf.best_params_)   
    
    # Predict for test data
    pred = clf.predict(X_test)
    scores = (clf.score(X_test,y_test))
    return scores, pred


# create cnn model to call for grid search
def baseline_model(nSNP):
    def create_conv_NN():
        nStride=3  # stride between convolutions
        nFilter=64 # filters
        model_cnn = Sequential()
    
        # add convolutional layer with l1 and l2 regularization
        model_cnn.add(Conv1D(nFilter, kernel_size=3, strides=nStride, input_shape=(nSNP,1), kernel_regularizer='l1_l2'))
        model_cnn.add(Conv1D(nFilter, kernel_size=3, activation='relu'))
        # dropout added for regularization
        model_cnn.add(Dropout(0.2))
        # add pooling layer: takes maximum of two consecutive values
        model_cnn.add(MaxPooling1D(pool_size=2))
        # Solutions above are linearized to accommodate a standard layer
        model_cnn.add(Flatten())
        model_cnn.add(Dense(64))
        # activation layer
        model_cnn.add(Activation('relu'))
        model_cnn.add(Dense(32))
        model_cnn.add(Activation('linear'))
        model_cnn.add(Dense(1)) 
        
        # Model Compiling 
        model_cnn.compile(loss='mean_squared_error', optimizer='adam')
        return model_cnn
    return create_conv_NN

# Convolutional neural network
def cnn(X_train, y_train, X_test, y_test):
    # need this to match dimensions
    X2_train = np.expand_dims(X_train, axis=2) 
    X2_test = np.expand_dims(X_test, axis=2) 
    
    nSNP=X_train.shape[1] 
    early_stopping = EarlyStopping()
     
    #build cnn regressor
    cnn = KerasRegressor(build_fn=baseline_model(nSNP), verbose=1)
    #hyper parameters
    batch_size = [64, 128]
    epochs = [150, 200]
    
    #parameter space
    param_grid = dict(epochs=epochs, batch_size = batch_size)
    # Grid search with 5 fold cross validation
    grid = GridSearchCV(estimator=cnn, param_grid=param_grid, cv=5)
    #fit training data
    grid_result = grid.fit(X2_train, y_train, validation_split=0.2, callbacks=[early_stopping])
    grid_result = grid.fit(X2_train, y_train)
    #get best model
    best_params = grid_result.best_params_
    print(best_params)
    
    # predict for test set     
    pred = grid.predict(X2_test)
    score = (grid.score(X2_test,y_test))
    
    return score, pred
    
# calculate average of items in a list
def Average(lst): 
    return sum(lst) / len(lst)


def main():
    #start time
    start = time.time()
    
    data = pd.read_csv('Filtered_Data.csv', header=0)
    dataset=data.iloc[:,60:]
    target=data.iloc[:, 4]
    
    #preprocess data
    X_train, y_train, X_test, y_test = data_preprocess(dataset, target)
    
    #store all cor and accuracy values for mlp
    cor_mlp = []
    acc_mlp = []
    #store all cor and accuracy values for mlp
    cor_cnn = []
    acc_cnn = []
    
    #200 iterations
    for i in range(200):
        # Multilayer perceptron
        print('MLP:')
        scores_mlp, pred_mlp = multi_layer_perceptron(X_train, y_train, X_test, y_test)
        cor_m = cal_correlation(pred_mlp, y_test, target)
        print('Correlation for MLP ',cor_m)
        cor_mlp.append(cor_m)
        acc_mlp.append(scores_mlp)
    
        # CNN
        print('CNN:')
        scores_cnn, pred_cnn = cnn(X_train, y_train, X_test, y_test)
        cor_c = cal_correlation(pred_cnn, y_test, target)
        print('Correlation for CNN ',cor_c)
        cor_cnn.append(cor_c)
        acc_cnn.append(scores_mlp)
    
    #average values
    print("Average accuracy for MLP: ", Average(acc_mlp))
    print("Average correlation for MLP: ", Average(cor_mlp))
    print("Average accuracy for CNN: ", Average(acc_cnn))
    print("Average correlation for CNN: ", Average(cor_cnn))
    
    #end time
    end = time.time()
    print("Time elapsed: ",end - start)
    
main()