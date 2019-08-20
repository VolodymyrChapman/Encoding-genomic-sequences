import tensorflow as tf
from tensorflow import keras
import pickle

labelled = pickle.load(open(r'C:\Users\chapmanvl\Documents\VC 2019 projects\Influenza A segment 6 NA\even_scaled_signals_labelled.txt', 'rb'))

# splitting into test/training sets
from sklearn.model_selection import train_test_split

# Splitting data into X and Y
X = labelled.iloc[:, 0:1467]
Y = labelled['colour']

# Splitting into training, CV and test data
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2)
X_train,X_CV,y_train, y_CV = train_test_split(X_train, y_train, test_size = 0.25)

# Checking dims
print(X_train.shape, y_train.shape)
print(X_CV.shape, y_CV.shape)
print(X_test.shape, y_test.shape)

#model
model = tf.keras.models.Sequential([tf.keras.layers.Conv1D(64, 3, activation = 'relu', input_shape = (None, 1467)),
    tf.keras.MaxPooling1D(2),
    tf.keras.layers.Conv1D(64, 3, activation = 'relu'),
    tf.keras.MaxPooling1D(2),
    tf.keras.layers.Dense(36, activation = 'relu'),
    tf.keras.layers.Dropout(0.2),
    tf.keras.layers.Dense(6, activation = 'softmax')
])
model.summary()

#compile model
from tensorflow.keras.optimizers import RMSprop
model.compile(loss = 'categorical_crossentropy', optimizer = RMSprop(lr = 0.001), metrics = ['acc'])


test_loss = model.evaluate(X_test, y_test)