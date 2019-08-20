import tensorflow as tf
from tensorflow import keras
import pickle
from keras import models
from keras.models import Sequential
from keras.layers import Conv1D
from keras.layers import MaxPooling1D
from keras.layers import Dense
from keras.layers import Dropout

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
model = tf.keras.models.Sequential([
    tf.keras.layers.Conv1D(64, 3, activation = 'relu', input_shape = (1, 1467)),
    tf.keras.MaxPooling1D(2),
    tf.keras.layers.Conv1D(64, 3, activation = 'relu'),
    tf.keras.MaxPooling1D(2),
    tf.keras.layers.Dense(72, activation = 'relu'),
    tf.keras.layers.Dropout(0.2),
    tf.keras.layers.Dense(6, activation = 'softmax')
])
model.summary()

#compile model
from tensorflow.keras.optimizers import RMSprop
model.compile(loss = 'categorical_crossentropy', optimizer = RMSprop(lr = 0.001), metrics = ['acc'])

#fit model to training set and implement validation
history = model.fit(X_train, y_train, epochs = 10, validation_data= (X_CV, y_CV))

#plot training and val accuracy over epochs
import matplotlib.pyplot as plt
acc = history.history['acc']
val_acc = history.history['val_acc']
loss = history.history['loss']
val_loss = history.history['val_loss']

epochs = range(len(acc))

plt.plot(epochs, acc, 'r', label='Training accuracy')
plt.plot(epochs, val_acc, 'b', label='Validation accuracy')
plt.title('Training and validation accuracy')
plt.legend(loc=0)
plt.figure()

plt.show()

#evaluate on test data
test_loss = model.evaluate(X_test, y_test)