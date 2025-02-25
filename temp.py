
import numpy as np
import matplotlib.pyplot as plt
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten, Merge
from keras.layers import Convolution2D, MaxPooling2D
from keras.utils import np_utils

data = np.load("TrainingSet.npy")
#np.random.shuffle(data)
X1 = data[:,:-3]
X2 = data[:,-3]
Y = data[:,-1] / data[:,-2]
X1 = np.reshape(X1, (len(X1),))
X1 = np.vstack(X1)
X1 = X1.reshape(X1.shape[0], 4, 4)
X1 = np.stack([np.real(X1), np.imag(X1)], axis=3)
X1_train = X1[:90000]
X2_train = X2[:90000]
Y_train = Y[:90000]
X1_test = X1[90000:]
X2_test = X2[90000:]
Y_test = Y[90000:]


branch1 = Sequential()
branch1.add(Convolution2D(64, (4,4), activation="elu", input_shape=(4,4,2)))
branch1.add(Flatten())

branch2 = Sequential()
branch2.add(Dense(1, input_dim=1, activation='relu'))

model = Sequential()
model.add(Merge([branch1, branch2], mode = 'concat'))
model.add(Dense(32, activation='elu'))
model.add(Dense(16, activation='elu'))
model.add(Dense(8, activation='elu'))
model.add(Dense(4, activation='elu'))
model.add(Dense(2, activation='elu'))
model.add(Dense(1, activation='sigmoid'))
model.compile(loss='mean_absolute_error',
             optimizer='adam')

n_epochs = 250
h = model.fit([X1_train, X2_train], Y_train, batch_size=32, epochs=n_epochs, verbose=1, validation_split=0.2)
s = model.evaluate([X1_test, X2_test], Y_test, batch_size=20)

dom = np.arange(n_epochs)
plt.plot(dom, h.history['loss'])
plt.xlabel('Epochs')
plt.ylabel('loss')
plt.show()
print(s)

'''

model = Sequential()

model.add(Dense(64, input_dim=33, activation='elu'))
model.add(Dense(32, activation='elu'))
model.add(Dense(16, activation='elu'))
model.add(Dense(8, activation='elu'))
model.add(Dense(4, activation='elu'))
model.add(Dense(2, activation='elu'))
model.add(Dense(1, activation='sigmoid'))
model.compile(loss='mean_absolute_error',
             optimizer='adam')

n_epochs = 120
h = model.fit(XL_train, YL_train, batch_size=32, epochs=n_epochs, verbose=1, validation_split=0.2)
s = model.evaluate(XL_test, YL_test, batch_size=20)

dom = np.arange(n_epochs)
plt.plot(dom, h.history['loss'])
plt.xlabel('Epochs')
plt.ylabel('loss')
plt.show()
print(s)

'''