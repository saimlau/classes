import matplotlib.pyplot as plt
import numpy as np
import util

from lwr import LocallyWeightedLinearRegression


def main(tau_values, train_path, valid_path, test_path, pred_path):
    """Problem: Tune the bandwidth paramater tau for LWR.

    Args:
        tau_values: List of tau values to try.
        train_path: Path to CSV file containing training set.
        valid_path: Path to CSV file containing validation set.
        test_path: Path to CSV file containing test set.
        pred_path: Path to save predictions.
    """
    # Load training set
    x_train, y_train = util.load_dataset(train_path, add_intercept=True)
    x_valid, y_valid = util.load_dataset(valid_path, add_intercept=True)
    x_test, y_test = util.load_dataset(test_path, add_intercept=True)

    # *** START CODE HERE ***
    # Search tau_values for the best tau (lowest MSE on the validation set)
    # Fit a LWR model with the best tau value
    # Run on the test set to get the MSE value
    # Save predictions to pred_path
    # Plot data
    MSE_values = np.ones(len(tau_values))
    models = []
    for i in range(len(tau_values)):
        plt.close()
        models.append(LocallyWeightedLinearRegression(tau_values[i]))
        models[i].fit(x_train, y_train)
        y_til = models[i].predict(x_valid)
        MSE_values[i] = np.mean([(y_til[i]-y_valid[i])**2 for i in range(len(y_valid))])
        plt.plot(x_train[:,1],y_train, "bx", label="train")
        plt.plot(x_valid[:,1],y_til, "ro", label="valid")
        plt.legend(loc="upper left")
        plt.ylabel("y")
        plt.xlabel("x")
        plt.title(fr'$\tau = {tau_values[i]},\;MSE\approx {MSE_values[i]}$')
        plt.show()
    indx = min(range(len(MSE_values)), key=lambda i: MSE_values[i])
    print(indx)
    y_test_til = models[indx].predict(x_test)
    MSE_test = np.mean([(y_test_til[i]-y_test[i])**2 for i in range(len(y_test))])
    print(MSE_test)
    plt.plot(x_train[:,1],y_train, "bx", label="train")
    plt.plot(x_test[:,1],y_test_til, "ro", label="test")
    plt.legend(loc="upper left")
    plt.ylabel("y")
    plt.xlabel("x")
    plt.title(fr'$\tau = {tau_values[indx]},\;MSE\approx {MSE_test}$')
    plt.show()
    # *** END CODE HERE ***

if __name__ == '__main__':
    main(tau_values=[3e-2, 5e-2, 1e-1, 5e-1, 1e0, 1e1],
         train_path='./train.csv',
         valid_path='./valid.csv',
         test_path='./test.csv',
         pred_path='./pred.txt')
