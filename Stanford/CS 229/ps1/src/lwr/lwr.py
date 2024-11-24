import matplotlib.pyplot as plt
import numpy as np
import util


def main(tau, train_path, eval_path):
    """Problem: Locally weighted regression (LWR)

    Args:
        tau: Bandwidth parameter for LWR.
        train_path: Path to CSV file containing dataset for training.
        eval_path: Path to CSV file containing dataset for evaluation.
    """
    # Load training set
    x_train, y_train = util.load_dataset(train_path, add_intercept=True)
    x_eval, y_eval = util.load_dataset(eval_path, add_intercept=True)

    # *** START CODE HERE ***
    model = LocallyWeightedLinearRegression(tau)
    # Fit a LWR model
    model.fit(x_train, y_train)
    y_til = model.predict(x_eval)
    # Get MSE value on the validation set
    MSE = np.mean([(y_til[i]-y_eval[i])**2 for i in range(len(y_eval))])
    print(MSE)
    # Plot validation predictions on top of training set
    plt.plot(x_train[:,1],y_train, "bx", label="train")
    plt.plot(x_eval[:,1],y_til, "ro", label="valid")
    plt.legend(loc="upper left")
    plt.ylabel("y")
    plt.xlabel("x")
    plt.title(fr'$\tau = {tau},\;MSE\approx {MSE}$')
    plt.show()
    # No need to save predictions
    # *** END CODE HERE ***


class LocallyWeightedLinearRegression():
    """Locally Weighted Regression (LWR).

    Example usage:
        > clf = LocallyWeightedLinearRegression(tau)
        > clf.fit(x_train, y_train)
        > clf.predict(x_eval)
    """

    def __init__(self, tau):
        super(LocallyWeightedLinearRegression, self).__init__()
        self.tau = tau
        self.x = None
        self.y = None

    def fit(self, x, y):
        """Fit LWR by saving the training set.

        """
        # *** START CODE HERE ***
        n = len(x)
        self.x = np.array(x).reshape(n,-1)
        self.y = np.array(y).reshape(n,1)
        W = lambda x_input : np.diag([np.exp(-np.linalg.norm(self.x[i,:]-x_input)**2/2/self.tau**2) for i in range(n)])
        self.theta = lambda x_input : np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(self.x.T,W(x_input)),self.x)),np.matmul(self.x.T,W(x_input))),self.y)
        # *** END CODE HERE ***

    def predict(self, x):
        """Make predictions given inputs x.

        Args:
            x: Inputs of shape (n, m).

        Returns:
            Outputs of shape (n,).
        """
        # *** START CODE HERE ***
        n = len(x)
        x = np.array(x).reshape(n,-1)
        y_pred = [x[i,:].dot(self.theta(np.array(x[i,:]))) for i in range(n)]
        return y_pred
        # *** END CODE HERE ***

if __name__ == '__main__':
    main(tau=5e-1,
         train_path='./train.csv',
         eval_path='./valid.csv')
