import numpy as np
import util
import matplotlib.pyplot as plt

def main(lr, train_path, eval_path, save_path):
    """Problem: Poisson regression with gradient ascent.

    Args:
        lr: Learning rate for gradient ascent.
        train_path: Path to CSV file containing dataset for training.
        eval_path: Path to CSV file containing dataset for evaluation.
        save_path: Path to save predictions.
    """
    # Load training set
    x_train, y_train = util.load_dataset(train_path, add_intercept=True)
    x_valid, y_valid = util.load_dataset(eval_path, add_intercept=True)

    # *** START CODE HERE ***
    # Fit a Poisson Regression model
    # Run on the validation set, and use np.savetxt to save outputs to save_path
    model = PoissonRegression()
    model.fit(x_train,y_train)
    y_pred = model.predict(x_valid)
    plt.scatter(y_valid, y_pred)
    plt.xlabel('true counts')
    plt.ylabel('predicted counts')
    plt.xlim((0,26))
    plt.ylim((0,26))
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.savefig('Pos4d.png')
    plt.show()
    np.savetxt(save_path, y_pred)
    # *** END CODE HERE ***


class PoissonRegression:
    """Poisson Regression.

    Example usage:
        > clf = PoissonRegression(step_size=lr)
        > clf.fit(x_train, y_train)
        > clf.predict(x_eval)
    """

    def __init__(self, step_size=1e-5, max_iter=10000000, eps=1e-5,
                 theta_0=None, verbose=True):
        """
        Args:
            step_size: Step size for iterative solvers only.
            max_iter: Maximum number of iterations for the solver.
            eps: Threshold for determining convergence.
            theta_0: Initial guess for theta. If None, use the zero vector.
            verbose: Print loss values during training.
        """
        self.theta = theta_0
        self.step_size = step_size
        self.max_iter = max_iter
        self.eps = eps
        self.verbose = verbose

    def fit(self, x, y):
        """Run gradient ascent to maximize likelihood for Poisson regression.
        Update the parameter by step_size * (sum of the gradient over examples)

        Args:
            x: Training example inputs. Shape (n_examples, dim).
            y: Training example labels. Shape (n_examples,).
        """
        # *** START CODE HERE ***
        n = len(x)
        x = np.array(x); y = np.array(y).reshape(n,1)
        dim = len(x[0])
        self.theta = np.zeros((dim,1))
        old_theta = self.theta.copy()
        for i in range(self.max_iter):
            tem1 = np.exp(np.matmul(x,self.theta)).reshape(n,1)
            for j in range(dim):
                self.theta[j] += self.step_size*sum(x[:,j].reshape(n,1)*(y-tem1))[0]
            if np.linalg.norm(self.theta-old_theta)<self.eps:
                return
            old_theta = self.theta.copy()
        print("Max iterations limit reached.")
        # *** END CODE HERE ***

    def predict(self, x):
        """Make a prediction given inputs x.

        Args:
            x: Inputs of shape (n_examples, dim).

        Returns:
            Floating-point prediction for each input, shape (n_examples,).
        """
        # *** START CODE HERE ***
        y_pred = np.exp(np.matmul(x,self.theta))
        return y_pred.flatten()
        # *** END CODE HERE ***

if __name__ == '__main__':
    main(lr=1e-5,
        train_path='train.csv',
        eval_path='valid.csv',
        save_path='poisson_pred.txt')
