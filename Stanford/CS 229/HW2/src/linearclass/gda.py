import numpy as np
import util
import matplotlib.pyplot as plt


def plot_data(X,Y,num,theta):
    plt.figure()
    titles = ["Dataset 1, GDA", "Dataset 2, GDA"]
    markers = ["8","x"]
    colors = ["red","blue"]
    n, _ = X.shape
    for i in range(n):
        plt.scatter(X[i][0], X[i][1], marker=markers[int(Y[i])], c=colors[int(Y[i])])
    xs = np.linspace(np.min(X[:,0]),np.max(X[:,0]), 1000)
    plt.plot(xs, -(theta[1]*xs+theta[0])/theta[2])
    plt.title(titles[num])
    plt.xlabel("x_1")
    plt.ylabel("x_2")
    plt.ylim(np.min(X[:,1])-1,np.max(X[:,1])+1)
    plt.xlim(np.min(X[:,0])-1,np.max(X[:,0])+1)
    plt.show()

def main(train_path, valid_path, save_path, num):
    """Problem: Gaussian discriminant analysis (GDA)

    Args:
        train_path: Path to CSV file containing dataset for training.
        valid_path: Path to CSV file containing dataset for validation.
        save_path: Path to save predicted probabilities using np.savetxt().
    """
    # Load dataset
    x_train, y_train = util.load_dataset(train_path, add_intercept=False)
    x_val, y_val = util.load_dataset(valid_path, add_intercept=False)

    # *** START CODE HERE ***
    model = GDA()
    model.fit(x_train,y_train)
    y_pred = model.predict(x_val)
    print('GDA Accuracy: %.2f' % np.mean( (y_pred == 1) == (y_val == 1)))
    np.savetxt(save_path, y_pred)
    plot_data(x_val,y_val, num, model.theta)
    # *** END CODE HERE ***


class GDA:
    """Gaussian Discriminant Analysis.

    Example usage:
        > clf = GDA()
        > clf.fit(x_train, y_train)
        > clf.predict(x_eval)
    """
    def __init__(self, step_size=0.01, max_iter=10000, eps=1e-5,
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

    def fit(self, x:np.ndarray, y):
        """Fit a GDA model to training set given by x and y by updating
        self.theta.

        Args:
            x: Training example inputs. Shape (n_examples, dim).
            y: Training example labels. Shape (n_examples,).
        """
        # *** START CODE HERE ***
        n,m = x.shape
        phi = sum(y)/n
        mu = {}
        mu[0] = np.sum(x[y==0, :], axis=0)/(n-sum(y))
        mu[1] = np.sum(x[y==1, :], axis=0)/sum(y)
        Sig = np.zeros((m,m))
        for i in range(n):
            Sig += np.outer(x[i,:]-mu[int(y[i])],x[i,:]-mu[int(y[i])])
        Sig = Sig/n
        self.theta = np.zeros(m+1)
        sig_inv = np.linalg.inv(Sig)
        self.theta[1:] = (sig_inv@((mu[1]-mu[0]).reshape(m,1))).flatten()
        self.theta[0] = 1/2*(mu[0].dot(np.matmul(sig_inv,mu[0].reshape(m,1)))-mu[1].dot(np.matmul(sig_inv,mu[1].reshape(m,1))))-np.log(1-phi)+np.log(phi)
        # *** END CODE HERE ***

    def predict(self, x):
        """Make a prediction given new inputs x.

        Args:
            x: Inputs of shape (n_examples, dim).

        Returns:
            Outputs of shape (n_examples,).
        """
        # *** START CODE HERE ***
        p = np.array([1/(1+np.exp(-x[i,:].dot(self.theta[1:])-self.theta[0])) for i in range(len(x))])
        return np.array(p>=0.5).astype(int)
        # *** END CODE HERE

if __name__ == '__main__':
    main(train_path='ds1_train.csv',
         valid_path='ds1_valid.csv',
         save_path='gda_pred_1.txt', num=0)

    main(train_path='ds2_train.csv',
         valid_path='ds2_valid.csv',
         save_path='gda_pred_2.txt', num=1)
