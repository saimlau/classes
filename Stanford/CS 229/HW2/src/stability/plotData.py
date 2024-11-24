import util
import numpy as np
import matplotlib.pyplot as plt

def calc_grad(X, Y, theta):
    """Compute the gradient of the loss with respect to theta."""
    count, _ = X.shape

    probs = 1. / (1 + np.exp(-X.dot(theta)))
    grad = (Y - probs).dot(X)
    
    return grad

def logistic_regression(X, Y):
    """Train a logistic regression model."""
    theta = np.zeros(X.shape[1])
    learning_rate = 0.1

    i = 0
    # while True:
    for k in range(100000):
        i += 1
        prev_theta = theta
        grad = calc_grad(X, Y, theta)
        theta = theta + learning_rate * grad
        if i % 10000 == 0:
            print(np.linalg.norm(theta))
            print('Finished %d iterations' % i)
        if np.linalg.norm(prev_theta - theta) < 1e-15:
            print('Converged in %d iterations' % i)
            break
    return theta

def plot_data(X,Y,num,theta):
    titles = ["Dataset A","Dataset B"]
    markers = ["x","8"]
    colors = ["red","blue"]
    n, _ = X.shape
    plt.subplot(121+num)
    for i in range(n):
        plt.scatter(X[i][1], X[i][2], marker=markers[int(Y[i])], c=colors[int(Y[i])])
    xs = np.linspace(0,np.max(X[:,1]), 1000)
    plt.plot(xs, -(theta[1]*xs+theta[0])/theta[2])
    plt.title(titles[num])
    plt.xlabel("x_1")
    plt.ylabel("x_2")

def main():
    plt.figure()
    print('==== Pllotting data set A ====')
    Xa, Ya = util.load_csv('ds1_a.csv', add_intercept=True)
    theta = logistic_regression(Xa, Ya)
    plot_data(Xa, Ya, 0, theta)

    print('\n==== Plotting data set B ====')
    Xb, Yb = util.load_csv('ds1_b.csv', add_intercept=True)
    theta = logistic_regression(Xb, Yb)
    plot_data(Xb, Yb, 1, theta)
    plt.show()


if __name__ == '__main__':
    main()