import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d
def create_edge(A, node_pairs, rest_length = 1):
    """
    Function that creates edges between nodes 
    input:
        - A: adjacency matrix
        - node_pairs: np.array with nodes that have an edge 
        - rest_length: rest length of edges, assume that all pairs have the same rest length
    """
    A[node_pairs[:,0],node_pairs[:,1]] = rest_length
    A[node_pairs[:,1],node_pairs[:,0]] = rest_length

def plot_cable_net(X,A,ax):
    """ 
    function that plots a cable network
    input:
        - X: matrix with nodes and their positions
        - A: adjacency matrix
        - ax: figure where we plot the network 
    output:
        - plot of cable network
    """
    N = np.shape(X)[1]
    for i in range(N):
        for j in range(i,N):
            if(A[i,j] != 0): 
                x = [X[:,i][0], X[:,j][0]]
                y = [X[:,i][1], X[:,j][1]]
                z = [X[:,i][2], X[:,j][2]]
                ax.scatter(x, y, z,  c='r', marker='o')
                ax.plot(x,y,z, '--b', alpha = 0.5)
    ax.scatter(x, y, z,  c='r', marker='o', label = 'node')
    ax.plot(x,y,z,'--b', label = 'cable', alpha = 0.5)
    plt.legend(fontsize = 20)

def plot_bar_net(X,B,ax):
    """ 
    function that plots a bar network
    input:
        - X: matrix with nodes and their positions
        - B: adjacency matrix
        - ax: figure where we plot the network 
    output:
        - plot of calbe network
    """
    N = np.shape(X)[1]
    for i in range(N):
        for j in range(i,N):
            if(B[i,j] == 1): 
                x = [X[:,i][0], X[:,j][0]]
                y = [X[:,i][1], X[:,j][1]]
                z = [X[:,i][2], X[:,j][2]]
                ax.plot(x,y,z, c = 'g', linewidth = '2.0', alpha = 0.7)
    ax.plot(x,y,z, c = 'g', linewidth = '2.0', label = 'bar', alpha = 0.7)
    plt.legend()

def E_cable_elast_net(X, A, L, k):
    """
    Function that calculates the total energy of a graph network of cables 
    input:
        - X: matrix of nodes and their positions 
        - A: adjency matrix for the graph network 
        - L: matrix including all rest lengths in the graph network 
        - k: material parameter
    output: 
        - energy: potential energy of the graph network 
    """
    energy = 0 
    N = np.shape(A)[0]
    for i in range(N):
        for j in range(i, N):
            l_ij = L[i,j]
            if((A[i,j] != 0) and np.linalg.norm(X[:,i] - X[:,j]) > l_ij):
                energy += k/(2*l_ij**2) * (np.linalg.norm(X[:,i] - X[:,j]) - l_ij)**2
    return energy

def E_ext(X, mg):
    """
    Function that computes the potential energy of all external loads in the graph network 
    input:
        - X: matrix of nodes and their positions 
        - mg: weight of the external loads, assumed to be homogeneous 
    """
    return np.sum(X[2]) * mg

def derivative_cable(x_i, x_j, l_ij, k):
    """
    function that computes the derivative of the potential energy of a cable between nodes x_i and x_j with regard to x_i
    input:
        - x_i: first node and the node we compute the derivative of 
        - x_j: the second node in the edge
        - l_ij: resting length of the edge 
        - k: material parameter
    output:
        - derivative of potential enregy of a cable with regard to x_i
    """
    norm = np.linalg.norm(x_i - x_j)
    if(norm > l_ij):
        return k/(l_ij**2) * (norm - l_ij) * (x_i - x_j)/norm 
    else:
        return np.zeros_like(x_i -x_j)

def derivative_external_load(x_i,mg_i):
    """
    function that computes the derivative of potential energy of an external load
    input:
        - x_i: external load node 
        - mg_i: weight of external load
    output:
        - numpy array containing the derivative 
    """
    output = np.zeros_like(x_i)
    output[-1] = mg_i
    return output

def E_bar_elast_net(X, B, L, pg, c):
    """
    Function that calculates the total energy of a graph network of bars 
    
    input:
        - X: matrix of nodes and their positions 
        - B: adjency matrix for the graph network 
        - L: matrix including all rest lengths in the graph network 
        - pg: bar density times acceleration of gravity 
        - c: material parameter
    
    output: 
        - energy: potential energy of the graph network 
    """
    energy = 0 
    N = np.shape(X)[1]
    for i in range(N):
        for j in range(i, N):
            if(B[i,j] != 0):
                l_ij = L[i,j]
                energy += c/(2*l_ij**2)*(np.linalg.norm(X[:,i] - X[:,j]) - l_ij)**2 
                energy += pg *l_ij/2 * (X[-1,i] + X[-1,j])
    return  energy

def derivative_bar(x_i, x_j, l_ij, c):
    """
    function that computes the derivative of the potential energy of a bar between nodes x_i and x_j with regard to x_i
    
    input:
        - x_i: first node and the node we compute the derivative of 
        - x_j: the second node in the edge
        - l_ij: resting length of the edge 
        - c: material parameter
    
    output:
        - derivative of potential enregy of a bar with regard to x_i
    """
    norm = np.linalg.norm(x_i - x_j)
    return c/(l_ij**2) * (norm - l_ij) * (x_i - x_j) /norm 

def derivative_bar_grav(x_i, l_ij, pg):
    """
    function that computes the derivative of gravitational potential energy of a bar 
    input:
        - x_i: node in a bar
        _ l_ij, resting length of bar
        - pg: line density of the bar multiplied by gravitational accceleration on the earth's surface
    output:
        - numpy array containing the derivative 
    """
    output = np.zeros_like(x_i)
    output[-1] = l_ij * pg/2
    return output

def gradient_E_cable_net(X, A, L, fixed_points, mg, k):
    """
    Function that computes the gradient of potential energy for a cable network
    input: 
        - X: matrix of nodes and their positions 
        - A: adjadjacency matrix 
        - L: matrix with resting lengths
        - mg: eight of external loads 
        - k: material parameter
    outout:
        - gradient_E: gradient of potential energy for a cabel network
    """
    N = np.shape(X)[1]
    gradient_E = np.zeros_like(X)
    for i in range(N):
        if(i in fixed_points):
            continue
        gradient_E[:,i] += derivative_external_load(X[:,i], mg)
        for j in  range(N):
            if(A[i,j]):
                gradient_E[:,i] += derivative_cable(X[:,i], X[:,j], L[i,j],k)
    return gradient_E

def E(X, A, B, L, fixed_points, mg, pg, k, c):
    """
    Function that computes the energy of a network that has cables and bars 
    input:
        - X: matrix of nodes and their positions 
        - A: adjacency matrix for cables
        - B: adjacency matrix for bars
        - L: resting lengths 
        - mg: weight of external nodes 
        - pg: line density of bar multiplied by earths gravitational acceleration on the surface 
        - k: material paramter of cable
        - c: material paramter of bar
    output:
        - Total potential energy of the graph network
    """
    return E_bar_elast_net(X, B, L, pg, c) + E_cable_elast_net(X, A, L, k) + E_ext(X, mg)


def gradient_E(X, A, B, L, fixed_points, mg, pg, k, c):
    """
    Function that computes the gradient of potential energy for a tensegrity structure
    
    input:
        - X: matrix of nodes and their positions 
        - A: adjacency matrix for cables
        - B: adjacency matrix for bars
        - L: resting lengths 
        - mg: weight of external nodes 
        - pg: line density of bar multiplied by earths gravitational acceleration on the surface 
        - k: material paramter of cable
        - c: material paramter of bar
    
    output:
        - gradient_E: gradient of potential energy
    """
    N = np.shape(X)[1]
    gradient_E = np.zeros_like(X)
    for i in range(N):
        if(i in fixed_points):
            continue
        gradient_E[:,i] += derivative_external_load(X[:,i], mg)
        for j in range(N):
            if(A[i,j] != 0):
                gradient_E[:,i] += derivative_cable(X[:,i], X[:,j], L[i,j], k)
            if(B[i,j] != 0):
                gradient_E[:,i] += derivative_bar(X[:,i], X[:,j], L[i,j], c) 
                gradient_E[:,i] += derivative_bar_grav(X[:,i], L[i,j], pg)

    return gradient_E

def f(X):
    """ 
    constraint function
    input:
        - X: matrix of nodes and their positions 
    output:
        - f(X)
    """
    return (X[0]**2 + X[1]**2)/20.0 

def df(X):
    """
    derivatives of constraint function 
    """
    Y = np.zeros_like(X)
    Y[0] = X[0]/10.0 
    Y[1] = X[1]/10.0
    return Y

def constraint(X):
    """
    function that computes the constraint inequality function on each node in the network 
    input:
        - X: matrix of nodes and their positions 
    output:
        - constraint inequality function 
    """
    return np.maximum(-(X[2] - f(X)), np.zeros_like(X[2]))

def gradient_constraint(X):
    """
    function that computes the gradient of the constraint inequality function 
    input: 
        - X: matrix of nodes and their positions 
    output:
        - gradient of constraint inequality function 
    """
    Y = df(X)
    Y[2] = -np.ones_like(X[2])
    return Y

def quad_penalty(X, A, B, L, fixed_points, mg, pg, k, c, mu):
    """ 
    function that computes the quadratic penalty function
    input:
        - X: matrix of nodes and their positions 
        - A: adjacency matrix for cables
        - B: adjacency matrix for bars
        - L: resting lengths 
        - mg: weight of external nodes 
        - pg: line density of bar multiplied by earths gravitational acceleration on the surface 
        - k: material paramter of cable
        - c: material paramter of bar
    output:
        - quadratic penalty function
    """
    return E(X, A, B, L, fixed_points, mg, pg, k, c) + mu/2 * np.sum(constraint(X)**2)


def gradient_quad_penalty(X, A, B, L, fixed_points, mg, pg, k, c, mu):
    """ 
    computes the gradient of the quadrative penalty function for inequality constraints
    input:
        - X: matrix of nodes and their positions 
        - A: adjacency matrix for cables
        - B: adjacency matrix for bars
        - L: resting lengths 
        - mg: weight of external nodes 
        - pg: line density of bar multiplied by earths gravitational acceleration on the surface 
        - k: material paramter of cable
        - c: material paramter of bar
        - mu: quadrative penaltty parameter
    output:
        - gradient of the quadrative penalty function for inequality constraints
    """
    return gradient_E(X, A, B, L, fixed_points, mg, pg, k, c) + mu * gradient_constraint(X) * constraint(X) 


def armijo_step(f, gradient_f, X, p_k, *args, alpha_0 = 1, c = 1e-1, rho = 1/2):
    """
    function that computes a steplength for armijo condition 
    input:
        - f: function we optimize
        - gradient_f: gradient of f
        - X: matrix of nodes and their positions
        - p_k: search direction
        - *args: arguments used in f and gradient_f
        - alpha_0: first steplength 
        - c: paramter for computing the armijo condition 
        - rh0: parameter we reduce alpha with 
    output:
        alpha: steplength 
        iter: number of iterations
    """
    alpha  = alpha_0
    iter = 0
    D,N = np.shape(X)
    M = D*N
    gradient = gradient_f(X, *args)
    f_call = f(X, *args)
    dot_call = np.dot(np.reshape(p_k, M), np.reshape(gradient, M))
    while(f(X + alpha * p_k, *args) > (f_call + c * alpha * dot_call) and iter < 100):
        alpha = rho * alpha
        iter +=1  
    return alpha, iter

def weak_wolfe(f, gradient_f, X, p_k, *args, alpha_0  =  1, c_1 = 1e-4, c_2 = 0.9):
    """
    function that computes a steplength for the weak Wolfe conditions
    input:
        - f: function we optimize
        - gradient_f: gradient of f
        - X: matrix of nodes and their positions
        - p_k: search direction
        - *args: arguments used in f and gradient_f
        - alpha_0: first steplength 
        - c_1,c_2: paramter for computing the armijo condition 
    output:
        alpha: steplength 
        iter: number of iterations
    """
    alpha_max = np.inf
    alpha_min = 0
    alpha = alpha_0
    iter = 0 
    D,N = np.shape(X)
    M = D*N
    f_call = f(X, *args) 
    gradient = gradient_f(X, *args)
    dot_call = np.dot(np.reshape(gradient, M), np.reshape(p_k, M))
    armijo_cond = False
    wolfe_cond = False
    while(iter < 100):
        iter +=1   
        if(f(X + alpha * p_k, *args) > (f_call + c_1 * alpha * dot_call)):
            alpha_max = alpha
            alpha = (alpha_max + alpha_min)/2
            armijo_cond = False 
        else:
            armijo_cond = True
        
        if(-np.dot(np.reshape(p_k, M), np.reshape(gradient_f(X + alpha * p_k, *args),M)) > -(c_2 * dot_call)):
            alpha_min = alpha    
            if(alpha_max == np.inf):
                alpha = 2*alpha
            else:
                alpha = (alpha_max + alpha_min)/2
            wolfe_cond = False
        else:
            wolfe_cond = True 
        
        if(armijo_cond and wolfe_cond):
            return alpha, iter
    return alpha, iter

def strong_wolfe(f, gradient_f, X, p_k, *args, alpha_0  =  2, c_1 = 1e-3, c_2 = 0.9):
    alpha_max = alpha_0
    alpha_min = 0.0
    iter = 0 
    D,N = np.shape(X)
    M = D*N
    f_call = f(X, *args)
    gradient = gradient_f(X, *args)
    dot_call = np.inner(np.reshape(gradient, M), np.reshape(p_k, M))
    Armijo = (f(X + alpha_max * p_k, *args) <= (f_call + c_1 * alpha_max * dot_call))
    curvature_low = (np.dot(np.reshape(p_k, M), np.reshape(gradient_f(X + alpha_max * p_k, *args),M)) >= (c_2 * dot_call))
    curvature_high = (np.dot(np.reshape(p_k, M), np.reshape(gradient_f(X + alpha_max * p_k, *args),M)) <= -(c_2 * dot_call))
    iter = 0 
    tot_iter = 0
    while(iter < 100 and (Armijo and (not curvature_low))):
        iter +=1 
        alpha_min = alpha_max
        alpha_max *=2 
        Armijo = (f(X + alpha_max * p_k, *args) <= (f_call + c_1 * alpha_max * dot_call))
        curvature_low = (np.dot(np.reshape(p_k, M), np.reshape(gradient_f(X + alpha_max * p_k, *args),M)) >= (c_2 * dot_call))
        curvature_high = (np.dot(np.reshape(p_k, M), np.reshape(gradient_f(X + alpha_max * p_k, *args),M)) <= -(c_2 * dot_call))
    if(Armijo and curvature_high and curvature_low):
        return alpha_max, iter
    tot_iter += iter
    iter = 0 
    alpha = np.copy(alpha_max)
    while(iter < 100 and not(Armijo and curvature_high and curvature_low)):
        iter +=1 
        if(Armijo and not(curvature_low)):
            alpha_min = alpha 
        else:
            alpha_max = alpha 
        alpha = (alpha_max + alpha_min)/2

        Armijo = (f(X + alpha * p_k, *args) <= (f_call + c_1 * alpha * dot_call))
        curvature_low = (np.dot(np.reshape(p_k, M), np.reshape(gradient_f(X + alpha * p_k, *args),M)) >= (c_2 * dot_call))
        curvature_high= (np.dot(np.reshape(p_k, M), np.reshape(gradient_f(X + alpha * p_k, *args),M)) <= -(c_2 * dot_call))
    tot_iter += iter
    return alpha, tot_iter


def gradient_descent(f, gradient_f, X, *args, tol = 1e-8, alpha_0 = 1, armijo = False, weak_w = False, strong_w = False):
    """
    function that does gradient descent with either armijo, weak wolfe or strong wolfe conditions or constant step size 
    input:
        - E: function we optimize
        - gradient_E: gradient of E
        - X: matrix of nodes and their positions
        - *args: arguments used in E and gradient_E
        - alpha_0: first steplength 
        - armijo: bool operator that tells us if we are doing armijo search 
        - weak_w: bool operator that tells us if we are doing weak wolfe search 
        - strong_w: bool operator that tells us if we are doing strong wolfe search 
    output:
        X: optimized position of nodes 
    """
    steps = 0
    Q = 0 
    d, N = np.shape(X)
    M = d*N
    error = []
    while(np.linalg.norm(np.reshape(gradient_f(X, *args), M)) > tol and steps < 1000):
        steps +=1 
        p_k = -gradient_f(X, *args)
        if(armijo):
            p_k = p_k
            alpha, q = armijo_step(f, gradient_f, X, p_k, *args)
            Q += q
        elif(weak_w):
            p_k = p_k
            alpha,q = weak_wolfe(f, gradient_f, X, p_k, *args)
            Q += q
        elif(strong_w):
            p_k = p_k
            alpha, q  = strong_wolfe(f, gradient_f, X, p_k, *args) 
            Q += q
        else:
            alpha = 1
        error.append(np.linalg.norm(np.reshape(gradient_f(X, *args),M)))
        X += alpha * p_k

    print(f'number of gradient descent steps = {steps}, number of step size optimization = {Q}, norm of gradient = {np.linalg.norm(np.reshape(gradient_f(X, *args),M))}, Energy = {f(X, *args)}')
    return X, error


def BFGS(f, gradient_f, X, *args, tol = 1e-10,  armijo = False, weak_w = False, strong_w = False):
    """
    function that does quasi Newton (Broyden-Fletcher-Goldfarb-Shanno) algorithm with either armijo, weak wolfe or strong wolfe conditions or constant step size 
    
    input:
        - E: function we optimize
        - gradient_E: gradient of E
        - X: matrix of nodes and their positions
        - *args: arguments used in E and gradient_E
        - armijo: bool operator that tells us if we are doing armijo search 
        - weak_w: bool operator that tells us if we are doing weak wolfe search 
        - strong_w: bool operator that tells us if we are doing strong wolfe search 
    
    output:
        X: optimized position of nodes 
    """
    steps = 0
    Q = 0 
    d,N = np.shape(X)
    M = d*N
    H = np.identity(M)
    error = []
    while(np.linalg.norm(np.reshape(gradient_f(X, *args),M)) > tol and steps < 1000):
        p_k = - H @ np.reshape(gradient_f(X, *args), M)
        p_k = np.reshape(p_k, (d,N))
        if(armijo):
            step_size,q = armijo_step(f, gradient_f, X, p_k, *args)
        elif(weak_w):
            step_size,q = weak_wolfe(f, gradient_f, X, p_k, *args)
        elif(strong_w):
            step_size,q = strong_wolfe(f, gradient_f, X, p_k, *args)

        else:
            step_size,q = (1,0)
        s = step_size * p_k
        y = gradient_f(X + s, *args) - gradient_f(X, *args)
        y = np.reshape(y, M)
        s = np.reshape(s, M)
        if(np.dot(y,s) == 0):
            print(f'stepsize = {step_size}, ||p_k|| = {np.linalg.norm(np.reshape(p_k, M))}')
            print(f'||y|| = {np.linalg.norm(y)}, ||s|| = {np.linalg.norm(s)}')
            print(f'y.T s  = {np.dot(y,s)}')
            break
        alpha = 1/np.inner(s, y)
        z =  H @ y
        H += alpha * (1 + np.inner(y, z) * alpha) * np.outer(s, s)  - alpha * (np.outer(z, s) + np.outer(s, z))
        X += step_size*p_k
        steps +=1 
        Q +=q
        error.append(np.linalg.norm(np.reshape(gradient_f(X, *args),M)))
    print(f'number of BFGS steps = {steps}, number of step size optimization = {Q}, norm of gradient = {np.linalg.norm(np.reshape(gradient_f(X, *args),M))}, Energy = {f(X, *args)}')
    return X, error 