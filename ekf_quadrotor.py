import numpy as np
import matplotlib.pyplot as plt

class quadrotor:

    def __init__(self, m, l, g, Ixx, Iyy, Izz):
        self.m = m
        self.l = l
        self.g = g
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Izz = Izz
        self.a1 = (Iyy - Izz) / Ixx
        self.a2 = (Izz - Ixx) / Iyy
        self.a3 = (Ixx - Iyy) / Izz
        self.b1 = l / Ixx
        self.b2 = l / Iyy
        self.b3 = l / Izz
        self._x = np.zeros(12)

    def step(self, u, dt):
        ux = np.cos(self._x[0]) * np.sin(self._x[2]) * np.cos(self._x[4]) + np.sin(self._x[0]) * np.sin(self._x[4])
        uy = np.cos(self._x[0]) * np.sin(self._x[2]) * np.sin(self._x[4]) - np.sin(self._x[0]) * np.cos(self._x[4])
        self._x[0] += dt * self._x[1]
        self._x[1] += dt * (self._x[3] * self._x[5] * self.a1 + self.b1 * u[1])
        self._x[2] += dt * self._x[3]
        self._x[3] += dt * (self._x[1] * self._x[5] * self.a2 + self.b2 * u[2])
        self._x[4] += dt * self._x[5]
        self._x[5] += dt * (self._x[3] * self._x[5] * self.a3 + self.b3 * u[3])
        self._x[6] += dt * self._x[7]
        self._x[7] += dt * ((np.cos(self._x[0]) * np.cos(self._x[2]) * u[0]) / self.m - self.g)
        self._x[8] += dt * self._x[9]
        self._x[9] += dt * (ux / self.m * u[0])
        self._x[10] += dt * self._x[11]
        self._x[11] += dt * (uy / self.m * u[0])  

        #new_x = np.array(self.x)
        #return new_x
    
    @property
    def x(self):
        return self._x
    
    @x.setter
    def x(self, init_state):
        self._x = init_state



def ekf_quadrotor(quad = None, x = None,u = None,dt = None,y = None,P = None,Q = None,R = None): 
  
    # for nonlinear dynamic system:
    #           x_k+1 = f(x_k) + w_k
    #           y_k   = h(x_k) + v_k
    # where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
    #       v ~ N(0,R) meaning v is gaussian noise with covariance R
    # Inputs:   f: function handle for f(x)
    #           x: "a priori" state estimate
    #           P: "a priori" estimated state covariance
    #           h: fanction handle for h(x)
    #           y: current measurement
    #           Q: process noise covariance
    #           R: measurement noise covariance
    # Output:   x: "a posteriori" state estimate
    #           P: "a posteriori" state covariance

    x1,A = jaccsd(x,quad,u,dt)
    P = np.multiply(A, P, np.transpose(A)) + Q
    y1,H = jaccsd(x1,quad,u,dt)
    P12 = np.multiply(P, np.transpose(H))
    K = np.multiply(P12, np.linalg.inv(np.multiply(H, P12) + R))
    x_est = x1 + (np.dot(K, np.transpose(y - y1))).ravel()

    P -= np.multiply(K, np.transpose(P12))
    
    # R=chol(H*P12+R);            #Sk, Cholesky factorization
    # U=P12/R;                    #K=U/R'; Faster because of back substitution
    # x=x1+U*(R'\(y-y1));         #Back substitution to get state update
    # K = U*inv(R');
    # P=P-U*U';                   #Covariance update, U*U'=P12/R/R'*P12'=K*P12.
    return x_est, P
    
def jaccsd(x0, quad = None,u = None,dt = None): 
    # JACCSD Jacobian through complex step differentiation
    # [z J] = jaccsd(f,x)
    # z = f(x)
    # J = f'(x)
    n = np.asarray(x0).size
    m = np.asarray(x0).size
    A = np.zeros((m,n))
    slack = 0.0001
    init_st = np.array(quad.x)
    quad.x = np.array(x0)
    quad.step(u,dt)
    z = np.array(quad.x)

    for k in np.arange(0,n).reshape(-1):
        x1 = np.zeros(n)
        x1[k] = slack
        quad.x = x0 + x1
        quad.step(u,dt)
        z1 = np.array(quad.x)
        A[:,k] = np.divide((z1 - z), slack)
    quad.x = init_st
    return z, A
    

def main():

    x0 = np.array([np.pi/8, 0.8, np.pi/8, 0.3, np.pi/8, 0.1, 0.4, 0.2, 0.4, 0.3, 0.2, 0.4])
    u0 = np.array([100, 1, 1, 0])
    dt = 0.01
    P = np.multiply(x0,x0) + 0.1
    Q = np.zeros(12) + 1e-10
    R = np.ones(12)

    quad = quadrotor(0.65, 0.23, 9.81, 7.5e-3, 7.5e-3, 1.3e-2)
    quad.x = np.array([np.pi/8, 0.8, np.pi/8, 0.3, np.pi/8, 0.1, 0.4, 0.2, 0.4, 0.3, 0.2, 0.4])
    x_est_save = np.array(x0)
    x_cur_save = np.array(x0)
    t = [0]

    for i in range(1,1000):
        t.append(i*dt)
        states_cur = quad.x
        y = states_cur + np.random.randn(1,12)
        states_est, P = ekf_quadrotor(quad, np.array(states_cur), u0, dt, y, P, Q, R)
        x_est_save = np.append(x_est_save, states_est[:])
        x_cur_save = np.append(x_cur_save, states_cur[:])

        quad.step(u0, dt)
    x_cur_save = x_cur_save.reshape(1000,12)
    x_est_save = x_est_save.reshape(1000,12)

    plt.subplot(3, 1, 1)
    plt.plot(t, x_cur_save[:,8], t, x_est_save[:,8])
    plt.ylabel("Position X")
    plt.legend()
    plt.grid()

    plt.subplot(3, 1, 2)
    plt.plot(t, x_cur_save[:,10], t, x_est_save[:,10])
    plt.ylabel("Position y")
    plt.legend()
    plt.grid()

    plt.subplot(3, 1, 3)
    plt.plot(t, x_cur_save[:,6], t, x_est_save[:,6])
    plt.xlabel("Time [s]")
    plt.ylabel("Position z")
    plt.legend()
    plt.grid()

    plt.show()

main()
