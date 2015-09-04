import numpy as np

class gaussian_noise_diag:
	def __init__(self, mu, cov, sz):
		self.mu = mu
		self.cov = cov
		self.sz = sz
	def get_mat(self):
		return np.diagflat(np.random.normal(self.mu, self.cov, self.sz))

class gaussian_noise_vec:
	def __init__(self, mu, cov, shape):
		self.mu = mu
		self.cov = cov
		self.shape = shape
	def get_mat(self):
		r, c = self.shape
		sz = r*c
		return np.random.normal(self.mu, self.cov, sz).reshape(self.shape)

class Kalman_filter:
	def __init__(self, A, B, C, err_cov=0.01, R_cov=1, Q_cov=1):
		self.A = A
		self.B = B
		self.C = C
		self.err_cov = err_cov
		self.R_cov = R_cov
		self.Q_cov = Q_cov
		self.dump = 'off'
	def set_A(self, A):
		self.A = A
	def set_B(self, B):
		self.B = B
	def set_C(slef, C):
		self.C = C
	def set_cov(self, err_cov, R_cov, Q_cov):
		if err_cov != 'default':
			self.err_cov = err_cov
		if R_cov != 'default':
			self.R_cov = R_cov
		if Q_cov != 'default':
			self.Q_cov = Q_cov
	def dump_matrices(self,dump):
		self.dump = dump
	def do_filter(self, mu, Sig, z, u, EKF_mu_callback='default', EKF_z_callback='default'):
		from numpy.linalg import inv
		mu = np.atleast_2d(mu).T
		z = np.atleast_2d(z).T
		u = np.atleast_2d(u).T
		rA, cA = np.shape(self.A)
		R = gaussian_noise_diag(0, self.R_cov, rA)
		rC, cC = np.shape(self.C)
		Q = gaussian_noise_diag(0, self.Q_cov, rC)
		err = gaussian_noise_vec(0, self.err_cov, mu.shape)
		if EKF_mu_callback=='default':
			mu_pre = self.A*mu + self.B*u + err.get_mat()
		else:
			mu_pre = EKF_mu_callback(mu) + err.get_mat()
		Sig_pre = self.A*Sig*self.A.T + np.fabs(R.get_mat())
		K = Sig_pre*self.C.T*inv(self.C*Sig_pre*self.C.T + np.fabs(Q.get_mat()))
		if EKF_z_callback=='default':
			pre_z = self.C*mu_pre
		else:
			pre_z = EKF_z_callback(mu_pre)
		mu_post = mu_pre + K*(z - pre_z)
		Sig_post = (np.eye(cC) - K*self.C)*Sig_pre
		if self.dump == 'on':
			print '{} mu'.format(mu.T)
			print '{} Sig'.format(np.diag(Sig))
		#	print u.T
			print '{} z'.format(z.T)
			print '{} mu_pre_t'.format(mu_pre.T)
			print '{} Sig_pre'.format(np.diag(Sig_pre))
		#	print '{} C*Sig_pre*C.T + Q'.format(inv(C*Sig_pre*C.T + Q))
		#	print '{} Sig_pre*C'.format(Sig_pre*C.T)
			print '{} K'.format(K)
		#	print '{} R'.format(np.diag(R))
		#	print '{} Q'.format(np.diag(Q))
		#	print '{} K*(z - C*mu_pre)'.format((K*(z - C*mu_pre)).T)
			print '{} mu_post'.format(mu_post.T)
			print ''
		return (mu_post.T, Sig_post, mu_pre.T)
