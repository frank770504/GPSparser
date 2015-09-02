import sys
import time
import math

import datetime
import gpxpy
import gpxpy.gpx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class gps_point:
	def __init__(self, lat, lon, ele, time, sp):
		self.lat = lat
		self.lon = lon
		self.ele = ele
		self.time = time
		self.spd = sp

def D2R(d):
	r = d*math.pi/180
	return r

def R2D(r):
	d = r*180/math.pi
	return d

def get_gps_list(name):
	gpxf = open(name, 'r')
	gpx = gpxpy.parse(gpxf)
	gps_list = []
	for track in gpx.tracks:
		for segment in track.segments:
			for point in segment.points:
				d = point.time
				p = gps_point(D2R(point.latitude), \
					D2R(point.longitude), \
					point.elevation, \
					time.mktime(d.timetuple()) + \
						1e-6 * d.microsecond, point.speed)
				gps_list.append(p)
	return gps_list

def get_pos_nplist(gps_list):
	pos_nplist = np.array([])
	for gps in gps_list:
		temp = np.array([gps.lon, gps.lat, gps.ele])
		pos_nplist = np.vstack([pos_nplist, temp]) if pos_nplist.size else temp
	return pos_nplist

def cal_dist(pos1, pos2):
	EarthRadius = 6371000
	# pos = [lon, lat, ele]
	p1 = (pos2[0] - pos1[0]) * math.cos( 0.5 * (pos2[1] + pos1[1]) )
	p2 = (pos2[1] - pos1[1])
	if np.size(pos1) == 3:
		distance = (EarthRadius + 0.5 * math.fabs( pos2[2] - pos1[2]) ) *\
			math.sqrt( p1**2 + p2**2 )
	elif np.size(pos1) == 2:
		distance = (EarthRadius) * math.sqrt( p1**2 + p2**2 )
	return distance

def get_gps_dist(gps1, gps2):
	pos1 = [gps1.lon, gps1.lat, gps1.ele]
	pos2 = [gps2.lon, gps2.lat, gps2.ele]
	dist = cal_dist(pos1, pos2)
	return dist

def get_nplist_dist(pos_nplist):
	sz = np.size(np.shape(pos_nplist[0]))
	if sz == 1:
		pos1 = pos_nplist[0].tolist()
	elif sz == 2:
		pos1 = np.array(pos_nplist[0]).reshape(-1,).tolist()
	dist = 0
	for pos in pos_nplist:
		if sz == 1:
			pos2 = pos.tolist()
		elif sz == 2:
			pos2 = np.array(pos).reshape(-1,).tolist()
		dist = dist + cal_dist(pos1, pos2)
		pos1 = pos2
	return dist

def get_boundary(pos_nplist):
	xmax = np.amax(pos_nplist[:,0])
	xmin = np.amin(pos_nplist[:,0])
	gap = (xmax - xmin)/10
	xmax = xmax + gap
	xmin = xmin - gap
	ymax = np.amax(pos_nplist[:,1])
	ymin = np.amin(pos_nplist[:,1])
	gap = (ymax - ymin)/10
	ymax = ymax + gap
	ymin = ymin - gap
	return (xmax, xmin, ymax, ymin)

def reform_nplist(nplist):
	row, col = np.shape(nplist)
	mnnplist = np.array([])
	for c in range(0,col,1):
		nptemp = nplist[:,c]
		sz = np.size(nptemp)
		mn = np.ones((1,sz))*np.mean(nptemp)
		mnnplist = np.vstack([mnnplist, mn]) if mnnplist.size > 0 else mn
		nplist[:,c] = (nptemp - mn)*10**7
	mnnplist = mnnplist.T
	return nplist, mnnplist

def turnback_nplist(nplist, mn):
	nplist = nplist[1:, :]
	back = nplist/(10**7) + mn
	return back

import Kalman as kf

def moving_average(a, n=3) :
	ret = np.cumsum(a, dtype=float)
	ret[n:] = ret[n:] - ret[:-n]
	return ret[n - 1:] / n

f_names = sys.argv
f_names = f_names[1:]

for name in f_names:
	gps_list = get_gps_list(name)
	pos_nplist = get_pos_nplist(gps_list)
	print "distance: {}".format(get_nplist_dist(pos_nplist))

pos_nplist, mn = reform_nplist(pos_nplist[:,0:2])

# TODO split to big scale and merge

sz = np.size(pos_nplist, axis=0)

vx_list = pos_nplist[1:,0]-pos_nplist[0:sz-1,0]
vy_list = pos_nplist[1:,1]-pos_nplist[0:sz-1,1]

vx = np.mean(np.abs(vx_list))
vy = np.mean(np.abs(vy_list))

covx = np.cov(np.abs(vx_list))
covy = np.cov(np.abs(vy_list))

Sig = np.diagflat([covx, covy, covx, covy])

covavg = np.mean([covx, covy])

mu_nplist = np.hstack([pos_nplist[0,:], 0, 0])
pre_nplist = np.hstack([pos_nplist[0,:], 0, 0])
Sig_nplist = np.diag(Sig)

vx_list = np.insert(vx_list, 0, 0)
vy_list = np.insert(vy_list, 0, 0)

V = np.vstack([vx_list, vy_list]).T
pos_nplist = np.hstack([pos_nplist, V])

A = np.matrix('1 0 1 0;\
               0 1 0 1;\
               0 0 1 0;\
               0 0 0 1')

C = np.eye(4)

B = np.matrix('0 0;\
               0 0;\
               1 0;\
               0 1')

kalman = kf.Kalman_filter(A, B, C, 0.01, covx, covy)

for z in pos_nplist:
	u = np.zeros(2)
	mu = mu_nplist[-1,:] if mu_nplist.size>4 else mu_nplist
	Sig = Sig_nplist[-1,:] if mu_nplist.size>4 else Sig_nplist
	Sig = np.diagflat(Sig)
	mu_nxt, Sig_nxt, pre = kalman.do_filter(mu, Sig, z, u)
	mu_nplist = np.vstack([mu_nplist, mu_nxt])
	pre_nplist = np.vstack([pre_nplist, pre])
	Sig_nplist = np.vstack([Sig_nplist, np.diag(Sig_nxt)])

mu_back = turnback_nplist(mu_nplist[:,0:2], mn)
print "KF distance: {}".format(get_nplist_dist(mu_back))

plt.figure()
plt.plot(pos_nplist[:,0], pos_nplist[:,1], '-*', label='raw')
plt.plot(moving_average(pos_nplist[:,0], 20), moving_average(pos_nplist[:,1], 20), '-*', label='smth')
plt.plot(mu_nplist[:,0], mu_nplist[:,1], '-+', label='post')
plt.title('GPS Path')
plt.xlabel('lat')
plt.ylabel('lon')
plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)

# First set up the figure, the axis, and the plot element we want to animate
xmax, xmin, ymax, ymin = get_boundary(pos_nplist)

fig = plt.figure()
ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
line, = ax.plot([], [], '.')

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
	x = pos_nplist[0:i,0]
	y = pos_nplist[0:i,1]
	#x = np.linspace(0, 2, 1000)
	#y = np.sin(2 * np.pi * (x - 0.01 * i))
	line.set_data(x, y)
	return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,\
                               frames=sz, interval=20, blit=True)

# Set up formatting for the movie files
#anim.save('basic_animation.mp4', writer='mencoder', fps=30)

plt.show()
