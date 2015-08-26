import sys
import time
import datetime
import gpxpy
import gpxpy.gpx
import math
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
	distance = (EarthRadius + 0.5 * math.fabs( pos2[2] - pos1[2]) ) *\
		math.sqrt( p1**2 + p2**2 )
	#distance = (EarthRadius) * math.sqrt( p1**2 + p2**2 )
	return distance

def get_gps_dist(gps1, gps2):
	pos1 = [gps1.lon, gps1.lat, gps1.ele]
	pos2 = [gps2.lon, gps2.lat, gps2.ele]
	dist = cal_dist(pos1, pos2)
	return dist

def get_nplist_dist(pos_nplist):
	pos1 = pos_nplist[0].tolist()
	dist = 0
	for pos in pos_nplist:
		pos2 = pos.tolist()
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

f_names = sys.argv
f_names = f_names[1:]

for name in f_names:
	gps_list = get_gps_list(name)
	pos_nplist = get_pos_nplist(gps_list)
	print "distance: {}".format(get_nplist_dist(pos_nplist))

plt.figure()
plt.plot(pos_nplist[:,0], pos_nplist[:,1], '.', label='path')
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
                               frames=pos_nplist.size, interval=20, blit=True)

# Set up formatting for the movie files
#anim.save('basic_animation.mp4', writer='mencoder', fps=30)

plt.show()
