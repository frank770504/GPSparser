import sys
import time
import datetime
import gpxpy
import gpxpy.gpx
import math
import numpy as np
import matplotlib

#matplotlib.use('Agg')
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

def get_gpx(name):
	gpxf = open(name, 'r')
	gpx = gpxpy.parse(gpxf)
	gps_list = []
	for track in gpx.tracks:
		for segment in track.segments:
			for point in segment.points:
				d = point.time
				p = gps_point(D2R(point.latitude), D2R(point.longitude), point.elevation, time.mktime(d.timetuple()) + 1e-6 * d.microsecond, point.speed)
				gps_list.append(p)
	return gps_list


def get_dist(gps1, gps2):
	EarthRadius = 6371000
	lonp = gps2.lon
	lon1 = gps1.lon
	latp = gps2.lat
	lat1 = gps1.lat
	p1 = (lonp - lon1) * math.cos( 0.5 * (latp + lat1) )
	p2 = (latp - lat1)
	distance = (EarthRadius + 0.5 * math.fabs( gps2.ele - gps1.ele) ) * math.sqrt( p1*p1 + p2*p2 )
	#distance = (EarthRadius) * math.sqrt( p1*p1 + p2*p2 )
	return distance

f_names = sys.argv
f_num = len(f_names)
f_names = f_names[1:f_num]


for name in f_names:
	gps_list = get_gpx(name)
	dis_list = []
	pos_list = np.array([])
	zer_list = np.array([])
	gps_0 = gps_list[0]
	dis = 0
	for gps in gps_list[0:]:
		#print "{} {} {} {}".format(gps.lat, gps.lon, gps.ele, gps.time)
		dis_list.append(get_dist(gps_0, gps))
		dis = dis + get_dist(gps_0, gps)
		gps_0 = gps
		if gps.spd!=0:
			#x1 = R2D(gps.lat)
			#x2 = R2D(gps.lon)
			x1 = gps.lat
			x2 = gps.lon
			temp = np.array([x2, x1])
			pos_list = np.vstack([pos_list, temp]) if pos_list.size else temp
		else:
			#x1 = R2D(gps.lat)
			#x2 = R2D(gps.lon)
			x1 = gps.lat
			x2 = gps.lon
			temp = np.array([x2, x1])
			zer_list = np.vstack([zer_list, temp]) if zer_list.size else temp
	print ""
	print "distance: {}".format(dis)

#for d in dis_list:
#	print d

#for p in pos_list:
#	print p

plt.figure(figsize=(10,5))
plt.plot(pos_list[:,0], pos_list[:,1], '.', label='path')
plt.title('GPS Path')
plt.xlabel('lat')
plt.ylabel('lon')
plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
plt.figure(figsize=(10,5))
plt.plot(zer_list[:,0], zer_list[:,1], '.', label='still')
plt.title('GPS Path still')
plt.xlabel('lat')
plt.ylabel('lon')
plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)

# First set up the figure, the axis, and the plot element we want to animate
xmax = np.amax(pos_list[:,0])
xmin = np.amin(pos_list[:,0])
gap = (xmax - xmin)/10
xmax = xmax + gap
xmin = xmin - gap
ymax = np.amax(pos_list[:,1])
ymin = np.amin(pos_list[:,1])
gap = (ymax - ymin)/10
ymax = ymax + gap
ymin = ymin - gap

fig = plt.figure()
ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
line, = ax.plot([], [], '.')

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,
# animation function.  This is called sequentially
def animate(i):
	x = pos_list[0:i,0]
	y = pos_list[0:i,1]
	#x = np.linspace(0, 2, 1000)
	#y = np.sin(2 * np.pi * (x - 0.01 * i))
	line.set_data(x, y)
	return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=pos_list.size, interval=20, blit=True)

# Set up formatting for the movie files
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
#anim.save('basic_animation.mp4', writer=writer)

plt.show()
