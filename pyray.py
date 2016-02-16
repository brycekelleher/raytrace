class vec3:
	def __init__(self, x, y, z):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
	def __getitem__(self, index):
		return (self.x, self.y, self.z)[index]
	def __str__(self):
		return "(%f, %f, %f)" % (self.x, self.y, self.z)
	def __add__(self, other):
		return vec3(self.x + other.x, self.y + other.y, self.z + other.z)
	def __sub__(self, other):
		return vec3(self.x - other.x, self.y - other.y, self.z - other.z)
	def __rmul__(self, other):
		return vec3(other * self.x, other * self.y, other * self.z)

def add3(a, b):
	return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
def sub3(a, b):
	return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def dot3(a, b):
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])
def cross3(a, b):
	return ((a[1] * b[2]) - (a[2] * b[1]), (a[2] * b[0]) - (a[0] * b[2]), (a[0] * b[1]) - (a[1] * b[0]));
def scalartriple(a, b, c):
	return dot3(a, cross3(b, c));

# The basic concept here is that the sign of the triple product of three vectors will change 
# depending on their handedness. Using this fact we can test which side the 'd' direction vector is
# compared to vectors from the origin to the triangle vertices. td measures whether the direction
# vector sits in the plane of the triangle. tt measures the distance from the origin to the triangle
def triangle_intersect(o, d, a, b, c):
	oa = a - o
	ob = b - o
	oc = c - o

	n = vec3(0, 0, 1)
	td = dot3(n, d)
	tt = (1.0 / td) * dot3(n, oa)

	u = scalartriple(d, oc, ob)
	v = scalartriple(d, oa, oc)
	w = scalartriple(d, ob, oa)

	denom = 1.0 / (u + v + w)
	u *= denom
	v *= denom
	w *= denom

	return (u, v, w, tt)

va = vec3( 0,  0, 0)
vb = vec3(10,  0, 0)
vc = vec3( 0, 10, 0)
origin = vec3(2.5, 2.5, 0)
direction = vec3(0, 0, 1)
imagew = 256
imageh = 256

#test_intersection(origin, direction, va, vb, vc)
#def test_intersection(o, d, a, b, c):
#	u, v, w = triangle_intersect(o, d, a, b, c)
#	r = (u * a) + (v * b) + (w * c)
#	print "u", u, "v", v, "w", w, "r", r

def intersected(u, v, w):
	if u < 0.0 or u > 1.0:
		return False
	if v < 0.0 or v > 1.0:
		return False
	if w < 0.0 or w > 1.0:
		return False
	return True

def shade(pixel):
	# calculate a view ray
	#origin = vec3(pixel[0] / float(imagew), pixel[1] / float(imageh), 0.0)
	origin = vec3(10 * (pixel[0] / imagew), 10 * (pixel[1] / imageh), 10.0)
	direction = vec3(0.0, 0.0, -1.0) 

	#intersect the ray
	u, v, w, t = triangle_intersect(origin, direction, va, vb, vc)
	if not intersected(u, v, w):
		return (0.0, 0.0, 0.0)

	# caclulate intersection point and return color
	r = vec3(1, 0, 0)
	g = vec3(0, 1, 0)
	b = vec3(0, 0, 1)

	color = (u * r) + (v * g) + (w * b)
	#print "origin", origin, "direction", direction, "u", u, "v", v, "w", w, "t", t, "color", color
	return color

from itertools import *

def enumerate_pixels(callback):
	cols = range(imagew)
	rows = range(imageh)
	image = {}

	for y, x in product(rows, cols):
		 pixel = (float(x), float(y))
		 color = callback(pixel)
		 image[pixel] = color
	return image

def write_ppm(filename, image):
	fp = open(filename, "w")
	fp.write("P3\n%i %i\n255\n" % (imagew, imageh))

	cols = range(imagew)
	rows = range(imageh)
	for y in rows:
		for x in cols:
			# convert from image space into sample space here (invert y)
			color = image[(x, imageh - 1 - y)]

			# 'tonemap' the colors
			r = int(255 * color[0])
			g = int(255 * color[1])
			b = int(255 * color[2])

			# write the pixel values to the ppm file
			fp.write("%i " % r)
			fp.write("%i " % g)
			fp.write("%i " % b)
			fp.write("\n")

image = enumerate_pixels(shade)

write_ppm("raytrace.ppm", image)

