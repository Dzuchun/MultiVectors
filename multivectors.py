import array
from math import cos, sin

"""
This module defines multivectors (MV).
More technical info on those: https://en.wikipedia.org/wiki/Geometric_algebra.
Initially inspired by: https://youtu.be/60z_hpEAtD8.
"""


class MultiVector:
    """
    A MV. Consists of scalar, vector, pseudovector (PV) and pseudoscalar (PS).
    """
    
    def __init__(self, real=0.0, vecx=0.0, vecy=0.0, vecz=0.0, axix=0.0, axiy=0.0, axiz=0.0, imag=0.0):
        """
        Initialises new multivector. All components have default value of zero.
        :param real: scalar component.
        :param vecx: vector x component.
        :param vecy: vector y component.
        :param vecz: vector z component.
        :param axix: PV x component.
        :param axiy: PV y component.
        :param axiz: PV z component.
        :param imag: PS component.
        """
        self.comp = array.array('d', [real, vecx, vecy, vecz, axix, axiy, axiz, imag])
    
    def __getitem__(self, item):
        return self.comp[item]
    
    def __setitem__(self, key, value):
        if isinstance(key, int):
            self.comp[key] = value
    
    def __neg__(self):
        """
        :return: the opposite MV, has all components negated.
        """
        return MultiVector(*[-c for c in self])
    
    def __add__(self, other):
        """
        :param other: MV to add.
        :return: MV that is sum of self and other.
        """
        return MultiVector(*[c1 + c2 for c1, c2 in zip(self, other)])
    
    def __sub__(self, other):
        """
        :param other: MV to subtract.
        :return: MV that is self - other.
        """
        return MultiVector(*[c1 - c2 for c1, c2 in zip(self, other)])
    
    def __mul__(self, other):
        """
        Geometric product of two MVs or a scalar multiplication. Geometric product is a quite complicated operation, so
        make sure that you understood underlying math before using it.
        :param other: MV or a scalar to multiply by.
        :return: MV that is a result of operation.
        """
        if isinstance(other, MultiVector):
            # multiplication by another MV
            return MultiVector(*[
                sum([
                    sum([
                        self[j] * other[k] * MultiVector.METRIC[j][k][i] for k in range(8)
                    ]) for j in range(8)
                ]) for i in range(8)
            ])
        elif isinstance(other, int) or isinstance(other, float):
            # multiplication by a scalar
            return MultiVector(*[c * other for c in self])
        else:
            return None
    
    def __rmul__(self, other):
        if isinstance(other, MultiVector):
            return other.__mul__(self)
        elif isinstance(other, int) or isinstance(other, float):
            return self.__mul__(other)
    
    def __str__(self):
        return fr'''
MultiVector:
[real={self[0]},
vec=({self[1]}; {self[2]}; {self[3]}),
axi=({self[4]}; {self[5]}; {self[6]}),
imag={self[7]}
]'''
    
    @classmethod
    def linear_combination(cls, system, comps):
        """
        Can be used to express new MV in terms of predefined MV system. Make sure you pass as much components, as there are MVs in the system.
        :param system: MV system.
        :param comps: new MV components.
        :return: new that is a linear combination of system MVs.
        """
        if len(system) == len(comps):
            return sum([c * v for c, v in zip(comps, system)], MultiVector.ZERO)
        else:
            return None
        
    @classmethod
    def get_vector(cls, x=0.0, y=0.0, z=0.0):
        """
        Used to create MV that only represents normal vector. Can be thought of as
        MultiVector.linear_combination(MultiVector.VECTOR_BASIS, (x, y, z)).
        Note, that this is not an actual implementation.
        :param x: vector x component.
        :param y: vector y component.
        :param z: vector z component.
        :return: MV with only vector components set.
        """
        return MultiVector(vecx=x, vecy=y, vecz=z)
    
    @classmethod
    def get_axis(cls, x=0.0, y=0.0, z=0.0):
        """
        Used to create MV that only represents PV (axis/axial vector). Can be thought of as
        MultiVector.linear_combination(MultiVector.AXIS_BASIS, (x, y, z))
        Note, that this is not an actual implementation.
        :param x: axis x component.
        :param y: axis y component.
        :param z: axis z component.
        :return: MV with only PV components set.
        """
        return MultiVector(axix=x, axiy=y, axiz=z)
    
    @classmethod
    def rotor(cls, axis, theta):
        return MultiVector.rotor2(axis, cos(theta), sin(theta))
    
    @classmethod
    def rotor2(cls, axis, theta_cos, theta_sin):
        return MultiVector(real=theta_cos, axix=axis[4] * theta_sin, axiy=axis[5] * theta_sin, axiz=axis[6] * theta_sin)
    
    def rotate(self, axis, theta):
        """
        Used to get self rotated around an axis. Note, that this rotation will include scaling,
        if axis vector is not normalised. If you need to perform this operation a lot, consider used
        rotate2 of rotate3 instead.
        :param axis: axis to rotate around.
        :param theta: angle to rotate by.
        :return: new MV representing result of the rotation.
        """
        return self.rotate2(axis, cos(theta / 2), sin(theta / 2))
    
    def rotate2(self, axis, cos2, sin2):
        """
        More in-depth version of rotate.
        """
        return self.rotate3(MultiVector.rotor2(axis, cos2, sin2), MultiVector.rotor2(axis, cos2, -sin2))
    
    def rotate3(self, l_rot, r_rot):
        """
        The deepest version of rotate.
        """
        return (l_rot * self) * r_rot
    
    def scal_part(self):
        return self[0]
    
    def vec_part(self):
        return self[1], self[2], self[3]
    
    def axi_part(self):
        return self[4], self[5], self[6]
    
    def imag_part(self):
        return self[7]


MultiVector.ZERO = MultiVector()
ONE = MultiVector.ONE = MultiVector(real=1.0)
VECX = MultiVector.VECX = MultiVector(vecx=1.0)
VECY = MultiVector.VECY = MultiVector(vecy=1.0)
VECZ = MultiVector.VECZ = MultiVector(vecz=1.0)
AXIX = MultiVector.AXIX = MultiVector(axix=1.0)
AXIY = MultiVector.AXIY = MultiVector(axiy=1.0)
AXIZ = MultiVector.AXIZ = MultiVector(axiz=1.0)
IMAG = MultiVector.IMAG = MultiVector(imag=1.0)
METRIC = [
    [ONE, VECX, VECY, VECZ, AXIX, AXIY, AXIZ, IMAG],
    [VECX, ONE, AXIX, -AXIZ, VECY, IMAG, -VECZ, AXIY],
    [VECY, -AXIX, ONE, AXIY, -VECX, VECZ, IMAG, AXIZ],
    [VECZ, AXIZ, -AXIY, ONE, IMAG, -VECY, VECX, AXIX],
    [AXIX, -VECY, VECX, IMAG, -ONE, -AXIZ, AXIY, -VECZ],
    [AXIY, IMAG, -VECZ, VECY, AXIZ, -ONE, -AXIX, -VECX],
    [AXIZ, VECZ, IMAG, -VECX, -AXIY, AXIX, -ONE, -VECY],
    [IMAG, AXIY, AXIZ, AXIX, -VECZ, -VECX, -VECY, -ONE]
]
MultiVector.METRIC = [[res.comp for res in row] for row in METRIC]
MultiVector.VECTOR_BASIS = [VECX, VECY, VECZ]
MultiVector.AXIS_BASIS = [AXIX, AXIY, AXIZ]

del ONE, VECX, VECY, VECZ, AXIX, AXIY, AXIZ, IMAG, METRIC

if __name__ == "__main__":
    import colorsys
    import math
    import random
    
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.axis([-200, 200, -200, 200])
    plt.axis("equal")
    
    # init points
    points = [
        MultiVector.get_vector(100 * math.pow(random.randrange(200) / 100, 2 / 5) * (random.randrange(0, 2) - 0.5),
                               100 * math.pow(random.randrange(200) / 100, 2 / 5) * (random.randrange(0, 2) - 0.5),
                               100 * math.pow(random.randrange(200) / 100, 2 / 5) * (random.randrange(0, 2) - 0.5))
        for i in
        range(50)]
    axis = MultiVector.get_axis(10, -20, 45)
    axis *= 1 / math.hypot(*axis)
    
    circles = []
    for p in points:
        c = plt.Circle((p[1], p[2]), radius=0.6)
        c.initial_point = p
        circles.append(c)
        ax.add_artist(c)
    
        
    def init():
        return animate(0)
    
    
    def animate(frames):
        theta = 0.1 * frames
        cos2 = math.cos(theta / 2)
        sin2 = math.sin(theta / 2)
        l_rot = MultiVector.rotor2(axis, cos2, sin2)
        r_rot = MultiVector.rotor2(axis, cos2, -sin2)
        
        for c in circles:
            p = c.initial_point.rotate3(l_rot, r_rot)
            c.set_center((p[1], p[2]))
            c.set_color(colorsys.hsv_to_rgb(frames / 200, 1.0, math.pow((200 + p[3])/400, 2)))
        
        return circles
    
    
    total_frames = 600
    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=total_frames, blit=True)
    # ani.save('animation.mp4', fps=60, extra_args=['-vcodec', 'libx264'])
    plt.show()
