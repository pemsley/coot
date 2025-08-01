import math

def rotation_z(theta):
    """Generates a rotation matrix about the z-axis."""
    return [math.cos(theta), -math.sin(theta), 0,
            math.sin(theta), math.cos(theta), 0,
            0, 0, 1]

def rotation_x(theta):
    """Generates a rotation matrix about the x-axis."""
    return [1, 0, 0,
            0, math.cos(theta), -math.sin(theta),
            0, math.sin(theta), math.cos(theta)]

def rotation_y(theta):
    """Generates a rotation matrix about the y-axis."""
    return [math.cos(theta), 0, math.sin(theta),
            0, 1, 0,
            -math.sin(theta), 0, math.cos(theta)]

def rotation_c2_xy(theta):
    """Generates a C2 rotation matrix in the xy-plane."""
    return [math.cos(2 * theta), math.sin(2 * theta), 0,
            math.sin(2 * theta), -math.cos(2 * theta), 0,
            0, 0, -1]

def rotation_c3_111():
    """Generates a C3 rotation matrix along the (1,1,1) axis."""
    return [0, 0, 1, 1, 0, 0, 0, 1, 0]

def rotation_c3_minus1_minus1_1():
    """Generates a C3 rotation matrix along the (-1,-1,1) axis."""
    return [0, 1, 0, 0, 0, 1, 1, 0, 0]

def rotation_c3_1_minus1_minus1():
    """Generates a C3 rotation matrix along the (1,-1,-1) axis."""
    return [0, 1, 0, 1, 0, 0, 0, 0, 1]

def rotation_c3_minus1_1_minus1():
    """Generates a C3 rotation matrix along the (-1,1,-1) axis."""
    return [0, 0, 1, 0, 1, 0, 1, 0, 0]

def rotation_c3_minus1_1_1():
    """Generates a C3 rotation matrix along the (-1,1,-1) axis."""
    return [0, -1, 0, 0, 0, 1, -1, 0, 0]

def rotation_c3_1_minus1_1():
    """Generates a C3 rotation matrix along the (-1,1,-1) axis."""
    return [0, -1, 0, 0, 0, -1, 1, 0, 0]

def rotation_c3_minus1_minus1_minus1():
    """Generates a C3 rotation matrix along the (-1,-1,-1) axis."""
    return [0, 1, 0, 0, 0, 1, 1, 0, 0]

def multiply_matrices(matrix1, matrix2):
    """Multiplies two 3x3 matrices represented as lists."""
    result = []
    for i in range(3):  # Iterate through rows of the result
        for j in range(3):  # Iterate through columns of the result
            element = 0
            for k in range(3):  # Iterate through elements for dot product
                element += matrix1[i * 3 + k] * matrix2[k * 3 + j]
            result.append(element)
    return result

def rotation_c3_1_1_minus1():
    """Generates a C3 rotation matrix along the (11-1) axis."""
    return [1, 0, 0, 0, 1, 0, 0, 0, 1]


def rotation_c2_icosahedral(phi):
    """
    Generates a C2 rotation matrix for the icosahedral group.

    Args:
        phi: The angle around the z-axis to rotate the C2 axis.

    Returns:
        A 3x3 rotation matrix as a flattened list.
    """
    c_phi = math.cos(phi)
    s_phi = math.sin(phi)

    # Rotation matrix for C2 about an axis in the xy plane.
    # The axis is rotated by phi around the z-axis.
    rotation_matrix = [
        c_phi * c_phi - s_phi * s_phi, 2 * c_phi * s_phi, 0,
        2 * c_phi * s_phi, s_phi * s_phi - c_phi * c_phi, 0,
        0, 0, -1
    ]

    return rotation_matrix



def rotation_c3_icosahedral(phi, theta):
    """
    Generates a C3 rotation matrix for the icosahedral group.

    Args:
        phi: The azimuthal angle (around the z-axis).
        theta: The polar angle (from the z-axis).

    Returns:
        A 3x3 rotation matrix as a flattened list.
    """
    c_phi = math.cos(phi)
    s_phi = math.sin(phi)
    c_theta = math.cos(theta)
    s_theta = math.sin(theta)

    # Rotation matrix for a C3 rotation about an axis defined by (theta, phi)
    # The axis is defined in spherical coordinates.

    a = 1/math.sqrt(3)

    u_x = s_theta * c_phi;
    u_y = s_theta * s_phi;
    u_z = c_theta;

    rotation_matrix = [
        c_phi*c_phi*(1-a)+a, c_phi*s_phi*(1-a) - u_z*a*math.sqrt(3), c_phi*c_theta*(1-a) + u_y*a*math.sqrt(3),
        s_phi*c_phi*(1-a) + u_z*a*math.sqrt(3), s_phi*s_phi*(1-a)+a, s_phi*c_theta*(1-a) - u_x*a*math.sqrt(3),
        c_theta*c_phi*(1-a) - u_y*a*math.sqrt(3), c_theta*s_phi*(1-a) + u_x*a*math.sqrt(3), c_theta*c_theta*(1-a)+a,
    ]

    return rotation_matrix

def rotation_c5_icosahedral(phi, theta):
    """
    Generates a C5 rotation matrix for the icosahedral group.

    Args:
        phi: The azimuthal angle (around the z-axis).
        theta: The polar angle (from the z-axis).

    Returns:
        A 3x3 rotation matrix as a flattened list.
    """
    c_phi = math.cos(phi)
    s_phi = math.sin(phi)
    c_theta = math.cos(theta)
    s_theta = math.sin(theta)

    # Rotation matrix for a C5 rotation about an axis defined by (theta, phi)
    # The axis is defined in spherical coordinates.

    angle = 2 * math.pi / 5  # 72 degrees

    u_x = s_theta * c_phi
    u_y = s_theta * s_phi
    u_z = c_theta

    c = math.cos(angle)
    s = math.sin(angle)
    one_minus_c = 1 - c

    rotation_matrix = [
        u_x * u_x * one_minus_c + c, u_x * u_y * one_minus_c - u_z * s, u_x * u_z * one_minus_c + u_y * s,
        u_y * u_x * one_minus_c + u_z * s, u_y * u_y * one_minus_c + c, u_y * u_z * one_minus_c - u_x * s,
        u_z * u_x * one_minus_c - u_y * s, u_z * u_y * one_minus_c + u_x * s, u_z * u_z * one_minus_c + c,
    ]

    return rotation_matrix

def generate_point_group_matrices(point_group_symbol):
    """Generates rotation matrices for a given point group."""

    def flatten(m33):
        return [item for sublist in m33 for item in sublist]

    def zeroify(m0):
        m1 = [   0 if math.fabs(x)     < 0.00000001 else x for x in m0]
        m2 = [ 0.5 if math.fabs(x-0.5) < 0.00000001 else x for x in m1]
        m3 = [-0.5 if math.fabs(x+0.5) < 0.00000001 else x for x in m2]
        return m3

    matrices = {"E": [1, 0, 0, 0, 1, 0, 0, 0, 1]}  # Identity

    point_groups = {
        "C2": {
            "C2(z)": rotation_z(math.pi),
        },
        "C3": {
            "C3(z)":   rotation_z(2 * math.pi / 3),
            "C3(z)^2": rotation_z(4 * math.pi / 3),
        },
        "C4": {
            "C4(z)":   rotation_z(math.pi / 2),
            "C4(z)^2": rotation_z(math.pi),
            "C4(z)^3": rotation_z(3 * math.pi / 2),
        },
        "C5": {
            "C5(z)":   rotation_z(2 * math.pi / 5),
            "C5(z)^2": rotation_z(4 * math.pi / 5),
            "C5(z)^3": rotation_z(6 * math.pi / 5),
            "C5(z)^4": rotation_z(8 * math.pi / 5),
        },
        "C6": {
            "C6(z)":   rotation_z(math.pi / 3),
            "C6(z)^2": rotation_z(2 * math.pi / 3),
            "C6(z)^3": rotation_z(math.pi),
            "C6(z)^4": rotation_z(4 * math.pi / 3),
            "C6(z)^5": rotation_z(5 * math.pi / 3),
        },
        "D2": {
            "C2(z)": rotation_z(math.pi),
            "C2(x)": [1, 0, 0, 0, -1, 0, 0, 0, -1],
            "C2(y)": [-1, 0, 0, 0, 1, 0, 0, 0, -1],
        },
        "D3": {
            "C3(z)":   rotation_z(2 * math.pi / 3),
            "C3(z)^2": rotation_z(4 * math.pi / 3),
            "C2(1)":   rotation_c2_xy(0),
            "C2(2)":   rotation_c2_xy(math.pi / 3),
            "C2(3)":   rotation_c2_xy(-math.pi / 3),
        },
        "D4": {
            "C4(z)":   rotation_z(math.pi / 2),
            "C4(z)^2": rotation_z(math.pi),
            "C4(z)^3": rotation_z(3 * math.pi / 2),
            "C2(1)":   rotation_c2_xy(0),
            "C2(2)":   rotation_c2_xy(math.pi / 4),
            "C2(3)":   rotation_c2_xy(3 * math.pi / 4),
            "C2(4)":   rotation_c2_xy(5 * math.pi / 4),
        },
        "D5": {
            "C5(z)":   rotation_z(2 * math.pi / 5),
            "C5(z)^2": rotation_z(4 * math.pi / 5),
            "C5(z)^3": rotation_z(6 * math.pi / 5),
            "C5(z)^4": rotation_z(8 * math.pi / 5),
            "C2(1)": rotation_c2_xy(0),
            "C2(2)": rotation_c2_xy(math.pi / 5),
            "C2(3)": rotation_c2_xy(3 * math.pi / 5),
            "C2(4)": rotation_c2_xy(5 * math.pi / 5),
            "C2(5)": rotation_c2_xy(7 * math.pi / 5),
        },
        "D6": {
            "C6(z)": rotation_z(math.pi / 3),
            "C6(z)^2": rotation_z(2 * math.pi / 3),
            "C6(z)^3": rotation_z(math.pi),
            "C6(z)^4": rotation_z(4 * math.pi / 3),
            "C6(z)^5": rotation_z(5 * math.pi / 3),
            "C2(1)": rotation_c2_xy(0),
            "C2(2)": rotation_c2_xy(math.pi / 6),
            "C2(3)": rotation_c2_xy(3 * math.pi / 6),
            "C2(4)": rotation_c2_xy(5 * math.pi / 6),
            "C2(5)": rotation_c2_xy(7 * math.pi / 6),
            "C2(6)": rotation_c2_xy(9 * math.pi / 6),
        },
        "T": {
            "C2(x)": rotation_x(math.pi),
            "C2(y)": rotation_y(math.pi),
            "C2(z)": rotation_z(math.pi),
            "C3(111)": rotation_c3_111(),
            "C3(-1-11)": rotation_c3_minus1_minus1_1(),
            "C3(1-1-1)": rotation_c3_1_minus1_minus1(),
            "C3(-11-1)": rotation_c3_minus1_1_minus1(),
            "C3(11-1)": rotation_c3_1_1_minus1(),
            "C3(-1-1-1)": rotation_c3_minus1_minus1_minus1(),
            "C3(-111)": rotation_c3_minus1_1_1(),
            "C3(1-11)": rotation_c3_1_minus1_1(),
        },
        "Th": {
            "i": [-1, 0, 0, 0, -1, 0, 0, 0, -1],
            "S6(111)": multiply_matrices(rotation_c3_111(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(-1-11)": multiply_matrices(rotation_c3_minus1_minus1_1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(1-1-1)": multiply_matrices(rotation_c3_1_minus1_minus1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(-11-1)": multiply_matrices(rotation_c3_minus1_1_minus1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(11-1)": multiply_matrices(rotation_c3_1_1_minus1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(-1-1-1)": multiply_matrices(rotation_c3_minus1_minus1_minus1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(-111)": multiply_matrices(rotation_c3_minus1_1_1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(1-11)": multiply_matrices(rotation_c3_1_minus1_1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
        },
        "Td": {
            "C2(x)": rotation_x(math.pi),
            "C2(y)": rotation_y(math.pi),
            "C2(z)": rotation_z(math.pi),
            "C3(111)": rotation_c3_111(),
            "C3(-1-11)": rotation_c3_minus1_minus1_1(),
            "C3(1-1-1)": rotation_c3_1_minus1_minus1(),
            "C3(-11-1)": rotation_c3_minus1_1_minus1(),
            "C3(11-1)": rotation_c3_1_1_minus1(),
            "C3(-1-1-1)": rotation_c3_minus1_minus1_minus1(),
            "C3(-111)": rotation_c3_minus1_1_1(),
            "C3(1-11)": rotation_c3_1_minus1_1(),
            "S4(x)": multiply_matrices(rotation_x(math.pi/2), flatten([[-1,0,0],[0,1,0],[0,0,1]])),
            "S4(y)": multiply_matrices(rotation_y(math.pi/2), flatten([[1,0,0],[0,-1,0],[0,0,1]])),
            "S4(z)": multiply_matrices(rotation_z(math.pi/2), flatten([[1,0,0],[0,1,0],[0,0,-1]])),
            "S4(x)^-1": multiply_matrices(rotation_x(3*math.pi/2), flatten([[-1,0,0],[0,1,0],[0,0,1]])),
            "S4(y)^-1": multiply_matrices(rotation_y(3*math.pi/2), flatten([[1,0,0],[0,-1,0],[0,0,1]])),
            "S4(z)^-1": multiply_matrices(rotation_z(3*math.pi/2), flatten([[1,0,0],[0,1,0],[0,0,-1]])),
            "sigma_d(xy)": flatten([[0,1,0],[1,0,0],[0,0,1]]),
            "sigma_d(xz)": flatten([[0,0,1],[0,1,0],[1,0,0]]),
            "sigma_d(yz)": flatten([[1,0,0],[0,0,1],[0,1,0]]),
            "sigma_d(x-y)": flatten([[0,-1,0],[-1,0,0],[0,0,1]]),
            "sigma_d(x-z)": flatten([[0,0,-1],[0,1,0],[-1,0,0]]),
            "sigma_d(y-z)": flatten([[1,0,0],[0,0,-1],[0,-1,0]]),
        },
        "O": {
            "C4(x)": rotation_x(math.pi/2),
            "C4(x)^-1": rotation_x(3*math.pi/2),
            "C4(y)": rotation_y(math.pi/2),
            "C4(y)^-1": rotation_y(3*math.pi/2),
            "C4(z)": rotation_z(math.pi/2),
            "C4(z)^-1": rotation_z(3*math.pi/2),
            "C2(x)": rotation_x(math.pi),
            "C2(y)": rotation_y(math.pi),
            "C2(z)": rotation_z(math.pi),
            "C3(111)": rotation_c3_111(),
            "C3(-1-11)": rotation_c3_minus1_minus1_1(),
            "C3(1-1-1)": rotation_c3_1_minus1_minus1(),
            "C3(-11-1)": rotation_c3_minus1_1_minus1(),
            "C3(11-1)": rotation_c3_1_1_minus1(),
            "C3(-1-1-1)": rotation_c3_minus1_minus1_minus1(),
            "C3(-111)": rotation_c3_minus1_1_1(),
            "C3(1-11)": rotation_c3_1_minus1_1(),
            "C2(xy)": flatten([[0,1,0],[1,0,0],[0,0,1]]),
            "C2(xz)": flatten([[0,0,1],[0,1,0],[1,0,0]]),
            "C2(yz)": flatten([[1,0,0],[0,0,1],[0,1,0]]),
            "C2(x-y)": flatten([[0,-1,0],[-1,0,0],[0,0,1]]),
            "C2(x-z)": flatten([[0,0,-1],[0,1,0],[-1,0,0]]),
            "C2(y-z)": flatten([[1,0,0],[0,0,-1],[0,-1,0]]),
        },
        "Oh": {
            "i": [-1, 0, 0, 0, -1, 0, 0, 0, -1],
            "C4(x)": rotation_x(math.pi/2),
            "C4(x)^-1": rotation_x(3*math.pi/2),
            "C4(y)": rotation_y(math.pi/2),
            "C4(y)^-1": rotation_y(3*math.pi/2),
            "C4(z)": rotation_z(math.pi/2),
            "C4(z)^-1": rotation_z(3*math.pi/2),
            "C2(x)": rotation_x(math.pi),
            "C2(y)": rotation_y(math.pi),
            "C2(z)": rotation_z(math.pi),
            "C3(111)": rotation_c3_111(),
            "C3(-1-11)": rotation_c3_minus1_minus1_1(),
            "C3(1-1-1)": rotation_c3_1_minus1_minus1(),
            "C3(-11-1)": rotation_c3_minus1_1_minus1(),
            "C3(11-1)": rotation_c3_1_1_minus1(),
            "C3(-1-1-1)": rotation_c3_minus1_minus1_minus1(),
            "C3(-111)": rotation_c3_minus1_1_1(),
            "C3(1-11)": rotation_c3_1_minus1_1(),
            "C2(xy)": flatten([[0,1,0],[1,0,0],[0,0,1]]),
            "C2(xz)": flatten([[0,0,1],[0,1,0],[1,0,0]]),
            "C2(yz)": flatten([[1,0,0],[0,0,1],[0,1,0]]),
            "C2(x-y)": flatten([[0,-1,0],[-1,0,0],[0,0,1]]),
            "C2(x-z)": flatten([[0,0,-1],[0,1,0],[-1,0,0]]),
            "C2(y-z)": flatten([[1,0,0],[0,0,-1],[0,-1,0]]),
            "sigma_h": flatten([[1,0,0],[0,1,0],[0,0,-1]]),
            "sigma_v(x)": flatten([[-1,0,0],[0,1,0],[0,0,1]]),
            "sigma_v(y)": flatten([[1,0,0],[0,-1,0],[0,0,1]]),
            "sigma_v(z)": flatten([[1,0,0],[0,1,0],[0,0,-1]]),
            "sigma_d(xy)": flatten([[0,1,0],[1,0,0],[0,0,1]]),
            "sigma_d(xz)": flatten([[0,0,1],[0,1,0],[1,0,0]]),
            "sigma_d(yz)": flatten([[1,0,0],[0,0,1],[0,1,0]]),
            "sigma_d(x-y)": flatten([[0,-1,0],[-1,0,0],[0,0,1]]),
            "sigma_d(x-z)": flatten([[0,0,-1],[0,1,0],[-1,0,0]]),
            "sigma_d(y-z)": flatten([[1,0,0],[0,0,-1],[0,-1,0]]),
            "S4(x)": multiply_matrices(rotation_x(math.pi/2), flatten([[-1,0,0],[0,1,0],[0,0,1]])),
            "S4(y)": multiply_matrices(rotation_y(math.pi/2), flatten([[1,0,0],[0,-1,0],[0,0,1]])),
            "S4(z)": multiply_matrices(rotation_z(math.pi/2), flatten([[1,0,0],[0,1,0],[0,0,-1]])),
            "S4(x)^-1": multiply_matrices(rotation_x(3*math.pi/2), flatten([[-1,0,0],[0,1,0],[0,0,1]])),
            "S4(y)^-1": multiply_matrices(rotation_y(3*math.pi/2), flatten([[1,0,0],[0,-1,0],[0,0,1]])),
            "S4(z)^-1": multiply_matrices(rotation_z(3*math.pi/2), flatten([[1,0,0],[0,1,0],[0,0,-1]])),
            "S6(111)": multiply_matrices(rotation_c3_111(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(-1-11)": multiply_matrices(rotation_c3_minus1_minus1_1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(1-1-1)": multiply_matrices(rotation_c3_1_minus1_minus1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(-11-1)": multiply_matrices(rotation_c3_minus1_1_minus1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(11-1)": multiply_matrices(rotation_c3_1_1_minus1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(-1-1-1)": multiply_matrices(rotation_c3_minus1_minus1_minus1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(-111)": multiply_matrices(rotation_c3_minus1_1_1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(1-11)": multiply_matrices(rotation_c3_1_minus1_1(), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
        },
        "I": {
            "C2(1)": rotation_c2_icosahedral(0),
            "C2(2)": rotation_c2_icosahedral(math.pi/5),
            "C2(3)": rotation_c2_icosahedral(2*math.pi/5),
            "C2(4)": rotation_c2_icosahedral(3*math.pi/5),
            "C2(5)": rotation_c2_icosahedral(4*math.pi/5),
            "C2(6)": rotation_c2_icosahedral(math.pi),
            "C2(7)": rotation_c2_icosahedral(6*math.pi/5),
            "C2(8)": rotation_c2_icosahedral(7*math.pi/5),
            "C2(9)": rotation_c2_icosahedral(8*math.pi/5),
            "C2(10)": rotation_c2_icosahedral(9*math.pi/5),
            "C3(1)": rotation_c3_icosahedral(0, 0),
            "C3(2)": rotation_c3_icosahedral(2*math.pi/5, 0),
            "C3(3)": rotation_c3_icosahedral(4*math.pi/5, 0),
            "C3(4)": rotation_c3_icosahedral(6*math.pi/5, 0),
            "C3(5)": rotation_c3_icosahedral(8*math.pi/5, 0),
            "C3(6)": rotation_c3_icosahedral(0, math.acos(1/math.sqrt(5))),
            "C3(7)": rotation_c3_icosahedral(2*math.pi/5, math.acos(1/math.sqrt(5))),
            "C3(8)": rotation_c3_icosahedral(4*math.pi/5, math.acos(1/math.sqrt(5))),
            "C3(9)": rotation_c3_icosahedral(6*math.pi/5, math.acos(1/math.sqrt(5))),
            "C3(10)": rotation_c3_icosahedral(8*math.pi/5, math.acos(1/math.sqrt(5))),
            "C3(11)": rotation_c3_icosahedral(0, -math.acos(1/math.sqrt(5))),
            "C3(12)": rotation_c3_icosahedral(2*math.pi/5, -math.acos(1/math.sqrt(5))),
            "C3(13)": rotation_c3_icosahedral(4*math.pi/5, -math.acos(1/math.sqrt(5))),
            "C3(14)": rotation_c3_icosahedral(6*math.pi/5, -math.acos(1/math.sqrt(5))),
            "C3(15)": rotation_c3_icosahedral(8*math.pi/5, -math.acos(1/math.sqrt(5))),
            "C5(1)": rotation_c5_icosahedral(0, 0),
            "C5(2)": rotation_c5_icosahedral(2*math.pi/5, 0),
            "C5(3)": rotation_c5_icosahedral(4*math.pi/5, 0),
            "C5(4)": rotation_c5_icosahedral(6*math.pi/5, 0),
            "C5(5)": rotation_c5_icosahedral(8*math.pi/5, 0),
            "C5(6)": rotation_c5_icosahedral(0, math.acos(math.sqrt(5)/3)),
            "C5(7)": rotation_c5_icosahedral(2*math.pi/5, math.acos(math.sqrt(5)/3)),
            "C5(8)": rotation_c5_icosahedral(4*math.pi/5, math.acos(math.sqrt(5)/3)),
            "C5(9)": rotation_c5_icosahedral(6*math.pi/5, math.acos(math.sqrt(5)/3)),
            "C5(10)": rotation_c5_icosahedral(8*math.pi/5, math.acos(math.sqrt(5)/3)),
            "C5(11)": rotation_c5_icosahedral(0, -math.acos(math.sqrt(5)/3)),
            "C5(12)": rotation_c5_icosahedral(2*math.pi/5, -math.acos(math.sqrt(5)/3)),
            "C5(13)": rotation_c5_icosahedral(4*math.pi/5, -math.acos(math.sqrt(5)/3)),
            "C5(14)": rotation_c5_icosahedral(6*math.pi/5, -math.acos(math.sqrt(5)/3)),
            "C5(15)": rotation_c5_icosahedral(8*math.pi/5, -math.acos(math.sqrt(5)/3)),
        },
        "Ih-broken": {
            "i": [-1, 0, 0, 0, -1, 0, 0, 0, -1],
            "C2(1)": rotation_c2_icosahedral(0),
            "C2(2)": rotation_c2_icosahedral(math.pi/5),
            "C2(3)": rotation_c2_icosahedral(2*math.pi/5),
            "C2(4)": rotation_c2_icosahedral(3*math.pi/5),
            "C2(5)": rotation_c2_icosahedral(4*math.pi/5),
            "C2(6)": rotation_c2_icosahedral(math.pi),
            "C2(7)": rotation_c2_icosahedral(6*math.pi/5),
            "C2(8)": rotation_c2_icosahedral(7*math.pi/5),
            "C2(9)": rotation_c2_icosahedral(8*math.pi/5),
            "C2(10)": rotation_c2_icosahedral(9*math.pi/5),
            "C3(1)": rotation_c3_icosahedral(0, 0),
            "C3(2)": rotation_c3_icosahedral(2*math.pi/5, 0),
            "C3(3)": rotation_c3_icosahedral(4*math.pi/5, 0),
            "C3(4)": rotation_c3_icosahedral(6*math.pi/5, 0),
            "C3(5)": rotation_c3_icosahedral(8*math.pi/5, 0),
            "C3(6)": rotation_c3_icosahedral(0, math.acos(1/math.sqrt(5))),
            "C3(7)": rotation_c3_icosahedral(2*math.pi/5, math.acos(1/math.sqrt(5))),
            "C3(8)": rotation_c3_icosahedral(4*math.pi/5, math.acos(1/math.sqrt(5))),
            "C3(9)": rotation_c3_icosahedral(6*math.pi/5, math.acos(1/math.sqrt(5))),
            "C3(10)": rotation_c3_icosahedral(8*math.pi/5, math.acos(1/math.sqrt(5))),
            "C3(11)": rotation_c3_icosahedral(0, -math.acos(1/math.sqrt(5))),
            "C3(12)": rotation_c3_icosahedral(2*math.pi/5, -math.acos(1/math.sqrt(5))),
            "C3(13)": rotation_c3_icosahedral(4*math.pi/5, -math.acos(1/math.sqrt(5))),
            "C3(14)": rotation_c3_icosahedral(6*math.pi/5, -math.acos(1/math.sqrt(5))),
            "C3(15)": rotation_c3_icosahedral(8*math.pi/5, -math.acos(1/math.sqrt(5))),
            "C5(1)": rotation_c5_icosahedral(0, 0),
            "C5(2)": rotation_c5_icosahedral(2*math.pi/5, 0),
            "C5(3)": rotation_c5_icosahedral(4*math.pi/5, 0),
            "C5(4)": rotation_c5_icosahedral(6*math.pi/5, 0),
            "C5(5)": rotation_c5_icosahedral(8*math.pi/5, 0),
            "C5(6)": rotation_c5_icosahedral(0, math.acos(math.sqrt(5)/3)),
            "C5(7)": rotation_c5_icosahedral(2*math.pi/5, math.acos(math.sqrt(5)/3)),
            "C5(8)": rotation_c5_icosahedral(4*math.pi/5, math.acos(math.sqrt(5)/3)),
            "C5(9)": rotation_c5_icosahedral(6*math.pi/5, math.acos(math.sqrt(5)/3)),
            "C5(10)": rotation_c5_icosahedral(8*math.pi/5, math.acos(math.sqrt(5)/3)),
            "C5(11)": rotation_c5_icosahedral(0, -math.acos(math.sqrt(5)/3)),
            "C5(12)": rotation_c5_icosahedral(2*math.pi/5, -math.acos(math.sqrt(5)/3)),
            "C5(13)": rotation_c5_icosahedral(4*math.pi/5, -math.acos(math.sqrt(5)/3)),
            "C5(14)": rotation_c5_icosahedral(6*math.pi/5, -math.acos(math.sqrt(5)/3)),
            "C5(15)": rotation_c5_icosahedral(8*math.pi/5, -math.acos(math.sqrt(5)/3)),
            "sigma_h(1)": multiply_matrices(rotation_c2_icosahedral(0), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(2)": multiply_matrices(rotation_c2_icosahedral(math.pi/5), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(3)": multiply_matrices(rotation_c2_icosahedral(2*math.pi/5), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(4)": multiply_matrices(rotation_c2_icosahedral(3*math.pi/5), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(5)": multiply_matrices(rotation_c2_icosahedral(4*math.pi/5), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(6)": multiply_matrices(rotation_c2_icosahedral(math.pi), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(7)": multiply_matrices(rotation_c2_icosahedral(6*math.pi/5), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(8)": multiply_matrices(rotation_c2_icosahedral(7*math.pi/5), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(9)": multiply_matrices(rotation_c2_icosahedral(8*math.pi/5), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "sigma_h(10)": multiply_matrices(rotation_c2_icosahedral(9*math.pi/5), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(1)": multiply_matrices(rotation_c3_icosahedral(0, 0), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(2)": multiply_matrices(rotation_c3_icosahedral(2*math.pi/5, 0), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(3)": multiply_matrices(rotation_c3_icosahedral(4*math.pi/5, 0), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
            "S6(4)": multiply_matrices(rotation_c3_icosahedral(6*math.pi/5, 0), [-1, 0, 0, 0, -1, 0, 0, 0, -1]),
        },
    }

    if point_group_symbol in point_groups:
        matrices.update(point_groups[point_group_symbol])
        zm = {}
        for key in matrices:
            zm[key] = zeroify(matrices[key])
        return zm
    else:
        return "Point group not supported"


if __name__ == "__main__":
    print(generate_point_group_matrices("C2"))
    print(generate_point_group_matrices("C3"))
    print(generate_point_group_matrices("C4"))
    print(generate_point_group_matrices("C5"))
    print(generate_point_group_matrices("D2"))
    print(generate_point_group_matrices("D3"))
    print(generate_point_group_matrices("D4"))

# https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22logical_operator%22%3A%22and%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22group%22%2C%22nodes%22%3A%5B%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22text%22%2C%22parameters%22%3A%7B%22attribute%22%3A%22rcsb_struct_symmetry.symbol%22%2C%22operator%22%3A%22exact_match%22%2C%22negation%22%3Afalse%2C%22value%22%3A%22D4%22%7D%7D%5D%2C%22logical_operator%22%3A%22and%22%7D%5D%2C%22label%22%3A%22text%22%7D%5D%7D%2C%22return_type%22%3A%22entry%22%2C%22request_options%22%3A%7B%22paginate%22%3A%7B%22start%22%3A0%2C%22rows%22%3A25%7D%2C%22results_content_type%22%3A%5B%22experimental%22%5D%2C%22sort%22%3A%5B%7B%22sort_by%22%3A%22score%22%2C%22direction%22%3A%22desc%22%7D%5D%2C%22scoring_strategy%22%3A%22combined%22%7D%2C%22request_info%22%3A%7B%22query_id%22%3A%221efe9d9f2c283fc40e6f0d8559eaf48f%22%7D%7D
