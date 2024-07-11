import shapely
import numpy as np
from shapely import Point, Polygon, LineString
from shapely.geometry.polygon import orient

import matplotlib.pyplot as plt

from cartagen4py.utils.geometry.angle import angle_2_pts
from cartagen4py.utils.geometry.polygon import enclosing_rectangle

from cartagen4py.utils.debug import plot_debug

def rectangle_transformation(polygon, factor):
    """
    Transform a polygon to a rectangle.
    """
    mbr = polygon.minimum_rotated_rectangle
    return shapely.affinity.scale(mbr, xfact=factor, origin=mbr.centroid)

def recursive_regression(polygon, sigma):
    """
    Regularize a polygon using recursive linear regression.

    Parameters
    ----------
    polygon : Polygon
        The polygon to regularize.
    sigma : float
        The standard deviation threshold above which
        the recursion continues.
    """
    # Calculate the standard deviation of the side from the regression line
    def __get_sigma(side, regression):
        # Storage for the sum of the squared distances from the regression line
        square_total = []
        start, end = np.array(regression[0]), np.array(regression[-1])
        # Loop through side vertexes 
        for vertex in side:
            v = np.array(vertex)
            # Calculate the squared distance to the regression line
            square_total.append((np.cross(end - start, v - start) / np.linalg.norm(end - start))**2)
        # Calculate the standard deviation of the vertexes
        # as the suqrae root of the mean of the sum of the
        # squared distances from the regression line
        return np.sqrt(np.mean(square_total))

    # Find the closest vertex from the mbr border
    def __closest_vertex(corner, vertices):
        distance = None
        closest = None
        pe = Point(corner)
        for i, vertex in enumerate(vertices):
            d = shapely.distance(pe, Point(vertex))
            if distance is None:
                distance, closest = d, i
            else:
                if d < distance:
                    distance, closest = d, i
        return closest

    if shapely.is_ccw(polygon.boundary):
        polygon.reverse()

    # Calculate the minimum rotated rectangle
    mbr = enclosing_rectangle(polygon, mode='input')

    # Retrieve first and last vertex
    mbr_coords = list(mbr.boundary.coords)[:-1]

    angle = angle_2_pts(Point(mbr_coords[0]), Point(mbr_coords[-1]))

    # Rotate the polygon to reach horizontality
    rotated_mbr = shapely.affinity.rotate(mbr, -angle, origin=mbr_coords[0], use_radians=True)
    rotated_polygon = shapely.affinity.rotate(polygon, -angle, origin=mbr_coords[0], use_radians=True)

    # Get the coordinates of the four corner points
    mbr_coords = list(rotated_mbr.boundary.coords)
    extent = mbr_coords[:-1]

    # Get rotated polygon vertexes
    coords = list(rotated_polygon.boundary.coords)[:-1]

    # Storage for the corners of the polygons
    corners = []
    # Loop through the mbr corners and store
    # the closest point of the polygon from the mbr corners
    for e in extent:
        corners.append(__closest_vertex(e, coords))

    # Storage for the sides of the polygon
    sides = []

    # Add a list of vertex coordinates as
    # each of the four sides of the polygon
    side = []
    for i, vertex in enumerate(coords):
        side.append(vertex)
        if i in corners:
            sides.append(side)
            side = [vertex]
    sides.append(side)

    # Remove the last side and merge it with the first
    # Work when starting at 0 or more.
    last = sides.pop()
    first = sides.pop(0)
    sides.insert(0, last + first)
    
    # Loop through each polygon side
    for i, side in enumerate(sides):
        # Retrieve the associated mbr corners
        j = i + 1
        if j == len(extent):
            j = 0

        # Get the mbr corners coordinates
        mbr_side = [ extent[i], extent[j] ]
        # Get the angle of the mbr side in degree
        side_angle = round(np.rad2deg(angle_2_pts(Point(mbr_side[0]), Point(mbr_side[1]))))

        # Store the orientation of the rectangle side
        vertical = False
        if side_angle == 90 or side_angle == -90:
            vertical = True
          
        # Calculate the mean x or y of the vertexes of
        # the side depending on the orientation
        if vertical:
            mean = np.mean([ v[0] for v in side ])
        else:
            mean = np.mean([ v[1] for v in side ])

        # Define the regression line start and end point as
        # the mbr corners with the mean x or y as its new x or y value.
        if vertical:
            reg = [(mean, mbr_side[0][1]), (mean, mbr_side[1][1])]
        else:
            reg = [(mbr_side[0][0], mean), (mbr_side[1][0], mean)]

        regression = shapely.LineString(reg)

        # Calculate the standard deviation of the vertexes
        sig = __get_sigma(side, reg)
        
        # Continue if the sigma is above the threshold
        if sig > sigma:
            # Make sure there is more than 3 vertexes in the polygon side
            if len(side) > 3:
                # Storage for the divided parts of the side by the regression line.
                subsides = []
                subside = []
                for index in range(0, len(side) - 1):
                    # Add the vertex to the subside
                    subside.append(side[index])
                    # Create the segment between this point and the next
                    segment = shapely.LineString([side[index], side[index + 1]])
                    # If segment intersects the regression line
                    if shapely.intersects(segment, regression):
                        intersection = shapely.intersection(segment, regression)
                        # subside.append(intersection.coords[0])
                        subsides.append(subside)
                        # subside = [intersection.coords[0]]
                        subside = []

                subside.append(side[-1])
                subsides.append(subside)

                plot_debug(
                    rotated_polygon, regression,
                    [ Point(reg[0]), Point(reg[1]) ],
                    [ Point(s) for s in side ]
                )

                # for index, subside in enumerate(subsides):
                #     if len(subside) > 2:
                #         # Calculate the minimum rotated rectange
                #         submbr, subangle = __get_mbr(Polygon(subside + [subside[0]]))

                #         # Calculate the mean y of the vertex of the subside
                #         submean = np.mean([ s[1] for s in subside ])

                #         if not shapely.is_ccw(submbr.boundary):
                #             submbr = submbr.reverse()

                #         submbr_coords = submbr.boundary.coords

                #         # Here, the subside is above the regression line
                #         v1, v2 = None, None

                #         if submean > mean:
                #             i1, i2 = 0, 1
                #         else:
                #             i1, i2 = 1, 0

                #         if index == 0:
                #             v1, v2 = subside[0], subside[__closest_vertex(submbr_coords[i1], subside)]
                #         elif index == (len(subsides) - 1):
                #             v1, v2 = subside[__closest_vertex(submbr_coords[i2], subside)], subside[-1]
                #         else:
                #             v1, v2 = subside[__closest_vertex(submbr_coords[i2], subside)],  subside[__closest_vertex(submbr_coords[i1], subside)]   

                        # plot_debug(new, shapely.LineString(reg), submbr, Point(v1))
                
                break

    return polygon

def feature_edge_reconstruction(polygon):
    """
    Regularize a polygon using feature edge reconstruction.
    """