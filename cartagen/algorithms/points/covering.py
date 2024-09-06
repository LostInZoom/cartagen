import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import LineString, Polygon, Point, MultiPoint
from shapely.ops import triangulate, unary_union
from shapely import contains
import matplotlib.pyplot as plt
from cartagen.utils.geometry.segment import get_segment_list

def hull_delaunay(points, length):
    """
    Create the hull of a set of points using the Delaunay triangulation.

    This algorithm was described by Duckham *et al.* :footcite:p:`duckham:2008` and
    is also described in Galton & Duckham :footcite:p:`galton:2006`
    It first computes the Delaunay triangulation of the points, and then removes iteratively
    the boundary edges that are longer than a given parameter.

    Parameters
    ----------
    points : list of Point
        The set of points.
    length : float
        The length in meters below which the Delaunay triangles edges are removed.

    Returns
    -------
    Polygon

    Warning
    -------
    This algorithm cannot create multiple polygons, unlike :func:`hull_swinging_arm`.
    Using a length too low can produce an invalid geometry.

    See Also
    --------
    hull_swinging_arm :
        Create a polygon from a set of points using the swinging arm algorithm.

    References
    ----------
    .. footbibliography::

    Examples
    --------
    >>> points = [ Point(1, 1), Point(6, 6), Point(1, 6), Point(6, 5), Point(2, 4), Point(2, 1), Point(1, 4) ]
    >>> delaunay_hull(points, 2.0)
    <POLYGON ((1 1, 1 4, 1 6, 6 6, 6 5, 2 4, 2 1, 1 1))>
    """
    multipoint = MultiPoint(points)
    triangles = triangulate(multipoint)
    hull = unary_union(triangles)

    removed = []
    shrinkable = True
    while(shrinkable):
        hull_ext = hull.exterior
        longest = -1
        max_length = length
        i = 0
        for triangle in triangles:
            if(i in removed):
                i += 1
                continue

            # first check: if the triangle is inside the hull, pass to next
            if(triangle.intersects(hull_ext)==False):
                i += 1
                continue

            # regularity check: if the three points of the triangle are on the boundary, removing this triangle would create an invalid geometry
            if(hull_ext.intersects(Point(triangle.exterior.coords[0])) and hull_ext.intersects(Point(triangle.exterior.coords[1])) 
               and hull_ext.intersects(Point(triangle.exterior.coords[2]))):
                i += 1
                continue
            
            # check that we do not create a multipolygon
            if(hull.difference(triangle).geom_type == 'MultiPolygon'):
                i += 1
                continue

            segments = get_segment_list(triangle)
            for segment in segments:
                segment_length = segment.length()
                if segment_length > max_length:
                    max_length = segment_length
                    longest = i
            i += 1
        
        if longest == -1:
            break
        
        hull = hull.difference(triangles[longest])
        removed.append(longest)

    return hull

"""
Created on Mon Jun 24 19:03:29 2024

@author: vpech
"""

def hull_swinging_arm(points, length, direction='ccw'):
    """
    Create the hull of a set of points using the swinging arm algorithm.

    This algorithm is described in Galton & Duckham. :footcite:p:`galton:2006`
    Starting at the most northen point, iteratively create the border of a polygon
    by adding the first point intersecting the line (the arm) rotating around the
    point.

    Parameters
    ----------
    points : list of Point
        The set of points.
    length : int
        The length of the line (or arm) that will rotate around the point.
    direction : str, optional
        The direction of the arm's rotation. Can be 'ccw' or 'cw'.

    Returns
    -------
    list of Polygon

    See Also
    --------
    hull_delaunay :
        Create the hull of a set of points using the Delaunay triangulation.

    Notes
    -----
    This algorithm can create multiple polygons, unlike :func:`hull_delaunay`

    References
    ----------
    .. footbibliography::
    
    Examples
    -------
    >>> points = [ Point(1, 1), Point(6, 6), Point(1, 6), Point(6, 5), Point(2, 4), Point(2, 1), Point(1, 4) ]           
    >>> swinging_arm(points, 5.0)
    [<POLYGON ((1 4, 1 1, 2 1, 2 4, 1 4))>]
    """

    ############################ Initialization ############################
    def init_gdf(points):
        """
        Initialize the geodataframe of points for the SwingingArm algorithm.
        
        Add two columns: ‘available’ and ‘notVisited’. Then, sort the points by 
        their x-coordinates in ascending order.

        Parameters
        ----------
        points : GeodatFrame
        
        Returns
        -------
        GeodatFrame
            Point data initialized.
        """
        points = points.drop_duplicates(subset=['geometry'])
        points["available"] = True
        points["notVisited"] = True
        return points.iloc[points.geometry.x.argsort().values].copy()

    def ymax_point(points):
        """
        Find the point with the maximal y-coordinate.

        Parameters
        ----------
        points : GeoDataframe

        Returns
        -------
        Point
            Algorithm starting point.
        """
        ind_max = points.geometry.y.nlargest(1).index[0]
        return points.loc[ind_max, 'geometry']

    def line(pt1, pt2):
        """
        Create a LineString object given two Point objects.

        Parameters
        ----------
        pt1 : Point
        pt2 : Point
        
        Returns
        -------
        LineString
        """
        return LineString([[pt1.x, pt1.y], [pt2.x, pt2.y]])

    def draw_swingingArm(points, polygons):
        """
        Draw the result of the swinging arm algorithm.

        Parameters
        ----------
        polygons : List

        Returns
        -------
        None.

        """
        if polygons:
            fig, ax1 = plt.subplots()
            points.plot(ax=ax1, color='dodgerblue')
            gs = gpd.GeoSeries(polygons)
            gpd.GeoDataFrame(geometry=gs).boundary.plot(ax=ax1, color='firebrick')
            print(polygons)
            print("Polygons drawn:", len(polygons))
        else:
            print("Nothing has been drawn. Arm too short...")
            
    ############################ Angle calculations ############################

    def angle(pt, pt0, pt00, direction):
        """
        Calculate the angle form by the tree given point in a given direction.

        Parameters
        ----------
        pt : Point
            Last point.
        pt0 : Point
            Center Point.
        pt00 : Point
            First point.
        direction : String
            Angle calculation direction. The default is 'ccw'.

        Returns
        -------
        Float
        """
        def angle_three_pts(point1, point2, point3):
            angle1 = angle_two_pts(point2, point1)
            angle2 = angle_two_pts(point2, point3)
            return (angle2 - angle1)%(2*np.pi)

        def angle_two_pts(point1, point2):
            x = point2.x - point1.x
            y = point2.y - point1.y
            return np.arctan2(y, x)
        
        angle_deg = angle_three_pts(pt00, pt0, pt)*180/np.pi
        if direction == 'cw':
            angle_deg = 360-angle_deg
        if angle_deg == 0:
            return np.nan
        return angle_deg

    def angle_min(ptIn, pt0, points):
        """
        Find the point with the minimal angle with the last segment drawn.

        Parameters
        ----------
        ptIn : GeoDataframe
            Points inside the buffer of the point considered
        pt0 : Point
            The point for which we need to find the next point.
        points : GeoDataframe
            Points of the original Geodataframe.

        Returns
        -------
        tmp_pt : Point
            The point with the minimal angle.
        ind_min : int
            indice of the minimal point in the ptIn GeoDataframe.
        original_index : int
            indice of the minimal point in the original GeoDataframe.
        """
        ind_min = ptIn.angle.idxmin()
        tmp_pt = ptIn.loc[ind_min, 'geometry']
        mini = ptIn.angle.min()
        gdf_occ = gpd.GeoDataFrame(geometry=ptIn.loc[ptIn.angle == mini, "geometry"])

        if len(gdf_occ):
            gdf_occ['l'] = [shapely.length(line(pt, pt0)) for pt in gdf_occ.geometry]
            ind_min = gdf_occ.l.idxmin()
            tmp_pt = gdf_occ.loc[ind_min, 'geometry']

        original_index = points[points.geometry == tmp_pt].index[0]
        return tmp_pt, ind_min, original_index

    ############################ Check next point ############################

    def valid_point(tmp_pt, pt0, lines):
        """
        Check if a point doesn't form a line that intersects with other lines.

        Parameters
        ----------
        tmp_pt : Point
            Supposed next point.
        pt0 : Point
            Previous point.
        lines : LineString List
            LineString already drawn.

        Returns
        -------
        bool
            if the supposed next point formed a line that intersects with the other
            lines and that is not the starting point return False.
        """
        lines = lines[:-1]
        li = line(tmp_pt, pt0)
        for e in lines:
            if e.intersects(li) and lines[0].coords[0] != tmp_pt.coords[0]:
                return False
        return True

    ############################ Calculate angle ############################

    def next_point(points, r, lines, direction, pt0, pt00=None):
        """
        Find a supposed next point.

        Parameters
        ----------
        points : GeoDataframe
            The original dataframe of point.
        r : float
            lenght of the line segment anchored to the starting point.
        lines : LineString List
            LineString that form the current polygon.
        pt0 : Point
            Current point.
        pt00 : Starting point, optional
            Which forms the previous line. The default is None.
        direction : string, optional
            Direction of the angle calculation. The default is 'cw'.

        Returns
        -------
        tmp_pt : Point
            Supposed next point.
        """
        points_in = [row[0] for row in zip(points.geometry, points.available, points.notVisited) 
                if row[0].within(pt0.buffer(r)) and row[0] != pt0 and row[0] != pt00 and row[1] and row[2]]

        ptIn = gpd.GeoDataFrame(geometry=gpd.GeoSeries(points_in))
        tmp_pt = None
        
        if not(ptIn.empty):
            if not(pt00):
                ptIn["angle"] = [angle(pt, pt0, Point(pt0.x, pt0.y+r), direction) for pt in ptIn.geometry]
                tmp_pt, ind, ind_ori = angle_min(ptIn, pt0, points)
            else:
                ptIn["angle"] = [angle(pt, pt0, pt00, direction) for pt in ptIn.geometry]
                tmp_pt, ind, ind_ori = angle_min(ptIn, pt0, points)

                while(not(valid_point(tmp_pt, pt0, lines))):
                    points.loc[ind_ori, "notVisited"] = False
                    ptIn = ptIn.drop(ind)
                    if not(ptIn.empty):
                        tmp_pt, ind, ind_ori = angle_min(ptIn, pt0, points)
                    else:
                        return tmp_pt
            points.loc[ind_ori, "available"] = False
        return tmp_pt

    def swingSimpleLoop(points, r, direction, poly=None):
        """
        Loop implementation to draw a polygon from given points. 

        Parameters
        ----------
        points : GeoDataframe
            The original dataframe of point of the current polygon.
        r : int
            Lenght of the line segment anchored to the starting point.
        poly : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        Godataframe
            Next geodataframe of points.
        Polygon
            Polygon drawn with the given geodataframe.
        """
        points = init_gdf(points)
        pt_dep0 = ymax_point(points)
        lines = []
        pt_next = next_point(points, r, lines, direction, pt_dep0)

        if not(pt_next):
            ind = points[points.geometry == pt_dep0].index[0]
            return points.drop(ind), poly
        
        lines.append(line(pt_dep0, pt_next))
        pt_dep = pt_dep0
        pt_prev = pt_next

        while(pt_next.coords[0] != pt_dep0.coords[0]):
            pt_next = next_point(points, r, lines, direction, pt_prev, pt_dep)
            if not(pt_next):
                ind = points[points.geometry == pt_dep0].index[0]
                return points.drop(ind), poly
            pt_dep = pt_prev
            pt_prev = pt_next
            lines.append(line(pt_dep, pt_prev))

        polygon = Polygon([list(l.coords[0]) for l in lines])
        points['inside'] = contains(polygon, points.geometry)
        new_points = points.loc[points["available"] ^ points["inside"], "geometry"]
        new_points = gpd.GeoDataFrame(geometry=new_points)
        return new_points, polygon

    new_points, poly = swingSimpleLoop(gpd.GeoDataFrame([ { 'geometry': p } for p in points ]), length, direction)
    polygons = []

    while(not(new_points.empty)):
        if isinstance(poly, Polygon):
            polygons.append(poly)
        new_points, poly = swingSimpleLoop(new_points, length, direction, poly)
    if isinstance(poly, Polygon):
        polygons.append(poly)

    return polygons