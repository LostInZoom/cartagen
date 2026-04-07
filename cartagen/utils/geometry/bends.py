import numpy as np
from shapely.geometry import Point, LineString, Polygon
from typing import List, Set, Optional, Tuple

from cartagen.utils.lines.smoothing.gaussian import smooth_gaussian
from cartagen.utils.geometry.angle import get_curvature
from cartagen.utils.geometry.line import get_line_middle_point, get_bend_side


class Bend:
    """
    A bend in a line (e.g. in a mountain road).
    
    The extent of the bend has to be computed in some way.
    A bend is defined between two inflection points.
    """
    
    _counter = 0
    
    def __init__(self, geometry: LineString, summit_sigma: float = 15.0):
        """
        Initialize a Bend object.
        
        Parameters
        ----------
        geometry : LineString
            The geometry of the bend.
        summit_sigma : float, optional
            Sigma value for smoothing the line when the summit is computed.
            Default is 15.0.
        """
        self.bend = geometry
        self.summit_sigma = summit_sigma
        self.summit = None
        self.summit_index = -1
        self.head = None
        
        # Assign unique ID
        self.id = Bend._counter
        Bend._counter += 1
    
    def get_geom(self) -> LineString:
        """
        Get the geometry of the bend.
                
        Returns
        -------
        LineString
        """
        return self.bend

    def get_bend_side(self) -> str:
        """
        Get the side of the interior of the bend regarding the line direction.
        
        The algorithm comes from the AGENT project.
        
        Returns
        -------
        str
            'left' or 'right' indicating the side of the bend.
        """
        return get_bend_side(self.bend)
    
    def close_bend(self) -> Polygon:
        """
        Close the bend by joining both extremities.
        
        Returns
        -------
        Polygon
            A polygon formed by closing the bend.
        """
        coords = list(self.bend.coords)
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        return Polygon(coords)
    
    def get_length(self) -> float:
        """
        Measure bend length, i.e. the length of the line.
        
        Returns
        -------
        float
            The length of the bend.
        """
        return self.bend.length
    
    def get_width(self) -> float:
        """
        Measure the width of the bend.
        
        The width is the length of the base segment drawn between 
        the inflection points that create the bend.
        
        Returns
        -------
        float
            The width of the bend.
        """
        coords = list(self.bend.coords)
        start = Point(coords[0])
        end = Point(coords[-1])
        return start.distance(end)
    
    def get_base_middle(self) -> Point:
        """
        Get the middle point of the base segment.
        
        Returns
        -------
        Point
            The middle point of the base segment.
        """
        coords = list(self.bend.coords)
        start = Point(coords[0])
        end = Point(coords[-1])
        return Point((start.x + end.x) / 2, (start.y + end.y) / 2)
    
    def get_bend_summit(self) -> Point:
        """
        Get the summit point of the bend.
        
        The summit is the point with maximum curvature after smoothing.
        
        Returns
        -------
        Point
            The summit of the bend.
        """
        if self.summit is not None:
            return self.summit
        
        coords = list(self.bend.coords)
        
        # Handle edge case: too few points
        if len(coords) < 3:
            mid_idx = len(coords) // 2
            self.summit = Point(coords[mid_idx])
            self.summit_index = mid_idx
            return self.summit
        
        # Smooth the line with the specified sigma
        smooth_line = smooth_gaussian(
            self.bend, 
            sigma=self.summit_sigma, 
            densify=True
        )
        
        smooth_coords = list(smooth_line.coords)
        
        # Handle edge case: smoothed line too short
        if len(smooth_coords) < 3:
            mid_idx = len(coords) // 2
            self.summit = Point(coords[mid_idx])
            self.summit_index = mid_idx
            return self.summit
        
        # Find point with maximum curvature
        max_curvature = 0.0
        smooth_summit = None
        
        for i in range(1, len(smooth_coords) - 1):
            p1 = Point(smooth_coords[i - 1])
            p2 = Point(smooth_coords[i])
            p3 = Point(smooth_coords[i + 1])
            
            curvature = get_curvature(p1, p2, p3)
            
            if curvature >= max_curvature:
                max_curvature = curvature
                smooth_summit = p2
        
        # If no summit found (all points aligned), use middle point
        if smooth_summit is None:
            mid_idx = len(coords) // 2
            self.summit = Point(coords[mid_idx])
            self.summit_index = mid_idx
            return self.summit
        
        # Find nearest vertex in original line to smooth summit
        min_distance = float('inf')
        for i, coord in enumerate(coords):
            point = Point(coord)
            distance = point.distance(smooth_summit)
            if distance < min_distance:
                min_distance = distance
                self.summit = point
                self.summit_index = i
        
        return self.summit
    
    def get_height(self) -> float:
        """
        Measure bend height.
        
        The height is the length of the segment between the summit 
        of the bend and the middle of the base segment.
        
        Returns
        -------
        float
            The height of the bend.
        """
        base_middle = self.get_base_middle()
        summit = self.get_bend_summit()
        return summit.distance(base_middle)
    
    def get_symmetry(self) -> float:
        """
        Compute a measure of bend symmetry.
        
        If D1 is the curvilinear length from bend start to bend summit 
        and D2 the curvilinear length from bend summit to bend end,
        symmetry is Min(D1,D2) / Max(D1,D2).
        
        Returns
        -------
        float
            A value between 0 and 1, where 1 is perfectly symmetric.
        """
        if self.summit is None:
            self.get_bend_summit()
        
        coords = list(self.bend.coords)
        
        # Calculate D1: distance from start to summit
        d1 = 0.0
        for i in range(self.summit_index):
            p1 = Point(coords[i])
            p2 = Point(coords[i + 1])
            d1 += p1.distance(p2)
        
        # Calculate D2: distance from summit to end
        d2 = 0.0
        for i in range(self.summit_index, len(coords) - 1):
            p1 = Point(coords[i])
            p2 = Point(coords[i + 1])
            d2 += p1.distance(p2)
        
        if max(d1, d2) == 0:
            return 1.0
        
        return min(d1, d2) / max(d1, d2)
    
    def get_shape_measure(self) -> float:
        """
        Compute a shape measure of the bend.
        
        Ported from PlaGe and F. Lecordix schematisation algorithm.
        
        Returns
        -------
        float
            The shape measure.
        """
        return (4.0 * self.get_height() + self.get_width()) / 5.0
    
    def get_size_measure(self) -> float:
        """
        Compute a size measure of the bend.
        
        Ported from PlaGe and F. Lecordix schematisation algorithm.
        
        Returns
        -------
        float
            The size measure.
        """
        height = self.get_height()
        width = self.get_width()
        length = self.get_length()
        return (4.0 * height * width + length * width) / 5.0
    
    def get_orientation(self) -> float:
        """
        Compute the orientation of the bend.
        
        The orientation is the angle of the segment between the middle 
        of the base segment and the summit.
        
        Returns
        -------
        float
            The orientation between 0 and 2*pi radians.
        """
        base_middle = self.get_base_middle()
        summit = self.get_bend_summit()
        
        # Calculate orientation
        dx = summit.x - base_middle.x
        dy = summit.y - base_middle.y
        
        orientation = np.arctan2(dy, dx)
        
        # Normalize to [0, 2*pi]
        if orientation < 0:
            orientation += 2 * np.pi
        
        return orientation
    
    def __hash__(self):
        return self.id
    
    def __eq__(self, other):
        if not isinstance(other, Bend):
            return False
        return self.id == other.id


class BendSeries:
    """
    A BendSeries feature is a part of a linear section composed 
    of a series of sinuous bends.
    
    The series is decomposed into individual bends using inflection points.
    """
    
    def __init__(self, geometry: LineString, sigma: float = 15.0, sample: int = None):
        """
        Initialize a BendSeries object.
        
        Parameters
        ----------
        geometry : LineString
            The geometry of the bend series.
        sigma_smoothing : float, optional
            Sigma value for Gaussian smoothing. Default is 15.0.
        """
        self.geometry = geometry
        self.sigma = sigma
        self.sample = sample
        self.inflection_pts = []
        self.bends = []
        
        # Compute the bend series
        self._compute_bend_series()
    
    def _compute_bend_series(self):
        """
        Compute the bend series by identifying inflection points and bends.
        """
        # Smooth the line
        smooth_line = smooth_gaussian(
            self.geometry, 
            sigma=self.sigma,
            sample=self.sample,
            densify=True
        )
        
        # Compute bend sequences (consecutive vertices with same turn direction)
        sequences = self._compute_bend_sequences(smooth_line)
        
        # Filter micro-inflections
        filtered_sequences = self._filter_bend_sequences(sequences, 1)
        
        # Compute inflection points
        self._compute_inflection_points(filtered_sequences)
        
        # Compute bends from inflection points
        self._compute_bends_from_inflection_pts()
    
    def _compute_bend_sequences(self, smooth_line: LineString) -> List[int]:
        """
        Identify inflection points in a polyline in a basic way.
        
        An inflection point is identified at each change of direction 
        of consecutive angles. Returns a list containing the number of 
        consecutive vertices (a sequence) with the same angle direction.
        
        Parameters
        ----------
        smooth_line : LineString
            The smoothed line to analyze.
        
        Returns
        -------
        List[int]
            List of sequence lengths (number of consecutive vertices 
            with same turn direction).
        """
        sequences = []
        coords = list(smooth_line.coords)
        
        # If the polyline has fewer than three points, there are no angles
        if len(coords) < 3:
            return sequences
        
        nb_vertices = 1
        previous_direction = None
        
        for i in range(len(coords) - 2):
            p1 = Point(coords[i])
            p2 = Point(coords[i + 1])
            p3 = Point(coords[i + 2])
            
            # Calculate angle using three points
            v1 = np.array([p2.x - p1.x, p2.y - p1.y])
            v2 = np.array([p3.x - p2.x, p3.y - p2.y])
            
            # Cross product to determine turn direction
            cross = v1[0] * v2[1] - v1[1] * v2[0]
            
            # Determine the direction of the angle
            current_direction = "left" if cross > 0 else "right"
            
            if i > 0:
                # Identify a change of direction
                if current_direction != previous_direction:
                    sequences.append(nb_vertices)
                    nb_vertices = 1
                else:
                    # No change of direction
                    nb_vertices += 1
            
            previous_direction = current_direction
        
        # Add the last sequence
        if nb_vertices > 0:
            sequences.append(nb_vertices)
        
        return sequences
    
    def _filter_bend_sequences(
        self, 
        bend_sequence: List[int], 
        filtering: int
    ) -> List[int]:
        """
        Filter out all micro-inflections in a polyline.
        
        Parameters
        ----------
        bend_sequence : List[int]
            The raw sequence of consecutive vertices.
        filtering : int
            The threshold for filtering (typically 1).
        
        Returns
        -------
        List[int]
            The filtered sequence.
        """
        filtered_sequence = bend_sequence.copy()
        
        if len(filtered_sequence) <= 1:
            return filtered_sequence
        
        # Handle case with only 2 sequences
        if len(filtered_sequence) == 2:
            if filtered_sequence[0] == filtering:
                filtered_vertex_nb = filtered_sequence[1] + 1
                return [filtered_vertex_nb]
            elif filtered_sequence[1] == filtering:
                filtered_vertex_nb = filtered_sequence[0] + 1
                return [filtered_vertex_nb]
            return filtered_sequence
        
        # Handle case with more than 2 sequences
        i = 0
        while i < len(filtered_sequence):
            if filtered_sequence[i] == filtering:
                # Case: remove micro-inflection at the beginning
                if i == 0:
                    filtered_vertex_nb = filtered_sequence[1] + 1
                    filtered_sequence.pop(0)
                    filtered_sequence[0] = filtered_vertex_nb
                    i = -1
                
                # Case: remove micro-inflection at the end
                elif i == len(filtered_sequence) - 1:
                    filtered_vertex_nb = filtered_sequence[-2] + 1
                    filtered_sequence.pop()
                    filtered_sequence[-1] = filtered_vertex_nb
                
                # Case: remove micro-inflection in the middle
                else:
                    filtered_vertex_nb = (
                        filtered_sequence[i - 1] + 
                        filtered_sequence[i + 1] + 1
                    )
                    filtered_sequence.pop(i - 1)
                    filtered_sequence.pop(i - 1)
                    filtered_sequence.insert(i - 1, filtered_vertex_nb)
                    i = -1
            
            i += 1
        
        return filtered_sequence
    
    def _compute_inflection_points(self, filtered_sequence: List[int]):
        """
        Determine inflection points from the input linestring and 
        the list of sequences of consecutive vertices in the same direction.
        
        Parameters
        ----------
        filtered_sequence : List[int]
            The filtered sequence of consecutive vertices.
        """
        self.inflection_pts = []
        coords = list(self.geometry.coords)
        
        # Add the first point
        self.inflection_pts.append(Point(coords[0]))
        
        # Determine inflection points on the original polyline
        position = 0
        for sequence_length in filtered_sequence:
            position += sequence_length
            if position < len(coords):
                self.inflection_pts.append(Point(coords[position]))
        
        # Add the last point
        self.inflection_pts.append(Point(coords[-1]))
    
    def _compute_bends_from_inflection_pts(self):
        """
        Compute the bends of the series from inflection points.
        
        A bend is the part of line between two consecutive inflection points.
        """
        self.bends = []
        coords = list(self.geometry.coords)
        
        # Find indices of inflection points in original coordinates
        inflection_indices = []
        for inflection_pt in self.inflection_pts:
            min_dist = float('inf')
            closest_idx = 0
            for i, coord in enumerate(coords):
                dist = inflection_pt.distance(Point(coord))
                if dist < min_dist:
                    min_dist = dist
                    closest_idx = i
            inflection_indices.append(closest_idx)
        
        # Create bends between consecutive inflection points
        for i in range(len(inflection_indices) - 1):
            start_idx = inflection_indices[i]
            end_idx = inflection_indices[i + 1]
            
            if start_idx == end_idx:
                continue
            
            # Extract sub-line between inflection points
            bend_coords = coords[start_idx:end_idx + 1]
            if len(bend_coords) >= 2:
                bend_line = LineString(bend_coords)
                self.bends.append(Bend(bend_line))
    
    def get_all_summits(self) -> Set[Point]:
        """
        Get all summits of the bend series.
        
        Returns
        -------
        Set[Point]
            Set of all bend summits.
        """
        summits = set()
        for bend in self.bends:
            summits.add(bend.get_bend_summit())
        return summits
    
    def get_inflection_pts(self) -> List[Point]:
        """
        Get the inflection points of the bend series.
        
        Returns
        -------
        List[Point]
            List of inflection points.
        """
        return self.inflection_pts
    
    def get_bends(self) -> List[Bend]:
        """
        Get the bends of the series.
        
        Returns
        -------
        List[Bend]
            List of Bend objects.
        """
        return self.bends