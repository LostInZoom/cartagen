from enum import Enum
from typing import List, Tuple, Optional
import numpy as np
from shapely.geometry import LineString, Point
from shapely.ops import substring

def coalescence_splitting(line, width, legibility_factor=1.7, join_style='round', quad_segs=16, debug=False):
    """
    Splits a line into parts when coalescence is detected.

    This algorithm proposed by Mustière :footcite:p:`mustiere:2001-a` :footcite:p:`mustiere:2001-b`
    subdivides the provided line into multiple parts with either
    no coalescence ('none'), or coalescence on one side ('left' or 'right'), or coalescence on both sides ('both'). 

    Parameters
    ----------
    line : LineString
        The line to detect the pastiness from.
    width : float
        The width of the offset used to detect the coalescence.
    legibility_factor : float, optional
        Detected coalescence chunks that are separated by a distance below
        legibility_factor * tolerance are grouped together.
    join_style : str, optional
        The type of caps at the start and end of the provided line.
        Possible values are 'round' or 'flat'.
    quad_segs : int, optional
        The number of point allowed per circle quadrant
        when interpolating points using round method.

    Returns
    -------
    Tuple[List[LineString], List[str]]
        Tuple of two -> A line and its coalescence ('none', 'left', 'right', 'both')

    References
    ----------
    .. footbibliography::
    """
    
    # Validate input
    if line is None or line.is_empty:
        return [], []
    
    if not isinstance(line, LineString):
        raise ValueError("Input must be a Shapely LineString")
    
    coords = np.array(line.coords)
    if len(coords) < 3:
        # Not enough points for coalescence detection
        return [line], ['none']
    
    # Compute offset curves on both sides
    offset_left = _offset_curve(line, width, join_style, quad_segs)
    offset_right = _offset_curve(line, -width, join_style, quad_segs)
    
    # Minimum distance for coalescence detection
    min_distance = width * legibility_factor
    
    # Initialize result containers
    sections = []
    coalescence_types = []
    
    # Track current section being built
    positions_for_current_line = []
    current_coalescence_type = 'none'
    
    # Debug information
    vectors = [] if debug else None
    vectors_boolean = [] if debug else None
    
    # Process each point of the line
    for i, coord in enumerate(coords):
        point = Point(coord)
        coord_array = np.array(coord)
        coalescence_detected = False
        
        # Find nearest points on offset curves
        nearest_left = _nearest_point_on_line(offset_left, point)
        nearest_right = _nearest_point_on_line(offset_right, point)
        
        # Check left side coalescence
        if nearest_left is None or _distance(nearest_left, coord_array) >= min_distance:
            coalescence_detected = True
            
            if debug and nearest_left is not None:
                vectors.append((nearest_left, coord_array))
                vectors_boolean.append(True)
            
            if current_coalescence_type == 'none':
                # Start new section with left coalescence
                if positions_for_current_line:
                    _add_new_section(sections, coalescence_types, positions_for_current_line, coord, current_coalescence_type)
                    positions_for_current_line = []
                current_coalescence_type = 'left'
            elif current_coalescence_type == 'right':
                # Both sides have coalescence
                current_coalescence_type = 'both'
        else:
            if debug and nearest_left is not None:
                vectors.append((nearest_left, coord_array))
                vectors_boolean.append(False)
        
        # Check right side coalescence
        if nearest_right is None or _distance(nearest_right, coord_array) >= min_distance:
            coalescence_detected = True
            
            if debug and nearest_right is not None:
                vectors.append((nearest_right, coord_array))
                vectors_boolean.append(True)
            
            if current_coalescence_type == 'none':
                # Start new section with right coalescence
                if positions_for_current_line:
                    _add_new_section(sections, coalescence_types, positions_for_current_line, coord, current_coalescence_type)
                    positions_for_current_line = []
                current_coalescence_type = 'right'
            elif current_coalescence_type == 'left':
                # Both sides have coalescence
                current_coalescence_type = 'both'
        else:
            if debug and nearest_right is not None:
                vectors.append((nearest_right, coord_array))
                vectors_boolean.append(False)
        
        # If no coalescence detected but we were in a coalescence section
        if not coalescence_detected:
            if current_coalescence_type != 'none':
                _add_new_section(sections, coalescence_types, positions_for_current_line, coord, current_coalescence_type)
                positions_for_current_line = []
                current_coalescence_type = 'none'
        
        positions_for_current_line.append(coord)
    
    # Add final section
    if positions_for_current_line:
        sections.append(LineString(positions_for_current_line))
        coalescence_types.append(current_coalescence_type)
    
    # Merge connected coalesced geometries
    sections, coalescence_types = _merge_coalesced_sections(
        sections, coalescence_types, min_distance
    )
    
    if debug:
        return sections, coalescence_types, vectors, vectors_boolean
    
    return sections, coalescence_types

def _offset_curve(line, distance, join_style, quad_segs):
    """
    Compute offset curve on one side of a line.
    
    Args:
        line: Input LineString
        distance: Offset distance (positive for left, negative for right)
    
    Returns:
        Offset LineString or None if operation fails
    """
    try:
        offset = line.offset_curve(distance, join_style=join_style, quad_segs=quad_segs)
        if offset.is_empty:
            return None
        # Ensure we return a LineString
        if offset.geom_type == 'LineString':
            return offset
        elif offset.geom_type == 'MultiLineString':
            # Return the longest segment if multiple
            return max(offset.geoms, key=lambda x: x.length)
        return None
    except Exception:
        return None


def _nearest_point_on_line(line: Optional[LineString], point: Point) -> Optional[np.ndarray]:
    """
    Find the nearest point on a line to a given point.
    
    Args:
        line: LineString to search on
        point: Reference point
    
    Returns:
        Coordinates of nearest point as numpy array or None
    """
    if line is None or line.is_empty:
        return None
    
    nearest = line.interpolate(line.project(point))
    return np.array([nearest.x, nearest.y])


def _distance(point1: np.ndarray, point2: np.ndarray) -> float:
    """Calculate Euclidean distance between two points."""
    return np.linalg.norm(point1 - point2)

def _add_new_section(
    sections: List[LineString],
    coalescence_types: List[str],
    positions: List[Tuple[float, float]],
    last_position: Tuple[float, float],
    coalescence_type: str
) -> None:
    """
    Add a new section to the results.
    
    Args:
        sections: List to append the new section to
        coalescence_types: List to append the coalescence type to
        positions: List of coordinates for the section
        last_position: Last position to add before creating the section
        coalescence_type: Type of coalescence for this section
    """
    if positions:
        positions_with_last = positions + [last_position]
        new_line = LineString(positions_with_last)
        sections.append(new_line)
        coalescence_types.append(coalescence_type)


def _merge_coalesced_sections(
    sections: List[LineString],
    coalescence_types: List[str],
    min_distance: float
) -> Tuple[List[LineString], List[str]]:
    """
    Merge connected coalesced geometries.
    
    If two coalesced sections are separated by a very short non-coalesced
    section, merge them together.
    
    Args:
        sections: List of LineString sections
        coalescence_types: List of coalescence types
        min_distance: Minimum distance threshold
    
    Returns:
        Merged sections and types
    """
    if len(sections) < 3:
        return sections, coalescence_types
    
    i = 2
    while i < len(sections):
        prev_prev_line = sections[i - 2]
        prev_prev_type = coalescence_types[i - 2]
        
        prev_line = sections[i - 1]
        prev_type = coalescence_types[i - 1]
        
        current_line = sections[i]
        current_type = coalescence_types[i]
        
        # Check if we should merge
        if (prev_prev_type != 'none' and
            current_type != 'none'):
            
            end_point = Point(prev_prev_line.coords[-1])
            start_point = Point(current_line.coords[0])
            
            if end_point.distance(start_point) < min_distance:
                # Merge the three sections
                all_coords = (list(prev_prev_line.coords)[:-1] + 
                            list(prev_line.coords)[:-1] + 
                            list(current_line.coords))
                
                merged_line = LineString(all_coords)
                
                # Determine merged type
                merged_type = prev_prev_type
                if (prev_prev_type != current_type or
                    (prev_type != 'none' and (prev_type != prev_prev_type or prev_type != current_type))):
                    merged_type = 'both'
                
                # Replace the three sections with the merged one
                sections[i - 2] = merged_line
                coalescence_types[i - 2] = merged_type
                
                # Remove the two following sections
                del sections[i - 1:i + 1]
                del coalescence_types[i - 1:i + 1]
                
                # Skip ahead
                i += 1
                continue
        
        i += 1
    
    return sections, coalescence_types