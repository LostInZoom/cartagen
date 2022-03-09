from shapely.geometry import Point
from evaluation.constraint_satisfaction import ConstraintSatisfaction

class ConstraintMonitor:

    def get_geometry(self) -> Point:
        """Get the point geometry locating the constraint monitor."""
        pass

    def get_importance(self) -> int:
        """Get the importance of the constraint monitor."""
        pass

    def get_satisfaction(self) -> ConstraintSatisfaction:
        """Get the current satisfaction of the constraint."""
        pass

    def get_subject(self):
        """Get the subject of the constraint monitor: it can be a geometry, or a group of geometries."""
        pass