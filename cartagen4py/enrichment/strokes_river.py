# -*- coding: utf-8 -*-
import math

import numpy as np
from shapely.geometry import LineString, Point, MultiLineString
from cartagen4py.utils.geometry import angle_between_2lines
from shapely import ops
import geopandas as gpd
import networkx as nx

