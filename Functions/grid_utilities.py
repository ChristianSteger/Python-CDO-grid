# Load modules
import numpy as np
from shapely.geometry import Polygon


###############################################################################

def gridcoord(x_cent, y_cent):
    """Compute edge coordinates.

    Compute edge coordinates from grid cell centre coordinates

    Parameters
    ----------
    x_cent: array_like
        Array (1-dimensional) with x-coordinates of grid centres [arbitrary]
    y_cent: array_like
        Array (1-dimensional) with y-coordinates of grid centres [arbitrary]

    Returns
    -------
    x_edge: array_like
        Array (1-dimensional) with x-coordinates of grid edges [arbitrary]
    y_edge: array_like
        Array (1-dimensional) with y-coordinates of grid edges [arbitrary]

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    # Check input arguments
    if len(x_cent.shape) != 1 or len(y_cent.shape) != 1:
        raise TypeError("number of dimensions of input arrays is not 1")
    if (np.any(np.diff(np.sign(np.diff(x_cent))) != 0) or
            np.any(np.diff(np.sign(np.diff(y_cent))) != 0)):
        sys.exit("input arrays are not monotonically in- or decreasing")

    # Compute grid spacing if not provided
    dx = np.diff(x_cent).mean()
    dy = np.diff(y_cent).mean()

    # Compute grid coordinates
    x_edge = np.hstack((x_cent[0] - (dx / 2.),
                        x_cent[:-1] + np.diff(x_cent) / 2.,
                        x_cent[-1] + (dx / 2.))).astype(x_cent.dtype)
    y_edge = np.hstack((y_cent[0] - (dy / 2.),
                        y_cent[:-1] + np.diff(y_cent) / 2.,
                        y_cent[-1] + (dy / 2.))).astype(y_cent.dtype)

    return x_edge, y_edge


###############################################################################

def gridframe(x_edge, y_edge, offset=0):
    """Compute frame around a grid.

    Compute frame around a grid with a certain offset from the outer boundary.

    Parameters
    ----------
    x_edge: array_like
        Array (1 or 2-dimensional) with x-coordinates of grid edges [arbitrary]
    y_edge: array_like
        Array (1 or 2-dimensional) with y-coordinates of grid edges [arbitrary]
    offset: int
        offset value from outer boundary [-]

    Returns
    -------
    x_frame: array_like
        Array (1-dimensional) with x-coordinates of frame [arbitrary]
    y_frame: array_like
        Array (1-dimensional) with y-coordinates of frame [arbitrary]

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    # Check input arguments
    if (x_edge.ndim != y_edge.ndim) or (y_edge.ndim > 2):
        raise TypeError("number of dimensions of input arrays is unequal"
                        + " or larger than 2")
    if y_edge.ndim == 1:
        if not ((offset * 2 + 2) <= np.minimum(len(x_edge), len(y_edge))):
            raise TypeError("offset value too large")
    else:
        if x_edge.shape != y_edge.shape:
            raise TypeError("dimension lengths of x_edge and y_edge are "
                            + "unequal")
        if not ((offset * 2 + 2) <= min(x_edge.shape)):
            raise TypeError("offset value too large")

    # Compute frame from 1-dimensional coordinate arrays
    if x_edge.ndim == 1:
        if offset > 0:
            x_edge = x_edge[offset:-offset]
            y_edge = y_edge[offset:-offset]
        x_frame = np.concatenate((x_edge[:-1],
                                  np.repeat(x_edge[-1], len(y_edge) - 1),
                                  x_edge[::-1][:-1],
                                  np.repeat(x_edge[0], len(y_edge) - 1)))
        y_frame = np.concatenate((np.repeat(y_edge[0], len(x_edge) - 1),
                                  y_edge[:-1],
                                  np.repeat(y_edge[-1], len(x_edge) - 1),
                                  y_edge[::-1][:-1]))

    # Compute frame from 2-dimensional coordinate arrays
    else:
        if offset > 0:
            x_edge = x_edge[offset:-offset, offset:-offset]
            y_edge = y_edge[offset:-offset, offset:-offset]
        x_frame = np.hstack((x_edge[0, :],
                             x_edge[1:-1, -1],
                             x_edge[-1, :][::-1],
                             x_edge[1:-1, 0][::-1]))
        y_frame = np.hstack((y_edge[0, :],
                             y_edge[1:-1, -1],
                             y_edge[-1, :][::-1],
                             y_edge[1:-1, 0][::-1]))

    return x_frame, y_frame


###############################################################################

def areagridcell(x_edge, y_edge):
    """Compute area of grid cells.

    Compute area of grid cells on a plane.

    Parameters
    ----------
    x_edge: array_like
        Array (1 or 2-dimensional) with x-coordinates of grid edges [arbitrary]
    y_edge: array_like
        Array (1 or 2-dimensional) with y-coordinates of grid edges [arbitrary]

    Returns
    -------
    area_gridcell: array_like
        Array (2-dimensional) with areas of grid cells [arbitrary]

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    # Check input arguments
    if ((x_edge.ndim not in [1, 2]) or (y_edge.ndim not in [1, 2]) or
            (x_edge.ndim != y_edge.ndim)):
        raise TypeError("input arrays must be both either 1- "
                        + "or 2-dimensional")
    if x_edge.ndim == 2:
        if x_edge.shape != y_edge.shape:
            raise TypeError("shapes of input arrays are inconsistent")

    # Calculate areas for grid cells
    if x_edge.ndim == 1:  # regular grid
        area_gc = np.multiply(np.diff(x_edge).reshape(1, (len(x_edge) - 1)),
                              np.diff(y_edge).reshape((len(y_edge) - 1), 1))
    else:  # potentially irregular grid
        area_gc = np.empty(np.asarray(x_edge.shape) - 1)
        for i in range(0, area_gc.shape[0]):
            for j in range(0, area_gc.shape[1]):
                x_vert = np.array([x_edge[i, j],
                                   x_edge[(i + 1), j],
                                   x_edge[(i + 1), (j + 1)],
                                   x_edge[i, (j + 1)]])
                y_vert = np.array([y_edge[i, j],
                                   y_edge[(i + 1), j],
                                   y_edge[(i + 1), (j + 1)],
                                   y_edge[i, (j + 1)]])
                polygon = Polygon(zip(x_vert, y_vert))
                area_gc[i, j] = polygon.area

    return area_gc


###############################################################################

def gridpolygon(x_edge, y_edge, x_vert, y_vert, agg_gc=np.array([])):
    """Compute intersection areas of polygon and grid.

    Compute area fractions of grid cells that are inside a polygon on a plane.

    Parameters
    ----------
    x_edge: array_like
        Array (2-dimensional) with x-coordinates of grid edges [arbitrary]
    y_edge: array_like
        Array (2-dimensional) with y-coordinates of grid edges [arbitrary]
    x_vert: array_like
        Array (1-dimensional) with x-coordinates of polygon's vertices
        [arbitrary]
    y_vert: array_like
        Array (1-dimensional) with y-coordinates of polygon's vertices
        [arbitrary]
    agg_gc: array_like, optional
        Array with decreasing integers. The values determine the
        aggregation of grid cells for processing

    Returns
    -------
    area_frac: array_like
        Array (2-dimensional) with area fractions [arbitrary]

    Notes
    -----
    Author: Christian Steger (christian.steger@env.ethz.ch)"""

    # Check input arguments
    if len(x_edge.shape) != 2 or len(y_edge.shape) != 2:
        raise TypeError("number of dimensions of x_edge and/or y_edge is "
                        + "not equal 2")
    if x_edge.shape != y_edge.shape:
        raise TypeError("dimension lengths of x_edge and y_edge are unequal")
    if len(x_vert.shape) != 1 or len(y_vert.shape) != 1:
        raise TypeError("number of dimensions of x_vert or y_vert is not "
                        + "equal to 1")
    if len(x_vert) != len(y_vert):
        raise TypeError("length of x_vert and y_vert is unequal")
    if np.any(np.diff(agg_gc) >= 0) or np.any(agg_gc <= 1):
        raise TypeError("agg_gc must be in descending order and contain "
                        + "values greater than one")

    # Allocate arrays
    area_frac = np.zeros((x_edge.shape[0] - 1, x_edge.shape[1] - 1),
                         dtype=np.float32)
    mask_proc = np.zeros((x_edge.shape[0] - 1, x_edge.shape[1] - 1),
                         dtype=bool)

    # Intersect polygon with grid cells
    polygon = Polygon(zip(x_vert, y_vert))
    for i in agg_gc:
        for j in range(0, area_frac.shape[0], i):
            for k in range(0, area_frac.shape[1], i):
                if not np.all(mask_proc[j:(j + i), k:(k + i)]):

                    ind_x = np.minimum((j + 1 + i), area_frac.shape[0])
                    ind_y = np.minimum((k + 1 + i), area_frac.shape[1])
                    x_box = np.array([x_edge[j, k],
                                      x_edge[ind_x, k],
                                      x_edge[ind_x, ind_y],
                                      x_edge[j, ind_y]])
                    y_box = np.array([y_edge[j, k],
                                      y_edge[ind_x, k],
                                      y_edge[ind_x, ind_y],
                                      y_edge[j, ind_y]])
                    box_poly = Polygon(zip(x_box, y_box))
                    
                    if polygon.contains(box_poly):
                        # Polygon contains box
                        mask_proc[j:(j + i), k:(k + i)] = True
                        area_frac[j:(j + i), k:(k + i)] = 1.0
                    elif not polygon.intersects(box_poly):
                        # Polygon does not intersect with box
                        #  (-> entirely outside)
                        mask_proc[j:(j + i), k:(k + i)] = True

    # Process individual (remaining) grid cells
    ind_0, ind_1 = np.where(np.invert(mask_proc))

    for i in range(len(ind_0)):
        
        gc_x = np.array([x_edge[ind_0[i], ind_1[i]],
                         x_edge[(ind_0[i] + 1), ind_1[i]],
                         x_edge[(ind_0[i] + 1), (ind_1[i] + 1)],
                         x_edge[ind_0[i], (ind_1[i] + 1)]])
        gc_y = np.array([y_edge[ind_0[i], ind_1[i]],
                         y_edge[(ind_0[i] + 1), ind_1[i]],
                         y_edge[(ind_0[i] + 1), (ind_1[i] + 1)],
                         y_edge[ind_0[i], (ind_1[i] + 1)]])
        gc_poly = Polygon(zip(gc_x, gc_y))
        area_frac[ind_0[i], ind_1[i]] = \
            polygon.intersection(gc_poly).area / gc_poly.area

    return area_frac
