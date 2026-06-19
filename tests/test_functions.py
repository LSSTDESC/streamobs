import numpy as np

from streamobs.functions import (CubicSplineInterpolation,
                                 FileCubicSplineInterpolation,
                                 FileLinearDensityCubicSplineInterpolation,
                                 LinearDensityCubicSplineInterpolation)

SPREAD_NODES = np.array([-13.0, -7.875, -2.75, 2.375, 7.5])

SPREAD_NODE_VALUES = np.array(
    [0.67851425, 0.12068391, 0.19419976, 0.15820543, 0.12868499]
)

SPREAD_INTERP_ARRAY = np.array(
    [
        0.28385327,
        0.18494655,
        0.12147638,
        0.09911923,
        0.10862358,
        0.13679006,
        0.17041929,
        0.19631621,
        0.20546274,
        0.20080746,
        0.18742835,
        0.17040343,
        0.15479583,
        0.1437526,
        0.13663492,
        0.1323599,
        0.12984461,
        np.nan,
        np.nan,
        np.nan,
    ]
)
INTENSITY_NODES = np.array(
    [
        -13.0,
        -11.04761905,
        -10.07142857,
        -9.0952381,
        -6.16666667,
        -5.92261905,
        -5.19047619,
        -4.45833333,
        -4.21428571,
        -3.9702381,
        -0.30952381,
        0.66666667,
        1.64285714,
        3.5952381,
        4.57142857,
        5.54761905,
        7.5,
    ]
)

INTENSITY_NODE_VALUES = np.array(
    [
        3.57488822e-08,
        8.15443503e-02,
        7.43156465e-02,
        8.38736587e-02,
        1.29506112e-01,
        1.48465666e-01,
        1.51251953e-08,
        9.32309494e-02,
        1.91713194e-02,
        4.61099534e-03,
        4.17791354e-02,
        4.89747482e-02,
        1.55361917e-02,
        2.27116146e-02,
        2.45111417e-02,
        6.95831709e-02,
        1.03856772e-07,
    ]
)

LINEAR_DENSITY_INTERP_ARRAY = np.array(
    [
        0.05301498,
        0.03867171,
        0.01956495,
        0.018129,
        0.03598436,
        0.02941,
        0.00440678,
        0.01117077,
        0.01486688,
        0.01937192,
        0.02381247,
        0.00740761,
        0.00514529,
        0.00808078,
        0.01046202,
        0.02456511,
        0.01387564,
        np.nan,
        np.nan,
        np.nan,
    ]
)


def test_CubicSplineInterpolation():
    # using atlas width from patrick_2022.csv
    interp = CubicSplineInterpolation(SPREAD_NODES, SPREAD_NODE_VALUES)

    result = interp(np.linspace(-10, 10, 20))

    np.testing.assert_allclose(result, SPREAD_INTERP_ARRAY)


def test_FileCubicSplineInterpolation():
    # using atlas width from patrick_2022.csv
    filepath = "./data/patrick_2022_splines.csv"
    stream_name = "atlas"
    spline_type = "spread"

    interp = FileCubicSplineInterpolation(
        filepath, spline_type=spline_type, stream_name=stream_name
    )
    result = interp(np.linspace(-10, 10, 20))

    np.testing.assert_allclose(result, SPREAD_INTERP_ARRAY)


def test_LinearDensityCubicSplineInterpolation():
    # using atlas width from patrick_2022.csv

    interp = LinearDensityCubicSplineInterpolation(
        spread_nodes=SPREAD_NODES,
        spread_node_values=SPREAD_NODE_VALUES,
        intensity_nodes=INTENSITY_NODES,
        intensity_node_values=INTENSITY_NODE_VALUES,
    )
    result = interp(np.linspace(-10, 10, 20))
    np.testing.assert_allclose(result, LINEAR_DENSITY_INTERP_ARRAY, atol=1e-7)


def test_FileLinearDensityCubicSplineInterpolation():
    # using atlas width from patrick_2022.csv
    filename = "./data/patrick_2022_splines.csv"
    stream_name = "atlas"
    interp = FileLinearDensityCubicSplineInterpolation(
        filename=filename, stream_name=stream_name
    )

    result = interp(np.linspace(-10, 10, 20))
    np.testing.assert_allclose(result, LINEAR_DENSITY_INTERP_ARRAY, atol=1e-7)
