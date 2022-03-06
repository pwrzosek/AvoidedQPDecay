# 1 -> 16 sites 2D
# 2 -> 18 sites 2D
# 3 -> 20 sites 2D
# 4 -> 26 sites 2D

if g_latticeType == 1
    g_latticeTypeName = "2D_square"
    g_lattice = [
      1  2  3  0  5  6  7  4  9 10 11  8 13 14 15 12
      4  5  6  7  8  9 10 11 12 13 14 15  0  1  2  3
      3  0  1  2  7  4  5  6 11  8  9 10 15 12 13 14
     12 13 14 15  0  1  2  3  4  5  6  7  8  9 10 11
    ]
    g_secondNeighbours = [
      5  6  7  4  9 10 11  8 13 14 15 12  1  2  3  0
      7  4  5  6 11  8  9 10 15 12 13 14  3  0  1  2
     15 12 13 14  3  0  1  2  7  4  5  6 11  8  9 10
     13 14 15 12  1  2  3  0  5  6  7  4  9 10 11  8
    ]
    g_cartesian = [
      (0, 0) (1, 0) (2, 0) (3, 0) (0, 1) (1, 1) (2, 1) (3, 1) (0, 2) (1, 2) (2, 2) (3, 2) (0, 3) (1, 3) (2, 3) (3, 3)
    ]
    g_unitA = 4
    g_unitB = 4
    g_eX = (1, 0)
    g_eY = (0, 1)
    g_kPath = [
      1  2  3  7 11  6  1
    ]
    g_kNames = [
      "G = (0, 0)", "(Pi/2, 0)", "X = (Pi, 0)", "(Pi, Pi/2)", "M = (Pi, Pi)", "S = (Pi/2, Pi/2)", "G = (0, 0)"
    ]
elseif g_latticeType == 2
    g_latticeTypeName = "2D_square"
    g_lattice = [
      1  2  3  4  5  0  7  8  9 10 11  6 13 14 15 16 17 12
      6  7  8  9 10 11 12 13 14 15 16 17  3  4  5  0  1  2
      5  0  1  2  3  4 11  6  7  8  9 10 17 12 13 14 15 16
     15 16 17 12 13 14  0  1  2  3  4  5  6  7  8  9 10 11
    ]
    g_secondNeighbours = [
      7  8  9 10 11  6 13 14 15 16 17 12  4  5  0  1  2  3
     11  6  7  8  9 10 17 12 13 14 15 16  2  3  4  5  0  1
     14 15 16 17 12 13  5  0  1  2  3  4 11  6  7  8  9 10
     16 17 12 13 14 15  1  2  3  4  5  0  7  8  9 10 11  6
    ]
    g_cartesian = [
      (0, 0) (1, 0) (2, 0) (3, 0) (4, 0) (5, 0) (0, 1) (1, 1) (2, 1) (3, 1) (4, 1) (5, 1) (0, 2) (1, 2) (2, 2) (3, 2) (4, 2) (5, 2)
    ]
    g_unitA = 6
    g_unitB = 6
    g_eX = (1, 1)
    g_eY = (-1, 1)
    g_kPath = [
      1  8  9  4  3  2  1
    ]
    g_kNames = [
      # TODO
    ]
elseif g_latticeType == 3
    g_latticeTypeName = "2D_square"
    g_lattice = [
      1  2  3  8  5  6  7 12  9 10 11 16 13 14 15 18 17  0 19  4
      4  5  6  7  8  9 10 11 12 13 14 15 16 17  0  1 18 19  2  3
     17  0  1  2 19  4  5  6  3  8  9 10  7 12 13 14 11 16 15 18
     14 15 18 19  0  1  2  3  4  5  6  7  8  9 10 11 12 13 16 17
    ]
    g_secondNeighbours = [
      5  6  7 12  9 10 11 16 13 14 15 18 17  0  1  2 19  4  3  8
     19  4  5  6  3  8  9 10  7 12 13 14 11 16 17  0 15 18  1  2
     13 14 15 18 17  0  1  2 19  4  5  6  3  8  9 10  7 12 11 16
     15 18 19  4  1  2  3  8  5  6  7 12  9 10 11 16 13 14 17  0
    ]
    g_cartesian = [
      (0, 0) (1, 0) (2, 0) (3, 0) (0, 1) (1, 1) (2, 1) (3, 1) (0, 2) (1, 2) (2, 2) (3, 2) (0, 3) (1, 3) (2, 3) (3, 3) (0, 4) (1, 4) (0, 5) (1, 5)
    ]
    g_unitA = 10
    g_unitB = 10
    g_eX = (2, 1)
    g_eY = (-1, 2)
    g_kPath = [
      # TODO
    ]
    g_kNames = [
      # TODO
    ]
elseif g_latticeType == 4
    g_latticeTypeName = "2D_square"
    g_lattice = [
      1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25  0
      5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25  0  1  2  3  4
     25  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
     21 22 23 24 25  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    ]
    g_secondNeighbours = [
      6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25  0  1  2  3  4  5
      4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25  0  1  2  3
     20 21 22 23 24 25  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
     22 23 24 25  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
    ]
    g_cartesian = [
      (0, 0) (1, 0) (2, 0) (3, 0) (4, 0) (0, 1) (1, 1) (2, 1) (3, 1) (4, 1) (0, 2) (1, 2) (2, 2) (3, 2) (4, 2) (0, 3) (1, 3) (2, 3) (3, 3) (4, 3) (0, 4) (1, 4) (2, 4) (3, 4) (4, 4) (0, 5)
    ]
    g_unitA = 26
    g_unitB = 26
    g_eX = (5, 1)
    g_eY = (-1, 5)
    g_kPath = [
      # TODO
    ]
    g_kNames = [
      # TODO
    ]
# Rectangular lattice is possible
# elseif g_latticeType == 5
#     g_latticeTypeName = "2D_rect_4x6"
#     g_lattice = [
#       1  2  3  4  5  0  7  8  9 10 11  6 13 14 15 16 17 12 19 20 21 22 23 18
#       6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23  0  1  2  3  4  5
#       5  0  1  2  3  4 11  6  7  8  9 10 17 12 13 14 15 16 23 18 19 20 21 22
#      18 19 20 21 22 23  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
#     ]
#     g_cartesian = [
#     ]
end

g_bravais = [
  (r[1] .* g_eX) .+ (r[2] .* g_eY) for r in g_cartesian
] # coordinates of eX and eY are given in AB coordinates system -> for each r from XY we get AB coordinates of r
g_brillouin = [
  r .* (2pi / g_unitA, 2pi / g_unitB) for r in g_bravais
] # have to be divided by corresponding units

# only if all k needed
# g_kPath = [
#   k for k in 1 : length(g_brillouin)
# ]
# g_kNames = [
#   string(round.(g_brillouin[k] ./ pi, digits = 2), " Pi") for k in g_kPath
# ]
nothing

# Comment on XY and AB coordinates :
# Consider A and B - orthogonal vectors defining periodicity of the system (due to pbc) - dependent on chosen lattice AB may be rotated with respect to XY system
# g_unitA (g_unitB) - length of A (B) in rotated coordinates system (note: e.g. for 20 site lattice g_unitA = 10, but in XY the same vector has length sqrt(5))
# g_eX (g_eY) - versor in XY coordinates system (also smallest vector with integer coordinates in AB coordinates system)


# Order of directions:
# 1 ->
# 2 ^
# 3 <-
# 4 v
# Order of second neighbours -> the same as order of quandrants in cartesian coordinate system
