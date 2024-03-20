g_impurityScale = 1. * g_t

if g_latticeType == 1
  g_impurityShape = [
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  ];
elseif g_latticeType == 2
  g_impurityShape = [
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  ];
elseif g_latticeType == 3
  g_impurityShape = [
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  ];
elseif g_latticeType == 4
  g_impurityShape = [
    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  ];
end
nothing
