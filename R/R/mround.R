# MATLAB's round: halves away from zero (R's round() takes halves to even).
mround <- function(v) sign(v) * floor(abs(v) + 0.5)
